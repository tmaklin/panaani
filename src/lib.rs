use std::collections::HashMap;

use itertools::Itertools;
use rand::Rng;

pub mod build;
pub mod clust;
mod dist;

#[derive(Clone)]
pub struct PanaaniParams {
    pub batch_step: usize,
    pub ani_threshold: f32,
}

impl Default for PanaaniParams {
    fn default() -> PanaaniParams {
        PanaaniParams {
            batch_step: 50,
            ani_threshold: 0.97,
        }
    }
}

pub fn dereplicate_iter(
    old_clusters: Vec<(String, String)>,
    out_prefix: &String,
    skani_params: Option<dist::SkaniParams>,
    kodama_params: Option<clust::KodamaParams>,
    ggcat_params: Option<build::GGCATParams>,
) -> Vec<(String, String)> {
    println!("Calculating ANIs...");
    let fastx_files = old_clusters.iter().map(|x| x.1.clone()).unique().collect();
    let ani_result = dist::ani_from_fastx_files(
        &fastx_files,
        &skani_params.unwrap_or(dist::SkaniParams::default()),
    );

    println!("Building dendrogram...");
    let clusters = clust::single_linkage_cluster(
        &ani_result,
        fastx_files.len(),
        kodama_params.unwrap_or(clust::KodamaParams::default()),
    );

    let mut old_cluster_to_new_cluster: HashMap<String, usize> = HashMap::new();
    fastx_files
        .iter()
        .sorted()
        .zip(clusters.iter())
        .for_each(|x| {
            old_cluster_to_new_cluster.insert(x.0.clone(), x.1.clone());
        });
    let new_clusters: Vec<(String, String)> = old_clusters
        .iter()
        .map(|x| {
            (
                x.0.clone(),
                out_prefix.to_owned()
                    + &old_cluster_to_new_cluster.get(&x.1).unwrap().to_string()
                    + ".dbg.fasta",
            )
        })
        .collect();

    println!("Building pangenome graphs...");
    build::build_pangenome_representations(
        &new_clusters,
        &ggcat_params.unwrap_or(build::GGCATParams::default()),
    );

    return new_clusters;
}

pub fn dereplicate(
    seq_files: &Vec<String>,
    initial_clusters: &Vec<String>,
    dereplicate_params: Option<PanaaniParams>,
    skani_params: Option<dist::SkaniParams>,
    kodama_params: Option<clust::KodamaParams>,
    ggcat_params: Option<build::GGCATParams>,
) -> Vec<(String, String)> {
    let my_params = dereplicate_params.unwrap_or(PanaaniParams::default());

    let mut iter: usize = 0;
    let mut iter_inputs: Vec<(String, String)> = seq_files
        .iter()
        .cloned()
        .zip(initial_clusters.iter().cloned())
        .collect();

    let mut n_remaining: usize = seq_files.len();
    while (iter + 1) * my_params.batch_step < n_remaining {
        let mut rng = rand::thread_rng();

        // horrible hack to use random file names within each batch
        iter_inputs = iter_inputs
            .chunks((iter + 1) * my_params.batch_step)
            .map(|x| {
                dereplicate_iter(
                    Vec::from(x),
                    &(iter.to_string() + "_" + &(rng.gen::<u64>() as u64).to_string() + "-"),
                    skani_params.clone(),
                    kodama_params.clone(),
                    ggcat_params.clone(),
                )
            })
            .flatten()
            .collect();

        n_remaining = iter_inputs
            .iter()
            .map(|x| x.1.clone())
            .unique()
            .collect::<Vec<String>>()
            .len();
        iter += 1;
    }

    let final_clusters = dereplicate_iter(
        iter_inputs,
        &"panANI-".to_string(),
        skani_params,
        kodama_params,
        ggcat_params,
    );

    return final_clusters;
}
