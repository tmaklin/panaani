// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use std::collections::HashMap;

use itertools::Itertools;
use log::info;
use rand::Rng;

pub mod build;
pub mod clust;
pub mod dist;

#[derive(Clone)]
pub struct PanaaniParams {
    pub batch_step: usize,
    pub batch_step_strategy: String,
    pub max_iters: usize,
}

impl Default for PanaaniParams {
    fn default() -> PanaaniParams {
        PanaaniParams {
	    batch_step: 50,
	    batch_step_strategy: "linear".to_string(),
	    max_iters: 10,
	}
    }
}

fn match_clustering_results(fastx_files: &[String],
			    old_clusters: &[(String, String)],
			    res_clusters: &[usize],
			    out_prefix: &String)
-> Vec<(String, String)> {

    let mut old_cluster_to_new_cluster: HashMap<&String, usize> = HashMap::new();
    fastx_files
        .iter()
        .sorted()
        .zip(res_clusters.iter())
        .for_each(|x| {
            old_cluster_to_new_cluster.insert(x.0, x.1.clone());
        });

    let new_clusters: Vec<(String, String)> = old_clusters
        .iter()
        .map(|x| {
            (
                x.0.clone(),
                out_prefix.to_owned()
                    + &old_cluster_to_new_cluster.get(&x.1).unwrap_or_else(|| { panic!("A fasta/fastq failed skani sketching!\nCheck log for records containing the message: 'WARN - File <path> is not a valid fasta/fastq file'.") } ).to_string()
                    + ".dbg.fasta",
            )
        })
        .collect();

    return new_clusters;
}

pub fn dereplicate_iter(
    old_clusters: &[(String, String)],
    out_prefix: &String,
    skani_params: Option<dist::SkaniParams>,
    kodama_params: Option<clust::KodamaParams>,
    ggcat_params: Option<build::GGCATParams>,
) -> Vec<(String, String)> {
    info!("Calculating ANIs...");
    let fastx_files = old_clusters.iter().map(|x| x.1.clone()).unique().collect();
    let ani_result = dist::ani_from_fastx_files(
        &fastx_files,
        &skani_params.unwrap_or(dist::SkaniParams::default()),
    );

    info!("Building dendrogram...");
    let clusters = clust::single_linkage_cluster(
        &ani_result,
        kodama_params.unwrap_or(clust::KodamaParams::default()),
    );

    let mut new_clusters: Vec<(String, String)> = match_clustering_results(&fastx_files, old_clusters, &clusters, out_prefix);

    info!("Building pangenome graphs...");
    build::build_pangenome_representations(
        &mut new_clusters,
        &ggcat_params.unwrap_or(build::GGCATParams::default()),
    );

    return new_clusters;
}

pub fn dereplicate(
    seq_files: &[String],
    initial_clusters: &[String],
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

    let mut batch_size = my_params.batch_step;
    let mut n_remaining: usize = seq_files.len();
    while batch_size < n_remaining && iter < my_params.max_iters {
        let mut rng = rand::thread_rng();

        // horrible hack to use random file names within each batch
        iter_inputs = iter_inputs
            .chunks(batch_size)
            .map(|x| {
                dereplicate_iter(
                    &Vec::from(x),
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
	match my_params.batch_step_strategy.as_str() {
	    "linear" => batch_size += my_params.batch_step,
	    "double" => batch_size *= 2,
	    &_ => batch_size += my_params.batch_step
	}
    }

    let final_clusters = dereplicate_iter(
        &iter_inputs,
        &"panANI-".to_string(),
        skani_params,
        kodama_params,
        ggcat_params,
    );

    return final_clusters;
}
