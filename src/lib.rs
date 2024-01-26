use std::collections::HashMap;
use std::path::PathBuf;

use itertools::Itertools;
use rand::Rng;

use ggcat_api::{GGCATConfig,GGCATInstance};

mod build;
pub mod clust;
mod dist;

pub fn dereplicate_iter(old_clusters: Vec<(String, String)>, out_prefix: &String, instance: &ggcat_api::GGCATInstance) -> Vec<(String, String)> {
    println!("Calculating ANIs...");
    let fastx_files = old_clusters.iter().map(|x| x.1.clone()).unique().collect();
    let ani_result = dist::ani_from_fastx_files(&fastx_files, &dist::SkaniParams::default());

    println!("Building dendrogram...");
    let clusters = clust::single_linkage_cluster(&ani_result, fastx_files.len());

    let mut old_cluster_to_new_cluster: HashMap<String, usize> = HashMap::new();
    fastx_files.iter().sorted().zip(clusters.iter()).for_each(|x| { old_cluster_to_new_cluster.insert(x.0.clone(), x.1.clone()); });
    let new_clusters: Vec<(String, String)> = old_clusters
	.iter()
	.map(|x| (x.0.clone(), out_prefix.to_owned() + &old_cluster_to_new_cluster.get(&x.1).unwrap().to_string() + ".dbg.fasta"))
	.collect();

    println!("Building pangenome graphs...");
    build::build_pangenome_representations(&new_clusters, instance);

    return new_clusters;
}

pub fn dereplicate(seq_files: &Vec<String>, initial_clusters: &Vec<String>, batch_step: &usize) -> Vec<(String, String)> {
    let instance = GGCATInstance::create(GGCATConfig {
	temp_dir: Some(PathBuf::from("/tmp")),
	memory: 8.0,
	prefer_memory: true,
	total_threads_count: 4,
	intermediate_compression_level: None,
	stats_file: None,
    });

    let mut iter = 0;
    let mut iter_inputs: Vec<(String, String)> = seq_files.iter().cloned().zip(initial_clusters.iter().cloned()).collect();

    let mut n_remaining = seq_files.len(); while (iter + 1)*batch_step
    < n_remaining { let mut rng = rand::thread_rng();

	// horrible hack to use random file names within each batch
	iter_inputs = iter_inputs
	    .chunks((iter + 1)*batch_step)
	    .map(|x| dereplicate_iter(Vec::from(x), &(iter.to_string() + "_" + &(rng.gen::<u64>() as u64).to_string() + "-" ), instance))
	    .flatten()
	    .collect();

	n_remaining = iter_inputs.iter().map(|x| x.1.clone()).unique().collect::<Vec<String>>().len();
	iter += 1;
    }

    let final_clusters = dereplicate_iter(iter_inputs, &"panANI-".to_string(), instance);

    return final_clusters;
}
