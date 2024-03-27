// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use std::cmp::Ordering;
use std::collections::HashMap;

use itertools::Itertools;
use log::info;
use log::trace;
use rand::Rng;

pub mod build;
pub mod clust;
pub mod dist;

#[derive(Clone)]
pub struct PanaaniParams {
    pub batch_step: usize,
    pub batch_step_strategy: String,
    pub max_iters: usize,
    pub temp_dir: String,
    pub guided: bool,
    pub external_clustering: Option<Vec<String>>,
    pub initial_batches: Option<Vec<String>>,
}

impl Default for PanaaniParams {
    fn default() -> PanaaniParams {
        PanaaniParams {
	    batch_step: 50,
	    batch_step_strategy: "linear".to_string(),
	    max_iters: 10,
	    temp_dir: "./".to_string(),
	    guided: false,
	    external_clustering: None,
	    initial_batches: None,
        }
    }
}

pub fn match_clustering_results(
    fastx_files: &[String],
    old_clusters: &[String],
    hclust_res: &[usize],
    out_prefix: &String,
) -> Vec<String> {
    let mut old_cluster_to_new_cluster: HashMap<&String, usize> = HashMap::new();
    fastx_files
        .iter()
        .sorted()
        .zip(hclust_res.iter())
        .for_each(|x| {
            old_cluster_to_new_cluster.insert(x.0, x.1.clone());
        });

    let new_clusters: Vec<String> = old_clusters
        .iter()
        .map(|x| {
            out_prefix.to_owned()
                + &old_cluster_to_new_cluster.get(&x).unwrap_or_else(|| { panic!("A fasta/fastq failed skani sketching!\nCheck log for records containing the message: 'WARN - File <path> is not a valid fasta/fastq file'.") } ).to_string()
                + ".dbg.fasta"
        })
        .collect();

    return new_clusters;
}

pub fn assign_seqs(seqs: &[String], clusters: &[String]) -> HashMap::<String, Vec<String>> {
    // Create hashmap mapping each cluster name to the sequences assigned to it
    let mut cluster_contents: HashMap<String, Vec<String>> = HashMap::new();
    seqs
	.iter()
	.zip(clusters)
	.for_each(|x|
		  {
		      if !cluster_contents.contains_key(x.1) {
			  cluster_contents.insert(x.1.clone(), Vec::new());
		      }
		      cluster_contents.get_mut(x.1).unwrap().push(x.0.clone());
		  }
	);
    return cluster_contents;
}

pub fn dereplicate_iter(
    prev_assignments: &HashMap<String, Vec<String>>,
    out_prefix: &String,
    skani_params: &Option<dist::SkaniParams>,
    kodama_params: &Option<clust::KodamaParams>,
    ggcat_params: &Option<build::GGCATParams>,
) -> HashMap<String, Vec<String>> {
    let seq_files = prev_assignments.iter().map(|x| x.1.clone()).flatten().collect::<Vec<String>>();
    let old_clusters = prev_assignments.iter().map(|x| vec![x.0.clone(); x.1.len()]).flatten().collect::<Vec<String>>();

    info!("Calculating ANIs...");
    let fastx_files = old_clusters.iter().cloned().unique().collect();
    let ani_result = dist::ani_from_fastx_files(
        &fastx_files,
        skani_params,
    );

    info!("Building dendrogram...");
    let hclust_res = clust::single_linkage_cluster(
        &ani_result,
        kodama_params,
    );

    let mut new_clusters: Vec<String> = match_clustering_results(&fastx_files, &old_clusters, &hclust_res, out_prefix);
    let mut new_assignments = assign_seqs(&seq_files, &new_clusters);
    // Singleton clusters should have the same name as in the previous round
    // dumb hack
    seq_files
	.iter()
	.zip(new_clusters.iter_mut())
	.for_each(|x| {
	    if new_assignments.get(x.1).unwrap().len() == 1 {
		*x.1 = x.0.clone();
	    }
	});
    new_assignments = assign_seqs(&seq_files, &new_clusters);

    info!("Building pangenome graphs...");
    build::build_pangenome_representations(
	&new_assignments,
        ggcat_params,
    );

    return new_assignments;
}

fn guide_batching(seq_files: &[String], kodama_params: &Option<clust::KodamaParams>) -> Vec<String> {
    let guide_params = dist::SkaniParams {
        kmer_subsampling_rate: 2500,
        marker_compression_factor: 2500,
        clip_tails: true,
        ..Default::default()
    };

    let fastx_files: Vec<String> = seq_files.iter().cloned().collect();
    let ani_result = dist::ani_from_fastx_files(
        &fastx_files,
        &Some(guide_params),
    );
    let hclust_res = clust::single_linkage_cluster(
        &ani_result,
        kodama_params,
    );

    let res = fastx_files
	.iter()
	.zip(hclust_res)
        .sorted_by(|k1, k2| match k1.1.cmp(&k2.1) {
            Ordering::Equal => k1.0.cmp(&k2.0),
            other => other,
        })
	.map(|x| x.0.clone())
	.collect();
    return res;
}

pub fn dereplicate(
    seq_files: &[String],
    dereplicate_params: &Option<PanaaniParams>,
    skani_params: &Option<dist::SkaniParams>,
    kodama_params: &Option<clust::KodamaParams>,
    ggcat_params: &Option<build::GGCATParams>,
) -> Vec<(String, String)> {
    trace!("Dereplicate input contains {} sequences in {} clusters", seq_files.len(), seq_files.iter().unique().collect::<Vec<&String>>().len());
    let my_params = dereplicate_params.clone().unwrap_or(PanaaniParams::default());

    // Create hashmap mapping each cluster name to the sequences assigned to it
    let mut cluster_contents = assign_seqs(seq_files, &my_params.external_clustering.unwrap_or(seq_files.to_vec()));

    let mut iter: usize = 0;
    let mut batch_size = my_params.batch_step;
    let mut n_remaining: usize = cluster_contents.len();

    while batch_size < n_remaining && iter < my_params.max_iters {
	info!("Iteration {} processing {} sequences in batches of {}...", iter + 1, n_remaining, batch_size);
        let mut rng = rand::thread_rng();

	let batch_assignments: Vec<String> = if iter == 0 && my_params.initial_batches.is_some() {
	    my_params.initial_batches.as_ref().unwrap().clone()
	} else if my_params.guided {
	    let current_clusters: Vec<String> = cluster_contents.iter().map(|x| x.0.clone()).collect();
	    guide_batching(&current_clusters, kodama_params)
	} else {
	    cluster_contents.iter().map(|x| x.0.clone()).collect()
	};

	// horrible hack to use random file names within each batch
        let new_clusters: Vec<HashMap<String, Vec<String>>> = batch_assignments
            .chunks(batch_size)
            .map(|x| {
		let mut batch_inputs: HashMap<String, Vec<String>> = HashMap::new();
		x.iter().for_each(|y| { batch_inputs.insert(y.clone(), cluster_contents.get(y).unwrap().clone()); });
                dereplicate_iter(
		    &batch_inputs,
                    &(my_params.temp_dir.to_string() + "/" + &iter.to_string() + "_" + &(rng.gen::<u64>() as u64).to_string() + "-"),
                    skani_params,
                    kodama_params,
                    ggcat_params,
                )
            })
            .collect();

	cluster_contents = assign_seqs(&new_clusters.iter().map(|x| x.iter().map(|y| y.1.clone()).flatten()).flatten().collect::<Vec<String>>(),
				       &new_clusters.iter().map(|x| x.iter().map(|y| vec![y.0.clone(); y.1.len()]).flatten()).flatten().collect::<Vec<String>>());

	n_remaining = cluster_contents.len();
        iter += 1;
        match my_params.batch_step_strategy.as_str() {
            "linear" => batch_size += my_params.batch_step,
            "double" => batch_size *= 2,
            &_ => batch_size += my_params.batch_step,
        }

	// If n_remaining/batch_size == 1 increase batch size so that
	// the last chunk contains more than a single sequence.
	while n_remaining % batch_size == 1 {
	    batch_size += 1;
	}
    }
    info!("Final iteration processing {} sequences...", n_remaining);

    let final_clusters = dereplicate_iter(
	&cluster_contents,
        &"panANI-".to_string(),
        skani_params,
        kodama_params,
        ggcat_params,
    );

    return final_clusters
	.iter()
	.map(|x| x.1.iter().cloned().zip(vec![x.0.clone(); x.1.len()]).collect::<Vec<(String, String)>>())
	.flatten()
        .sorted_by(|k1, k2| match k1.1.cmp(&k2.1) {
            Ordering::Equal => k1.0.cmp(&k2.0),
            other => other,
        })
	.collect();
}
