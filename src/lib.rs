use std::path::PathBuf;
use std::sync::mpsc::channel;
use std::collections::HashMap;

use itertools::Itertools;
use rand::Rng;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;

use ggcat_api::{GGCATConfig,GGCATInstance};

pub fn compare_fastx_files(reference: &String, query: &String,
		       sketch_params: &skani::params::SketchParams,
		       cmd_params: &skani::params::CommandParams) -> skani::types::AniEstResult{
    let sketches = skani::file_io::fastx_to_sketches(&vec![reference.clone(), query.clone()], &sketch_params, true);
    let adjust_ani = skani::regression::get_model(sketch_params.c, false);
    let map_params = skani::chain::map_params_from_sketch(sketches.first().unwrap(), false, &cmd_params, &adjust_ani);
    return skani::chain::chain_seeds(sketches.first().unwrap(), sketches.last().unwrap(), map_params);
}

pub fn ani_from_fastx_files(fastx_files: &Vec<String>) -> Vec<skani::types::AniEstResult>{
    // skani parameters
    let m = 1000;
    let c = 30;
    let k = 15;
    let sketch_params = skani::params::SketchParams::new(m, c, k, false, false);
    let cmd_params = skani::params::CommandParams {
        screen: false,
        screen_val: 0.00,
        mode: skani::params::Mode::Dist,
        out_file_name: "".to_string(),
        ref_files: vec![],
        query_files: vec![],
        refs_are_sketch: false,
        queries_are_sketch: false,
        robust: false,
        median: false,
        sparse: false,
        full_matrix: false,
        max_results: 10000000,
        individual_contig_q: false,
        individual_contig_r: false,
        min_aligned_frac: 0.0,
        keep_refs: false,
        est_ci: false,
        learned_ani: false,
        detailed_out: false,
	rescue_small: false,
	distance: true,
    };

    let (sender, receiver) = channel();
    fastx_files.iter().cloned().combinations(2).par_bridge().for_each_with(sender, |s, pair| {
	let _ = s.send(compare_fastx_files(pair.first().unwrap(),
					   pair.last().unwrap(),
					   &sketch_params, &cmd_params));
    });
    let mut ani_result: Vec<skani::types::AniEstResult> = receiver.iter().collect();

    // Ensure output order is same regardless of parallelization
    ani_result.sort_by_key(|k| (k.ref_file.clone(), k.query_file.clone()));
    return ani_result;
}

pub fn cut_dendrogram(dendr: &kodama::Dendrogram<f32>, height: f32) -> Vec<usize> {
    let cutoff = 1.0 - height;
    let num_seqs = dendr.observations();
    let num_nodes = 2 * num_seqs - 1;

    let mut num_groups = 0;
    let mut membership = vec![None; num_nodes];

    for (cluster_index, step) in dendr.steps().iter().enumerate().rev() {
	let cluster = cluster_index + num_seqs;
	if step.dissimilarity <= cutoff {
	    if membership[cluster].is_none() {
		membership[cluster] = Some(num_groups);
		num_groups += 1;
	    }

	    membership[step.cluster1] = membership[cluster];
	    membership[step.cluster2] = membership[cluster];
	}
    }

    let mut groups = Vec::with_capacity(num_seqs);
    for group in membership.into_iter().take(num_seqs) {
	if let Some(group) = group {
	    groups.push(group);
	} else {
	    groups.push(num_groups);
	    num_groups += 1;
	}
    }

    return groups;
}

pub fn single_linkage_cluster(ani_result: &Vec<skani::types::AniEstResult>, num_seqs: usize) -> Vec<usize> {
    let mut ani: Vec<f32> = ani_result.into_iter().map(|x| if x.ani > 0.0 && x.ani < 1.0 && !x.ani.is_nan() { 1.0 - x.ani } else { 1.0 }).collect();
    let dend = kodama::linkage(&mut ani, num_seqs, kodama::Method::Single);

    return cut_dendrogram(&dend, 0.97);
}

pub fn single_linkage_cluster2(ani_result: &Vec<(String, String, f32, f32, f32)>, num_seqs: usize) -> Vec<usize> {
    let mut ani: Vec<f32> = ani_result.into_iter().map(|x| if x.2 > 0.0 && x.2 < 1.0 && !x.2.is_nan() { 1.0 - x.2 } else { 1.0 }).collect();
    let dend = kodama::linkage(&mut ani, num_seqs, kodama::Method::Single);

    return cut_dendrogram(&dend, 0.97);
}

pub fn build_pangenome_graph(inputs: Vec<ggcat_api::GeneralSequenceBlockData>,
			     input_seq_names: &Vec<String>,
			     prefix: &String, instance: &ggcat_api::GGCATInstance) {

    let kmer_size = 51;
    let min_multiplicity = 1;


    let graph_file = PathBuf::from(prefix.to_owned());
    instance.build_graph(
	inputs,
	graph_file,
	Some(input_seq_names),
	kmer_size,
	4,
	false,
	None,
	false,
	min_multiplicity,
	ggcat_api::ExtraElaboration::GreedyMatchtigs,
    );
}

pub fn open_ggcat_inputs(seq_files: &Vec<String>) -> Vec<ggcat_api::GeneralSequenceBlockData> {
    let mut ggcat_inputs: Vec<ggcat_api::GeneralSequenceBlockData> = Vec::new();
    for file in seq_files {
	ggcat_inputs.push(ggcat_api::GeneralSequenceBlockData::FASTA((
	    PathBuf::from(file),
	    None,
	)));
    }
    return ggcat_inputs;
}

pub fn build_pangenome_representations(seq_files: &Vec<(String, String)>,
				       instance: &ggcat_api::GGCATInstance) {

    let mut files_in_cluster: HashMap<String, Vec<String>> = HashMap::new();

    for val in seq_files.iter() {
	if files_in_cluster.contains_key(&val.1) {
	    files_in_cluster.get_mut(&val.1).unwrap().push(val.0.clone());
	} else {
	    files_in_cluster.insert(val.1.clone(), vec![val.0.clone()]);
	}
    }

    for (graph_name, files) in files_in_cluster {
	println!("Building graph {}...", graph_name);
	let mut ggcat_input_names: Vec<String> = Vec::new();
	for file in files {
	    ggcat_input_names.push(file);
	}

	let ggcat_inputs = open_ggcat_inputs(&ggcat_input_names);
	build_pangenome_graph(ggcat_inputs, &ggcat_input_names, &graph_name, instance);
    }

}

pub fn dereplicate_iter(old_clusters: Vec<(String, String)>, out_prefix: &String, instance: &ggcat_api::GGCATInstance) -> Vec<(String, String)> {
    println!("Calculating ANIs...");
    let fastx_files = old_clusters.iter().map(|x| x.1.clone()).unique().collect();
    let ani_result = ani_from_fastx_files(&fastx_files);

    println!("Building dendrogram...");
    let clusters = single_linkage_cluster(&ani_result, fastx_files.len());

    let mut old_cluster_to_new_cluster: HashMap<String, usize> = HashMap::new();
    fastx_files.iter().sorted().zip(clusters.iter()).for_each(|x| { old_cluster_to_new_cluster.insert(x.0.clone(), x.1.clone()); });
    let new_clusters: Vec<(String, String)> = old_clusters
	.iter()
	.map(|x| (x.0.clone(), out_prefix.to_owned() + &old_cluster_to_new_cluster.get(&x.1).unwrap().to_string() + ".dbg.fasta"))
	.collect();

    println!("Building pangenome graphs...");
    build_pangenome_representations(&new_clusters, instance);

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

    let mut n_remaining = seq_files.len();
    while (iter + 1)*batch_step < n_remaining {
	let mut rng = rand::thread_rng();

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
