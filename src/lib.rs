use std::path::PathBuf;
use std::sync::mpsc::channel;

use itertools::Itertools;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;

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
        min_aligned_frac: 0.15,
        keep_refs: false,
        est_ci: false,
        learned_ani: true,
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

pub fn cut_dendrogram(dendr: &kodama::Dendrogram<f32>, height: f32) -> Vec<Vec<usize>> {
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

    let mut seqs_by_group = vec![Vec::new(); num_groups];
    for (seq_index, group) in groups.iter().enumerate() {
	seqs_by_group[*group].push(seq_index);
    }

    return seqs_by_group;
}

pub fn single_linkage_cluster(ani_result: &Vec<skani::types::AniEstResult>, num_seqs: usize) -> Vec<Vec<usize>> {
    let mut ani: Vec<f32> = ani_result.into_iter().map(|x| if x.ani > 0.0 && x.ani < 1.0 && !x.ani.is_nan() { 1.0 - x.ani } else { 1.0 }).collect();
    let dend = kodama::linkage(&mut ani, num_seqs, kodama::Method::Single);

    return cut_dendrogram(&dend, 0.95);
}

pub fn single_linkage_cluster2(ani_result: &Vec<(String, String, f32, f32, f32)>, num_seqs: usize) -> Vec<Vec<usize>> {
    let mut ani: Vec<f32> = ani_result.into_iter().map(|x| if x.2 > 0.0 && x.2 < 1.0 && !x.2.is_nan() { 1.0 - x.2 } else { 1.0 }).collect();
    let dend = kodama::linkage(&mut ani, num_seqs, kodama::Method::Single);

    return cut_dendrogram(&dend, 0.95);
}

pub fn build_pangenome_graph(inputs: Vec<ggcat_api::GeneralSequenceBlockData>,
			     input_seq_names: &Vec<String>,
			     prefix: &String, instance: &ggcat_api::GGCATInstance) {
    let graph_file = PathBuf::from(prefix.to_owned() + ".dbg.fasta");
    instance.build_graph(
	inputs,
	graph_file,
	Some(input_seq_names),
	31,
	4,
	false,
	None,
	true,
	1,
	ggcat_api::ExtraElaboration::UnitigLinks,
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

pub fn build_pangenome_representations(seq_files: &Vec<String>, seq_clusters: &Vec<Vec<usize>>,
				       out_prefix: &String, instance: &ggcat_api::GGCATInstance) {
    let mut i = 0;
    for group in seq_clusters {
	println!("Building graph {}...", i.to_string() + ".dbg.fasta");
	let prefix = out_prefix.to_owned() + &i.to_string();

	let mut ggcat_input_names: Vec<String> = Vec::new();
	for seq in group {
	    ggcat_input_names.push(seq_files[seq.clone()].clone());
	}

	let ggcat_inputs = open_ggcat_inputs(&ggcat_input_names);
	build_pangenome_graph(ggcat_inputs, &ggcat_input_names, &prefix, instance);
	i = i + 1;
    }
}
pub fn dereplicate_iter(seq_files: &Vec<String>, instance: &ggcat_api::GGCATInstance) -> (Vec<usize>, usize) {
    println!("Calculating ANIs...");
    let ani_result = ani_from_fastx_files(seq_files);

    println!("Building dendrogram...");
    let seqs_by_group = single_linkage_cluster(&ani_result, seq_files.len());

    println!("Building pangenome graphs...");
    build_pangenome_representations(&seq_files, &seqs_by_group, &"./".to_string(), instance);

    let n_clusters = seqs_by_group.len();
    let clusters: Vec<usize> = seqs_by_group.into_iter().flatten().collect();

    return (clusters, n_clusters);
}
