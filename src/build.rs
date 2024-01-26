use std::collections::HashMap;
use std::path::PathBuf;

use ggcat_api::{GGCATConfig,GGCATInstance};

#[derive(Clone)]
pub struct GGCATParams {
    // k-mer sketching
    kmer_size: u32,
    kmer_min_multiplicity: u64,

    // Graph construction
    minimizer_length: Option<usize>,
    no_reverse_complement: bool,
    unitig_type: ggcat_api::ExtraElaboration,

    // Resources
    threads: u32,
    memory: u32,
    prefer_memory: bool,
    temp_dir_path: String,

    // Intermediate outputs
    intermediate_compression_level: Option<u32>,
    stats_file: Option<PathBuf>,
}

impl Default for GGCATParams {
    fn default() -> GGCATParams {
	GGCATParams {
	    kmer_size: 51,
	    kmer_min_multiplicity: 1,

	    minimizer_length: None,
	    no_reverse_complement: false,
	    unitig_type: ggcat_api::ExtraElaboration::GreedyMatchtigs,

	    threads: 1,
	    memory: 4,
	    prefer_memory: true,
	    temp_dir_path: "/tmp".to_string(),

	    intermediate_compression_level: None,
	    stats_file: None,
	}
    }
}

pub fn build_pangenome_graph(inputs: Vec<ggcat_api::GeneralSequenceBlockData>,
			     input_seq_names: &Vec<String>,
			     prefix: &String, params: &GGCATParams) {

    let instance = GGCATInstance::create(GGCATConfig {
	temp_dir: Some(PathBuf::from(params.temp_dir_path.clone())),
	memory: params.memory as f64,
	prefer_memory: params.prefer_memory,
	total_threads_count: params.threads as usize,
	intermediate_compression_level: params.intermediate_compression_level,
	stats_file: params.stats_file.clone(),
    });

    let graph_file = PathBuf::from(prefix.to_owned());
    instance.build_graph(
	inputs,
	graph_file,
	Some(input_seq_names),
	params.kmer_size as usize,
	params.threads as usize,
	params.no_reverse_complement,
	params.minimizer_length,
	false, // No colors
	params.kmer_min_multiplicity as usize,
	params.unitig_type,
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

pub fn build_pangenome_representations(seq_files: &Vec<(String, String)>, params: &GGCATParams) {
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
	build_pangenome_graph(ggcat_inputs, &ggcat_input_names, &graph_name, params);
    }

}
