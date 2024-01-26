use std::collections::HashMap;
use std::path::PathBuf;

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
