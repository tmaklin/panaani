// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use std::collections::HashMap;
use std::io::Read;
use std::path::PathBuf;

use log::debug;
use log::trace;

use ggcat_api::{GGCATConfig, GGCATInstance};

#[derive(Clone)]
pub struct GGCATParams {
    // k-mer sketching
    pub kmer_size: u32,
    pub kmer_min_multiplicity: u64,

    // Graph construction
    pub minimizer_length: Option<usize>,
    pub no_reverse_complement: bool,
    pub unitig_type: ggcat_api::ExtraElaboration,

    // Resources
    pub threads: u32,
    pub memory: u32,
    pub temp_dir_path: String,

    // Output
    pub out_prefix: String,

    // Intermediate outputs
    pub intermediate_compression_level: Option<u32>,
    pub stats_file: Option<PathBuf>,
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
            temp_dir_path: "/tmp".to_string(),

	    out_prefix: "".to_string(),

            intermediate_compression_level: None,
            stats_file: None,
        }
    }
}

fn build_pangenome_graph(input_seq_names: &[String], prefix: &String, instance: &GGCATInstance, params: &GGCATParams) {
    debug!("Building graph {} from {} sequences:", prefix, input_seq_names.len());
    input_seq_names.iter().for_each(|x| { debug!("\t{}", x) });

    let graph_file = PathBuf::from(params.out_prefix.clone() + prefix);
    let inputs: Vec<ggcat_api::GeneralSequenceBlockData> = input_seq_names
        .iter()
        .map(|x| ggcat_api::GeneralSequenceBlockData::FASTA((PathBuf::from(x), None)))
        .collect();

    let mut buf = gag::BufferRedirect::stdout().unwrap();
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
    let mut output = String::new();
    buf.read_to_string(&mut output).unwrap();
    drop(buf);
    for line in output.lines() {
	trace!("{}", line);
    }
}

pub fn build_pangenome_representations(
    seq_files: &[String],
    clusters: &mut [String],
    opt: &Option<GGCATParams>,
) {
    let params = opt.clone().unwrap_or(GGCATParams::default());
    let mut files_in_cluster: HashMap<String, Vec<String>> = HashMap::new();

    seq_files.iter().zip(clusters.iter()).for_each(|x| {
        if files_in_cluster.contains_key(x.1) {
            files_in_cluster.get_mut(x.1).unwrap().push(x.0.clone());
        } else {
            files_in_cluster.insert(x.1.clone(), vec![x.0.clone()]);
        }
    });

    let config = GGCATConfig {
        temp_dir: Some(PathBuf::from(params.temp_dir_path.clone())),
        memory: params.memory as f64,
        prefer_memory: false,
        total_threads_count: params.threads as usize,
        intermediate_compression_level: params.intermediate_compression_level,
        stats_file: params.stats_file.clone(),
    };

    let mut buf = gag::BufferRedirect::stdout().unwrap();
    let mut output = String::new();
    buf.read_to_string(&mut output).unwrap();
    let instance = ggcat_api::GGCATInstance::create(config);
    drop(buf);
    for line in output.lines() {
	debug!("{}", line);
    }

    files_in_cluster
        .iter()
	.filter(|x| x.1.len() > 1)
        .for_each(|x| build_pangenome_graph(x.1, x.0, &instance, &params));

    seq_files.iter().zip(clusters.iter_mut()).for_each(|x| {
        if files_in_cluster.get(x.1).unwrap().len() == 1 {
            *x.1 = x.0.clone();
        }
    });
}
