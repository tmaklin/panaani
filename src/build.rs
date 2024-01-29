// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use std::collections::HashMap;
use std::path::PathBuf;

use parallel_processor::enable_counters_logging;
use parallel_processor::memory_data_size::MemoryDataSize;
use parallel_processor::memory_fs::MemoryFs;
use parallel_processor::phase_times_monitor::PHASES_TIMES_MONITOR;
use std::cmp::max;
use std::fs::create_dir_all;
use std::sync::atomic::Ordering;
use std::time::Duration;
use log::info;
use log::warn;

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

            intermediate_compression_level: None,
            stats_file: None,
        }
    }
}

fn open_ggcat_inputs(seq_files: &[String]) -> Vec<ggcat_api::GeneralSequenceBlockData> {
    let ggcat_inputs: Vec<ggcat_api::GeneralSequenceBlockData> = seq_files
        .iter()
        .map(|x| ggcat_api::GeneralSequenceBlockData::FASTA((PathBuf::from(x), None)))
        .collect();
    return ggcat_inputs;
}

fn build_pangenome_graph(input_seq_names: &[String], prefix: &String, params: &GGCATParams) {
    info!("Building graph {} from {} sequences...", prefix, input_seq_names.len());

    let config = GGCATConfig {
        temp_dir: Some(PathBuf::from(params.temp_dir_path.clone())),
        memory: params.memory as f64,
        prefer_memory: false,
        total_threads_count: params.threads as usize,
        intermediate_compression_level: params.intermediate_compression_level,
        stats_file: params.stats_file.clone(),
    };

    // Increase the maximum allowed number of open files
    if let Err(err) = fdlimit::raise_fd_limit() {
        warn!("Failed to increase the maximum number of open files: {}", err);
    }

    ggcat_config::PREFER_MEMORY.store(false, Ordering::Relaxed);

    if let Some(temp_dir) = &config.temp_dir {
        create_dir_all(temp_dir).unwrap();
    }

    if let Some(stats_file) = &config.stats_file {
        enable_counters_logging(stats_file, Duration::from_millis(1000), |val| {
            val["phase"] = PHASES_TIMES_MONITOR.read().get_phase_desc().into();
        });
    }

    MemoryFs::init(
        MemoryDataSize::from_bytes(
            (config.memory * (MemoryDataSize::OCTET_GIBIOCTET_FACTOR as f64)) as usize,
        ),
        16 * config.total_threads_count,
        max(1, config.total_threads_count / 4),
        8192,
    );

    let instance = Some(Box::leak(Box::new(GGCATInstance(config)))).unwrap();

    let graph_file = PathBuf::from(prefix.to_owned());

    let inputs = open_ggcat_inputs(input_seq_names);
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

pub fn build_pangenome_representations(seq_files: &[(String, String)], params: &GGCATParams) {
    let mut files_in_cluster: HashMap<String, Vec<String>> = HashMap::new();

    seq_files.iter().for_each(|x| {
        if files_in_cluster.contains_key(&x.1) {
            files_in_cluster.get_mut(&x.1).unwrap().push(x.0.clone());
        } else {
            files_in_cluster.insert(x.1.clone(), vec![x.0.clone()]);
        }
    });

    files_in_cluster
        .iter()
        .for_each(|x| build_pangenome_graph(x.1, x.0, params));
}
