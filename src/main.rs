// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use clap::{Parser, Subcommand};
use itertools::Itertools;

use kodama::{Method, linkage};

use skani::chain;
use skani::file_io;
use skani::params::*;
use skani::regression;
use skani::types::AniEstResult;

use rayon;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use std::sync::mpsc::channel;

use ggcat_api::{
    GGCATConfig,
    GGCATInstance,
    ExtraElaboration,
};

use ggcat_api::GeneralSequenceBlockData::FASTA;
use std::path::PathBuf;

#[derive(Parser)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    // Print testing stuff
    Dist {
	#[arg(group = "input")]
	seq_files: Vec<String>,
    },
    Build {
	#[arg(group = "input")]
	seq_files: Vec<String>,
    },
    Cluster {
	#[arg(group = "input")]
	dist_file: String,
    },
    Dereplicate {
	#[arg(group = "input")]
	seq_files: Vec<String>,
    },
}

fn main() {
    println!("panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters");

    // skani parameters
    let m = 1000;
    let c = 30;
    let k = 15;
    let sketch_params = SketchParams::new(m, c, k, false, false);
    let cmd_params = CommandParams {
        screen: false,
        screen_val: 0.00,
        mode: Mode::Dist,
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

    let cli = Cli::parse();

    // Subcommands:
    match &cli.command {

	// Run the full pipeline
	Some(Commands::Dereplicate { seq_files }) => {
	    rayon::ThreadPoolBuilder::new()
		.num_threads(4)
		.thread_name(|i| format!("rayon-thread-{}", i))
		.build_global()
		.unwrap();

	    let sketches = file_io::fastx_to_sketches(seq_files, &sketch_params, true);
	    let adjust_ani = regression::get_model(sketch_params.c, false);
	    println!("Calculating ANIs...");

	    let (sender, receiver) = channel();

	    sketches.iter().combinations(2).par_bridge().for_each_with(sender, |s, pair| {
		let map_params = chain::map_params_from_sketch(pair.first().unwrap(), false, &cmd_params, &adjust_ani);
		let res = chain::chain_seeds(pair.first().unwrap(), pair.last().unwrap(), map_params);
		println!("{}\t{}\t{}\t{}\t{}", pair.first().unwrap().file_name, pair.last().unwrap().file_name,
			 res.ani, res.align_fraction_query, res.align_fraction_ref);
		let mut ret = 1.0;
		if res.ani > 0.0 && res.ani < 1.0 && !res.ani.is_nan() {
		    ret = 1.0 - res.ani;
		}
		s.send(ret);
	    });
	    let mut ani_result: Vec<f32> = receiver.iter().collect();

	    println!("{}", ani_result.len());
	    let num_seqs = seq_files.len();
	    let dend = linkage(&mut ani_result, num_seqs, Method::Single);

	    let cutoff = 1.0 - 0.95;
	    let mut num_groups = 0;
	    let num_nodes = 2 * num_seqs - 1;
	    let mut membership = vec![None; num_nodes];

	    println!("Building dendrogram...");
	    for (cluster_index, step) in dend.steps().iter().enumerate().rev() {
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

	    println!("Clustering...");
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

	    let mut i = 0;
	    println!("Building pangenome graphs...");
	    let instance = GGCATInstance::create(GGCATConfig {
		temp_dir: Some(PathBuf::from("/tmp")),
		memory: 2.0,
		prefer_memory: true,
		total_threads_count: 1,
		intermediate_compression_level: None,
		stats_file: None,
	    });
	    i = 0;
	    for group in seqs_by_group {
		println!("Building graph {}...", i.to_string() + ".dbg");
		let graph_file = PathBuf::from(i.to_string() + ".dbg");
		let mut ggcat_inputs: Vec<ggcat_api::GeneralSequenceBlockData> = Vec::new();
		for seq in group {
		    let file = &seq_files[seq];
		    println!("{}\t{}", file, i);
		    ggcat_inputs.push(ggcat_api::GeneralSequenceBlockData::FASTA((
			PathBuf::from(file),
			None,
		    )));
		}
		instance.build_graph(
		    ggcat_inputs,
		    graph_file,
		    Some(&seq_files),
		    31,
		    4,
		    false,
		    None,
		    true,
		    1,
		    ExtraElaboration::UnitigLinks,
		);
		i = i + 1;
	    }
	}

	// Calculate distances between some input fasta files
	Some(Commands::Dist { seq_files }) => {
	    let sketches = file_io::fastx_to_sketches(seq_files, &sketch_params, true);
	    let adjust_ani = regression::get_model(sketch_params.c, false);
	    let mut ani_result: Vec<f32> = Vec::new();
	    for pair in sketches.iter().combinations(2) {
		let map_params = chain::map_params_from_sketch(pair.first().unwrap(), false, &cmd_params, &adjust_ani);
		ani_result.push(1.0 - chain::chain_seeds(pair.first().unwrap(), pair.last().unwrap(), map_params).ani);
	    }
	}

	// Build a de Bruijn graph from some input fasta files
	Some(Commands::Build { seq_files,  }) => {
	    let graph_file = PathBuf::from("/tmp/dbg.fa");
	    let instance = GGCATInstance::create(GGCATConfig {
		temp_dir: Some(PathBuf::from("/tmp")),
		memory: 2.0,
		prefer_memory: true,
		total_threads_count: 4,
		intermediate_compression_level: None,
		stats_file: None,
	    });
	}

	// Cluster a distance matrix
	Some(Commands::Cluster { dist_file,  }) => {
	    let mut distances = vec![0.99, 0.94, 0.5];
	    let n_obs = 3;
	    let dend = linkage(&mut distances, n_obs, Method::Single);
	}
	None => {}
    }
}
