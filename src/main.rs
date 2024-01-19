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

use ggcat_api::{
    GGCATConfig,
    GGCATInstance,
};

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
    let c = 125;
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
	    // TODO implement
	}

	// Calculate distances between some input fasta files
	Some(Commands::Dist { seq_files }) => {
	    let sketches = file_io::fastx_to_sketches(seq_files, &sketch_params, true);
	    let adjust_ani = regression::get_model(sketch_params.c, false);
	    for pair in sketches.iter().combinations(2) {
		let map_params = chain::map_params_from_sketch(pair.first().unwrap(), false, &cmd_params, &adjust_ani);
		let ani_result = chain::chain_seeds(pair.first().unwrap(), pair.last().unwrap(), map_params);
		println!("{}", ani_result.ani);
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
