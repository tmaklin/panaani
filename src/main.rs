// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use std::collections::HashSet;

use clap::{Parser, Subcommand};
use itertools::Itertools;

pub mod build;
pub mod clust;
pub mod dist;

#[derive(Parser)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    // Print testing stuff
    Dereplicate {
	#[arg(group = "input", required = true)]
	seq_files: Vec<String>,

	// Resources
	#[arg(short = 't', long = "threads", default_value_t = 1)]
	threads: u32,

	#[arg(short = 'm', long = "memory", default_value_t = 4)]
	memory: u32,

	#[arg(long = "tmp-dir", required = false)]
	temp_dir_path: Option<String>,

	// Dereplicate parameters
	#[arg(short = 'b', long = "batch-step", default_value_t = 50)]
	batch_step: usize,
    },

    Dist {
	#[arg(group = "input")]
	seq_files: Vec<String>,
    },
    Build {
	#[arg(group = "input", required = true)]
	seq_files: Vec<String>,

    },
    Cluster {
	#[arg(group = "input")]
	dist_file: String,
    },
}

fn main() {
    println!("panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters");
    let cli = Cli::parse();

    // Subcommands:
    match &cli.command {

	// Run the full pipeline
	Some(Commands::Dereplicate { seq_files,
				     threads,
				     memory,
				     temp_dir_path,
				     batch_step,
	}) => {
	    rayon::ThreadPoolBuilder::new()
		.num_threads(*threads as usize)
		.thread_name(|i| format!("rayon-thread-{}", i))
		.build_global()
		.unwrap();

	    let params: panaani::PanaaniParams = panaani::PanaaniParams { batch_step: *batch_step };
	    let ggcat_params: panaani::build::GGCATParams = panaani::build::GGCATParams {
		temp_dir_path: temp_dir_path.clone().unwrap_or("./".to_string()),
		threads: *threads,
		memory: *memory,
		..Default::default() };

	    let clusters = panaani::dereplicate(&seq_files, &seq_files, Some(params), None, Some(ggcat_params));
	    let n_clusters = clusters.iter().map(|x| x.1.clone()).unique().collect::<Vec<String>>().len();

	    println!("Created {} clusters", n_clusters);
	    for cluster in clusters {
		println!("{}\t{}", cluster.0, cluster.1);
	    }
	}

	// Calculate distances between some input fasta files
	Some(Commands::Dist { seq_files }) => {
	    let results = dist::ani_from_fastx_files(seq_files, &dist::SkaniParams::default());
	    for res in results {
		println!("{}\t{}\t{}\t{}\t{}",
			 res.ref_file,
			 res.query_file,
			 res.ani,
			 res.align_fraction_ref,
			 res.align_fraction_query
			 );
	    }
	}

	// Build a de Bruijn graph from some input fasta files
	Some(Commands::Build { seq_files  }) => {
	    let ggcat_inputs = build::open_ggcat_inputs(seq_files);
	    build::build_pangenome_graph(ggcat_inputs, seq_files, &("out".to_string() + ".dbg.fasta"), &build::GGCATParams::default());
	}

	// Cluster distance data created with `skani dist` or `panaani dist`.
	Some(Commands::Cluster { dist_file,  }) => {
            let f = std::fs::File::open(dist_file).unwrap();
	    let mut reader = csv::ReaderBuilder::new()
		.delimiter(b'\t')
		.has_headers(false)
		.from_reader(f);

	    let mut seq_names: HashSet<String> = HashSet::new();
	    let mut res: Vec<(String, String, f32, f32, f32)> = Vec::new();
	    for line in reader.records().into_iter() {
		let record = line.unwrap();
		res.push((
		    record[0].to_string().clone(),
		    record[1].to_string().clone(),
		    record[2].parse().unwrap(),
		    record[3].parse().unwrap(),
		    record[4].parse().unwrap()
		));
		seq_names.insert(record[0].to_string());
		seq_names.insert(record[0].to_string());
	    }
	    res.sort_by_key(|k| (k.0.clone(), k.1.clone()));

	    clust::single_linkage_cluster2(&res, seq_names.len());
	}
	None => {}
    }
}
