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
	    let batch_step = 50;

	    let clusters = panaani::dereplicate(&seq_files, &seq_files, &batch_step, None, None);
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
