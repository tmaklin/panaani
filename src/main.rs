// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use std::path::PathBuf;
use std::collections::HashSet;

use clap::{Parser, Subcommand};

use ggcat_api::{GGCATConfig,GGCATInstance};

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
	    let instance = GGCATInstance::create(GGCATConfig {
		temp_dir: Some(PathBuf::from("/tmp")),
		memory: 8.0,
		prefer_memory: true,
		total_threads_count: 4,
		intermediate_compression_level: None,
		stats_file: None,
	    });

	    let batch_step = 50;
	    let mut iter = 0;
	    let mut iter_inputs = seq_files.clone();
	    let mut batch_size = batch_step;

	    while batch_size < iter_inputs.len() {
		println!("iter: {}\tProcessing inputs in batches of {} sequences", iter, batch_size);
		let mut rng = rand::thread_rng();

		// horrible hack to use random file names within each batch
		let res: Vec<(Vec<String>, Vec<String>, Vec<usize>, usize)> = iter_inputs
		    .chunks(batch_size)
		    .map(|x| panaani::dereplicate_iter(&Vec::from(x), &(iter.to_string() + "_" + &(rng.gen::<u64>() as u64).to_string() + "-" ), instance))
		    .collect();

		// payment for horrible hack
		let mut hash: HashSet<String> = HashSet::new();
		for i in 0..res.len() {
		    for j in 0..res[i].1.len() {
			hash.insert(res[i].1[j].clone());
		    }
		}
		iter_inputs = hash.into_iter().collect();

		batch_size += batch_step;
		iter += 1;
	    }
	    println!("iter: {}\tAll inputs fit in the same batch", iter);
	    let (seqs, pangenome_files, clusters, n_clusters) = panaani::dereplicate_iter(&iter_inputs, &"panANI-".to_string(), instance);
	    println!("Created {} clusters", n_clusters);
	}

	// Calculate distances between some input fasta files
	Some(Commands::Dist { seq_files }) => {
	    let results = panaani::ani_from_fastx_files(seq_files);
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
	    let ggcat_inputs = panaani::open_ggcat_inputs(seq_files);
	    let instance = GGCATInstance::create(GGCATConfig {
		temp_dir: Some(PathBuf::from("./")),
		memory: 2.0,
		prefer_memory: true,
		total_threads_count: 4,
		intermediate_compression_level: None,
		stats_file: None,
	    });

	    panaani::build_pangenome_graph(ggcat_inputs, seq_files, &"out".to_string(), instance);
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

	    panaani::single_linkage_cluster2(&res, seq_names.len());
	}
	None => {}
    }
}
