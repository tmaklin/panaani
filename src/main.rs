// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::HashSet;

use clap::Parser;
use itertools::Itertools;
use log::{info, Record, Level, Metadata};
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;

mod build;
mod cli;
mod clust;
mod dist;

struct Logger;

impl log::Log for Logger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            println!("{} - {}", record.level(), record.args());
        }
    }

    fn flush(&self) {}
}

fn init_log(log_max_level: usize) {
    stderrlog::new()
	.module(module_path!())
	.quiet(false)
	.verbosity(log_max_level)
	.timestamp(stderrlog::Timestamp::Off)
	.init()
	.unwrap();
}

fn init(threads: usize, log_max_level: usize) {
    init_log(log_max_level);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .thread_name(|i| format!("rayon-thread-{}", i))
        .build_global()
        .unwrap();
}

fn read_input_list(input_list_file: &String) -> Vec<String> {
    let f = std::fs::File::open(input_list_file).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(f);

    let mut seq_files: Vec<String> = Vec::new();
    reader.records().into_iter().for_each(|line| {
        let record = line.unwrap();
	seq_files.push(record[0].to_string().clone());
    });
    seq_files
}

fn read_seq_assignments(seq_files_in: &[String], seq_assignments_file: &String) -> Vec<(String, String)> {
    let f = std::fs::File::open(seq_assignments_file).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(f);

    // Read the cluster assignments into a HashMap to get the order correct
    let mut seq_assignments: HashMap<String, String> = HashMap::new();
    reader.records().into_iter().for_each(|line| {
        let record = line.unwrap();
	seq_assignments.insert(record[0].to_string(), record[1].to_string());
    });
    return seq_files_in
	.iter()
	.map(|x| (x.clone(), seq_assignments.get(x).unwrap_or_else(|| { panic!("Input sequence {} was not found in {}!", x, seq_assignments_file) }).clone()))
	.sorted_by(|k1, k2| match k1.1.cmp(&k2.1) {
	    Ordering::Equal => k1.0.cmp(&k2.0),
            other => other,
	})
	.collect::<Vec<(String, String)>>();
}

fn main() {
    let cli = cli::Cli::parse();

    // Subcommands:
    match &cli.command {
        // Run the full pipeline
        Some(cli::Commands::Dereplicate {
            seq_files,
            input_list,
            batch_step,
            linkage_method,
            skani_kmer_size,
            kmer_subsampling_rate,
            marker_compression_factor,
            rescue_small,
            clip_tails,
            median,
            adjust_ani,
            min_aligned_frac,
            ggcat_kmer_size,
            kmer_min_multiplicity,
            minimizer_length,
            no_reverse_complement,
            unitig_type,
            intermediate_compression_level,
            threads,
            memory,
            temp_dir_path,
            ani_threshold,
	    verbose,
	    max_iters,
	    batch_step_strategy,
	    out_prefix,
	    guided_batching,
	    external_clustering_file,
	    initial_batches_file,
        }) => {
	    init_log(if *verbose { 2 } else { 1 });

            let skani_params = panaani::dist::SkaniParams {
                kmer_size: *skani_kmer_size,
                kmer_subsampling_rate: *kmer_subsampling_rate,
                marker_compression_factor: *marker_compression_factor,
                rescue_small: *rescue_small,

                clip_tails: *clip_tails,
                median: *median,
                adjust_ani: *adjust_ani,

                min_aligned_frac: *min_aligned_frac,
		progress: *verbose,
                ..Default::default()
            };

            let kodama_params = panaani::clust::KodamaParams {
                cutoff: *ani_threshold,
                method: if linkage_method.is_some() {
                    match linkage_method.as_ref().unwrap().as_str() {
                        "single" => kodama::Method::Single,
                        "complete" => kodama::Method::Complete,
                        "average" => kodama::Method::Average,
                        "weighted" => kodama::Method::Weighted,
                        "ward" => kodama::Method::Ward,
                        "centroid" => kodama::Method::Centroid,
                        "median" => kodama::Method::Median,
                        &_ => kodama::Method::Single,
                    }
                } else {
                    kodama::Method::Single
                },
                ..Default::default()
            };

            let ggcat_params = panaani::build::GGCATParams {
                kmer_size: *ggcat_kmer_size,
                kmer_min_multiplicity: *kmer_min_multiplicity,
                minimizer_length: if minimizer_length.is_some() {
                    *minimizer_length
                } else {
                    None
                },
                no_reverse_complement: *no_reverse_complement,
		progress: *verbose,
                unitig_type: if unitig_type.is_some() {
                    match unitig_type.as_ref().unwrap().as_str() {
                        "greedymatchtigs" => ggcat_api::ExtraElaboration::GreedyMatchtigs,
                        "unitiglinks" => ggcat_api::ExtraElaboration::UnitigLinks,
                        "eulertigs" => ggcat_api::ExtraElaboration::Eulertigs,
                        "pathtigs" => ggcat_api::ExtraElaboration::Pathtigs,
                        &_ => ggcat_api::ExtraElaboration::GreedyMatchtigs,
                    }
                } else {
                    ggcat_api::ExtraElaboration::GreedyMatchtigs
                },
                intermediate_compression_level: if intermediate_compression_level.is_some() {
                    *intermediate_compression_level
                } else {
                    None
                },
                temp_dir_path: temp_dir_path.clone().unwrap_or("./".to_string()),
                threads: *threads,
                memory: *memory,
		out_prefix: out_prefix.clone().unwrap_or("".to_string()),
                ..Default::default()
            };

	    // TODO seq_files should be mutable by default to avoid cloning
	    let mut seq_files_in: Vec<String> = seq_files.clone();
	    if input_list.is_some() {
		seq_files_in.append(read_input_list(input_list.as_ref().unwrap()).as_mut());
	    }

            let params: panaani::PanaaniParams = panaani::PanaaniParams {
                batch_step: *batch_step,
                batch_step_strategy: batch_step_strategy.clone(),
                max_iters: *max_iters,
		temp_dir: temp_dir_path.clone().unwrap_or("/tmp".to_string()),
		guided: *guided_batching,
		external_clustering: if external_clustering_file.is_some() {
		    Some(read_seq_assignments(&seq_files_in, &external_clustering_file.as_ref().unwrap()).iter().map(|x| x.1.clone()).collect())
		} else {
		    None
		},
		initial_batches: if initial_batches_file.is_some() {
		    Some(read_seq_assignments(&seq_files_in, &initial_batches_file.as_ref().unwrap()).iter().map(|x| x.0.clone()).collect())
		} else {
		    None
		},
		..Default::default()
            };

	    panaani::build::init_ggcat(&Some(ggcat_params.clone()));

            let clusters = panaani::dereplicate(
                &seq_files_in,
                &Some(params),
                &Some(skani_params),
                &Some(kodama_params),
                &Some(ggcat_params),
            );
            let n_clusters = clusters.iter().map(|x| x.1.clone()).unique().collect::<Vec<String>>().len();

            info!("Created {} clusters", n_clusters);
            clusters
                .iter()
                .for_each(|x| println!("{}\t{}", x.0, x.1));
        }

        // Calculate distances between some input fasta files
        Some(cli::Commands::Dist {
            seq_files,
	    input_list,
            threads,
            skani_kmer_size,
            kmer_subsampling_rate,
            marker_compression_factor,
            rescue_small,
            clip_tails,
            median,
            adjust_ani,
            min_aligned_frac,
	    verbose
        }) => {
	    init(*threads as usize, if *verbose { 2 } else { 1 });

            let skani_params = dist::SkaniParams {
                kmer_size: *skani_kmer_size,
                kmer_subsampling_rate: *kmer_subsampling_rate,
                marker_compression_factor: *marker_compression_factor,
                rescue_small: *rescue_small,

                clip_tails: *clip_tails,
                median: *median,
                adjust_ani: *adjust_ani,

                min_aligned_frac: *min_aligned_frac,
		progress: *verbose,
                ..Default::default()
            };

	    // TODO seq_files should be mutable by default to avoid cloning
	    let mut seq_files_in: Vec<String> = seq_files.clone();
	    if input_list.is_some() {
		seq_files_in.append(read_input_list(input_list.as_ref().unwrap()).as_mut());
	    }

            let results = dist::ani_from_fastx_files(&seq_files_in, &Some(skani_params));
	    results.iter().for_each(|x| { println!("{}\t{}\t{}", x.0, x.1, x.2) });
        }

        // Build pangenome representations from input fasta files and their clusters
        Some(cli::Commands::Build {
            seq_files,
	    input_list,
            external_clustering_file,
	    target_cluster,
            threads,
            memory,
            temp_dir_path,
            ggcat_kmer_size,
            kmer_min_multiplicity,
            minimizer_length,
            no_reverse_complement,
            unitig_type,
            intermediate_compression_level,
	    verbose,
	    out_prefix,
        }) => {
	    init_log(if *verbose { 2 } else { 1 });

            let ggcat_params = panaani::build::GGCATParams {
                kmer_size: *ggcat_kmer_size,
                kmer_min_multiplicity: *kmer_min_multiplicity,
                minimizer_length: if minimizer_length.is_some() {
                    *minimizer_length
                } else {
                    None
                },
                no_reverse_complement: *no_reverse_complement,
		progress: *verbose,
                unitig_type: if unitig_type.is_some() {
                    match unitig_type.as_ref().unwrap().as_str() {
                        "greedymatchtigs" => ggcat_api::ExtraElaboration::GreedyMatchtigs,
                        "unitiglinks" => ggcat_api::ExtraElaboration::UnitigLinks,
                        "eulertigs" => ggcat_api::ExtraElaboration::Eulertigs,
                        "pathtigs" => ggcat_api::ExtraElaboration::Pathtigs,
                        &_ => ggcat_api::ExtraElaboration::GreedyMatchtigs,
                    }
                } else {
                    ggcat_api::ExtraElaboration::GreedyMatchtigs
                },
                intermediate_compression_level: if intermediate_compression_level.is_some() {
                    *intermediate_compression_level
                } else {
                    None
                },
                temp_dir_path: temp_dir_path.clone().unwrap_or("./".to_string()),
                threads: *threads,
                memory: *memory,
		out_prefix: out_prefix.clone().unwrap_or("".to_string()),
                ..Default::default()
            };

	    panaani::build::init_ggcat(&Some(ggcat_params.clone()));

	    // TODO seq_files should be mutable by default to avoid cloning
	    let mut seq_files_in: Vec<String> = seq_files.clone();
	    if input_list.is_some() {
		seq_files_in.append(read_input_list(input_list.as_ref().unwrap()).as_mut());
	    }

	    let external_clusters: Vec<(String, String)> = read_seq_assignments(&seq_files_in, &external_clustering_file.as_ref().unwrap());
	    let mut seq_to_cluster = panaani::assign_seqs(&external_clusters.iter().map(|x| x.0.clone()).collect::<Vec<String>>(),
							  &external_clusters.iter().map(|x| x.1.clone()).collect::<Vec<String>>());

	    if target_cluster.is_some() {
		let mut target_to_seqs: HashMap<String, Vec<String>> = HashMap::new();
		target_to_seqs.insert(target_cluster.as_ref().unwrap().clone(), seq_to_cluster.get(target_cluster.as_ref().unwrap()).unwrap().clone());
		seq_to_cluster = target_to_seqs;
	    }

            panaani::build::build_pangenome_representations(
		&seq_to_cluster,
                &Some(ggcat_params),
            );
        }

        // Cluster distance data created with `skani dist` or `panaani dist`.
        Some(cli::Commands::Cluster {
            dist_file,
            ani_threshold,
            linkage_method,
	    verbose,
	    out_prefix,
        }) => {
	    init(1, if *verbose { 2 } else { 1 });

            let kodama_params = clust::KodamaParams {
                cutoff: *ani_threshold,
                method: if linkage_method.is_some() {
                    match linkage_method.as_ref().unwrap().as_str() {
                        "single" => kodama::Method::Single,
                        "complete" => kodama::Method::Complete,
                        "average" => kodama::Method::Average,
                        "weighted" => kodama::Method::Weighted,
                        "ward" => kodama::Method::Ward,
                        "centroid" => kodama::Method::Centroid,
                        "median" => kodama::Method::Median,
                        &_ => kodama::Method::Single,
                    }
                } else {
                    kodama::Method::Single
                },
                ..Default::default()
            };

            let f = std::fs::File::open(dist_file).unwrap();
            let mut reader = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_reader(f);

            let mut seq_names: HashSet<String> = HashSet::new();
            let mut res: Vec<(String, String, f32)> = Vec::new();
            for line in reader.records().into_iter() {
                let record = line.unwrap();
                res.push((
                    record[0].to_string().clone(),
                    record[1].to_string().clone(),
                    record[2].parse::<f32>().unwrap(),
                ));
                seq_names.insert(record[0].to_string());
                seq_names.insert(record[1].to_string());
            }
	    res.sort_by(|k1, k2| match k1.0.cmp(&k2.0) {
		Ordering::Equal => k1.1.cmp(&k2.1),
		other => other,
            });

	    let old_clusters = seq_names.iter().map(|x| x).cloned().collect::<Vec<String>>();
            let hclust_res = clust::single_linkage_cluster(&res, &Some(kodama_params));

	    let prefix = out_prefix.clone().unwrap_or("".to_string()) + &"panANI-".to_string();
	    let new_clusters: &mut Vec<String> = &mut
		panaani::match_clustering_results(&old_clusters, &old_clusters, &hclust_res, &prefix);

	    let mut files_in_cluster: HashMap<String, Vec<String>> = HashMap::new();
	    seq_names.iter().zip(new_clusters.iter()).for_each(|x| {
		if files_in_cluster.contains_key(x.1) {
		    files_in_cluster.get_mut(x.1).unwrap().push(x.0.clone());
		} else {
		    files_in_cluster.insert(x.1.clone(), vec![x.0.clone()]);
		}
	    });
	    seq_names.iter().zip(new_clusters.iter_mut()).for_each(|x| {
		if files_in_cluster.get(x.1).unwrap().len() == 1 {
		    *x.1 = x.0.clone();
		}
	    });

	    old_clusters.iter().zip(new_clusters.iter()).for_each(|x| { println!("{}\t{}", x.0, x.1) } );
        }

        // Calculate distances between some input fasta files
        Some(cli::Commands::Assign {
            query_files,
	    query_files_list,
	    ref_files_list,
            threads,
	    verbose,
            skani_kmer_size,
            kmer_subsampling_rate,
            marker_compression_factor,
            rescue_small,
            clip_tails,
            median,
            adjust_ani,
            min_aligned_frac,
	    ani_threshold,
        }) => {
	    init(*threads as usize, if *verbose { 2 } else { 1 });

            let skani_params = dist::SkaniParams {
                kmer_size: *skani_kmer_size,
                kmer_subsampling_rate: *kmer_subsampling_rate,
                marker_compression_factor: *marker_compression_factor,
                rescue_small: *rescue_small,

                clip_tails: *clip_tails,
                median: *median,
                adjust_ani: *adjust_ani,

                min_aligned_frac: *min_aligned_frac,
		progress: *verbose,
                ..Default::default()
            };

	    let cmd_params = skani::params::CommandParams {
		screen: false,
		screen_val: 0.00,
		mode: skani::params::Mode::Dist,
		out_file_name: "".to_string(),
		ref_files: vec![],
		query_files: vec![],
		refs_are_sketch: false,
		queries_are_sketch: false,
		robust: skani_params.clip_tails,
		median: skani_params.median,
		sparse: false,
		full_matrix: false,
		max_results: 10000000,
		individual_contig_q: false,
		individual_contig_r: false,
		min_aligned_frac: 0.0,
		keep_refs: false,
		est_ci: skani_params.bootstrap_ci,
		learned_ani: skani_params.adjust_ani,
		detailed_out: false,
		rescue_small: skani_params.rescue_small,
		distance: true,
	    };
	    let adjust_ani = skani::regression::get_model(skani_params.kmer_subsampling_rate.into(), false);

	    let mut query_files_in: Vec<String> = query_files.clone();
	    if query_files_list.is_some() {
		query_files_in.append(read_input_list(query_files_list.as_ref().unwrap()).as_mut());
	    }

	    let mut ref_files_in: Vec<String> = Vec::new();
	    ref_files_in.append(read_input_list(ref_files_list.as_ref().unwrap()).as_mut());

	    let ref_db = dist::sketch_fastx_files(&ref_files_in, Some(skani::params::SketchParams::new(
		skani_params.marker_compression_factor as usize,
		skani_params.kmer_subsampling_rate as usize,
		skani_params.kmer_size as usize,
		false,
		false,
	    )));

	    let query_db = dist::sketch_fastx_files(&query_files_in, Some(skani::params::SketchParams::new(
		skani_params.marker_compression_factor as usize,
		skani_params.kmer_subsampling_rate as usize,
		skani_params.kmer_size as usize,
		false,
		false,
	    )));

	    let query_dists = ref_db
		.iter()
		.map(|r| { query_db
			   .par_iter()
			   .map(|q| {
			       (q.file_name.clone(),
				r.file_name.clone(),
				skani::chain::chain_seeds(
				    r,
				    q,
				    skani::chain::map_params_from_sketch(
					r,
					false,
					&cmd_params,
					&adjust_ani,
				    ),
				)
			       )
			   })
			   .collect::<Vec<(String, String, skani::types::AniEstResult)>>()
		})
		.flatten()
		.map(|x| {
		    (x.0,
		     x.1,
		     dist::filter_ani(x.2.ani, x.2.align_fraction_ref, x.2.align_fraction_query, skani_params.min_aligned_frac as f32, skani_params.min_aligned_frac as f32)
		    )
		})
		.collect::<Vec<(String, String, f32)>>();

	    // Check that all queries were assigned
	    let mut all_assigned = true;
	    let mut best_match: HashMap<String, (String, f32, bool)> = HashMap::new();
	    query_dists
		.iter()
		.for_each(|x| {
		    if !best_match.contains_key(&x.0) {
			best_match.insert(x.0.clone(), (x.1.clone(), x.2.clone(), false));
		    } else if x.2 > best_match.get(&x.0).unwrap().1 {
			let assigned_twice: bool = (best_match.get(&x.0).unwrap().1 > *ani_threshold && x.2 > *ani_threshold) || best_match.get(&x.0).unwrap().2;
			*best_match.get_mut(&x.0).unwrap() = (x.1.clone(), x.2.clone(), assigned_twice);
		    }
		});

	    let mut all_unambiguous = true;
	    best_match
		.iter()
		.for_each(|x| { all_assigned &= x.1.1 > *ani_threshold; all_unambiguous &= !x.1.2 });

	    if all_assigned && all_unambiguous {
		info!("Assigned {}/{} queries unambiguously to reference database (ANI threshold {})", query_db.len(), query_db.len(), ani_threshold);
		best_match
		    .iter()
		    .for_each(|x| { println!("{}\t{}", x.0, x.1.0); });
	    } else if all_unambiguous {
		let n_assigned: usize = best_match.iter().filter(|x| x.1.1 > *ani_threshold).count();
		info!("Assigned {}/{} queries unambiguously to reference database (ANI threshold {})", n_assigned, query_db.len(), ani_threshold);
		info!("{}/{} queries could not be assigned to any reference", query_db.len() - n_assigned,  query_db.len());
		best_match
		    .iter()
		    .for_each(|x| { if x.1.1 > *ani_threshold { println!("{}\t{}", x.0, x.1.0); } else { println!("{}\t{}", x.0, "new_cluster"); } });
	    } else {
		let n_assigned: usize = best_match.iter().filter(|x| x.1.1 > *ani_threshold).count();
		let n_ambiguous: usize = best_match.iter().filter(|x| x.1.2).count();
		info!("Assigned {}/{} queries unambiguously to reference database (ANI threshold {})", n_assigned - n_ambiguous, query_db.len(), ani_threshold);
		info!("{}/{} queries could not be assigned to any reference", query_db.len() - n_assigned,  query_db.len());
		info!("{}/{} queries were assigned to multiple references", n_ambiguous, query_db.len());
		best_match
		    .iter()
		    .for_each(|x| { if x.1.1 > *ani_threshold && !x.1.2 { println!("{}\t{}", x.0, x.1.0); } else if x.1.1 > *ani_threshold && x.1.2 { println!("{}\t{}", x.0, "ambiguous"); } else { println!("{}\t{}", x.0, "new_cluster"); } });
	    }
	}
        None => {}
    }
}
