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
use std::io::Read;

use clap::Parser;
use itertools::Itertools;
use log::{info, trace, Record, Level, LevelFilter, Metadata};

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

fn init_ggcat(opt: &Option<build::GGCATParams>) {
    // GGCAT API force initializes rayon::ThreadPool using build_global
    // so chaining skani -> kodama -> ggcat requires calling the GGCAT
    // API before running skani to get parallelism working correctly.
    let params = opt.clone().unwrap_or(build::GGCATParams::default());
    let config = ggcat_api::GGCATConfig {
        temp_dir: Some(std::path::PathBuf::from(params.temp_dir_path.clone())),
        memory: params.memory as f64,
        prefer_memory: false,
        total_threads_count: params.threads as usize,
        intermediate_compression_level: params.intermediate_compression_level,
        stats_file: params.stats_file.clone(),
    };

    // GGCATInstance is static in the API and can be retrieved by calling
    // GGCATInstance::create again, so no need to return the value.
    let mut buf = gag::BufferRedirect::stdout().unwrap();
    let _ = ggcat_api::GGCATInstance::create(config);
    let mut output = String::new();
    buf.read_to_string(&mut output).unwrap();
    drop(buf);
    for line in output.lines() {
	trace!("{}", line);
    }
}

fn init_ggcat2(opt: &Option<panaani::build::GGCATParams>) {
    // GGCAT API force initializes rayon::ThreadPool using build_global
    // so chaining skani -> kodama -> ggcat requires calling the GGCAT
    // API before running skani to get parallelism working correctly.
    let params = opt.clone().unwrap_or(panaani::build::GGCATParams::default());
    let config = ggcat_api::GGCATConfig {
        temp_dir: Some(std::path::PathBuf::from(params.temp_dir_path.clone())),
        memory: params.memory as f64,
        prefer_memory: false,
        total_threads_count: params.threads as usize,
        intermediate_compression_level: params.intermediate_compression_level,
        stats_file: params.stats_file.clone(),
    };

    // GGCATInstance is static in the API and can be retrieved by calling
    // GGCATInstance::create again, so no need to return the value.
    let mut buf = gag::BufferRedirect::stdout().unwrap();
    let _ = ggcat_api::GGCATInstance::create(config);
    let mut output = String::new();
    buf.read_to_string(&mut output).unwrap();
    drop(buf);
    for line in output.lines() {
	trace!("{}", line);
    }
}

fn init_log(log: &'static Logger, log_max_level: LevelFilter) {
    let _ = log::set_logger(log).map(|()| log::set_max_level(log_max_level ));
}

fn init(threads: usize, log: &'static Logger, log_max_level: LevelFilter) {
    init_log(log, log_max_level);
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

fn main() {
    static LOG: Logger = Logger;
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
        }) => {
	    init_log(&LOG, if *verbose { LevelFilter::Info } else { LevelFilter::Warn });

            let params: panaani::PanaaniParams = panaani::PanaaniParams {
                batch_step: *batch_step,
                batch_step_strategy: batch_step_strategy.clone(),
                max_iters: *max_iters,
		temp_dir: temp_dir_path.clone().unwrap_or("/tmp".to_string()),
            };

            let skani_params = panaani::dist::SkaniParams {
                kmer_size: *skani_kmer_size,
                kmer_subsampling_rate: *kmer_subsampling_rate,
                marker_compression_factor: *marker_compression_factor,
                rescue_small: *rescue_small,

                clip_tails: *clip_tails,
                median: *median,
                adjust_ani: *adjust_ani,

                min_aligned_frac: *min_aligned_frac,
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

	    init_ggcat2(&Some(ggcat_params.clone()));

            let clusters = panaani::dereplicate(
                &seq_files_in,
                &seq_files_in,
                &Some(params),
                &Some(skani_params),
                &Some(kodama_params),
                &Some(ggcat_params),
            );
            let n_clusters = clusters.iter().unique().collect::<Vec<&String>>().len();

            info!("Created {} clusters", n_clusters);
            seq_files
                .iter()
                .zip(clusters.iter())
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
	    init(*threads as usize, &LOG, if *verbose { LevelFilter::Info } else { LevelFilter::Warn });

            let skani_params = dist::SkaniParams {
                kmer_size: *skani_kmer_size,
                kmer_subsampling_rate: *kmer_subsampling_rate,
                marker_compression_factor: *marker_compression_factor,
                rescue_small: *rescue_small,

                clip_tails: *clip_tails,
                median: *median,
                adjust_ani: *adjust_ani,

                min_aligned_frac: *min_aligned_frac,
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
            external_clusters,
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
	    init_log(&LOG, if *verbose { LevelFilter::Info } else { LevelFilter::Warn });

            let ggcat_params = build::GGCATParams {
                kmer_size: *ggcat_kmer_size,
                kmer_min_multiplicity: *kmer_min_multiplicity,
                minimizer_length: if minimizer_length.is_some() {
                    *minimizer_length
                } else {
                    None
                },
                no_reverse_complement: *no_reverse_complement,
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

	    init_ggcat(&Some(ggcat_params.clone()));

	    let clusters: &mut Vec<String> = &mut Vec::new();

	    // TODO seq_files should be mutable by default to avoid cloning
	    let mut seq_files_in: Vec<String> = seq_files.clone();
	    if input_list.is_some() {
		seq_files_in.append(read_input_list(input_list.as_ref().unwrap()).as_mut());
	    }

	    if external_clusters.is_some() {
		let f = std::fs::File::open(external_clusters.clone().unwrap()).unwrap();
		let mut reader = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .has_headers(false)
                    .from_reader(f);

		let mut seq_to_cluster: HashMap<String, String> = HashMap::new();
		reader.records().into_iter().for_each(|line| {
                    let record = line.unwrap();
		    seq_to_cluster.insert(record[0].to_string().clone(), record[1].to_string().clone());
		});

		if target_cluster.is_some() {
		    seq_files_in.iter().for_each(|seq| {
			let cluster = seq_to_cluster.get(seq).unwrap().clone();
			clusters.push(if cluster == target_cluster.clone().unwrap() { seq_to_cluster.get(seq).unwrap().clone() } else { seq.clone() });
		    });
		} else {
		    seq_files_in.iter().for_each(|seq| {
			clusters.push(seq_to_cluster.get(seq).unwrap().clone());
		    });
		}
	    } else {
		seq_files_in.iter().for_each(|seq| clusters.push(seq.clone()));
	    }

            build::build_pangenome_representations(
                &seq_files_in,
		clusters,
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
	    init(1 as usize, &LOG, if *verbose { LevelFilter::Info } else { LevelFilter::Warn });

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
        None => {}
    }
}
