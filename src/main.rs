// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use std::collections::HashSet;

use clap::Parser;
use itertools::Itertools;
use log::{info, Record, Level, LevelFilter, Metadata};

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

fn main() {
    static LOG: Logger = Logger;
    let cli = cli::Cli::parse();

    // Subcommands:
    match &cli.command {
        // Run the full pipeline
        Some(cli::Commands::Dereplicate {
            seq_files,
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
        }) => {
	    let _ = log::set_logger(&LOG).map(|()| log::set_max_level(if *verbose { LevelFilter::Info } else { LevelFilter::Warn } ));
            rayon::ThreadPoolBuilder::new()
                .num_threads(*threads as usize)
                .thread_name(|i| format!("rayon-thread-{}", i))
                .build_global()
                .unwrap();

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
                ..Default::default()
            };

            let clusters = panaani::dereplicate(
                &seq_files,
                &seq_files,
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
	    let _ = log::set_logger(&LOG).map(|()| log::set_max_level(if *verbose { LevelFilter::Info } else { LevelFilter::Warn } ));
            rayon::ThreadPoolBuilder::new()
                .num_threads(*threads as usize)
                .thread_name(|i| format!("rayon-thread-{}", i))
                .build_global()
                .unwrap();

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

            let results = dist::ani_from_fastx_files(seq_files, &Some(skani_params));

            for res in results {
                println!(
                    "{}\t{}\t{}",
                    res.0,
                    res.1,
                    res.2,
                );
            }
        }

        // Build pangenome representations from input fasta files and their clusters
        Some(cli::Commands::Build {
            seq_files,
            external_clusters,
            threads,
            memory,
            temp_dir_path,
            ggcat_kmer_size,
            kmer_min_multiplicity,
            minimizer_length,
            no_reverse_complement,
            unitig_type,
            intermediate_compression_level,
	    verbose
        }) => {
	    let _ = log::set_logger(&LOG).map(|()| log::set_max_level(if *verbose { LevelFilter::Info } else { LevelFilter::Warn } ));
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
                ..Default::default()
            };

            build::build_pangenome_representations(
                &seq_files
                    .iter()
                    .cloned()
                    .zip(external_clusters.iter().cloned())
                    .collect::<Vec<(String, String)>>(),
                &ggcat_params,
            );
        }

        // Cluster distance data created with `skani dist` or `panaani dist`.
        Some(cli::Commands::Cluster {
            dist_file,
            ani_threshold,
            linkage_method,
	    verbose
        }) => {
	    let _ = log::set_logger(&LOG).map(|()| log::set_max_level(if *verbose { LevelFilter::Info } else { LevelFilter::Warn } ));
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
                    record[2].parse().unwrap(),
                ));
                seq_names.insert(record[0].to_string());
                seq_names.insert(record[0].to_string());
            }
            res.sort_by_key(|k| (k.0.clone(), k.1.clone()));

            clust::single_linkage_cluster(&res, &Some(kodama_params));
        }
        None => {}
    }
}
