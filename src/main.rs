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

mod build;
mod clust;
mod dist;

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
        // Input files
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

        // Clustering parameters
        #[arg(long = "ani-threshold", default_value_t = 0.97)]
        ani_threshold: f32,

        #[arg(long = "linkage-method", required = false)]
        linkage_method: Option<String>,

        // ANI estimation parameters
        #[arg(long = "skani-kmer-size", default_value_t = 15)]
        skani_kmer_size: u8,

        #[arg(long = "kmer-subsampling-rate", default_value_t = 30)]
        kmer_subsampling_rate: u16,

        #[arg(long = "marker-compression-factor", default_value_t = 1000)]
        marker_compression_factor: u16,

        #[arg(long = "rescue-small", default_value_t = false)]
        rescue_small: bool,

        #[arg(long = "clip-tails", default_value_t = false)]
        clip_tails: bool,

        #[arg(long = "median", default_value_t = false)]
        median: bool,

        #[arg(long = "adjust-ani", default_value_t = false)]
        adjust_ani: bool,

        #[arg(long = "min-af", default_value_t = 0.)]
        min_aligned_frac: f64,

        // de Bruijn graph construction parameters
        #[arg(long = "ggcat-kmer-size", default_value_t = 51)]
        ggcat_kmer_size: u32,

        #[arg(long = "min-kmer-count", default_value_t = 1)]
        kmer_min_multiplicity: u64,

        #[arg(long = "minimzer-length", required = false)]
        minimizer_length: Option<usize>,

        #[arg(long = "no-rc", default_value_t = false)]
        no_reverse_complement: bool,

        #[arg(long = "unitig-type", required = false)]
        unitig_type: Option<String>,

        #[arg(long = "prefer-memory", default_value_t = true)]
        prefer_memory: bool,

        #[arg(long = "intermediate-compression", required = false)]
        intermediate_compression_level: Option<u32>,
    },

    Dist {
        // Input files
        #[arg(group = "input", required = true)]
        seq_files: Vec<String>,

        // Resources
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        threads: u32,

        // ANI estimation parameters
        #[arg(long = "skani-kmer-size", default_value_t = 15)]
        skani_kmer_size: u8,

        #[arg(long = "kmer-subsampling-rate", default_value_t = 30)]
        kmer_subsampling_rate: u16,

        #[arg(long = "marker-compression-factor", default_value_t = 1000)]
        marker_compression_factor: u16,

        #[arg(long = "rescue-small", default_value_t = false)]
        rescue_small: bool,

        #[arg(long = "clip-tails", default_value_t = false)]
        clip_tails: bool,

        #[arg(long = "median", default_value_t = false)]
        median: bool,

        #[arg(long = "adjust-ani", default_value_t = false)]
        adjust_ani: bool,

        #[arg(long = "min-af", default_value_t = 0.)]
        min_aligned_frac: f64,
    },
    Build {
        // Input files
        #[arg(group = "input", required = true)]
        seq_files: Vec<String>,

	#[arg(long = "external-clustering", required = true)]
        external_clusters: Vec<String>,

        // Resources
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        threads: u32,

        #[arg(short = 'm', long = "memory", default_value_t = 4)]
        memory: u32,

        #[arg(long = "tmp-dir", required = false)]
        temp_dir_path: Option<String>,

        // de Bruijn graph construction parameters
        #[arg(long = "ggcat-kmer-size", default_value_t = 51)]
        ggcat_kmer_size: u32,

        #[arg(long = "min-kmer-count", default_value_t = 1)]
        kmer_min_multiplicity: u64,

        #[arg(long = "minimzer-length", required = false)]
        minimizer_length: Option<usize>,

        #[arg(long = "no-rc", default_value_t = false)]
        no_reverse_complement: bool,

        #[arg(long = "unitig-type", required = false)]
        unitig_type: Option<String>,

        #[arg(long = "prefer-memory", default_value_t = true)]
        prefer_memory: bool,

        #[arg(long = "intermediate-compression", required = false)]
        intermediate_compression_level: Option<u32>,
    },
    Cluster {
        #[arg(group = "input")]
        dist_file: String,

        // Clustering parameters
        #[arg(long = "ani-threshold", default_value_t = 0.97)]
        ani_threshold: f32,

        #[arg(long = "linkage-method", required = false)]
        linkage_method: Option<String>,
    },
}

fn main() {
    println!("panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters");
    let cli = Cli::parse();

    // Subcommands:
    match &cli.command {
        // Run the full pipeline
        Some(Commands::Dereplicate {
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
            prefer_memory,
            intermediate_compression_level,
            threads,
            memory,
            temp_dir_path,
            ani_threshold,
        }) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(*threads as usize)
                .thread_name(|i| format!("rayon-thread-{}", i))
                .build_global()
                .unwrap();

            let params: panaani::PanaaniParams = panaani::PanaaniParams {
                batch_step: *batch_step,
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
                prefer_memory: *prefer_memory,
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
                Some(params),
                Some(skani_params),
                Some(kodama_params),
                Some(ggcat_params),
            );
            let n_clusters = clusters
                .iter()
                .map(|x| x.1.clone())
                .unique()
                .collect::<Vec<String>>()
                .len();

            println!("Created {} clusters", n_clusters);
            for cluster in clusters {
                println!("{}\t{}", cluster.0, cluster.1);
            }
        }

        // Calculate distances between some input fasta files
        Some(Commands::Dist {
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
        }) => {
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

            let results = dist::ani_from_fastx_files(seq_files, &skani_params);

            for res in results {
                println!(
                    "{}\t{}\t{}\t{}\t{}",
                    res.ref_file,
                    res.query_file,
                    res.ani,
                    res.align_fraction_ref,
                    res.align_fraction_query
                );
            }
        }

        // Build pangenome representations from input fasta files and their clusters
        Some(Commands::Build {
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
            prefer_memory,
            intermediate_compression_level,
        }) => {
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
                prefer_memory: *prefer_memory,
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
                &seq_files.iter().cloned().zip(external_clusters.iter().cloned()).collect(),
                &ggcat_params,
            );
        }

        // Cluster distance data created with `skani dist` or `panaani dist`.
        Some(Commands::Cluster {
            dist_file,
            ani_threshold,
            linkage_method,
        }) => {
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
            let mut res: Vec<(String, String, f32, f32, f32)> = Vec::new();
            for line in reader.records().into_iter() {
                let record = line.unwrap();
                res.push((
                    record[0].to_string().clone(),
                    record[1].to_string().clone(),
                    record[2].parse().unwrap(),
                    record[3].parse().unwrap(),
                    record[4].parse().unwrap(),
                ));
                seq_names.insert(record[0].to_string());
                seq_names.insert(record[0].to_string());
            }
            res.sort_by_key(|k| (k.0.clone(), k.1.clone()));

            clust::single_linkage_cluster(&res, kodama_params);
        }
        None => {}
    }
}
