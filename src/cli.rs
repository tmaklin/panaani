// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(version)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
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
        #[arg(
            short = 'b',
            long = "batch-step",
            default_value_t = 50,
            help_heading = "Dereplication"
        )]
        batch_step: usize,

        #[arg(
            long = "batch-step-strategy",
            default_value = "double",
            help_heading = "Dereplication"
        )]
        batch_step_strategy: String,

        #[arg(
            long = "max-iters",
            default_value_t = 10,
            help_heading = "Dereplication"
        )]
        max_iters: usize,

        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,

        // ANI estimation parameters
        #[arg(
            long = "skani-kmer-size",
            default_value_t = 15,
            help_heading = "ANI estimation"
        )]
        skani_kmer_size: u8,

        #[arg(
            long = "kmer-subsampling-rate",
            default_value_t = 30,
            help_heading = "ANI estimation"
        )]
        kmer_subsampling_rate: u16,

        #[arg(
            long = "marker-compression-factor",
            default_value_t = 1000,
            help_heading = "ANI estimation"
        )]
        marker_compression_factor: u16,

        #[arg(
            long = "min-af",
            default_value_t = 0.15,
            help_heading = "ANI estimation"
        )]
        min_aligned_frac: f64,

        #[arg(
            long = "rescue-small",
            default_value_t = false,
            help_heading = "ANI estimation"
        )]
        rescue_small: bool,

        #[arg(
            long = "clip-tails",
            default_value_t = false,
            help_heading = "ANI estimation"
        )]
        clip_tails: bool,

        #[arg(
            long = "median",
            default_value_t = false,
            help_heading = "ANI estimation"
        )]
        median: bool,

        #[arg(
            long = "adjust-ani",
            default_value_t = false,
            help_heading = "ANI estimation"
        )]
        adjust_ani: bool,

        // Clustering parameters
        #[arg(
            long = "ani-threshold",
            default_value_t = 0.97,
            help_heading = "ANI clustering"
        )]
        ani_threshold: f32,

        #[arg(
            long = "linkage-method",
            required = false,
            help_heading = "ANI clustering"
        )]
        linkage_method: Option<String>,

        // de Bruijn graph construction parameters
        #[arg(
            long = "ggcat-kmer-size",
            default_value_t = 51,
            help_heading = "Pangenome construction"
        )]
        ggcat_kmer_size: u32,

        #[arg(
            long = "min-kmer-count",
            default_value_t = 1,
            help_heading = "Pangenome construction"
        )]
        kmer_min_multiplicity: u64,

        #[arg(
            long = "minimzer-length",
            required = false,
            help_heading = "Pangenome construction"
        )]
        minimizer_length: Option<usize>,

        #[arg(
            long = "no-rc",
            default_value_t = false,
            help_heading = "Pangenome construction"
        )]
        no_reverse_complement: bool,

        #[arg(
            long = "unitig-type",
            required = false,
            help_heading = "Pangenome construction"
        )]
        unitig_type: Option<String>,

        #[arg(
            long = "intermediate-compression",
            required = false,
            help_heading = "Pangenome construction"
        )]
        intermediate_compression_level: Option<u32>,
    },

    Dist {
        // Input files
        #[arg(group = "input", required = true)]
        seq_files: Vec<String>,

        // Resources
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        threads: u32,

        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,

        // ANI estimation parameters
        #[arg(
            long = "skani-kmer-size",
            default_value_t = 15,
            help_heading = "ANI estimation"
        )]
        skani_kmer_size: u8,

        #[arg(
            long = "kmer-subsampling-rate",
            default_value_t = 30,
            help_heading = "ANI estimation"
        )]
        kmer_subsampling_rate: u16,

        #[arg(
            long = "marker-compression-factor",
            default_value_t = 1000,
            help_heading = "ANI estimation"
        )]
        marker_compression_factor: u16,

        #[arg(
            long = "min-af",
            default_value_t = 0.15,
            help_heading = "ANI estimation"
        )]
        min_aligned_frac: f64,

        #[arg(
            long = "rescue-small",
            default_value_t = false,
            help_heading = "ANI estimation"
        )]
        rescue_small: bool,

        #[arg(
            long = "clip-tails",
            default_value_t = false,
            help_heading = "ANI estimation"
        )]
        clip_tails: bool,

        #[arg(
            long = "median",
            default_value_t = false,
            help_heading = "ANI estimation"
        )]
        median: bool,

        #[arg(
            long = "adjust-ani",
            default_value_t = false,
            help_heading = "ANI estimation"
        )]
        adjust_ani: bool,
    },
    Build {
        // Input files
        #[arg(group = "input", required = true)]
        seq_files: Vec<String>,

        #[arg(long = "external-clustering", required = true, help_heading = "Input")]
        external_clusters: Option<String>,

	#[arg(long = "target", required = false, help_heading = "Input")]
        target_cluster: Option<String>,

        // Resources
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        threads: u32,

        #[arg(short = 'm', long = "memory", default_value_t = 4)]
        memory: u32,

        #[arg(long = "tmp-dir", required = false)]
        temp_dir_path: Option<String>,

        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,

        // de Bruijn graph construction parameters
        #[arg(
            long = "ggcat-kmer-size",
            default_value_t = 51,
            help_heading = "Pangenome construction"
        )]
        ggcat_kmer_size: u32,

        #[arg(
            long = "min-kmer-count",
            default_value_t = 1,
            help_heading = "Pangenome construction"
        )]
        kmer_min_multiplicity: u64,

        #[arg(
            long = "minimzer-length",
            required = false,
            help_heading = "Pangenome construction"
        )]
        minimizer_length: Option<usize>,

        #[arg(
            long = "no-rc",
            default_value_t = false,
            help_heading = "Pangenome construction"
        )]
        no_reverse_complement: bool,

        #[arg(
            long = "unitig-type",
            required = false,
            help_heading = "Pangenome construction"
        )]
        unitig_type: Option<String>,

        #[arg(
            long = "intermediate-compression",
            required = false,
            help_heading = "Pangenome construction"
        )]
        intermediate_compression_level: Option<u32>,
    },
    Cluster {
        #[arg(group = "input")]
        dist_file: String,

        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,

        // Clustering parameters
        #[arg(
            long = "ani-threshold",
            default_value_t = 0.97,
            help_heading = "ANI estimation"
        )]
        ani_threshold: f32,

        #[arg(
            long = "linkage-method",
            required = false,
            help_heading = "ANI estimation"
        )]
        linkage_method: Option<String>,
    },
}
