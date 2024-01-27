// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
use std::sync::mpsc::channel;

use itertools::Itertools;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;

#[derive(Clone)]
pub struct SkaniParams {
    // k-mer sketching
    pub kmer_size: u8,
    pub kmer_subsampling_rate: u16,
    pub marker_compression_factor: u16,
    pub rescue_small: bool,

    // ANI estimation
    pub clip_tails: bool,
    pub median: bool,
    pub adjust_ani: bool,

    // Results reporting
    pub min_aligned_frac: f64,
    pub bootstrap_ci: bool,
}

impl Default for SkaniParams {
    fn default() -> SkaniParams {
        SkaniParams {
            kmer_size: 15,
            kmer_subsampling_rate: 30,
            marker_compression_factor: 1000,
            rescue_small: false,

            clip_tails: false,
            median: false,
            adjust_ani: false,

            min_aligned_frac: 0.0,
            bootstrap_ci: false,
        }
    }
}

pub fn ani_from_fastx_files(
    fastx_files: &Vec<String>,
    skani_params: &SkaniParams,
) -> Vec<skani::types::AniEstResult> {

    let sketch_params = skani::params::SketchParams::new(
        skani_params.marker_compression_factor as usize,
        skani_params.kmer_subsampling_rate as usize,
        skani_params.kmer_size as usize,
        false,
        false,
    );
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
        min_aligned_frac: skani_params.min_aligned_frac,
        keep_refs: false,
        est_ci: skani_params.bootstrap_ci,
        learned_ani: skani_params.adjust_ani,
        detailed_out: false,
        rescue_small: skani_params.rescue_small,
        distance: true,
    };

    let sketches = skani::file_io::fastx_to_sketches(
        &fastx_files,
        &sketch_params,
        true,
    );
    let adjust_ani = skani::regression::get_model(sketch_params.c, false);

    let (sender, receiver) = channel();
    sketches
        .iter()
        .cloned()
        .combinations(2)
        .par_bridge()
        .for_each_with(sender, |s, pair| {
            let _ = s.send(skani::chain::chain_seeds(
		pair.first().unwrap(),
		pair.last().unwrap(),
		skani::chain::map_params_from_sketch(
		    pair.first().unwrap(),
		    false,
		    &cmd_params,
		    &adjust_ani)
            ));
        });

    let mut ani_result: Vec<skani::types::AniEstResult> = receiver.iter().collect();

    // Ensure output order is same regardless of parallelization
    ani_result.sort_by_key(|k| (k.ref_file.clone(), k.query_file.clone()));
    return ani_result;
}
