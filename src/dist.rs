use std::sync::mpsc::channel;

use itertools::Itertools;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;

pub fn compare_fastx_files(reference: &String, query: &String,
		       sketch_params: &skani::params::SketchParams,
		       cmd_params: &skani::params::CommandParams) -> skani::types::AniEstResult{
    let sketches = skani::file_io::fastx_to_sketches(&vec![reference.clone(), query.clone()], &sketch_params, true);
    let adjust_ani = skani::regression::get_model(sketch_params.c, false);
    let map_params = skani::chain::map_params_from_sketch(sketches.first().unwrap(), false, &cmd_params, &adjust_ani);
    return skani::chain::chain_seeds(sketches.first().unwrap(), sketches.last().unwrap(), map_params);
}

pub fn ani_from_fastx_files(fastx_files: &Vec<String>) -> Vec<skani::types::AniEstResult>{
    // skani parameters
    let m = 1000;
    let c = 30;
    let k = 15;
    let sketch_params = skani::params::SketchParams::new(m, c, k, false, false);
    let cmd_params = skani::params::CommandParams {
        screen: false,
        screen_val: 0.00,
        mode: skani::params::Mode::Dist,
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
        min_aligned_frac: 0.0,
        keep_refs: false,
        est_ci: false,
        learned_ani: false,
        detailed_out: false,
	rescue_small: false,
	distance: true,
    };

    let (sender, receiver) = channel();
    fastx_files.iter().cloned().combinations(2).par_bridge().for_each_with(sender, |s, pair| {
	let _ = s.send(compare_fastx_files(pair.first().unwrap(),
					   pair.last().unwrap(),
					   &sketch_params, &cmd_params));
    });
    let mut ani_result: Vec<skani::types::AniEstResult> = receiver.iter().collect();

    // Ensure output order is same regardless of parallelization
    ani_result.sort_by_key(|k| (k.ref_file.clone(), k.query_file.clone()));
    return ani_result;
}
