#[derive(Clone)]
pub struct KodamaParams {
    // Hierarchical clustering
    pub method: kodama::Method,
    pub cutoff: f32,
}

impl Default for KodamaParams {
    fn default() -> KodamaParams {
        KodamaParams {
            method: kodama::Method::Single,
            cutoff: 0.97,
        }
    }
}

pub fn cut_dendrogram(dendr: &kodama::Dendrogram<f32>, height: f32) -> Vec<usize> {
    let cutoff = 1.0 - height;
    let num_seqs = dendr.observations();
    let num_nodes = 2 * num_seqs - 1;

    let mut num_groups = 0;
    let mut membership = vec![None; num_nodes];

    for (cluster_index, step) in dendr.steps().iter().enumerate().rev() {
        let cluster = cluster_index + num_seqs;
        if step.dissimilarity <= cutoff {
            if membership[cluster].is_none() {
                membership[cluster] = Some(num_groups);
                num_groups += 1;
            }

            membership[step.cluster1] = membership[cluster];
            membership[step.cluster2] = membership[cluster];
        }
    }

    let mut groups = Vec::with_capacity(num_seqs);
    for group in membership.into_iter().take(num_seqs) {
        if let Some(group) = group {
            groups.push(group);
        } else {
            groups.push(num_groups);
            num_groups += 1;
        }
    }

    return groups;
}

pub fn single_linkage_cluster(
    ani_result: &Vec<skani::types::AniEstResult>,
    num_seqs: usize,
    params: KodamaParams,
) -> Vec<usize> {
    let mut ani: Vec<f32> = ani_result
        .into_iter()
        .map(|x| {
            if x.ani > 0.0 && x.ani < 1.0 && !x.ani.is_nan() {
                1.0 - x.ani
            } else {
                1.0
            }
        })
        .collect();
    let dend = kodama::linkage(&mut ani, num_seqs, params.method);

    return cut_dendrogram(&dend, params.cutoff);
}

pub fn single_linkage_cluster2(
    ani_result: &Vec<(String, String, f32, f32, f32)>,
    num_seqs: usize,
    params: KodamaParams,
) -> Vec<usize> {
    let mut ani: Vec<f32> = ani_result
        .into_iter()
        .map(|x| {
            if x.2 > 0.0 && x.2 < 1.0 && !x.2.is_nan() {
                1.0 - x.2
            } else {
                1.0
            }
        })
        .collect();
    let dend = kodama::linkage(&mut ani, num_seqs, params.method);

    return cut_dendrogram(&dend, params.cutoff);
}
