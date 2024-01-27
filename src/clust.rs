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

fn cut_dendrogram(dendr: &kodama::Dendrogram<f32>, height: f32) -> Vec<usize> {
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
    ani_result: &Vec<(String, String, f32, f32, f32)>,
    params: KodamaParams,
) -> Vec<usize> {

    let mut flattened_similarity_matrix: Vec<f32> = ani_result.into_iter().map(|x| x.2).collect();
    let num_seqs = (0.5*(f64::sqrt((8*flattened_similarity_matrix.len() + 1) as f64) + 1.0)).round() as usize;
    let dend = kodama::linkage(&mut flattened_similarity_matrix, num_seqs, params.method);

    return cut_dendrogram(&dend, params.cutoff);
}
