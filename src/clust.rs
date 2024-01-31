// panaani: Pangenome-aware dereplication of bacterial genomes into ANI clusters
//
// Copyright (c) Tommi Mäklin <tommi 'at' maklin.fi>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
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
    ani_result: &Vec<(String, String, f32)>,
    opt: &Option<KodamaParams>,
) -> Vec<usize> {

    let params = opt.clone().unwrap_or(KodamaParams::default());
    let mut flattened_similarity_matrix: Vec<f32> = ani_result.into_iter().map(|x| 1.0 - x.2).collect();
    let num_seqs = (0.5*(f64::sqrt((8*flattened_similarity_matrix.len() + 1) as f64) + 1.0)).round() as usize;
    let dend = kodama::linkage(&mut flattened_similarity_matrix, num_seqs, params.method);

    return cut_dendrogram(&dend, params.cutoff);
}
