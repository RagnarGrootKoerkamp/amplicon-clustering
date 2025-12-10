use std::{
    path::PathBuf,
    sync::{atomic::AtomicUsize, Mutex, RwLock},
};

use clap::Parser;
use log::info;
use needletail;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sassy::profiles::Iupac;

#[derive(Debug, clap::Parser)]
struct Args {
    input: PathBuf,
    #[clap(short, long)]
    threshold: Option<usize>,
    #[clap(short, long)]
    alpha: f32,
    #[clap(short, long)]
    relative: Option<f32>,
}

/// The first element is the representative.
struct Cluster {
    root: Vec<u8>,
    /// (Sequence, Dist)
    seqs: Mutex<Vec<(Vec<u8>, pa_types::Cost)>>,
}

impl Cluster {
    fn new(seq: Vec<u8>) -> Self {
        Self {
            root: seq.clone(),
            seqs: Mutex::new(vec![(seq, 0)]),
        }
    }
    fn representative(&self) -> &[u8] {
        &self.root
    }
    fn push(&self, seq: Vec<u8>, dist: pa_types::Cost) {
        self.seqs.lock().unwrap().push((seq, dist));
    }
}

fn main() {
    env_logger::Builder::from_default_env()
        .format_timestamp_millis()
        .init();

    let args = Args::parse();
    let mut reader = needletail::parse_fastx_file(&args.input).unwrap();
    let mut reads: Vec<Vec<u8>> = vec![];
    while let Some(record) = reader.next() {
        let record = record.unwrap();
        reads.push(record.seq().into_owned());
    }

    let num_reads = reads.len();
    info!("Number of reads: {}", num_reads);
    let total_len: usize = reads.iter().map(|r| r.len()).sum();
    let avg_len = total_len as f64 / num_reads as f64;
    info!("Average read length: {:.1}", avg_len);

    let global_clusters: RwLock<Vec<Cluster>> = RwLock::new(vec![]);

    // A*PA2
    // let params = astarpa2::AstarPa2Params::full();
    // let mut aligner = params.make_aligner(false);

    // sassy
    let mut searcher = sassy::Searcher::<Iupac>::new_rc_with_overhang(args.alpha);

    let done = AtomicUsize::new(0);

    let start = std::time::Instant::now();

    let threshold = |read: &Vec<u8>| -> usize {
        let threshold = if let Some(r) = args.relative {
            (r * read.len() as f32) as usize
        } else {
            args.threshold
                .expect("Either --threshold or --relative must be specified")
        };
        threshold
    };

    reads
        .into_par_iter()
        .for_each_with(searcher.clone(), |searcher, read: Vec<u8>| {
            let threshold = threshold(&read);

            'clusters_changed: loop {
                let clusters = global_clusters.read().unwrap();

                let mut best = (i32::MAX, 0);
                for (i, cluster) in clusters.iter().enumerate() {
                    let repr = cluster.representative();
                    // let cost = aligner.align(&read, repr).0;

                    let matches = searcher.search(&read, repr, threshold);
                    let best_match = matches.iter().min_by_key(|m| m.cost);
                    if let Some(best_m) = best_match {
                        let cost = best_m.cost;
                        best = best.min((cost, i));
                    }
                }
                if best.0 < threshold as i32 {
                    let j = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let trhp = (j + 1) as f64 / start.elapsed().as_secs_f64();
                    eprintln!(
                        "{j:>6}/{num_reads} {trhp:>5.0}/s Best cost: {:>4} (len {:>4}) => append to {:>3}",
                        best.0,
                        read.len(),
                        best.1
                    );
                    clusters[best.1].push(read, best.0);
                    break;
                } else {
                    let old_len = clusters.len();
                    drop(clusters);
                    let mut clusters = global_clusters.write().unwrap();
                    if clusters.len() != old_len {
                        // clusters changed, retry
                        continue 'clusters_changed;
                    }

                    let j = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let trhp = (j + 1) as f64 / start.elapsed().as_secs_f64();
                    eprintln!(
                        "{j:>6}/{num_reads} {trhp:>5.0}/s Best cost: {:>4} (len {:>4}) => start new cluster {:>3}",
                        best.0,
                        read.len(),
                        clusters.len()
                    );
                    clusters.push(Cluster::new(read));
                    break;
                }
            }
        });

    let duration = start.elapsed();
    eprintln!(
        "Clustering completed in {:.1} seconds",
        duration.as_secs_f64()
    );
    eprintln!("Cluster sizes");
    let clusters = global_clusters.into_inner().unwrap();
    for (i, cluster) in clusters.iter().enumerate() {
        let seqs = cluster.seqs.lock().unwrap();
        let mut lens = seqs.iter().map(|(s, _)| s.len()).collect::<Vec<usize>>();
        lens.sort_unstable();
        eprintln!(
            "Cluster {:>3}: {:>4} seqs with lens {lens:?}",
            i,
            seqs.len()
        );
    }

    // For each cluster, compute the distance from each seq to all others, and
    // take the one with lowest median.
    for (ci, cluster) in clusters.iter().enumerate() {
        let seqs = cluster.seqs.lock().unwrap();
        let mut distances: Vec<Vec<i32>> = vec![vec![0; seqs.len()]; seqs.len()];
        eprintln!("Cluster {ci}:");
        let mut best_med = i32::MAX;
        for (j, (seq_j, _)) in seqs.iter().enumerate() {
            for (k, (seq_k, _)) in seqs.iter().enumerate() {
                let threshold = threshold(&seq_k);
                let matches = searcher.search(seq_k, seq_j, threshold);
                let best_match = matches.iter().min_by_key(|m| m.cost);
                let dist = best_match.map_or(i32::MAX, |m| m.cost);
                distances[j][k] = dist;
            }
            distances[j].sort_unstable();
            let med = distances[j][distances[j].len() / 2];
            eprint!(" {med}");
            best_med = best_med.min(med);
        }
        eprintln!(" => {best_med}");
    }
}
