use std::{
    cmp::Reverse,
    collections::HashMap,
    path::PathBuf,
    sync::{atomic::AtomicUsize, Mutex, RwLock},
};

use clap::Parser;
use log::info;
use needletail;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
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

    /// Discard too short reads
    #[clap(long, default_value_t = 425)]
    min_len: usize,
    /// Discard too long reads
    #[clap(long, default_value_t = 550)]
    max_len: usize,

    /// Discard reads in too smal clusters
    #[clap(long, default_value_t = 5)]
    min_cluster_size: usize,
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

    // Print histogram of read lengths
    {
        const B: usize = 25;
        let mut lens = HashMap::new();
        for read in &reads {
            *lens.entry(read.len() / B).or_default() += 1;
        }
        let mut lens_vec: Vec<(usize, usize)> = lens.into_iter().collect();
        lens_vec.sort_unstable();
        eprintln!("Read length histogram:");
        for (bucket, count) in lens_vec {
            eprintln!(" {:>4}..{:>4}: {}", bucket * B, (bucket + 1) * B, count);
        }
    }

    // Filter reads by length
    {
        eprintln!("Filtering for length {}..{}", args.min_len, args.max_len);
        reads.retain(|r| args.min_len <= r.len() && r.len() <= args.max_len);
    }

    // Sort reads by decreasing length
    reads.sort_unstable_by_key(|r| Reverse(r.len()));

    let global_clusters: RwLock<Vec<Cluster>> = RwLock::new(vec![]);

    // A*PA2
    // let params = astarpa2::AstarPa2Params::full();
    // let mut aligner = params.make_aligner(false);

    // sassy
    let mut searcher = sassy::Searcher::<Iupac>::new_rc_with_overhang(args.alpha);

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

    let num_reads = reads.len();
    let reads = Mutex::new(reads);
    let done = AtomicUsize::new(0);

    (0..num_reads)
        .into_par_iter()
        .for_each_with(searcher.clone(), |searcher, _| {
            let j = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let trhp = (j + 1) as f64 / start.elapsed().as_secs_f64();
            let read = std::mem::take(&mut reads.lock().unwrap()[j]);
            assert!(!read.is_empty());

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

    let mut clusters = global_clusters.into_inner().unwrap();
    {
        eprintln!("Cluster sizes");
        // Sort clusters by size.
        clusters.sort_by_cached_key(|c| Reverse(c.seqs.lock().unwrap().len()));
        // retain only sufficiently large clusters
        clusters.retain(|c| c.seqs.lock().unwrap().len() >= args.min_cluster_size);
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
    }

    // For each cluster, compute the distance from each seq to all others, and
    // take the one with lowest median.
    for (ci, cluster) in clusters.iter().enumerate() {
        let seqs = cluster.seqs.lock().unwrap();
        eprintln!("Cluster {ci}:");
        let mut best_med = (i32::MAX, 0);
        for (j, (seq_j, _)) in seqs.iter().enumerate() {
            let mut dists = seqs
                .par_iter()
                .map_with(searcher.clone(), |searcher, (seq_k, _)| {
                    let threshold = threshold(&seq_k);
                    let matches = searcher.search(seq_k, seq_j, threshold);
                    let best_match = matches.iter().min_by_key(|m| m.cost);
                    let dist = best_match.map_or(i32::MAX, |m| m.cost);
                    dist
                })
                .collect::<Vec<i32>>();
            dists.sort_unstable();
            let med = dists[dists.len() / 2];
            eprint!(" {med}");
            best_med = best_med.min((med, j));
        }
        let (best_med_dist, best_idx) = best_med;
        eprintln!(" => best {best_med_dist}");
    }
}
