use std::{
    cmp::Reverse,
    collections::HashMap,
    path::PathBuf,
    sync::{
        atomic::{AtomicUsize, Ordering::Relaxed},
        Mutex,
    },
};

use clap::Parser;
use log::info;
use needletail;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use sassy::profiles::Iupac;

#[derive(Debug, clap::Parser)]
struct Args {
    input: PathBuf,
    #[clap(short, long)]
    threshold: Option<usize>,
    #[clap(short, long, default_value_t = 0.2)]
    alpha: f32,
    #[clap(short, long, default_value = "0.15")]
    relative: Option<f32>,

    #[clap(short = 'j', long)]
    threads: Option<usize>,

    /// Discard too short reads
    #[clap(long, default_value_t = 450)]
    min_len: usize,
    /// Discard too long reads
    #[clap(long, default_value_t = 550)]
    max_len: usize,

    /// Discard reads in too smal clusters
    #[clap(long, default_value_t = 5)]
    min_cluster_size: usize,
}

/// The first element is the representative.
struct Cluster<'r> {
    /// (Sequence, Dist)
    seqs: Vec<(&'r [u8], pa_types::Cost)>,
}

impl<'r> Cluster<'r> {
    fn new(seq: &'r [u8]) -> Self {
        Self {
            seqs: vec![(seq, 0)],
        }
    }
    fn push(&mut self, seq: &'r [u8], dist: pa_types::Cost) {
        self.seqs.push((seq, dist));
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

    let global_clusters: Vec<Mutex<Cluster>> = (0..5000)
        .map(|_| Mutex::new(Cluster::new(&reads[0])))
        .collect();

    // sassy
    let searcher = sassy::Searcher::<Iupac>::new_rc_with_overhang(args.alpha)
        .only_best_match()
        .without_trace();

    let start = std::time::Instant::now();

    let threshold = |read: &[u8]| -> usize {
        let threshold = if let Some(r) = args.relative {
            (r * read.len() as f32) as usize
        } else {
            args.threshold
                .expect("Either --threshold or --relative must be specified")
        };
        threshold
    };

    let num_reads = reads.len();
    let idx = AtomicUsize::new(0);

    eprintln!("Create channel");
    let (tx, _rx) = tokio::sync::broadcast::channel::<Vec<u8>>(10000);
    drop(_rx);

    std::thread::scope(|scope| {
        for tid in 0..args.threads.unwrap_or(6) {
            let mut searcher = searcher.clone();
            let global_clusters = &global_clusters;
            let mut local_roots = vec![];
            let reads = &reads;
            let idx = &idx;
            let tx = tx.clone();
            let mut rx = tx.subscribe();
            scope.spawn(move || {
                eprintln!("Spawned thread {tid}");
                loop {
                    let idx = idx.fetch_add(1, Relaxed);
                    let trhp = (idx + 1) as f64 / start.elapsed().as_secs_f64();
                    let Some(read) = reads.get(idx) else {break;};
                    assert!(!read.is_empty());
                    let threshold = threshold(&read);

                    'clusters_changed: loop {
                        while let Ok(new_root) = rx.try_recv() {
                            // eprintln!("{tid:>4} received new root of len {}", new_root.len());
                            local_roots.push(new_root);
                        }
                        let local_roots: Vec<&[u8]> = local_roots.iter().map(|r| r.as_slice()).collect();


                        let mut best = (i32::MAX, 0);

                        let mut matches = searcher.search_texts(&read, &local_roots, threshold, Some(threshold/2));
                        matches.sort_by_key(|m| (m.1.cost, m.0));
                        let best_match = matches.get(0);
                        if let Some(best_m) = best_match {
                            let cost = best_m.1.cost;
                            best = best.min((cost, best_m.0));
                        }

                        if matches.len() > 1 {
                            eprintln!("{tid: >4} matches: {:?}", matches.iter().map(|m| m.1.cost).collect::<Vec<_>>().as_slice().iter().collect::<Vec<_>>());
                        }

                        if best.0 < threshold as i32 {
                            eprintln!(
                                "{tid:>4} {idx:>6}/{num_reads} {trhp:>5.0}/s Best cost: {:>4} (len {:>4}) => append to {:>3}",
                                best.0,
                                read.len(),
                                best.1
                            );
                            global_clusters[best.1].lock().unwrap().push(read, best.0);
                            break;
                        } else {
                            if !rx.is_empty() {
                                eprintln!("{tid:>4} {idx:>6} retry",);
                                continue 'clusters_changed;
                            }

                            eprintln!(
                                "{tid:>4} {idx:>6}/{num_reads} {trhp:>5.0}/s Best cost: {:>4} (len {:>4}) => start new cluster {}",
                                best.0,
                                read.len(),
                                local_roots.len()
                            );
                            tx.send(read.to_vec()).unwrap();
                            break;
                        }
                    }
                }
            });
        }
    });

    let duration = start.elapsed();
    eprintln!(
        "Clustering completed in {:.1} seconds",
        duration.as_secs_f64()
    );

    let mut clusters: Vec<_> = global_clusters
        .into_iter()
        .map(|c| c.into_inner().unwrap())
        .collect();
    {
        eprintln!("Cluster sizes");
        // Sort clusters by size.
        clusters.sort_by_cached_key(|c| Reverse(c.seqs.len()));
        // retain only sufficiently large clusters
        clusters.retain(|c| c.seqs.len() >= args.min_cluster_size);
        for (i, cluster) in clusters.iter().enumerate() {
            let seqs = &cluster.seqs;
            let mut lens = seqs.iter().map(|(s, _)| s.len()).collect::<Vec<usize>>();
            lens.sort_unstable();
            eprintln!(
                "Cluster {:>3}: {:>4} seqs with lens {lens:?}",
                i,
                seqs.len()
            );
        }
    }

    let mut trace_searcher = searcher.clone().with_trace();

    // For each cluster, compute the distance from each seq to all others, and
    // take the one with lowest median.
    for (ci, cluster) in clusters.iter().enumerate() {
        let seqs = &cluster.seqs;
        eprint!("Cluster {ci}:");
        let mut best_med = (i32::MAX, 0);
        for (j, (seq_j, _)) in seqs.iter().enumerate() {
            let mut dists = seqs
                .par_iter()
                .map_with(searcher.clone(), |searcher, (seq_k, _)| {
                    let threshold = threshold(seq_j);
                    let matches = searcher.search(seq_k, seq_j, threshold);
                    let best_match = matches.iter().min_by_key(|m| m.cost);
                    let dist = best_match.map_or(i32::MAX, |m| m.cost);
                    dist
                })
                .collect::<Vec<i32>>();
            dists.sort_unstable();
            let med = dists[dists.len() / 2];
            // eprint!(" {med}");
            best_med = best_med.min((med, j));
        }
        let (best_med_dist, best_idx) = best_med;
        eprintln!(" => best median dist {best_med_dist}");
        let root = &seqs[best_idx].0;
        let threshold = threshold(&root) + 20;
        // print alignment of root to all others
        let mut best_matches = vec![];
        for (seq, _) in seqs.iter() {
            let matches = trace_searcher.search(root, seq, threshold);
            if let Some(best_match) = matches.into_iter().min_by_key(|m| m.cost) {
                best_matches.push((seq, best_match));
            }
        }
        best_matches.sort_unstable_by_key(|(_, m)| m.cost);
        for (seq, m) in best_matches {
            let pretty = m.pretty_print(
                None,
                root,
                seq,
                sassy::pretty_print::PrettyPrintDirection::Pattern,
                10,
                sassy::pretty_print::PrettyPrintStyle::Compact,
            );
            eprintln!("{pretty}");
        }
    }
}
