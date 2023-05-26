use std::time::Duration;

use a0kzg::{
    kzg::{PowTauG1Affine, PowTauG2Affine},
    Prover, Scalar, TrustedTau, UntrustedProver, UntrustedVerifier, Verifier,
};
use bls12_381::{G1Affine, G2Affine};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ff::Field;
use once_cell::sync::Lazy;

static TAU: Lazy<TrustedTau> = Lazy::new(|| TrustedTau(Scalar::random(rand::thread_rng())));

static SIZES: &[(usize, &[usize])] = &[
    (32, &[1]),
    (64, &[1]),
    (128, &[1]),
    (256, &[1]),
    (512, &[1]),
];

static POINTS: Lazy<Vec<(Scalar, Scalar)>> = Lazy::new(|| {
    (0..1024)
        .map(|_| {
            (
                Scalar::random(rand::thread_rng()),
                Scalar::random(rand::thread_rng()),
            )
        })
        .collect::<Vec<_>>()
});

static G1_AFFINE: Lazy<Vec<G1Affine>> =
    Lazy::new(|| TAU.g1_affine_iter().take(POINTS.len() - 1).collect());
static G2_AFFINE: Lazy<Vec<G2Affine>> =
    Lazy::new(|| TAU.g2_affine_iter().take(POINTS.len() - 1).collect());

fn zkg(c: &mut Criterion) {
    let affine = (
        PowTauG1Affine(|| G1_AFFINE.iter()),
        PowTauG2Affine(|| G2_AFFINE.iter()),
    );

    let prover = UntrustedProver::new(affine.clone());
    let verifier = UntrustedVerifier::new(affine.clone());

    let mut group = c.benchmark_group("zkg affine untrusted");
    for &(set_size, proof_sizes) in SIZES.iter() {
        let id = BenchmarkId::from_parameter(format!("commit {set_size}"));
        group.throughput(Throughput::Elements(set_size as u64));
        group.bench_with_input(id, &POINTS[0..set_size], |b, points| {
            b.iter(|| prover.poly_commitment_from_set(points))
        });

        let v = prover.poly_commitment_from_set(&POINTS[0..set_size]);
        let (p, c) = v;
        let c = c.into();

        for &proof_size in proof_sizes {
            let id = BenchmarkId::from_parameter(format!("prove {set_size}:{proof_size}"));
            group.throughput(Throughput::Elements(proof_size as u64));
            group.bench_with_input(id, &(&p, &POINTS[0..proof_size]), |b, &(poly, points)| {
                b.iter(|| prover.prove(poly.clone(), points))
            });
            let pi = prover.prove(p.clone(), &POINTS[0..proof_size]);
            let pi = pi.into();

            let id = BenchmarkId::from_parameter(format!("verify {set_size}:{proof_size}"));
            group.throughput(Throughput::Elements(proof_size as u64));
            group.bench_with_input(
                id,
                &(&c, &POINTS[0..proof_size]),
                |b, &(commitment, points)| b.iter(|| verifier.verify(commitment, points, &pi)),
            );
            assert!(verifier.verify(&c, &POINTS[0..proof_size], &pi));
        }
    }
    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(60).measurement_time(Duration::from_secs(20));
    targets = zkg
}
criterion_main!(benches);
