"use strict";

// Standard normal via Box-Muller. Not seeded to match numpy bit-for-bit --
// each Run produces a fresh stochastic realization, same spirit as the
// Pyodide site (numpy there isn't seed-matched to anything external either).
function randn() {
  let u = 0,
    v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

function randnVec(n) {
  return Array.from({ length: n }, randn);
}

function randVec(n) {
  return Array.from({ length: n }, Math.random);
}

// Sample from N(mean, cov) given mean vector and covariance matrix (not its
// Cholesky factor -- factorization happens here so callers can pass K
// straight from a kernel function).
function mvnSample(mean, cov) {
  const L = window.linalg.cholesky(cov);
  const z = randnVec(mean.length);
  const Lz = window.linalg.matVec(L, z);
  return mean.map((m, i) => m + Lz[i]);
}

// Laplace(0, scale) via inverse-CDF sampling.
function laplaceSample(scale = 1) {
  const u = Math.random() - 0.5;
  return -scale * Math.sign(u) * Math.log(1 - 2 * Math.abs(u));
}

// Standard normal truncated to [lo, hi], via simple rejection sampling.
// Fine for the moderate, not-too-narrow ranges these figures use (e.g.
// [-1, 1]) -- not efficient for extreme/narrow truncation ranges.
function truncNormSample(lo, hi) {
  for (let i = 0; i < 1000; i++) {
    const z = randn();
    if (z >= lo && z <= hi) return z;
  }
  return Math.max(lo, Math.min(hi, 0)); // fallback, practically unreachable
}

window.rnd = { randn, randnVec, randVec, mvnSample, laplaceSample, truncNormSample };
