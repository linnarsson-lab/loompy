# Original algorithm was published by John Marioni and colleagues as EmptyDrops (Lun, A. et al. Distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data.)

# This implementation is based on the code in cellranger v3.0 by 10x Genomics

# Copyright 2018 10X Genomics, Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import logging

import numpy as np
import scipy.sparse as sparse
import scipy.stats as sp_stats

# Simple Good-Turing estimator.
# Based on S implementation in
# William A. Gale & Geoffrey Sampson (1995) Good-turing frequency estimation without tears,
# Journal of Quantitative Linguistics, 2:3, 217-237, DOI: 10.1080/09296179508590051


class SimpleGoodTuringError(Exception):
	pass


def _averaging_transform(r, nr):
	d = np.concatenate((np.ones(1, dtype=int), np.diff(r)))
	dr = np.concatenate((
		0.5 * (d[1:] + d[0:-1]),
		np.array((d[-1],), dtype=float),
	))
	return nr.astype(float) / dr


def _rstest(r, coef):
	return r * np.power(1 + 1 / r, 1 + coef)


def simple_good_turing(xr, xnr):
	"""Make a Simple Good-Turing estimate of the frequencies.

	Args:
	xr (np.array(int)): Non-zero item frequencies
	xnr (np.array(int)): Non-zero frequencies of frequencies
	Returns:
	(rstar (np.array(float)), p0 (float)):
		rstar: The adjusted non-zero frequencies
		p0: The total probability of unobserved items
	"""

	xr = xr.astype(float)
	xnr = xnr.astype(float)

	xN = np.sum(xr * xnr)

	# Get Linear Good-Turing estimate
	xnrz = _averaging_transform(xr, xnr)
	slope, intercept, _, _, _ = sp_stats.linregress(np.log(xr), np.log(xnrz))

	if slope > -1:
		raise SimpleGoodTuringError("The log-log slope is > -1 (%d); the SGT estimator is not applicable to these data." % slope)

	xrst = _rstest(xr, slope)
	xrstrel = xrst / xr

	# Get traditional Good-Turing estimate
	xrtry = xr == np.concatenate((xr[1:] - 1, np.zeros(1)))
	xrstarel = np.zeros(len(xr))
	xrstarel[xrtry] = (xr[xrtry] + 1) / xr[xrtry] * np.concatenate((xnr[1:], np.zeros(1)))[xrtry] / xnr[xrtry]

	# Determine when to switch from GT to LGT estimates
	tursd = np.ones(len(xr))
	for i in range(len(xr)):
		if xrtry[i]:
			tursd[i] = float(i + 2) / xnr[i] * np.sqrt(xnr[i + 1] * (1 + xnr[i + 1] / xnr[i]))

	xrstcmbrel = np.zeros(len(xr))
	useturing = True
	for r in range(len(xr)):
		if not useturing:
			xrstcmbrel[r] = xrstrel[r]
		else:
			if np.abs(xrstrel[r] - xrstarel[r]) * (1 + r) / tursd[r] > 1.65:
				xrstcmbrel[r] = xrstarel[r]
			else:
				useturing = False
				xrstcmbrel[r] = xrstrel[r]

	# Renormalize the probabilities for observed objects
	sumpraw = np.sum(xrstcmbrel * xr * xnr / xN)

	xrstcmbrel = xrstcmbrel * (1 - xnr[0] / xN) / sumpraw
	p0 = xnr[0] / xN

	return (xr * xrstcmbrel, p0)


def sgt_proportions(frequencies):
	"""Use Simple Good-Turing estimate to adjust for unobserved items

	Args:
	frequencies (np.array(int)): Nonzero frequencies of items
	Returns:
	(pstar (np.array(float)), p0 (float)):
		pstar: The adjusted non-zero proportions
		p0: The total probability of unobserved items
	"""
	if len(frequencies) == 0:
		raise ValueError("Input frequency vector is empty")
	if np.count_nonzero(frequencies) != len(frequencies):
		raise ValueError("Frequencies must be greater than zero")

	freqfreqs = np.bincount(frequencies.astype(np.int64))
	assert freqfreqs[0] == 0
	use_freqs = np.flatnonzero(freqfreqs)

	if len(use_freqs) < 10:
		raise SimpleGoodTuringError("Too few non-zero frequency items (%d). Aborting SGT." % len(use_freqs))

	rstar, p0 = simple_good_turing(use_freqs, freqfreqs[use_freqs])

	# rstar contains the smoothed frequencies.
	# Map each original frequency r to its smoothed rstar.
	rstar_dict = dict(zip(use_freqs, rstar))

	rstar_sum = np.sum(freqfreqs[use_freqs] * rstar)
	rstar_i = np.fromiter((rstar_dict[f] for f in frequencies), dtype=float, count=len(frequencies))
	pstar = (1 - p0) * (rstar_i / rstar_sum)

	assert np.isclose(p0 + np.sum(pstar), 1)
	return (pstar, p0)


def adjust_pvalue_bh(p):
	""" Multiple testing correction of p-values using the Benjamini-Hochberg procedure """
	descending = np.argsort(p)[::-1]
	# q = p * N / k where p = p-value, N = # tests, k = p-value rank
	scale = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(scale * p[descending]))

	# Return to original order
	return q[np.argsort(descending)]


def eval_multinomial_loglikelihoods(matrix, profile_p, max_mem_gb=0.1):
	"""Compute the multinomial log PMF for many barcodes
	Args:
	matrix (scipy.sparse.csc_matrix): Matrix of UMI counts (feature x barcode)
	profile_p (np.ndarray(float)): Multinomial probability vector
	max_mem_gb (float): Try to bound memory usage.
	Returns:
	log_likelihoods (np.ndarray(float)): Log-likelihood for each barcode
	"""
	gb_per_bc = float(matrix.shape[0] * matrix.dtype.itemsize) / (1024**3)
	bcs_per_chunk = max(1, int(round(max_mem_gb / gb_per_bc)))
	num_bcs = matrix.shape[1]

	loglk = np.zeros(num_bcs)

	for chunk_start in range(0, num_bcs, bcs_per_chunk):
		chunk = slice(chunk_start, chunk_start + bcs_per_chunk)
		matrix_chunk = matrix[:, chunk].transpose().toarray()
		n = matrix_chunk.sum(1)
		loglk[chunk] = sp_stats.multinomial.logpmf(matrix_chunk, n, p=profile_p)
	return loglk


def simulate_multinomial_loglikelihoods(profile_p, umis_per_bc, num_sims=1000, jump=1000, n_sample_feature_block=1000000, verbose=False):
	"""Simulate draws from a multinomial distribution for various values of N.
	Uses the approximation from Lun et al. ( https://www.biorxiv.org/content/biorxiv/early/2018/04/04/234872.full.pdf )
	Args:
	profile_p (np.ndarray(float)): Probability of observing each feature.
	umis_per_bc (np.ndarray(int)): UMI counts per barcode (multinomial N).
	num_sims (int): Number of simulations per distinct N value.
	jump (int): Vectorize the sampling if the gap between two distinct Ns exceeds this.
	n_sample_feature_block (int): Vectorize this many feature samplings at a time.
	Returns:
	(distinct_ns (np.ndarray(int)), log_likelihoods (np.ndarray(float)):
	distinct_ns is an array containing the distinct N values that were simulated.
	log_likelihoods is a len(distinct_ns) x num_sims matrix containing the
		simulated log likelihoods.
	"""
	distinct_n = np.flatnonzero(np.bincount(umis_per_bc.astype(np.int64)))

	loglk = np.zeros((len(distinct_n), num_sims), dtype=float)

	sampled_features = np.random.choice(len(profile_p), size=n_sample_feature_block, p=profile_p, replace=True)
	k = 0

	log_profile_p = np.log(profile_p)

	for sim_idx in range(num_sims):
		curr_counts = np.ravel(sp_stats.multinomial.rvs(distinct_n[0], profile_p, size=1))

		curr_loglk = sp_stats.multinomial.logpmf(curr_counts, distinct_n[0], p=profile_p)

		loglk[0, sim_idx] = curr_loglk

		for i in range(1, len(distinct_n)):
			step = distinct_n[i] - distinct_n[i - 1]
			if step >= jump:
				# Instead of iterating for each n, sample the intermediate ns all at once
				curr_counts += np.ravel(sp_stats.multinomial.rvs(step, profile_p, size=1))
				curr_loglk = sp_stats.multinomial.logpmf(curr_counts, distinct_n[i], p=profile_p)
				assert not np.isnan(curr_loglk)
			else:
				# Iteratively sample between the two distinct values of n
				for n in range(distinct_n[i - 1] + 1, distinct_n[i] + 1):
					j = sampled_features[k]
					k += 1
					if k >= n_sample_feature_block:
						# Amortize this operation
						sampled_features = np.random.choice(len(profile_p), size=n_sample_feature_block, p=profile_p, replace=True)
						k = 0
					curr_counts[j] += 1
					curr_loglk += log_profile_p[j] + np.log(float(n) / curr_counts[j])

			loglk[i, sim_idx] = curr_loglk

	return distinct_n, loglk


def compute_ambient_pvalues(umis_per_bc, obs_loglk, sim_n, sim_loglk):
	"""Compute p-values for observed multinomial log-likelihoods
	Args:
	umis_per_bc (nd.array(int)): UMI counts per barcode
	obs_loglk (nd.array(float)): Observed log-likelihoods of each barcode deriving from an ambient profile
	sim_n (nd.array(int)): Multinomial N for simulated log-likelihoods
	sim_loglk (nd.array(float)): Simulated log-likelihoods of shape (len(sim_n), num_simulations)
	Returns:
	pvalues (nd.array(float)): p-values
	"""
	assert len(umis_per_bc) == len(obs_loglk)
	assert sim_loglk.shape[0] == len(sim_n)

	# Find the index of the simulated N for each barcode
	sim_n_idx = np.searchsorted(sim_n, umis_per_bc)
	num_sims = sim_loglk.shape[1]

	num_barcodes = len(umis_per_bc)

	pvalues = np.zeros(num_barcodes)

	for i in range(num_barcodes):
		num_lower_loglk = np.sum(sim_loglk[sim_n_idx[i], :] < obs_loglk[i])
		pvalues[i] = float(1 + num_lower_loglk) / (1 + num_sims)
	return pvalues


def estimate_profile_sgt(matrix, barcode_indices, nz_feat):
	""" Estimate a gene expression profile by Simple Good Turing.
	Args:
	raw_mat (sparse matrix): Sparse matrix of all counts
	barcode_indices (np.array(int)): Barcode indices to use
	nz_feat (np.array(int)): Indices of features that are non-zero at least once
	Returns:
	profile (np.array(float)): Estimated probabilities of length len(nz_feat).
	"""
	# Initial profile estimate
	prof_mat = matrix[:, barcode_indices]

	profile = np.ravel(prof_mat[nz_feat, :].sum(axis=1))
	zero_feat = np.flatnonzero(profile == 0)

	# Simple Good Turing estimate
	p_smoothed, p0 = sgt_proportions(profile[np.flatnonzero(profile)])

	# Distribute p0 equally among the zero elements.
	p0_i = p0 / len(zero_feat)

	profile_p = np.repeat(p0_i, len(nz_feat))
	profile_p[np.flatnonzero(profile)] = p_smoothed

	assert np.isclose(profile_p.sum(), 1.0)
	return profile_p


# Construct a background expression profile from barcodes with <= T UMIs
def est_background_profile_sgt(matrix, use_bcs):
	""" Estimate a gene expression profile on a given subset of barcodes.
		Use Good-Turing to smooth the estimated profile.
	Args:
	matrix (scipy.sparse.csc_matrix): Sparse matrix of all counts
	use_bcs (np.array(int)): Indices of barcodes to use (col indices into matrix)
	Returns:
	profile (use_features, np.array(float)): Estimated probabilities of length use_features.
	"""
	# Use features that are nonzero anywhere in the data
	use_feats = np.flatnonzero(np.asarray(matrix.sum(1)))

	# Estimate background profile
	bg_profile_p = estimate_profile_sgt(matrix, use_bcs, use_feats)

	return (use_feats, bg_profile_p)


# Sten Linnarsson's version (Aug 2019)
def call_cells(matrix: sparse.csr_matrix, expected_n_cells: int = 5000) -> np.ndarray:
	"""
	Determine likely true cells among the barcodes by contrasting with the ambient RNA profile

	Args:
		matrix: 			expression matrix
		expected_n_cells:	expected number of true cells in the sample

	Returns:
		calls:	vector of bools indicating true cell barcodes
	"""
	n_barcodes = matrix.shape[1]
	expected_n_cells = min(expected_n_cells, n_barcodes // 5)
	total_umis = np.array(matrix.sum(axis=0))[0]  # total UMIs per barcode
	# upper limit of UMIs for barcodes considered ambient, calculated as greatest UMI count after removing twice the expected number of cells
	max_ambient_umis = np.percentile(total_umis, 100 * (n_barcodes - expected_n_cells * 2) / n_barcodes)
	# median number of UMIs among the top expected_n_cells barcodes
	median_initial_umis = np.median(total_umis[total_umis > np.percentile(total_umis, 100 * (n_barcodes - expected_n_cells) / n_barcodes)])
	min_cell_umis = int(max(500, median_initial_umis * 0.1))  # 10% of median, but at least 500 UMIs

	# Ambient RNA beads, covering the range 20 to max_amient_umis
	ambient_bcs = (total_umis < max_ambient_umis) & (total_umis > 20)
	if ambient_bcs.sum() == 0:
		# No beads were ambient, because cells had very low UMIs
		logging.warning("No ambient RNA beads were found; maybe sample had too few cells?")
		return max_ambient_umis, np.ones_like(total_umis)
	try:
		eval_features, ambient_profile_p = est_background_profile_sgt(matrix, ambient_bcs)
	except SimpleGoodTuringError as e:
		logging.error(e)
		return max_ambient_umis, np.ones_like(total_umis)

	# Evaluate candidate barcodes
	eval_bcs = total_umis > min_cell_umis
	eval_mat = matrix[eval_features, :][:, eval_bcs]

	# Compute observed log-likelihood of barcodes being generated from ambient RNA
	obs_loglk = eval_multinomial_loglikelihoods(eval_mat, ambient_profile_p)

	# Simulate log likelihoods
	distinct_ns, sim_loglk = simulate_multinomial_loglikelihoods(ambient_profile_p, total_umis[eval_bcs], num_sims=1000, verbose=True)

	# Compute p-values
	pvalues = compute_ambient_pvalues(total_umis[eval_bcs], obs_loglk, distinct_ns, sim_loglk)
	pvalues_adj = adjust_pvalue_bh(pvalues)
	pvalues_adj_all = np.ones_like(total_umis)
	pvalues_adj_all[eval_bcs] = pvalues_adj
	return max_ambient_umis, pvalues_adj_all
