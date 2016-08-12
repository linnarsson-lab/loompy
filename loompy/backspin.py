#!/usr/bin/env python

# Copyright (c) 2015, Amit Zeisel, Gioele La Manno and Sten Linnarsson
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# This .py file can be used as a library or a command-line version of BackSPIN, 
# This version of BackSPIN was implemented by Gioele La Manno.
# The BackSPIN biclustering algorithm was developed by Amit Zeisel and is described
# in Zeisel et al. Cell types in the mouse cortex and hippocampus revealed by 
# single-cell RNA-seq Science 2015 (PMID: 25700174, doi: 10.1126/science.aaa1934). 
#
# Building using pyinstaller:
# pyinstaller -F backSPIN.py -n backspin-mac-64-bit
#

from __future__ import division
from numpy import *
import sys
import os
import logging

class Results:
		pass

def calc_loccenter(x, lin_log_flag):
	M,N = x.shape
	if N==1 and M>1:
		x = x.T
	M,N = x.shape
	loc_center = zeros(M)
	min_x = x.min(1)
	x = x - min_x[:,newaxis]
	for i in range(M):
		ind = where(x[i,:]>0)[0]
		if len(ind) != 0:
			if lin_log_flag == 1:
				w = x[i,ind]/sum(x[i,ind], 0)
			else:
				w = (2**x[i,ind])/sum(2**x[i,ind], 0)
			loc_center[i] = sum(w*ind, 0)       
		else:
			loc_center[i] = 0

	return loc_center

def _calc_weights_matrix(mat_size, wid):
	'''Calculate Weight Matrix
	Parameters
	----------
	mat_size: int
		dimension of the distance matrix
	wid: int
		parameter that controls the width of the neighbourood
	Returns
	-------
	weights_mat: 2-D array
		the weights matrix to multiply with the distance matrix

	'''
	#calculate square distance from the diagonal
	sqd = (arange(1,mat_size+1)[newaxis,:] - arange(1,mat_size+1)[:,newaxis])**2
	#make the distance relative to the mat_size
	norm_sqd = sqd/wid
	#evaluate a normal pdf
	weights_mat = exp(-norm_sqd/mat_size)
	#avoid useless precision that would slow down the matrix multiplication
	weights_mat -= 1e-6
	weights_mat[weights_mat<0] = 0
	#normalize row and column sum
	weights_mat /= sum(weights_mat,0)[newaxis,:]
	weights_mat /= sum(weights_mat,1)[:, newaxis]
	#fix asimmetries
	weights_mat = (weights_mat + weights_mat.T) / 2.
	return weights_mat


def _sort_neighbourhood( dist_matrix, wid ):
	'''Perform a single iteration of SPIN
	Parameters
	----------
	dist_matrix: 2-D array
		distance matrix
	wid: int
		parameter that controls the width of the neighbourood
	Returns
	-------
	sorted_ind: 1-D array
		indexes that order the matrix

	'''
	assert wid > 0, 'Parameter wid < 0 is not allowed'
	mat_size = dist_matrix.shape[0]
	#assert mat_size>2, 'Matrix is too small to be sorted'
	weights_mat = _calc_weights_matrix(mat_size, wid)
	#Calculate the dot product (can be very slow for big mat_size)
	mismatch_score = dot(dist_matrix, weights_mat)
	energy, target_permutation = mismatch_score.min(1), mismatch_score.argmin(1)
	max_energy = max(energy)
	#Avoid points that have the same target_permutation value
	sort_score = target_permutation - 0.1 * sign( (mat_size/2 - target_permutation) ) * energy/max_energy
	#sort_score = target_permutation - 0.1 * sign( 1-2*(int(1000*energy/max_energy) % 2) ) * energy/max_energy # Alternative
	# Sorting the matrix
	sorted_ind = sort_score.argsort(0)[::-1]
	return sorted_ind


def sort_mat_by_neighborhood(dist_matrix, wid, times):
	'''Perform several iterations of SPIN using a fixed wid parameter
	Parameters
	----------
	dist_matrix: 2-D array
		distance matrix
	wid: int
		parameter that controls the width of the neighbourood
	times: int
		number of repetitions
	Returns
	-------
	indexes: 1-D array
		indexes that order the matrix

	'''
	# original indexes
	indexes = arange(dist_matrix.shape[0])
	for i in range(times):
		#sort the sitance matrix according the previous iteration
		tmpmat = dist_matrix[indexes,:] 
		tmpmat = tmpmat[:,indexes]
		sorted_ind = _sort_neighbourhood(tmpmat, wid);
		#resort the original indexes
		indexes = indexes[sorted_ind]
	return indexes


def _generate_widlist(data, axis=1, step=0.6):
	'''Generate a list of wid parameters to execute sort_mat_by_neighborhood
	Parameters
	----------
	data: 2-D array
		the data matrix
	axis: int
		the axis to take in consideration
	step: float
		the increment between two successive wid parameters
	Returns
	-------
	wid_list: list of int
		list of wid parameters to run SPIN

	'''
	max_wid = data.shape[axis]*0.6
	new_wid = 1
	wid_list = []
	while new_wid < (1+step)*max_wid:
		wid_list.append( new_wid )
		new_wid = int(ceil( new_wid + new_wid*(step) +1))
	return wid_list[::-1]



def SPIN(dt, widlist=[10,1], iters=30, axis='both'):
	"""Run the original SPIN algorithm
	Parameters
	----------
	dt: 2-D array
		the data matrix
	widlist: float or list of int
		If float is passed, it is used as step parameted of _generate_widlist, 
		and widlist is generated to run SPIN.
		If list is passed it is used directly to run SPIN.
	iters: int
		number of repetitions for every wid in widlist
	axis: int
		the axis to take in consideration (must be 0, 1 or 'both')
	step: float
		the increment between two successive wid parameters
	Returns
	-------
	indexes: 1-D array (if axis in [0,1]) or tuple of 1-D array (if axis = 'both')
		indexes that sort the data matrix
	Notes
	-----
	Typical usage
	sorted_dt0 = SPIN(dt, iters=30, axis=0)
	sorted_dt1 = SPIN(dt, iters=30, axis=1)
	dt = dt[sorted_dt0,:]
	dt = dt[:,sorted_dt1]
	"""
	IXc = arange(dt.shape[1])
	IXr = arange(dt.shape[0])
	assert axis in ['both', 0,1], 'axis must be 0, 1 or \'both\' '
	#Sort both axis
	if axis == 'both':
		CCc = 1 - corrcoef(dt.T) 
		CCr = 1 - corrcoef(dt) 
		if type(widlist) != list:
			widlist_r = _generate_widlist(dt, axis=0, step=widlist)
			widlist_c = _generate_widlist(dt, axis=1, step=widlist)
		logging.debug("Sorting genes")
		for wid in widlist_r:
			logging.debug("Neighbourhood = " + str(wid))
			INDr = sort_mat_by_neighborhood(CCr, wid, iters)
			CCr = CCr[INDr,:][:,INDr]
			IXr = IXr[INDr]
		logging.debug("Sorting cells")
		for wid in widlist_c:
			logging.debug("Neighbourhood = " + str(wid))
			INDc = sort_mat_by_neighborhood(CCc, wid, iters)
			CCc = CCc[:,INDc][INDc,:]
			IXc= IXc[INDc]
		return IXr, IXc
	#Sort rows
	elif axis == 0:
		CCr = 1 - corrcoef(dt)
		if type(widlist) != list:
			widlist = _generate_widlist(dt, axis=0, step=widlist)
		logging.debug("Sorting genes")
		for wid in widlist:
			logging.debug("Neighbourhood = " + str(wid))
			INDr = sort_mat_by_neighborhood(CCr, wid, iters)
			CCr = CCr[INDr,:][:,INDr]
			IXr = IXr[INDr]
		return IXr
	#Sort columns
	elif axis == 1:
		CCc = 1 - corrcoef(dt.T)
		if type(widlist) != list:
			widlist = _generate_widlist(dt, axis=1, step=widlist)
		logging.debug("Sorting cells")
		for wid in widlist:
			logging.debug("Neighbourhood = " + str(wid))
			INDc = sort_mat_by_neighborhood(CCc, wid, iters)
			CCc = CCc[:,INDc][INDc,:]
			IXc = IXc[INDc]
		return IXc


def backSPIN(data, numLevels=2, first_run_iters=10, first_run_step=0.05, runs_iters=8 ,runs_step=0.25,\
	split_limit_g=2, split_limit_c=2, stop_const = 1.15, low_thrs=0.2):
	'''Run the backSPIN algorithm
	Parameters
	----------
	data: 2-D array
		the data matrix, rows should be genes and columns single cells/samples
	numLevels: int
		the number of splits that will be tried
	first_run_iters: float
		the iterations of the preparatory SPIN
	first_run_step: float
		the step parameter passed to _generate_widlist for the preparatory SPIN
	runs_iters: int
		the iterations parameter passed to the _divide_to_2and_resort.
		influences all the SPIN iterations except the first
	runs_step: float
		the step parameter passed to the _divide_to_2and_resort.
		influences all the SPIN iterations except the first
	wid: float
		the wid of every iteration of the splitting and resorting
	split_limit_g: int
		If the number of specific genes in a subgroup is smaller than this number
		 splitting of that subgrup is not allowed
	split_limit_c: int
		If the number cells in a subgroup is smaller than this number splitting of
		that subgrup is not allowed
	stop_const: float
		minimum score that a breaking point has to reach to be suitable for splitting
	low_thrs: float
		genes with average lower than this threshold are assigned to either of the 
		splitting group reling on genes that are higly correlated with them

	Returns
	-------
	results: Result object
		The results object contain the following attributes
		genes_order: 1-D array
			indexes (a permutation) sorting the genes 
		cells_order: 1-D array
			indexes (a permutation) sorting the cells 
		genes_gr_level: 2-D array
			for each depth level contains the cluster indexes for each gene
		cells_gr_level:
			for each depth level contains the cluster indexes for each cell
		cells_gr_level_sc:
			score of the splitting
		genes_bor_level:
			the border index between gene clusters
		cells_bor_level:
			the border index between cell clusters

	Notes
	-----
	Typical usage
	
	'''
	assert numLevels>0, '0 is not an available depth for backSPIN, use SPIN instead'
	#initialize some varaibles
	genes_bor_level = [[] for i in range(numLevels)] 
	cells_bor_level = [[] for i in range(numLevels)] 
	N,M = data.shape
	genes_order = arange(N)
	cells_order = arange(M)
	genes_gr_level = zeros((N,numLevels+1))
	cells_gr_level = zeros((M,numLevels+1))
	cells_gr_level_sc = zeros((M,numLevels+1))

	# Do a Preparatory SPIN on cells
	logging.debug("Preparatory SPIN")
	ix1 = SPIN(data, widlist=_generate_widlist(data, axis=1, step=first_run_step), iters=first_run_iters, axis=1)
	cells_order = cells_order[ix1]

	#For every level of depth DO:
	for i in range(numLevels): 
		k=0 # initialize group id counter
		# For every group generated at the parent level DO:
		for j in range( len( set(cells_gr_level[:,i]) ) ): 
			# Extract the a data matrix of the genes at that level
			g_settmp = nonzero(genes_gr_level[:,i]==j)[0] #indexes of genes in the level j
			c_settmp = nonzero(cells_gr_level[:,i]==j)[0] #indexes of cells in the level j
			datatmp = data[ ix_(genes_order[g_settmp], cells_order[c_settmp]) ]
			# If we are not below the splitting limit for both genes and cells DO:
			if (len(g_settmp)>split_limit_g) & (len(c_settmp)>split_limit_c): 
				# Split and SPINsort the two halves
				if i == numLevels-1:
					divided = _divide_to_2and_resort(datatmp, wid=runs_step, iters_spin=runs_iters,\
						stop_const=stop_const, low_thrs=low_thrs, sort_genes=True)
				else:
					divided = _divide_to_2and_resort(datatmp, wid=runs_step, iters_spin=runs_iters,\
						stop_const=stop_const, low_thrs=low_thrs, sort_genes=False)
				# _divide_to_2and_resort retruns an empty array in gr2 if the splitting condition was not satisfied
				if divided:
					sorted_data_resort1, genes_resort1, cells_resort1,\
					gr1, gr2, genesgr1, genesgr2, score1, score2 = divided
					# Resort from the previous level
					genes_order[g_settmp] = genes_order[g_settmp[genes_resort1]]
					cells_order[c_settmp] = cells_order[c_settmp[cells_resort1]]
					# Assign a numerical identifier to the groups
					genes_gr_level[g_settmp[genesgr1],i+1] = k
					genes_gr_level[g_settmp[genesgr2],i+1] = k+1
					cells_gr_level[c_settmp[gr1],i+1] = k
					cells_gr_level[c_settmp[gr2],i+1] = k+1
					# Not really clear what sc is
					cells_gr_level_sc[c_settmp[gr1],i+1] = score1
					cells_gr_level_sc[c_settmp[gr2],i+1] = score2
					# Augment the counter of 2 becouse two groups were generated from one
					k = k+2
				else:
					# The split is not convenient, keep everithing the same
					genes_gr_level[g_settmp,i+1] = k
					# if it is the deepest level: perform gene sorting
					if i == numLevels-1:
						if (datatmp.shape[0] > 2 )and (datatmp.shape[1] > 2):
							genes_resort1 = SPIN(datatmp, widlist=runs_step, iters=runs_iters, axis=0)
							genes_order[g_settmp] = genes_order[g_settmp[genes_resort1]]
					cells_gr_level[c_settmp,i+1] = k
					cells_gr_level_sc[c_settmp,i+1] = cells_gr_level_sc[c_settmp,i]
					# Augment of 1 becouse no new group was generated
					k = k+1
			else:
				# Below the splitting limit: the split is not convenient, keep everithing the same
				genes_gr_level[g_settmp,i+1] = k
				cells_gr_level[c_settmp,i+1] = k
				cells_gr_level_sc[c_settmp,i+1] = cells_gr_level_sc[c_settmp,i]
				# Augment of 1 becouse no new group was generated
				k = k+1
		
		# Find boundaries
		genes_bor_level[i] = r_[0, nonzero(diff(genes_gr_level[:,i+1])>0)[0]+1, data.shape[0] ]
		cells_bor_level[i] = r_[0, nonzero(diff(cells_gr_level[:,i+1])>0)[0]+1, data.shape[1] ]

	#dataout_sorted = data[ ix_(genes_order,cells_order) ]

	results = Results()
	results.genes_order = genes_order
	results.cells_order = cells_order
	results.genes_gr_level = genes_gr_level
	results.cells_gr_level = cells_gr_level
	results.cells_gr_level_sc = cells_gr_level_sc
	results.genes_bor_level = genes_bor_level
	results.cells_bor_level = cells_bor_level

	return results
	
	

def _divide_to_2and_resort(sorted_data, wid, iters_spin=8, stop_const = 1.15, low_thrs=0.2 , sort_genes=True):
	'''Core function of backSPIN: split the datamatrix in two and resort the two halves

	Parameters
	----------
	sorted_data: 2-D array
		the data matrix, rows should be genes and columns single cells/samples
	wid: float
		wid parameter to give to widlist parameter of th SPIN fucntion
	stop_const: float
		minimum score that a breaking point has to reach to be suitable for splitting
	low_thrs: float
		if the difference between the average expression of two groups is lower than threshold the algorythm 
		uses higly correlated gens to assign the gene to one of the two groups

	Returns
	-------
	'''
	
	# Calculate correlation matrix for cells and genes
	Rcells = corrcoef(sorted_data.T)
	Rgenes = corrcoef(sorted_data)
	# Look for the optimal breaking point
	N = Rcells.shape[0]
	score = zeros(N)
	for i in range(2,N-2):
		if i == 2:
			tmp1 = sum( Rcells[:i,:i] )
			tmp2 = sum( Rcells[i:,i:] )
			score[i] = (tmp1+tmp2) / float(i**2 + (N-i)**2)
		else:
			tmp1 += sum(Rcells[i-1,:i]) + sum(Rcells[:i-1,i-1]);
			tmp2 -= sum(Rcells[i-1:,i-1]) + sum(Rcells[i-1,i:]);
			score[i] = (tmp1+tmp2) / float(i**2 + (N-i)**2)
	
	breakp1 = argmax(score)
	score1 = Rcells[:breakp1,:breakp1]
	score1 = triu(score1)
	score1 = mean( score1[score1 != 0] )
	score2 = Rcells[breakp1:,breakp1:]
	score2 = triu(score2)
	score2 = mean( score2[score2 != 0] )
	avg_tot = triu(Rcells)
	avg_tot = mean( avg_tot[avg_tot != 0] )

	# If it is convenient to break
	if (max([score1,score2])/avg_tot) > stop_const:
		# Divide in two groups
		gr1 = arange(N)[:breakp1]
		gr2 = arange(N)[breakp1:]
		# and assign the genes into the two groups
		mean_gr1 = sorted_data[:,gr1].mean(1)
		mean_gr2 = sorted_data[:,gr2].mean(1)
		concat_loccenter_gr1 = c_[ calc_loccenter(sorted_data[:,gr1], 2), calc_loccenter(sorted_data[:,gr1][...,::-1], 2) ]
		concat_loccenter_gr2 = c_[ calc_loccenter(sorted_data[:,gr2], 2), calc_loccenter(sorted_data[:,gr2][...,::-1], 2) ]
		center_gr1, flip_flag1 = concat_loccenter_gr1.min(1), concat_loccenter_gr1.argmin(1)
		center_gr2, flip_flag2 = concat_loccenter_gr2.max(1), concat_loccenter_gr2.argmax(1)
		sorted_data_tmp = array( sorted_data )
		sorted_data_tmp[ix_(flip_flag1==1,gr1)] = sorted_data[ix_(flip_flag1==1,gr1)][...,::-1]
		sorted_data_tmp[ix_(flip_flag2==1,gr2)] = sorted_data[ix_(flip_flag2==1,gr2)][...,::-1]
		loc_center = calc_loccenter(sorted_data_tmp, 2)

		imax = zeros(loc_center.shape)
		imax[loc_center<=breakp1] = 1
		imax[loc_center>breakp1] = 2

		genesgr1 = where(imax==1)[0]
		genesgr2 = where(imax==2)[0]
		if size(genesgr1) == 0:
			IN = argmax(mean_gr1)
			genesgr1 = array([IN])
			genesgr2 = setdiff1d(genesgr2, IN)
		elif size(genesgr2) == 0:
			IN = argmax(mean_gr2)
			genesgr2 = array([IN])
			genesgr1 = setdiff1d(genesgr1, IN)
		
		logging.debug("Splitting %s in %s and %s" % (str(sorted_data.shape),str((genesgr1.shape[0],gr1.shape[0])),str((genesgr2.shape[0],gr2.shape[0]))))

		# Data of group1
		datagr1 = sorted_data[ix_(genesgr1,gr1)]
		# zero center
		datagr1 = datagr1 - datagr1.mean(1)[:,newaxis]
		# Resort group1
		if min( datagr1.shape ) > 1:
			if sort_genes:
				genesorder1,cellorder1 = SPIN(datagr1, widlist=wid, iters=iters_spin, axis='both')
			else:
				cellorder1 = SPIN(datagr1, widlist=wid, iters=iters_spin, axis=1)
				genesorder1 = arange(datagr1.shape[0])
		elif len(genesgr1) == 1:
			genesorder1 = 0
			cellorder1 = argsort( datagr1[0,:] )
		elif len(gr1) == 1:
			cellorder1 = 0
			genesorder1 = argsort( datagr1[:,0] )

		# Data of group2
		datagr2 = sorted_data[ix_(genesgr2,gr2)]
		# zero center
		datagr2 = datagr2 - datagr2.mean(1)[:,newaxis]
		# Resort group2
		if min( datagr2.shape )>1:
			if sort_genes:
				genesorder2, cellorder2 = SPIN(datagr2, widlist=wid, iters=iters_spin, axis='both')
			else:
				cellorder2 = SPIN(datagr2, widlist=wid, iters=iters_spin, axis=1)
				genesorder2 = arange(datagr2.shape[0])
		elif len(genesgr2) == 1:
			genesorder2 = 0
			cellorder2 = argsort(datagr2[0,:])
		elif len(gr2) == 1:
			cellorder2 = 0
			genesorder2 = argsort(datagr2[:,0])
		
		# contcatenate cells and genes indexes
		genes_resort1 = r_[genesgr1[genesorder1], genesgr2[genesorder2] ]
		cells_resort1 = r_[gr1[cellorder1], gr2[cellorder2] ]
		genesgr1 = arange(len(genesgr1))
		genesgr2 = arange(len(genesgr1), len(sorted_data[:,0]))
		# resort
		sorted_data_resort1 = sorted_data[ix_(genes_resort1,cells_resort1)]

		return sorted_data_resort1, genes_resort1, cells_resort1, gr1, gr2, genesgr1, genesgr2, score1, score2

	else:
		logging.debug("Low splitting score was: %.4f" % (max([score1,score2])/avg_tot))
		return False

	
def loom_backspin(ds):
	numLevels=2 # -d
	first_run_iters=10 # -t
	first_run_step=0.1 # -s
	runs_iters=8 # -T
	runs_step=0.3 # -S
	split_limit_g=2 # -g
	split_limit_c=2 # -c
	stop_const = 1.15 # -k
	low_thrs=0.2 # -r

	logging.debug("  numLevels: " + str(numLevels))
	logging.debug("  first_run_iters: " + str(first_run_iters))
	logging.debug("  first_run_step: " + str(first_run_step))
	logging.debug("  runs_iters: " + str(runs_iters))
	logging.debug("  runs_step: " + str(runs_step))
	logging.debug("  split_limit_g: " + str(split_limit_g))
	logging.debug("  split_limit_c: " + str(split_limit_c))
	logging.debug("  stop_const: " + str(stop_const))
	logging.debug("  low_thrs: " + str(low_thrs))

	if ds.row_attrs.has_key("_Excluded"):
		data = ds[(1 - ds._Excluded).astype("bool"), :]
	else:
		logging.warn("No genes have been excluded, running BackSPIN on entire matrix")
		data = ds[:,:]

	logging.info("Using %i genes and %i cells" % data.shape)
	results = backSPIN(data, numLevels, first_run_iters, first_run_step, runs_iters, runs_step,\
		split_limit_g, split_limit_c, stop_const, low_thrs)

	logging.debug("Writing result to file")

	for level, groups in enumerate( results.genes_gr_level.T ):
		temp = zeros((ds.shape[0],))
		temp[where(ds._Excluded == 0)] = array([int(el) for el in groups])
		ds.set_attr('BackSPIN_level_%i_group' % level, temp, axis=0)
	for level, groups in enumerate( results.cells_gr_level.T ):
		ds.set_attr('BackSPIN_level_%i_group' % level, array([int(el) for el in groups]), axis=1)

	logging.info("BackSPIN all done")

