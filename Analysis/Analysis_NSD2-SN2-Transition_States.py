from __future__ import print_function
import mdtraj as mdt
import os, shutil
import re
import numpy as np
from numpy import array
import collections
from collections import defaultdict
import pandas as pd
import itertools
import warnings

warnings.filterwarnings("ignore", message="top= kwargs ignored since this file parser does not support it")

Protein = 'NSD2'
sim_time = '100ns'
peptide = 'H3K36'

number_replicates = 1
count = 1

while (count <= number_replicates):
 traj = mdt.load('/path-to-trajectory/production_{}_{}_{}_{}.h5'.format(Protein, peptide, sim_time, count))
 
 topology=traj.topology                  #here the topolgy of the trajectory is made
 table, bonds = topology.to_dataframe()  #here the atoms for the SN2 criteria are selected
 df = pd.DataFrame(table)
 
 df_search = df[(df['chainID'] == 1) & (df['resName'] == 'LYS') & (df['resSeq'] == 36) & (df['name'] == 'CE')]
 K36_CE = df_search.index.item()
 df_search = df[(df['chainID'] == 1) & (df['resName'] == 'LYS') & (df['resSeq'] == 36) & (df['name'] == 'NZ')]
 K36_NZ = df_search.index.item()
 df_search = df[(df['chainID'] == 2) & (df['resName'] == 'SAM') & (df['resSeq'] == 1804) & (df['name'] == 'N')]
 SAM_N = df_search.index.item()
 df_search = df[(df['chainID'] == 2) & (df['resName'] == 'SAM') & (df['resSeq'] == 1804) & (df['name'] == 'SD')]
 SAM_SD = df_search.index.item()
 df_search = df[(df['chainID'] == 2) & (df['resName'] == 'SAM') & (df['resSeq'] == 1804) & (df['name'] == 'CE')]
 SAM_CE = df_search.index.item()
 
 NSD2_dist = [[K36_NZ,SAM_CE],[K36_NZ,SAM_N]]
 NSD2_109 = [[K36_CE, K36_NZ, SAM_CE]]
 NSD2_linear = [[SAM_SD, SAM_CE, K36_NZ]]
 
 dist_NAC = mdt.compute_distances(traj,NSD2_dist, periodic=True, opt=True)
 angle_109 = mdt.compute_angles(traj, array(NSD2_109), periodic=True, opt=True)
 angle_linear = mdt.compute_angles(traj, array(NSD2_linear), periodic=True, opt=True)

 #write all distances in a list but only from first atom pair
 dist_NAC_list = np.ndarray.tolist(dist_NAC)
 list_dist_NAC_merged = list(itertools.chain.from_iterable(dist_NAC_list))
 NAC = list_dist_NAC_merged[::2]

 #tranform all dihedrals (BB, SC) from rad in degrees and write them in a list
 angle_109=np.rad2deg(angle_109)
 angle_linear=np.rad2deg(angle_linear)
 
 angle_109_list=np.ndarray.tolist(angle_109)
 angle_linear_list=np.ndarray.tolist(angle_linear)
 angle_109_list_merged = list(itertools.chain.from_iterable(angle_109_list))
 angle_linear_list_merged = list(itertools.chain.from_iterable(angle_linear_list))
 
 #criteria lists in array
 criteria_array = np.column_stack([NAC,  angle_linear_list_merged, angle_109_list_merged])
 df = pd.DataFrame(criteria_array)
 df_sliced = df[(df[0] <0.4) & (df[1] >150) & (df[1] <210) & (df[2] >79) & (df[2] <139)] #+-30°, 4.0 A distance NZ-CE, 150°-210° angle_linear, 79°-139° angle_109
 

 indices = df_sliced.index.array
 indices_df = pd.DataFrame(indices)
 temp = indices_df.values.tolist()
 INT = list(itertools.chain.from_iterable(temp))
 
 print('number of SN2 frames in {} {} in {}:'.format(Protein, peptide,count), len(INT)) #shows you the sumber of SN2 TS like conformations
 
 count = count+1