#!/usr/bin/env python
#@file: combiningfiles.py
#@author: Rituparna Samanta (rsamant2@jhu.edu)
#@brief: This combines multiple energy landscape files to a single one.

import sys
import os
import glob
import pandas as pd
import numpy as np

def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return(unique_list)

def combiningallfiles(folder, peptide_list, partitions):
    start_index = [ f'{float(ii):.6f}' for ii in partitions[:-1] ]
    end_index = [ f'{float(ii):.6f}' for ii in partitions[1:] ]

    for peptide_name in peptide_list:

        labelin1 = folder + '/' + peptide_name + '/' + peptide_name + '_combined.dat'
        labelout1 = folder + '/' + peptide_name + '/' + peptide_name + '_franklin2019_landscape.dat'
        print("Processing data for", labelin1)

        all_tables = []
        for si, ei in zip(start_index, end_index):
            label= folder+'/'+peptide_name+'/'+peptide_name+'_franklin2019_'+ei+'_'+si+'_landscape.dat'
            all_tables.append( pd.read_csv(label) )

        df = pd.concat(all_tables)
        df.to_csv(labelin1, index=False)
        df_read = pd.read_csv(labelin1, delimiter=" ")
        X = unique(df_read['zcoord'])
        Y = unique(df_read['angle'])

        mindata = []

        for i in range(len(X)):

            for j in range(len(Y)):


                newarr = df_read[df_read['angle']==Y[j]]
                secondarr = newarr[newarr['zcoord']==X[i]]

                arr = np.array(secondarr['total_score'])

                if len(arr) == 0:
                    print("Skipping zcoord", X[i], "and angle", Y[j], "for", labelin1, "due to missing data")
                    continue

                minpos = np.argmin(arr)

                mindata.append(secondarr.iloc[minpos,0:23])
        mindatadf = pd.DataFrame(mindata)
        mindatadf.to_csv(labelout1, sep=" ", index=False)

        del df, df_read
        del X, Y

   # fclose(fileID);




