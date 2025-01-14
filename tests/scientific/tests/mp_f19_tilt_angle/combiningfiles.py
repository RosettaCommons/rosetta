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

def combiningallfiles(folder, peptide_list):

    for element_num in range(len(peptide_list)):
        start_index = ['-60.000000','-40.000000','-20.000000','0.000000','20.000000', '40.000000']
        end_index = ['-40.000000','-20.000000', '0.000000', '20.000000','40.000000', '60.000000']

        peptide_name = peptide_list[element_num]
        labelin1 = folder + '/' + peptide_name + '/' + peptide_name + '_combined.dat'
        labelout1 = folder + '/' + peptide_name + '/' + peptide_name + '_franklin2019_landscape.dat'
        print(labelin1)
        print(labelout1)


        label0= folder+'/'+peptide_name+'/'+peptide_name+'_franklin2019_'+str(end_index[0])+'_'+str(start_index[0])+'_landscape.dat'
        df0 = pd.read_csv(label0)
        label1= folder+'/'+peptide_name+'/'+peptide_name+'_franklin2019_'+str(end_index[1])+'_'+str(start_index[1])+'_landscape.dat'
        df1 = pd.read_csv(label1)
        label2= folder+'/'+peptide_name+'/'+peptide_name+'_franklin2019_'+str(end_index[2])+'_'+str(start_index[2])+'_landscape.dat'
        df2 = pd.read_csv(label2)
        label3= folder+'/'+peptide_name+'/'+peptide_name+'_franklin2019_'+str(end_index[3])+'_'+str(start_index[3])+'_landscape.dat'
        df3 = pd.read_csv(label3)
        label4= folder+'/'+peptide_name+'/'+peptide_name+'_franklin2019_'+str(end_index[4])+'_'+str(start_index[4])+'_landscape.dat'
        df4 = pd.read_csv(label4)
        label5= folder+'/'+peptide_name+'/'+peptide_name+'_franklin2019_'+str(end_index[5])+'_'+str(start_index[5])+'_landscape.dat'
        df5 = pd.read_csv(label5)



        df = pd.concat([df0,df1,df2,df3,df4,df5])
        df.to_csv(labelin1, index=False)
        df_read = pd.read_csv(labelin1, delimiter=" ")
        X = unique(df_read['zcoord'])
        Y = unique(df_read['angle'])

        mindata = []

        try:
            for i in range(len(X)):

                for j in range(len(Y)):


                    newarr = df_read[df_read['angle']==Y[j]]
                    secondarr = newarr[newarr['zcoord']==X[i]]

                    arr = np.array(secondarr['total_score'])


                    minpos = np.argmin(arr)

                    mindata.append(secondarr.iloc[minpos,0:23])
        except:
            print("ERROR When processing data in file", labelin1)
            raise

        mindatadf = pd.DataFrame(mindata)
        mindatadf.to_csv(labelout1, sep=" ", index=False)

        del df0
        del df1
        del df2
        del df3, df4, df5
        del df, df_read
        del X, Y

   # fclose(fileID);




