This integration test tests repacking a symmetric pose with strong
electron density constraints.


The ccp4 file was generated with this command

  phenix.model_vs_data <pdb_file.pdb> <structure_factors.ent> map="2mFo-DFc"

and output below. Note the final highest resolution is
1.9A. Frank Dimaio suggests using grid spacing ~45% of the
resolution. So the grid spacing should be .85A but for
integration test, use the larger grid spacing of 3A to make it faster.

This integration test uses 1xu1FH_D.pdb from the make_symmdef_file
integration test to conserve space. As a consistency check, the
generated symmetry definition file,
'../../tests/make_symmdef_file/1xu1FH_D.symm' should match
'inputs/1xu1FH_D.symm'

/////// phenix.model_vs_data output ////////////

pdb /home/momeara/scr/data/pdb/xu/pdb1xu1.ent
sf /home/momeara/scr/data/structure_factors/xu/r1xu1sf.ent
begin:  /home/momeara/scr/data/pdb/xu/pdb1xu1.ent
miller array size: 44519
miller array size: 44519
  Unit cell:        (59.339, 91.839, 102.268, 90, 90, 90)
  Space group:      P 21 21 21 (No. 19) number of symmetry operations: 4
  Number of models: 1
  Model #1:
    Number of residues in alternative conformations: 0
    Residue content:
      element    : 1
      water      : 130
      amino_acid : 524
    Atoms:
      content: 4280
        Ni count:    1 occupancy sum: 1.00
        C count: 2619 occupancy sum: 2616.00
        S count:   33 occupancy sum: 33.00
        O count:  873 occupancy sum: 873.00
        N count:  754 occupancy sum: 751.00
      ADP (min,max,mean):
        all           (4280 atoms): 12.3   88.4   30.2  
        side chains   (2053 atoms): 12.5   88.4   33.7  
        main chains   (2096 atoms): 12.3   71.4   26.5  
        macromolecule (4149 atoms): 12.3   88.4   30.1  
        ligands       (1 atoms): 18.8   18.8   18.8  
        solvent       (130 atoms): 18.6   51.9   33.5  
      mean bonded (Bi-Bj) : 3.25
      occupancies (min,max,mean)       : 0.00 1.00 1.00
      number_of_anisotropic            : 4150   
      number_of_non_positive_definite  : 0
    Stereochemistry statistics (mean, max, count) - overall:
      bonds            :   0.0087   0.1008 4254
      angles           :   1.3139   9.6412 5730
      dihedrals        :  16.4309  88.0722 1559
      chirality        :   0.0893   0.3084 628
      planarity        :   0.0044   0.0233 737
      non-bonded (min) :   2.2193
    Stereochemistry statistics (mean, max, count) - macromolecule:
      bonds            :   0.0087   0.1008 4254
      angles           :   1.3139   9.6412 5730
      dihedrals        :  16.4309  88.0722 1559
      chirality        :   0.0893   0.3084 628
      planarity        :   0.0044   0.0233 737
      non-bonded (min) :   2.3433
    Stereochemistry statistics (mean, max, count) - ligands:
      bonds            :   0.0000   0.0000 0
      angles           :   0.0000   0.0000 0
      dihedrals        :   0.0000   0.0000 0
      chirality        :   0.0000   0.0000 0
      planarity        :   0.0000   0.0000 0
      non-bonded (min) :   0.0000
    Stereochemistry statistics - solvent:
      non-bonded (min) :   2.4925
    Molprobity statistics:
      Ramachandran plot, number of:
        outliers : 0     (0.00  %)
        allowed  : 10    (1.95  %)
        favored  : 502   (98.05 %)
      Rotamer outliers        : 10 (2.18 %)
      Cbeta deviations >0.25A : 0
      All-atom clashscore     : 4.75 (steric overlaps >0.4A per 1000 atoms)
  Data:
    data_label                           : r1xu1sf,_refln.F_meas_au,_refln.F_meas_sigma_au 
    high_resolution                      : 1.90  
    low_resolution                       : 35.69  
    completeness_in_range                : 1.00   
    completeness(d_min-inf)              : 1.00   
    completeness(6A-inf)                 : 0.91   
    wilson_b                             : 24.2   
    number_of_reflections                : 44515    
    number_of_reflections(non-anomalous) : 44515    
    test_set_size                        : 0.0999   
    test_flag_value                      : 0 
    number_of_Fobs_outliers              : 4        
    twinned                              : False 
    anomalous_flag                       : False 
  Model_vs_Data:
    r_work(re-computed)                : 0.1803 
    r_free(re-computed)                : 0.2102 
    bulk_solvent_(k_sol,b_sol)         : 0.40 60.00   
    overall_anisotropic_scale_(b_cart) : 1.20 -1.24 0.04 -0.00 0.00 0.00 
    solvent_content_estimated_via_mask : 43.9 %
    mFo-DFc map: positive and negative peak numbers:
      >  3 sigma:  507
      >  6 sigma:  7
      >  9 sigma:  1
      < -3 sigma:  264
      < -6 sigma:  5
      < -9 sigma:  0
  Information extracted from PDB file header:
    program_name    : REFMAC
    year            : 4
    r_work          : 0.167
    r_free          : 0.203
    high_resolution : 1.9
    low_resolution  : 30.0
    sigma_cutoff    : 0.0
    matthews_coeff  : 2.1
    solvent_cont    : 40.0 %
    TLS             : True (number of groups: 6)
  After applying resolution and sigma cutoffs:
    n_refl_cutoff : 44514
    r_work_cutoff : 0.1803
    r_free_cutoff : 0.2102
  Statistics in resolution bins:
    SIGMAA vs Resolution:
     resolution(A)  sigmaa
            35.811   0.941
             8.173   0.941
             5.698   0.944
             4.544   0.957
             3.855   0.958
             3.391   0.945
             3.057   0.946
             2.805   0.941
             2.609   0.940
             2.453   0.942
             2.328   0.943
             2.225   0.947
             2.142   0.938
             2.074   0.954
             2.020   0.954
             1.976   0.951
             1.943   0.951
             1.918   0.945
             1.903   0.945
             1.895   0.947
 Bin  Resolution Range  Compl.    Nwork  Rwork   Scale     <Fobs>
   1 35.6974 -  8.1589    0.72      465 0.2277  0.9965    546.198
   2  8.1589 -  6.4886    0.92      548 0.2185  0.9684    492.555
   3  6.4886 -  5.6721    0.87      517 0.2263  0.9307    498.759
   4  5.6721 -  5.1552    0.89      515 0.2008  0.9648    562.053
   5  5.1552 -  4.7866    0.91      528 0.1709  0.9888    676.739
   6  4.7866 -  4.5050    0.88      517 0.1657  1.0169    682.474
   7  4.5050 -  4.2798    0.90      501 0.1461  1.0277    705.744
   8  4.2798 -  4.0937    0.89      517 0.1667  1.0406    653.365
   9  4.0937 -  3.9363    0.87      479 0.1683  1.0010    588.191
  10  3.9363 -  3.8007    0.87      513 0.1754  1.0240    582.195
  11  3.8007 -  3.6820    0.98      552 0.1825  1.0216    526.675
  12  3.6820 -  3.5768    0.85      479 0.1628  1.0171    512.757
  13  3.5768 -  3.4827    0.84      469 0.1691  1.0123    503.798
  14  3.4827 -  3.3978    1.00      567 0.1691  1.0210    449.483
  15  3.3978 -  3.3206    0.84      474 0.1735  1.0099    435.406
  16  3.3206 -  3.2500    1.00      575 0.1881  1.0056    393.467
  17  3.2500 -  3.1850    0.79      432 0.1769  1.0270    382.777
  18  3.1850 -  3.1250    1.00      542 0.1817  0.9908    370.553
  19  3.1250 -  3.0692    0.75      443 0.1820  0.9813    340.370
  20  3.0692 -  3.0172    1.00      559 0.1857  0.9967    321.643
  21  3.0172 -  2.9686    0.79      435 0.1825  0.9912    313.964
  22  2.9686 -  2.9229    1.00      578 0.2116  0.9860    282.766
  23  2.9229 -  2.8799    0.74      391 0.1859  0.9835    271.032
  24  2.8799 -  2.8394    1.00      564 0.2154  0.9827    264.009
  25  2.8394 -  2.8010    0.84      483 0.2068  0.9817    248.052
  26  2.8010 -  2.7647    0.92      517 0.2036  0.9782    235.404
  27  2.7647 -  2.7301    1.00      548 0.2050  0.9663    229.476
  28  2.7301 -  2.6972    0.72      390 0.1940  0.9805    233.086
  29  2.6972 -  2.6659    1.00      573 0.2075  0.9693    224.233
  30  2.6659 -  2.6359    0.99      548 0.1852  0.9670    218.167
  31  2.6359 -  2.6073    0.72      411 0.1845  0.9683    222.158
  32  2.6073 -  2.5799    1.00      549 0.1893  0.9599    203.006
  33  2.5799 -  2.5535    0.94      534 0.1908  0.9568    202.213
  34  2.5535 -  2.5283    0.76      430 0.1954  0.9629    201.080
  35  2.5283 -  2.5040    1.00      542 0.1978  0.9727    191.269
  36  2.5040 -  2.4806    0.97      529 0.1927  0.9549    186.399
  37  2.4806 -  2.4580    0.69      384 0.1930  0.9589    182.424
  38  2.4580 -  2.4363    1.00      570 0.1904  0.9700    178.681
  39  2.4363 -  2.4153    1.00      523 0.1623  1.0101    190.959
  40  2.4153 -  2.3950    0.65      383 0.1598  0.9729    182.452
  41  2.3950 -  2.3754    1.00      529 0.1758  0.9877    181.926
  42  2.3754 -  2.3564    1.00      574 0.1902  0.9921    176.247
  43  2.3564 -  2.3380    0.78      418 0.1837  0.9977    175.455
  44  2.3380 -  2.3201    0.86      487 0.1773  1.0039    175.217
  45  2.3201 -  2.3028    1.00      546 0.1745  1.0064    167.929
  46  2.3028 -  2.2860    1.00      549 0.1710  1.0125    167.244
  47  2.2860 -  2.2697    0.65      349 0.1742  1.0011    159.205
  48  2.2697 -  2.2538    1.00      592 0.1636  1.0208    171.048
  49  2.2538 -  2.2384    1.00      540 0.1722  1.0294    164.322
  50  2.2384 -  2.2234    0.99      545 0.1819  1.0308    158.630
  51  2.2234 -  2.2087    0.60      328 0.1605  1.0151    160.867
  52  2.2087 -  2.1945    1.00      544 0.1821  0.9967    143.500
  53  2.1945 -  2.1806    1.00      559 0.1618  1.0148    155.727
  54  2.1806 -  2.1671    0.99      560 0.1644  1.0273    149.309
  55  2.1671 -  2.1538    0.59      327 0.1606  1.0109    157.832
  56  2.1538 -  2.1410    1.00      530 0.1724  0.9952    140.134
  57  2.1410 -  2.1284    1.00      580 0.1718  1.0387    142.970
  58  2.1284 -  2.1161    1.00      519 0.1820  1.0275    133.995
  59  2.1161 -  2.1040    0.62      361 0.1498  1.0208    142.911
  60  2.1040 -  2.0923    0.94      518 0.1641  1.0026    127.932
  61  2.0923 -  2.0808    1.00      535 0.1517  1.0163    132.766
  62  2.0808 -  2.0696    1.00      557 0.1606  1.0615    136.788
  63  2.0696 -  2.0585    0.91      488 0.1565  1.0288    139.390
  64  2.0585 -  2.0478    0.62      345 0.1683  1.0079    124.714
  65  2.0478 -  2.0372    1.00      556 0.1563  1.0192    123.833
  66  2.0372 -  2.0269    1.00      533 0.1526  1.0235    125.338
  67  2.0269 -  2.0167    1.00      576 0.1594  1.0153    123.371
  68  2.0167 -  2.0068    0.96      533 0.1590  1.0201    120.441
  69  2.0068 -  1.9971    0.90      461 0.1653  1.0107    116.735
  70  1.9971 -  1.9875    0.90      506 0.1560  1.0107    116.930
  71  1.9875 -  1.9781    0.90      510 0.1495  0.9930    115.429
  72  1.9781 -  1.9689    0.88      467 0.1587  1.0007    107.108
  73  1.9689 -  1.9599    0.91      522 0.1605  0.9992    103.532
  74  1.9599 -  1.9510    0.91      493 0.1600  1.0016     98.184
  75  1.9510 -  1.9423    0.90      513 0.1692  0.9913     96.951
  76  1.9423 -  1.9338    0.89      468 0.1728  0.9686     96.263
  77  1.9338 -  1.9254    0.91      504 0.1705  0.9933     97.129
  78  1.9254 -  1.9171    0.91      477 0.1645  0.9803     95.672
  79  1.9171 -  1.9090    0.90      502 0.1703  0.9552     89.589
  80  1.9090 -  1.9010    0.87      491 0.1910  0.9226     83.575

done:  /home/momeara/scr/data/pdb/xu/pdb1xu1.ent
