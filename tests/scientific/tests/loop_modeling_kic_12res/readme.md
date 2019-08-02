## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year
Phuong T. Nguyen, tranphuonguns@gmail.com, PI: Vladimir Yarov-Yarovoy, Jun 2019

## PURPOSE OF THE TEST
#### What does the benchmark test and why?
This scientific benchmark tests loop modeling protocol using different algorithms:
- classical methods: Cyclic coordinate descent (ccd) and kinematic loop closure (kic)
- more recent kinematic-based methods: Next-generation kic (ngk) and kic with fragments

## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?
A starndard 12-residue benchmark of 7 proteins, ranging from 100-400 amino acids taken from the Kortemme's lab full benchmark set: https://guybrush.ucsf.edu/benchmarks/benchmarks/loop_modeling
Input proteins were preminimized with Rosetta, native coordinates of loop regions were removed along with sidechains of nearby residues within 10A. Non-native conformations of these loops were then added. Reference native structures are also provided for analysis.

## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?
The Next-generation kic and kic with fragments protocols run LoopModelerMover through RosettaScripts interface. This Mover currently only supports kinematic-based sampling. This was used in the Kortemme lab's benchmark sets, but with much older Rosetta version and scoring function. Here ref2015 scoring function was used. Note that, despite simulation configuration was called as config="kic", this mover appears to invoke ramping on repulsive and rama by default. From my understanding, this is NGK method and being called "kic" in this Mover. RMSDMetric was used to calculate difference in generated loop conformations compared to native. No superposition was performed as the protocol does not move other parts of protein.

The classic protocols Cyclic coordinate descent (ccd) and kinematic loop closure (kic) run through loopmodel application. Rmsd calculation is intergrated in the application.

The ccd and kic with fragments used fragments insertion while classical kic and ngk do not.
## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?
total_score and rmsd (loop_rms) to the native were used as metrics. These values are plotted with expection of forming a funnel shape. However, in actual modeling process, we don't know the native conformation, which makes the score vs rmsd plot is not very useful. We rely mostly on the far-perfect Rosetta scoring function to make decision. Here we chose 10% cutoff for score vs rmsd to describe ultimate pass/fail. The test will pass if we can find any structures in the top 10% scoring models that has rmsd to the native less than a cutoff value, which is 1A in this case for all structures.

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset. 
Due to time-comsuming nature of kinematic sampling and high numbers of nstructs needed to somewhat qualitatively describe efficiency of the protocol, initial local runs severly cut down the numbers of sampling cycles in all methods. Rmsd cutoff of 1A resulted in a significant numbers of failures. However, we still see funnel-like shape in some of proteins. The structure 2tgi has positive scores, probably because of not properly preminimized before loop modeling. However, it does not appear to severely affect the loop modeling accuracy.
This will be updated with more details after full-test runs on server. 

## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results. 
All methods need both ends of loop regions to perform kinematic sampling, thus cannot model loop regions at termini.
Cyclic coordinate descent (ccd) is the only method can do the job, with condition not using kinematic sampling in the refinement stage. "refine_ccd" needs to be used, instead of "refine_kic".
Phuong T. Nguyen (tigerous) is not an expert in Rosetta loop modeling. 

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?
The result is decent, despite of the severely cut-down numbers of cycles in initial local tests.
The dataset is not big and diverse enough to draw a full conclusion, but maybe good enough to describe a working algorithm.