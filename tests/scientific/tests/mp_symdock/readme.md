## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year
The benchmark was set up by Julia Koehler Leman (julia.koehler.leman@gmail.com) in March 2019.
The research was performed in Jeff Gray's lab, before Julia Koehler Leman moved to the lab of Richard Bonneau. 

## PURPOSE OF THE TEST
#### What does the benchmark test and why?
We test Rosetta's ability to dock symmetric membrane proteins, i.e. how well the symmetry framework and membrane framework RosettaMP interface with the docking protocol. We basically check how well Rosetta can recapitulate native subunit interfaces in symmetric membrane proteins, how much the symmetry framework samples in the membrane and how well the membrane scorefunction scores the decoys. 

## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?
We test on 5 proteins that were part of the original benchmark set for testing RosettaMP described in the paper (Alford, Koehler Leman et. al. PlosCompBio 2015). This benchmark set tests C2, C3, C4 and C5 symmetries. 

The input files were downloaded from the OPM database (1BL8 didn't have the correct bioassembly in PDBTM), cleaned, spanfiles were created with mp_span_from_pdb, then I used mp_range_relax to refine the structure. The outer helices in 1BL8 often came off the rest of the structure and started flopping around, the closest RMSD structures had high scores due to clashes, but I took one of those anyway as the starting structure. 

## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?

The protocol is described in (Alford, Koehler Leman et. al. PlosCompBio 2015) with command lines in the paper supplement. However, some flags had to be adjusted to make this run properly. 

Note that for -in:file:native I am using the relaxed native but the symmetry framework starts with the asymmetric unit as the input structure. Therefore, the RMSD will never be super close to the native (even though with the adjusted not-optimal flags the RMSDs are actually amazing compared to what I see in the figures in the paper supplement). Not-optimal flags means that to get decoys in a reasonable amount of time, I had so use some flags that should not be used (like one for docking prepack for instance). I can imagine that some code has changed for low-res and high-res filtering steps as well, the docking:ppk flag makes sure the models pass these filters which they otherwise would not. 

The protocol is the MPSymDock protocol, which runs the SymDock protocol underneath, which is based on the DockingProtocol. So we have a low-res docking step before it does high-res docking. Because of this, some decoys wouldn't pass the low-res stage, therefore the score file can have a mixture of centroid and high-res decoys with differing headers. That's the reason why the plotting script goes by index and not by column tag.

Runtimes are <100s per decoys: 100s x 5 proteins x 1000 decoys = 500,000 CPU secs which is <140 CPU hours. Note that the runtimes for these are with the docking:ppk flag and will likely increase (and possibly substantially) once better flags are found. 

## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?

We use interface score I_sc vs. RMSD to the native. The native here is the relaxed OPM structure, not the symmetrized version created by the make_symdef_file script. 

Cutoffs were defined by looking at the I_sc vs. RMSD distribution of the first run. The I_sc vs RMSD plots look substantially better than in the paper, not quite sure why. Differences are inclusion of the '-score:weights mpframework_docking_fa_2015.wts' flag and the '-docking::dock_ppk true'. Cutoffs need to be redone once the protocol has been optimized and better flags are found. 

Note that the I_sc in a normal range would be <0, maybe in the range between -20 to 0, not in negative hundreds as we sometimes see here. 

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset. 

The baseline are relaxed structures from the OPM database. 

1afo - C2 symmetry of single TM helix
1bl8 - C4 symmetry of potassium channel
2mpn - C2 symmetry of 2TMspan helix - looks like a handshake
2oar - C5 symmetry of a channel
2uuh - C3 symmetry of leukotriene synthase

The RMSDs that are sampled are pretty good but the interface scores are sometimes a little too good and it needs to be figured out where this is coming from. Note that this is a proof-of-concept protocol that could and should be improved code-wise and properly benchmarked. 

## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results. 

The original benchmarking was done by Rebecca Alford and this benchmark set, even though the same proteins, was set up by JKLeman, so there are slight differences. For instance, while relaxing the structures, I noticed that mp_relax moved 1 protein out of the membrane. Using mp_range_relax the 1BL8 structure of the potassium channel had the outer helices flopping around, meaning they didn't keep their interface with the other helices. This makes me wonder how stable these structures are in the first place and should be played around with more!!!

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?

The protocol should certainly be improved code-wise. This protocol is a proof-of-concept so far, needs to be benchmarked on a larger dataset with more diverse structures and symmetries. I'd like to see some beta-barrels in the benchmark set. Quality measures are standard, so they are fine. 

In an earlier version of this protocol the funnel plots looked pretty great, but again, these are docking_ppk structures and had horrible total scores due to clashes. I also noticed that the models (without docking_ppk flag) were too close and tight and therefore have horrible scores compared to the natives. So adjusting the protocol to move the structures slightly further apart (perpendicular to the symmetry axis) should give better scoring models that are closer to the native without needing the docking_ppk flag.

