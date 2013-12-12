#!/usr/bin/env python

# [2/devel] chef: tail output.txt protocols.looprelax: ===  Getting Statistics 
# protocols.looprelax: ===
# protocols.looprelax: ===
# protocols::checkpoint: Deleting checkpoints of Loopbuild
# protocols.loop_build.LoopBuildMover: loop_rms 0.498187
# protocols.loop_build.LoopBuildMover: total_energy -674.763
# protocols.loop_build.LoopBuildMover: chainbreak 0.0321136
# protocols.jd2.JobDistributor: 1srp_0001 reported success in 1238 seconds
# protocols.jd2.JobDistributor: no more batches to process... 
# protocols.jd2.JobDistributor: 1 jobs considered, 1 jobs attempted in 1239 seconds
# 
# [2/devel] chef: cat score.sc 
# SEQUENCE: 
# SCORE: total_score chainbreak  dslf_ca_dih dslf_cs_ang dslf_ss_dih dslf_ss_dst      fa_atr      fa_dun fa_intra_rep      fa_pair       fa_rep       fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb loop_rms      omega   p_aa_pp pro_close      rama       ref total_energy  description 
# SCORE:    -674.763       0.032       0.000       0.000       0.000       0.000   -1808.915     305.711        5.419      -36.300      396.301      965.860    -114.668    -145.631     -38.193     -59.127     0.498     6.096   -53.231     1.777   -12.340   -87.520      -674.763 1srp_0001

import os, glob
import argparse
from pylab import *

parser = argparse.ArgumentParser()
parser.add_argument('cluster')
arguments = parser.parse_args()

pattern = os.path.join(arguments.cluster, '*/output.txt')
outputs = glob.glob(pattern)

scores = []
rmsds = []

for output in outputs:
    with open(output) as file:
        for line in file:
            if line.startswith('protocols.loop_build.LoopBuildMover: loop_rms'):
                rmsd = float(line.split()[2])
                rmsds.append(rmsd)
            if line.startswith('protocols.loop_build.LoopBuildMover: total_energy'):
                score = float(line.split()[2])
                scores.append(score)

print scores
print rmsds

title(arguments.cluster)
plot(rmsds, scores, 'o')
xlabel('RMSD')
ylabel('Score')

if not os.fork():
    show()
