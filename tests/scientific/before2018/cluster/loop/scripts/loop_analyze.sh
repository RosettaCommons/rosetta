#!/bin/bash
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

proteins="154l  1ctm  1dts  1ede  1hfc  1msc  1pbe  1prn  1rro  1srp  1thg  1tib  1xif  2cyp  2exo  2rn2  2sns  3cla  3hsc  4enl
1arp  1cyo  1eco  1ezm  1ivd  1onc  1pmy  1rcf  1scs  1tca  1thw  1tml  2cpl  2ebn  2pgd  2sil  2tgi  3cox  451c  4i1b"

echo Name MedianLoopRMS 10%LoopRMS > ./../.results.log
for i in $(echo $proteins ); do
	./../scripts/sd $i.out looprms > scorerms
	rms50=$(cat scorerms | ./../scripts/stat.py | ./../scripts/col 9 )
	rms5=$(cat  scorerms | ./../scripts/stat.py | ./../scripts/col 8 )
	echo $i $rms50 $rms5 | awk '{ printf "%5s %8.2f %8.2f\n",$1,$2,$3 }' >> ./../.results.log
done

rm -f scorerms
