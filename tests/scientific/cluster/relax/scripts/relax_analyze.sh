#!/bin/bash
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

proteins="1a19  1acf  1b3a  1bkr  1c8c  1cei  1dhn  1enh  1fkb  1hz6  1kpe  1nps  1ptq  1scj  1tif  1tul  1urn  1vie  1wit  2chf  4ubp \
1a32  1ail  1bgf  1bm8  1c9o  1cg5  1e6i  1ew4  1fna  1ig5  1lis  1opd  1r69  1shf  1tig  1ubi  1utg  1vls  256b  2ci2  5cro \
1a68  1aiu  1bk2  1bq9  1cc8  1ctf  1elw  1eyv  1gvp  1iib  1lou  1pgx  1rnb  1ten  1tit  1ugh  1vcc  1who  2acy  2vik"

echo Name DecoyRMS DecoyEnergy  NativeRMS  NativeEnergy > ./../.results.log
for i in $(echo $proteins ); do
	./../scripts/sd $i/$i.out rms score description > scorerms
	decoyrms=$(grep -v native scorerms | ./../scripts/stat.py | ./../scripts/col 9 )
	decoyenergy=$(grep -v native scorerms | ./../scripts/col 2 | ./../scripts/stat.py | ./../scripts/col 9 )
	nativerms=$(grep  native scorerms | ./../scripts/stat.py | ./../scripts/col 9 )
	nativeenergy=$(grep  native scorerms | ./../scripts/col 2 | ./../scripts/stat.py | ./../scripts/col 9 )
	echo $i $decoyrms $decoyenergy $nativerms $nativeenergy | awk '{ printf "%5s %8.2f %8.1f %8.2f %8.1f\n",$1,$2,$3,$4,$5 }' >> ./../.results.log
done

rm -f scorerms
