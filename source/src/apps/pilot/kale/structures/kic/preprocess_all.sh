#!/usr/bin/env sh

trap "" HUP

for pdb in *.loop; do
    pdb=${pdb%.*}

    if [ ! -e $pdb.pdb ]; then 
        echo "Preprocessing $pdb..."
        ./preprocess.sh $pdb &
    fi
done

wait
