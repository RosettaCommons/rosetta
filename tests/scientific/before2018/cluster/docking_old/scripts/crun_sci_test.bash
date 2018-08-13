if [ -z "$2" ]
then
    echo "Usage: crun.bash <pdb> <config>"
    echo "Sets up and launches a condor run"
    echo "<config> is the name of the condor config file"
    echo "without the '.config' extension"
    echo "The extension on the config file must be '.config'"
    exit
fi

pdb=$1

#Name of config file (not including extension)
#extension must be '.config'
config=$2

#Check for files before starting condor

if [ ! -f $pdb.pdb ]
then
    echo $pdb.pdb not found -- exiting!
    exit
fi

#Set up directories (replaces pdb_dir_maker.pl)
mkdir -p $pdb/outerr

exe=$pdb.$config.bash
cscript=$pdb.$config.con

#Get variables from the config file
source $config.config

#Set up output directory
mkdir -p $pdb/$prefix
pdbpath=../$pdb/$prefix/

#Set up command line
command="../scripts/rrun_mini.sh $compiler $pdb $prefix -out:prefix $prefix $protocol_flags -s ../$pdb.pdb -database $database -out::path:pdb $pdbpath -out:file:o $pdb"

#Set up an executable to run RosettaDock
echo \#\!/bin/bash > $exe

echo $command >> $exe
chmod +x $exe

#Make a relevant condor script
cp ../scripts/base.con $cscript

echo "Executable = ./$exe" >> $cscript
echo "pdb = $pdb" >> $cscript
echo "Queue $Njobs" >> $cscript

condor_submit $cscript


