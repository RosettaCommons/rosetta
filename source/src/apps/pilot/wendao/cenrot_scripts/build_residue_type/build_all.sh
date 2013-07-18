for aa in `cat aa.lst`
do
	echo $aa
	./build_restypes.py $aa > residue_types/${aa}.params
done
