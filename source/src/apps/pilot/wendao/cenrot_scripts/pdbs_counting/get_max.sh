for aa in `cat aa.lst`
do
	sort -n $aa.dat | tail -n 20 > $aa.last
	echo $aa
	#drop the last 19 value may caused by wrong bb str
	head -n1 $aa.last
done
