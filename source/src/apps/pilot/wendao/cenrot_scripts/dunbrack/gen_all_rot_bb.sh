cat nnrot.dat | while read aa
do
#	./kmeans_bb_dependent_rotlib.py $aa
	./kmeans_adaptive_kernel_density_bb_dependent_rotlib.py $aa
done
