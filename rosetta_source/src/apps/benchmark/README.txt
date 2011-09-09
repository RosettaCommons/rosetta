
Performance Benchmarks

The performance benchmarks in this folder are short pieces of code
that are run many times in order to accurately measure how long they
takes to run. In contrast, profile tests, located in
rosetta_tests/profile are whole protocols that are run on realistic
data to measure how long they take to run.

See the performance benchmarks section on RosettaTests

    http://rosettatests.graylab.jhu.edu/analyze



##################################
# To add a performance benchmark #
##################################

1) Create a file <mytest>.bench.hh
2) Add your file to rosetta_source/src/apps.src.settings
3) Subclass Benchmark and implement the 'setUp', 'run' and 'tearDown' functions.
4) Instantiate an instance of your class at the top of benchmark.cc with a name

############################################
#To run the performance benchmarks locally:#
############################################

1) Compile the benchmark application in release mode 
    
     cd rosetta_source
     ./scons.py bin/benchmark.<platform/compiler>release mode=release -j<n_cores>

2) Run the benchark application from this directory

     cd rosetta_source/src/apps/benchmark
     ../../../bin/benchmark.default.<platform/compiler>release \
       -database <path_to_database> \
       -benchmark_scale <multiple run lenght by this value> \
       [-run_one_benchmark <name of benchmark>

3) Look at output for results --or--
4) Look at rosetta_source/src/apps/benchmark/_performance_ for results
