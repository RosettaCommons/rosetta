
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

1) Create a file <TestName>.bench.hh
2) Subclass Benchmark and implement the 'setUp', 'run' and 'tearDown' functions.
3) Instantiate an instance of your class at the top of benchmark.cc with a name

############################################
#To run the performance benchmarks locally:#
############################################

1) Compile the benchmark application in release mode 
    
     cd source
     ./scons.py performance_benchmark mode=release -j<n_cores>

2) Run the benchark application from this directory

     cd rosetta_source/src/apps/benchmark
     ../../../bin/performance_benchmark.default.<platform/compiler>release \
       -database <path_to_database> \
       -benchmark_scale <multiple run length by this value> \
       [-run_one_benchmark <name of benchmark>]

3) Look at output for results --or--
4) Look at rosetta_source/src/apps/benchmark/performance/_performance_ for results
