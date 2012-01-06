FEATURES SCIENTIFIC BENCHMARK

(documentation updated june 2011, by Matthew O'Meara)


INTRODUCTION


The Features Scientific Benchmark measures the differences between
structural features coming from different sample sources, e.g. crystal
structures and your favorite prediction protocol.  For documentation
please see

        https://wiki.rosettacommons.org/index.php/FeaturesScientificBenchmark

This test is regularly run on the RosettaTests cluster.  Results from
the test are available at

        http://rosettatests.graylab.jhu.edu/analyze/features


HOW TO RUN FEATURE SCIENTIFIC BENCHMARK

  * Obtain sample sources you want to compare, for example
      - pdbs of native structures
      - structures generated using your favorite Rosetta protocol
      - a pre-generated features dataset from
           https://svn.rosettacommons.org/source/mini.data/tests/scientific/cluster/features/sample_sources
  
  * extract features from the sample sources into an SQLite3 features database
      - see e.g. sample_sources/LOCAL_EXAMPLE
  
  * Run the analysis scripts from this directory
      - ./compare_sample_sources.R --help

  * Look for results in build/
