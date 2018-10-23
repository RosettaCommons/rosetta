# How to set up a scientific benchmark test

## Scientific test directory structure
Several tests are located in the `Rosetta/main/tests directory`. The directory structure is the following:

* `Rosetta/main/tests/integration` contains integration tests
* `Rosetta/main/tests/benchmark` contains files required for the benchmark test server, i.e. the framework which runs the scientific tests, you might have to look them over when debugging
* `Rosetta/main/tests/scientific` contains the scientific tests
* `Rosetta/main/tests/scientific/tests` contains the implementations of the tests, one directory per test
* `Rosetta/main/tests/scientific/data` submodule that contains the input data if >5 MB per test
* `Rosetta/main/tests/scientific/tests/_template_` template directory with all necessary input files

## Where and how to add your test
1.	To set up a scientific benchmark, you will need to think about
    * what your input data are
        * you will want to have a diverse benchmark set that is realistic and not just structures that give you the best possible output
    * the command line you are running
        * make sure you are running either a baseline protocol (for instance de novo structure prediction, ddg_monomer, etc) that has been run by many people for many years or a state-of-the-art protocol: since these tests are computationally expensive, we need to ensure what we are testing makes sense
    * what the output data looks like
        * in many cases we analyze score files but you might want to look at something different
    * what your quality measures are
        * are the quality measures you are using reasonable? Are there better measures in existence? 
    * what cutoffs you want to define to make the test a pass/fail
        * you should have a general idea when you run the test locally but if you don’t know the specific cutoffs yet, don’t worry; you can set them later (see below)
        * however, you want to think about HOW you define a pass/fail. It is encouraged to define it with respect to cutoffs that you define, i.e. running a single, self-contained test defines the initial cutoffs. We don’t prohibit cutoffs that are defined over time, for example running it 10x and comparing to previous runs (similar to regression tests), but this behavior is not necessarily encouraged. 
    * the entire test should not run more than 1000-2000 CPU hours
        * the runtime depends on the number of structures you run the benchmark on, the number of output models, and the protocol runtime per output model
2.	`cd Rosetta/main/tests/scientific`
3.	run `git submodule update --init —recursive` to get the submodule containing the input data. You will now see a ` Rosetta/main/tests/scientific/data` directory
4.	if your input data is larger than 5 MB, create a directory in the data submodule directory where you drop your input files, otherwise you can keep them in the tests directory (see below)
    * commit your changes in that submodule
    * move out of the data directory (`cd ..`) and commit your changes again
    * [during that process: if you use zsh or check for `git status` you will notice when git complains about uncommitted files]
5.	`cd Rosetta/main/tests/scientific/tests`
6.	copy the contents of the template directory to create your own scientific test directory via `cp -r Rosetta/main/tests/scientific/tests/_template_ <my_awesome_test>`
7.	`cd <my_awesome_test>`: tests are set up by individual steps numbered sequentially so that they can be run individually without having to rerun the entire pipeline. Don’t run anything yet (we’ll get to that further below). 
8.	Put together your test, don’t run anything yet: edit each file in this directory wherever you find the `==> EDIT HERE` tags. 
    * if you are not running a RosettaScripts protocol, you can delete the fast_relax.xml file
    * `0.compile.py`: this script is for compilation, likely you won’t need to edit this file
    * `1.submit.py`: this contains the command line and the target proteins you are running. The debug option does not refer to the release vs debug run in Rosetta but rather to the debug version of running the scientific test for faster setup. 
    * `2.analyze.py`: this script analyzes the results from the score files or whatever files you are interested in. It reads the cutoffs file and compares this run’s results against them. It will write the `result.txt` file for easy reading. A few basic functions for data analysis are at the bottom of this script. If you need your specialized function, please add it there. 
    * `3.plot.py`: this script will plot the results via matplotlib into `plot_results.png` with subplots for each protein. It draws vertical and horizontal lines for the cutoffs. 
    * `9.finalize.py`: this script gathers all the information from the readme, the results plot and the results and creates a html page `index.html` that displays everything
    * add all relevant papers to the `citation` file
    * edit the `cutoffs` file, header line starts with `#` and has one protein per row, different measures in columns 
    * make sure you add your email to the `observers`
    * fill out all the sections in the `readme.md` in as much detail as you can. Keep in mind that after you have left your lab, someone else in the community should be able to understand, run and maintain your test! This is automatically linked to the Gollum wiki. 
    * add your test to the `run` function at the bottom of `Rosetta/main/tests/benchmark/scientific/command.py` to hook it into the test server framework
    * a couple of notes on the test:
        * we aim for plots to show the data, not just tables. Please include tables and/or the data in `json` format, but plots are highly encouraged as digesting data visually is faster and easier for us for debugging. 
        * we like `json` format because it’s easy to read the data directly into python dictionaries. 
        * a test should ideally be self-contained to obtain the cutoffs. As described above, try to avoid regression-test-like behavior where you compare to previous results or gather results from multiple runs over time. 
        * be consistent with the test name, keep a single name that is also the name of the test directory. It appears in multiple files and consistency makes debugging easier
9.	Install python3.6 or later if you haven’t already, anything earlier won’t work. Make sure you have an alias set where you an specify the python version, for instance `python3`. Every python script in the scientific tests will need the `python3` prefix to run them properly!
10.	Run your tests locally in debug mode:
    * `cd Rosetta/main/tests/benchmark`
    * `python3 benchmark.py --compiler <clang or else> --skip --debug scientific.<my_awesome_test>`
        * the `--skip` flag is to skip compilation, only recommended if you have an up-to-date version of master compiled in release mode (Sergey advises against skipping)
        * the `--debug` flag is to run in debug mode which is highly recommended for debugging (i.e. you create 2 decoys instead of 1000s)
    * this creates a directory `Rosetta/main/tests/benchmark/results/<os>.scientific.<my_awesome_test>` where it creates softlinks to the files in `Rosetta/main/tests/scientific/tests/<my_awesome_test>` and then it will likely crash in one way or another
    * for step-by-step debugging `cd Rosetta/main/tests/benchmark/results/<os>.scientific.<my_awesome_test>` and debug each script individually, starting from the lowest number, by running for instance `python3 1.submit.py`
    * note that several other files and folders are created in the process: 
        * `config.json` which contains the configuration settings
        * an `output` directory is created that contains the subdirectories for each protein
        * an `hpc-logs` directory is created that contains the Rosetta run logs. You might have to check them out to debug your run if it crashed in the Rosetta run step. 
        * each numbered script will also have an associated `.json` file that contains the variables you want to carry over into the next step
        * `9.finalize.output.json` contains all the variables and results saved
        * `plot_results.png` with the results
        * `index.html` with the gathered results, failures and details you have written up in the readme. While all the files are accessible on the test server later, this file is the results summary that people will look at
        * `output.results.json` will tell you whether the tests passed or failed
11.	Submit your branch for testing
    * once you are finished debugging locally, commit all of your changes to your branch
    * create a pull-request
    * run the test on the test server
        * since scientific tests require a huge amount of computational time, you might want to lower your `nstruct` for debugging your run on the test server. If you do that, don’t forget to increase it later once the tests run successfully
    * once the tests run as you want, merge your branch into `master`
        * The `scientific` branch is an extra branch that grabs the latest master version every few weeks to run all scientific tests on. **DO NOT MERGE YOUR BRANCH INTO THE SCIENTIFIC BRANCH!!!**
        * tell Sergey Lyskov (sergey.lyskov@gmail.com) that your test is ready to be continuously run on the scientific branch
12.	Celebrate! Congrats, you have added a new scientific test and contributed to Rosetta’s greatness. :D

# The science behind your test: Scientific test template
Please use this template to describe your scientific test in the `readme.md` as described above. Also check out the `fast_relax` test for ideas of what we are looking for. 

\## AUTHOR AND DATE
\#### Who set up the benchmark? Please add name, email, PI, month and year

\## PURPOSE OF THE TEST
\#### What does the benchmark test and why?

\## BENCHMARK DATASET
\#### How many proteins are in the set?
\#### What dataset are you using? Is it published? If yes, please add a citation.
\#### What are the input files? How were the they created?

\## PROTOCOL
\#### State and briefly describe the protocol.
\#### Is there a publication that describes the protocol?
\#### How many CPU hours does this benchmark take approximately?

\## PERFORMANCE METRICS
\#### What are the performance metrics used and why were they chosen?
\#### How do you define a pass/fail for this test?
\#### How were any cutoffs defined?

\## KEY RESULTS
\#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
\#### Describe outliers in the dataset. 

\## DEFINITIONS AND COMMENTS
\#### State anything you think is important for someone else to replicate your results. 

\## LIMITATIONS
\#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
\#### How could the benchmark be improved?
\#### What goals should be hit to make this a "good" benchmark?

