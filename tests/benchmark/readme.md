# BenchmarkServer test scripts and emulation API
-----
#### BenchmarkServer API
Python scripts in `main/tests/benchmark/tests/` is scripts used to run tests by Testing Server.
1. Each script treated as test suite and could run family of tests.
2. For example `build.py` provide `build.debug`, `build.release`, `build.header`, `build.ui` and others tests.
3. Each test name could have have any characters in they name (including dot `.`) except forward slash (`/`).
4. It is also possible to add test suite as directory instead of single file. That way test can have access to supplemental files in that directory and clearly 'own them' without polluting script directories with files specific to particular tests. To define test-suite as directory you will need to place script file into `.../suite-name/command.py` file.
5. In addition to `suite-name.py` files each test can also supply optional documentation in markup format by providing `suite-name.md` and/or `suite-name.<test-name>.md` files. If documentation files is present they will be rendered in Benchmark web site on test page. If test-suite defined as directory instead of file then documentation file should be placed into `.../suite-name/description.md` file.

-----
#### Emulation API
To facilitate development of a new test this directory provide `benchmark.py` script which is designed to emulate BenchmarkServer API and allow developers to run test locally. A few notes:

1. Minimal requirements for Python version is not Python-3.6, so you will need to have Python-3.6 or later installed on your local machine in order to run tests in emulation mode.
2. To run test use `cd main/tests/benchmark && python3.6 benchmark.py <test-suite>.<test-name>` (you can also use path to `test-script.py` instead of specifying `test-suite.test-name` pair).
3. `benchmark.py` accept various command line options which allow you to speficy test platform and resources available to run the test. You can also create `benchmark.linux.ini` and `benchmark.mac.ini` files to specify default platform (see `benchmark.ini.template` for example).
4. You can use `--skip-compile` and `--debug` command line flags to instruct test-suites to run test without compilation (reuse previously created binaries) and in debug mode (boolean `debug` is passed as an argument to test-suite so developers could use this flag to by-pass certain test steps if needed)
5. When test completed in emulation mode the results will became available in `main/tests/benchmark/results/test-suite.test-name/` directory.
6. To see full list of options available for `benchmark.py` run `benchmark.py --help` in terminal.
----
#### emulation API for step-by-step tests
Some tests might considerable amount of CPU time to complete. So during development of such tests it is impractical to always run all steps of the tests due the time and computation resources required to do this. To solve this we provide step-by-step API which allow us to split test into multiple Python scripts independent of each other. To run tests which use step-by-step API in emulation mode you will need to:
1. At first run test using normal emulation API (see above) by invoking `benchmark.py` script. This will create `main/tests/benchmark/results/test-suite.test-name/` directory and symlink all tests scripts and data files and create necessary config files.
2. After initial setup of step [1] is completed you can safely terminate test run by using `ctrl-c`.
3. To continue execution of test or to re-run any particular steps: cd into test result directory (ie `cd main/tests/benchmark/results/test-suite.test-name/` and initiate step execution by using `python-3.6 <step-script-name>.py`
