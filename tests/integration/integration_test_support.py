import filecmp
import difflib
from os import path
import os
import json

import pygments
import pygments.formatters
import pygments.lexers

class IntegrationConfig(object):
    def __init__(self, sub_test_parameters = ["database", "minidir", "bin"], common_patterns = None, ignore_files=["command.sh"]):
        """Configuration for integration analysis.

        sub_test_paramters - Replace test parameters values with placeholders during diff.        

        common_patterns - Equivalent patterns to be ignore in diff generation. Provided as dict of form:
                            { pattern_name : (ref_pattern, new_pattern) }
                 
                          Eg. { "branch_directory" : ("build_dir/branch_a", "build_dir/branch_b") }

        ignore_files - Files to ignore. Defaults to ["command.sh"]

        """

        self.sub_test_parameters = sub_test_parameters
        self.common_patterns = common_patterns
        self.ignore_files = ignore_files

class IntegrationRunParameters(object):
    def __init__(self, ref_parameters, new_parameters):
        self.ref_parameters = ref_parameters
        self.new_parameters = new_parameters

    @classmethod
    def from_run_directories(cls, ref_dir, new_dir):
        with open(path.join(ref_dir, "test_parameters.json")) as ref_file:
            ref_parameters = json.load(ref_file)

        with open(path.join(new_dir, "test_parameters.json")) as new_file:
            new_parameters = json.load(new_file)

        return IntegrationRunParameters(ref_parameters, new_parameters)
        
class IntegrationAnalysis(object):
    def __init__(self, ref_dir, new_dir, config = IntegrationConfig()):
        """Configure integration analysis.

        ref_dir - Reference result directory.
        new_dir - New result directory.
        config - IntegrationConfig config object.

        """

        self.ref_dir = ref_dir
        self.new_dir = new_dir
        self.config = config
        self.run_parameters = IntegrationRunParameters.from_run_directories(ref_dir, new_dir)
        
        self.perform_analysis()
    
    def perform_analysis(self):
        dir_comparer = filecmp.dircmp(self.ref_dir, self.new_dir, ignore=self.config.ignore_files)
        
        self.tests_removed = dir_comparer.left_list
        self.tests_added = dir_comparer.right_list
        
        self.test_results = {}
        
        for test in dir_comparer.common_dirs:
            self.test_results[test] = \
                IntegrationTestAnalysis(
                    path.join(self.ref_dir, test),
                    path.join(self.new_dir, test),
                    self.config,
                    self.run_parameters)
    
    def write_test_report(self, outdir):
        if not path.exists(outdir):
            os.makedirs(outdir)

        with open(path.join(outdir, "tests_failed.txt"), "w") as failed, open(path.join(outdir, "tests_passed.txt"), "w") as passed:
            for test, result in self.test_results.iteritems():
                if result.passed:
                    passed.write(test + "\n")
                else:
                    failed.write(test + "\n")
                    result.write_test_report(path.join(outdir, test))

class IntegrationTestAnalysis(object):
    def __init__(self, ref_dir, new_dir, config = IntegrationConfig(), run_parameters=None):
        """Result analysis for a single integration test result.

        ref_dir - Reference result directory.
        new_dir - New result directory.
        config - IntegrationConfifg config object.
        """
        
        self.ref_dir = ref_dir
        self.new_dir = new_dir
        self.config = config
        self.run_parameters = run_parameters
        
        self.perform_analysis()

    def read_ref_result_file(self, result_file):
        """Filter ref result file."""
        if self.config.sub_test_parameters:
            result = result_file.read()
            for param, value in self.run_parameters.ref_parameters.items():
                if param in self.config.sub_test_parameters:
                    result = result.replace(value, param)

            return result.split("\n")
        else:
            return result_file.readlines()

    def read_new_result_file(self, result_file):
        """Filter new result file."""
        if self.config.sub_test_parameters:
            result = result_file.read()
            for param, value in self.run_parameters.new_parameters.items():
                if param in self.config.sub_test_parameters:
                    result = result.replace(value, param)

            return result.split("\n")
        else:
            return result_file.readlines()
    
    def perform_analysis(self):
        dir_comparer = filecmp.dircmp(self.ref_dir, self.new_dir, ignore=self.config.ignore_files)
        
        candidate_diffs = []
        candidate_diffs.extend(dir_comparer.diff_files)
        for subdir, subdir_diff in dir_comparer.subdirs.iteritems():
            candidate_diffs.extend(path.join(subdir, f) for f in subdir_diff.diff_files)
        
        self.diff_files = [target_file for target_file in candidate_diffs if self.perform_file_diff(target_file)]

    def perform_file_diff(self, target_file):
        with open(path.join(self.ref_dir, target_file)) as ref_file, open(path.join(self.new_dir, target_file)) as new_file:
            return list(difflib.unified_diff(
                self.read_ref_result_file(ref_file), self.read_new_result_file(new_file),
                path.join(self.ref_dir, target_file), path.join(self.new_dir, target_file),
                n=3))
                
    @property
    def passed(self):
        return not self.diff_files
    
    def write_test_report(self, outdir):
        for f in self.diff_files:
            output_basename = path.join(outdir, f)
            if path.dirname(output_basename) and not path.exists(path.dirname(output_basename)):
                os.makedirs(path.dirname(output_basename))

            #Write raw diffs
            with open(output_basename + ".diff", "w") as diffout:
                diffout.writelines(( l + "\n" for l in self.perform_file_diff(f)))

            #Write html highlighted diffs
            with open(output_basename + ".diff") as diffin, open(output_basename + ".html", "w") as htmlout:
                pygments.highlight(
                        "".join(diffin.readlines()),
                        pygments.lexers.get_lexer_by_name("diff"),
                        pygments.formatters.get_formatter_by_name("html", full=True),
                        htmlout)

            #Write ref and new file versions
            if not path.exists(path.dirname(path.join(outdir, "ref", f))):
                os.makedirs(path.dirname(path.join(outdir, "ref", f)))
            with open(path.join(self.ref_dir, f), "r") as infile, open(path.join(outdir, "ref", f), "w") as outfile:
                outfile.writelines(self.read_ref_result_file(infile))

            if not path.exists(path.dirname(path.join(outdir, "new", f))):
                os.makedirs(path.dirname(path.join(outdir, "new", f)))
            with open(path.join(self.new_dir, f), "r") as infile, open(path.join(outdir, "new", f), "w") as outfile:
                outfile.writelines(self.read_new_result_file(infile))
