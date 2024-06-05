#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   scientific/protein_data_bank_diagnostic.py
## @brief  Benchmark script for testing how Rosetta can load PBD's from www.rcsb.org
## @author Sergey Lyskov

import os, json, enum, time, shutil
from collections import namedtuple, OrderedDict

script_name = os.path.abspath(__file__)  # keep this line above imp.load_source(...) line below, beacuse later change value of __file__ variable

_api_version_ = '1.1'

_number_of_jobs_ = 1024

# A bit of Python magic here, what we trying to say is this: from ../../__init__ import *, but init is calculated from file location
import importlib.util, sys
importlib.util.spec_from_file_location(__name__, '/'.join(__file__.split('/')[:-3]) +  '/__init__.py').loader.exec_module(sys.modules[__name__])

Job = namedtuple('Job', 'name pdbs path rosetta_dir command_line')  #  hpc_job_id
#Job.__new__.__defaults__ = (None, )


class TestMode(enum.Enum):
    def __str__(self): return f'{self.name}'  # (self.__class__.__name__, self.name)

    fast = enum.auto()
    full = enum.auto()

    cif = enum.auto()


class InputType(enum.Enum):
    def __str__(self): return f'{self.name}'  # (self.__class__.__name__, self.name)

    pdb = enum.auto()
    cif = enum.auto()


_timeouts_ = {
    TestMode.fast : 8,
    TestMode.full : 32,

    TestMode.cif  : 8,
}


def download_all_input_files(kind, config):
    ''' Download all PDB files from www.rcsb.org and store them at given prefix (from config).
        Return a map PDB-code --> absolute path to PDB file
    '''
    prefix = os.path.abspath(config['prefix'] + f'/rcsb/{kind}')
    if not os.path.isdir(prefix): os.makedirs(prefix)

    # taken from https://www.rcsb.org/pages/download/ftp
    # Rsync only the PDB format coordinates  /pub/pdb/data/structures/divided/pdb (Aproximately 20 GB)
    #${RSYNC} -rlpt -v -z --delete --port=$PORT ${SERVER}/data/structures/divided/pdb/ $MIRRORDIR > $LOGFILE 2>/dev/null
    #SERVER=rsync.rcsb.org::ftp                                   # RCSB PDB server name
    #PORT=33444                                                    # port RCSB PDB server is using

    suffix = { InputType.pdb: 'pdb', InputType.cif : 'mmCIF'} [kind]

    if config['debug']: print('WARNING: Debug mode, skipping PBD downloading...')
    else: execute(f"Downloading {str(kind).upper()}s from www.rcsb.org...", f'cd {prefix} && rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/{suffix}/ {prefix}')

    pattern = { InputType.pdb: 'pdb*.ent.gz', InputType.cif : '*.cif.gz'} [kind]
    files = execute(f"Getting list of available {str(kind).upper()}s...", f"cd {prefix} && find -name '{pattern}'", return_='output', silence_output=True).split()

    #print(f'files:{files}')

    pdbs = {}
    for p in files: # expecting ./34/pdb134d.ent.gz or ./34/134d.cif.gz
        if   kind == InputType.pdb: code = p.split('/')[-1][len('pdb'):-len('.ent.gz')]
        elif kind == InputType.cif: code = p.split('/')[-1][:4]

        pdbs[code] = prefix + p[1:]  # removing leading '.'

    return pdbs


_index_html_template_ = '''\
<html>
<H2>Scientific test: Protein Data Bank diagnostic</H2>
<h3>## AUTHOR AND DATE
</h3>
<p>Steven M Lewis, smlewi@gmail.com, Cyrus Biotechnology, Oct 2018</p>
<p>Sergey Lyskov, sergey.lyskov@gmail.com, Gray lab, Johns Hopkins University, Oct 2018</p>
<p>
</p>

<h3>FAILURES</h3>
<p>
The PDB diagnostic was run on <b>{len_pdbs:,}</b> PDB files. <b>{len_passed:,}</b> files passed the test and <b>{len_failed:,}</b> failed!
</p>

<head>
    <meta charset="utf-8"/>
    <title>Protein Data Bank Diagnostic test results</title>

  <style>
    fixed {{background-color: #eee; white-space: pre-wrap; font-family: Monaco, 'Liberation Mono', Courier, monospace; font-size:12px; }}
  </style>
</head>
<body>

<h3>## PURPOSE OF THE TEST</h3>
<p>
This test tries to load every PDB file in the PDB database and classifies the failures that occur.
The command line below shows what was done; broadly all versions of this test examine load-time problems and more expensive versions
(<fixed>-PDB_diagnostic::skip_pack_and_min&nbsp;false</fixed>) also check for errors during scoring, packing, and minimization.
</p>

<p><i>"Hunting down these bugs is the most fun thing you can do on a Thursday morning"</i> - Andy Watkins, probably.</p>

<p>An individual PDB passes or fails this test based on whether it errors out or completes the diagnostic.  The test as a whole passes or fails based on a "reference results" system, like an expected result in a unit test.  About 700 PDBs and 1200 CIFs fail at the time of this writing; the purpose of the test is to document the failures and watch for new ones, so PDBs failing in an expected manner does not constitute an overall test failure.  The test will fail if PDBs pass or fail <b>UNEXPECTEDLY</b>, where the expectation is defined by the reference results (see below).</p>

<p>If you find that this page is telling you the test failed because there are <b>FEWER</b> errors: <b>GREAT</b>!  You fixed some bugs!  You can update the reference results following the instructions at the bottom of the page.</p>

<p>If you find that there are <b>MORE</b> failures than expected, especially non-timeout failures, consider it a warning that recent code changes may have introduced bugs into the PDB reading machinery. If you don't know what's going on: post to Slack or devel. <b>DO NOT</b> just update the reference results in this case.</p>

<p>If you want to know more about what this test does - pester Steven Lewis to write proper documentation.</p>

<br/>
<p>Command line used: <fixed>{command_line}</fixed>
<br/>
<p>
    <b>{len_failed:,}</b> total PDBs failed with the following error codes:
    <div style="white-space: pre-wrap; font-family: Monaco, 'Liberation Mono', Courier, monospace; font-size:12px;">{failed}</div>
</p>
<br/>
<p>
    {explanation}
    <div style="white-space: pre-wrap; font-family: Monaco, 'Liberation Mono', Courier, monospace; font-size:12px;">{explanations}</div>
</p>
{note}
<br/>
<p>
 To update reference results please copy the files below into the main repository:<br/>
 &nbsp;&nbsp;&nbsp;&nbsp;<a href="reference-results.{mode}.new.json">reference-results.{mode}.new.json</a> → <a href="https://github.com/RosettaCommons/rosetta">「Rosetta repository」</a> as <a href="https://github.com/RosettaCommons/rosetta/tree/main/tests/benchmark/tests/scientific/protein_data_bank_diagnostic/reference-results.{mode}.json">tests/benchmark/tests/scientific/protein_data_bank_diagnostic/reference-results.{mode}.json</a><br/>
 &nbsp;&nbsp;&nbsp;&nbsp;<a href="blocklist.{mode}.new.json">blocklist.{mode}.new.json</a> → <a href="https://github.com/RosettaCommons/rosetta">「Rosetta repository」</a> as <a href="https://github.com/RosettaCommons/rosetta/tree/main/tests/benchmark/tests/scientific/protein_data_bank_diagnostic/reference-results.json">tests/benchmark/tests/scientific/protein_data_bank_diagnostic/blocklist.json</a> (note the cif/pdb/fast mode is ignored: the blocklist is segfaults and huge PDBs and is shared)<br/>
</p>
</body></html>
'''

_code_html_template_ = '''\
<html><head>
<title>Protein Data Bank Diagnostic test results for PDB failed with {code} error</title>
<pre>{logs}</pre>
</body></html>
'''
#<pre style="white-space: pre-wrap; word-wrap: break-word;">{command_line}</pre></p>


def get_blocklist(rosetta_dir, mode):
    with open(f'{rosetta_dir}/tests/benchmark/tests/scientific/protein_data_bank_diagnostic/blocklist.{mode}.json') as f: blocklist = json.load(f)
    return blocklist


def remove_blocklisted_pdbs(pdbs, rosetta_dir, mode):
    for p in get_blocklist(rosetta_dir, mode)['ignore']: pdbs.pop(p, None)


def protein_data_bank_diagnostic(mode, rosetta_dir, working_dir, platform, config, hpc_driver, verbose, debug):
    #exit_code, output = execute('Checking if execute_through_pty is working...', 'echo "start" && sleep 10 && echo "one" && sleep 10 && sleep 10 && sleep 10 && echo "almost finish!" && sleep 10', return_='tuple')
    #print(f'Result: exit_code:{exit_code}, output:{output}')
    # for i in range(1000):
    #     exit_code, output = execute('Checking if execute_through_pty is working...', 'echo "start" && sleep 0.1', return_='tuple')


    res, output, build_command_line = build_rosetta(rosetta_dir, platform, config)
    with open(working_dir + '/build-log.txt', 'wb') as f: f.write( to_bytes(output) )

    if res: return { _StateKey_ : _S_build_failed_,  _ResultsKey_ : {},
                     _LogKey_ : 'Building rosetta failed!\n{}\n{}\n'.format(build_command_line, output) }
    else:
        extension = calculate_extension(platform)

        command_line = f'{rosetta_dir}/source/test/timelimit.py {_timeouts_[mode]} {rosetta_dir}/source/bin/PDB_diagnostic.{extension} -no_color -out:file:score_only /dev/null -jd2::delete_old_poses true -ignore_unrecognized_res false -load_PDB_components true -packing::pack_missing_sidechains false -packing::repack_only true -s {{input_file}}'
        command_line += ' -ignore_zero_occupancy false'
        command_line += ' -in:file:obey_ENDMDL true'  # necessary to read PDB-sourced multimodel NMR files

        if   mode in [TestMode.fast, TestMode.cif]: command_line += ' -PDB_diagnostic::skip_pack_and_min true  -PDB_diagnostic::reading_only true'
        elif mode == TestMode.full: command_line += ' -PDB_diagnostic::skip_pack_and_min false -PDB_diagnostic::reading_only false'
        else: raise BenchmarkError(f'protein_data_bank_diagnostic: Unknown mode {mode}!')

        input_kind = { TestMode.fast : InputType.pdb, TestMode.full : InputType.pdb, TestMode.cif : InputType.cif } [mode]

        all_pdbs = download_all_input_files(input_kind, config)

        number_of_jobs = _number_of_jobs_ if mode == TestMode.fast else _number_of_jobs_ * 1

        if debug:
            _all_pdbs_ = all_pdbs
            number_of_jobs, all_pdbs = 2, { k : all_pdbs[k] for k in list(all_pdbs.keys())[10:10] } # 18

            all_pdbs.update( { k : _all_pdbs_[k] for k in '100d 101d 1bz9 2agd 2gtx 1bbi'.split() } )  # 1ehd 106d 145d 154d 159d 6enz 6eor
            all_pdbs['aaaa'] = all_pdbs['100d']
            #number_of_jobs , pdbs = 5, { k : pdbs[k] for k in list(pdbs.keys())[:1000] }
            #pdbs = { str(i): '_'+str(i)+'_' for i in range(10) }


        remove_blocklisted_pdbs(all_pdbs, rosetta_dir, mode)


        hpc_logs = f'{working_dir}/.hpc-logs';  os.makedirs(hpc_logs)

        #script_name = f'{rosetta_dir}/tests/benchmark/tests/scientific/protein_data_bank_diagnostic.py'

        jobs,  hpc_job_ids = [], []
        keys = list( all_pdbs.keys() );  keys.sort()
        slice = ( len(keys) + number_of_jobs - 1 ) // number_of_jobs
        for i in range(0, len(keys), slice):
            job_index = i // slice

            name = f'{job_index:04d}'

            path =  f'{working_dir}/{name}'
            os.makedirs(path)

            jobs.append( Job( name = name,
                              pdbs = { k : all_pdbs[k] for k in keys[i:i+slice] },
                              path = path, rosetta_dir = rosetta_dir,
                              command_line = command_line,
            ) )
            #print(f'job:{job_index} i:{i}, len:{len(jobs[-1].pdbs)} command_line:{jobs[-1].command_line}')

            job_json_file = f'{path}/{job_index:04d}.input.json'
            with open(job_json_file, 'w') as f: json.dump(jobs[-1]._asdict(), f, sort_keys=True, indent=2)

            hpc_job_ids.append( hpc_driver.execute(executable = sys.executable,
                                                   arguments = f'{script_name} {job_json_file}',
                                                   working_dir = path, log_dir = hpc_logs, time=24,
                                                   name=f'protein-data-bank-diagnostic-{name}', block=False)
            )

        if not debug:
            hpc_driver.wait_until_complete(hpc_job_ids)
            time.sleep(64)  # waiting for NFS caching

        results = dict(passed=[], failed=dict())
        for job in jobs:
            if not os.path.exists( f'{job.path}/{job.name}.output.json' ):
                raise RuntimeError(f'Output "output/{job.name}/{job.name}.output.json" for job {job.name} not found -- check output/.hpc-logs/.hpc.protein-data-bank-diagnostic-{job.name}.errors.0.log and output/.hpc-logs/.hpc.protein-data-bank-diagnostic-{job.name}.output.0.log in the file-list-view for more details.')
            with open(f'{job.path}/{job.name}.output.json') as f: r = json.load(f)

            for code, pdbs in r['failed'].items():
                results['failed'].setdefault(code, {})
                for pdb, path in pdbs.items(): results['failed'][code][pdb] = os.path.relpath(path, working_dir)

            results['passed'].extend(r['passed'])


        with open(f'{working_dir}/results.aggregated.json', 'w') as f: json.dump(results, f, sort_keys=True, indent=2)

        len_failed = len(all_pdbs) - len(results['passed'])


        # Generating new reference results
        reference_results = dict(passed=results['passed'], failed = { code : list( results['failed'][code].keys() ) for code in results['failed'] } )
        with open(f'{working_dir}/reference-results.{mode}.new.json', 'w') as f: json.dump(reference_results, f, sort_keys=True, indent=2)

        state, explanations = _S_passed_, ''
        with open(f'{rosetta_dir}/tests/benchmark/tests/scientific/protein_data_bank_diagnostic/reference-results.{mode}.json') as f: reference = json.load(f)
        for pdb in reference['passed']:
            if pdb not in results['passed']:
                category = None
                for c, pdbs in results['failed'].items():
                    if pdb in pdbs:
                        category, log_path = c, pdbs[pdb]
                        break

                if category:
                    state = _S_failed_
                    explanations += f'PDB <a href="https://www.rcsb.org/structure/{pdb.upper()}">{pdb.upper()}</a> was passing test before but now failed with <a href="pdbs.{category}.html">「{category}」</a> error! Its run-log could be found in <a href="{log_path}">「{pdb.upper()}」</a>\n'

        explanation = '' if state == _S_passed_ else '<p>Test marked as <b>FAILED</b> due to following errors:</p>'


        # Generating new blocklist based on PDB's in `exceed_timeout` and old blocklist
        blocklist = get_blocklist(rosetta_dir, mode)['ignore']
        for p in results['failed'].get('exceed_timeout', []): blocklist.append(p)
        with open(f'{working_dir}/blocklist.{mode}.new.json', 'w') as f: json.dump( dict(ignore = sorted( set(blocklist) ) ), f, sort_keys=True, indent=2)


        new_passed_pdbs = []
        for pdb in results['passed']:
            if pdb not in reference['passed']: new_passed_pdbs.append(pdb)  # f' {pdb.upper()}'

        note = f"<br/><p>NOTE: {len(new_passed_pdbs)} PDB's passed the tests but was not listed in reference results.</p>" if new_passed_pdbs else ''

        with open(f'{working_dir}/index.html', 'w') as f:

            failed = []
            for c in sorted( results['failed'].keys(), key=lambda k: -len(results['failed'][k]) ):
                failed += f'{len(results["failed"][c]):,}'.rjust(7) + f' <a href=pdbs.{c}.html>「{c}」</a>\n'

            f.write(
                _index_html_template_.format(
                    command_line = command_line,
                    len_pdbs = len(all_pdbs), len_passed = len(results['passed']), len_failed = len_failed,
                    failed = ''.join(failed), #failed = ''.join( ( f'<p>{len(results["failed"][c]):,6} <a href=pdbs.{c}.html>「{c}」</a></p>' for c in sorted( results['failed'].keys(), key=lambda k: -len(results['failed'][k]) ) ) )
                    explanation = explanation, explanations = explanations, mode=mode, note=note,
                )
            )

        for c in results['failed']:
            with open(f'{working_dir}/pdbs.{c}.html', 'w') as f:
                f.write(_code_html_template_.format(code = c, logs = ''.join( ( f'<a href={log}>{pdb.upper()}</a>\n' for pdb, log in results['failed'][c].items() ) ) ) )

        # Here we overwrite pass/fail status with a percentage bar.
        # The original design of this code was very strict, failure occurred if any
        # PDB failed that used to pass.  This is too fragile.
        # Now we impose that 5% of the total number of PDBs can fail before failing
        # the whole test.
        # We leave all preceding logic in place to make the detailed HTML report.
        # Hard failures like a compilation problem should have already returned by now
        # so resetting `state` is ok.
        threshold_percent = 0.05  # up to 5% failure allowed
        threshold = float(len(all_pdbs)) * threshold_percent
        if (len_failed > threshold):
            state = _S_failed_
        else:
            state = _S_passed_


        return { _StateKey_  : state, _ResultsKey_ : {},
                 _LogKey_ : f'Protein Data Bank Diagnostic is finished.\n{len(results["passed"]):,} PDB files tested - {len_failed:,} structures failed.\nPlease see HTML output for details.  The test passes if more than 95% of PDBs pass.',
                 _IgnoreKey_ : [ os.path.basename(hpc_logs) + '/' + f for f in os.listdir(hpc_logs)]
        }



def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    observers_file_name = f'{repository_root}/tests/benchmark/tests/scientific/protein_data_bank_diagnostic/observers.' + test
    if os.path.isfile(observers_file_name): shutil.copy(observers_file_name, f'{working_dir}/observers')

    if test in ['', 'fast']: return protein_data_bank_diagnostic(TestMode.fast, repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == "full":     return protein_data_bank_diagnostic(TestMode.full, repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    elif test == 'cif':      return protein_data_bank_diagnostic(TestMode.cif,  repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    else: raise BenchmarkError(f'Unknown scripts test: {test}!')



class PDB_Diagnostic_Codes(enum.Enum):
    def __str__(self): return f'{self.name}'  # (self.__class__.__name__, self.name)

    unknown               = enum.auto()

    #ace                   = enum.auto()

    sugar_variant                       = enum.auto()
    zero_atom_restype                   = enum.auto()
    unknown_atom_name                   = enum.auto()
    unrecognized_residue                = enum.auto()
    unrecognized_experimental_technique = enum.auto()
    unrecognized_variant                = enum.auto()
    unrecognized_compound               = enum.auto()
    multiple_disulfides                 = enum.auto()

    partly_recognized     = enum.auto()

    fill_missing_atoms    = enum.auto()
    rotno                 = enum.auto()

    zero_length_xyzVector = enum.auto()
    zero_nres             = enum.auto()

    exceed_timeout        = enum.auto()

    # assert_segfault       = enum.auto()  # we are lumping these in with unknown
    # misc_segfault         = enum.auto()

    reroot_disconnected   = enum.auto()
    atom_base_not_bonded  = enum.auto()
    first_base_not_p      = enum.auto()
    base_of_chi           = enum.auto()
    atom_base_not_found   = enum.auto()
    reset_icoor_root      = enum.auto()

    res_map_range         = enum.auto()

    unrecognized_stub_id  = enum.auto()

    alias_missing_atom    = enum.auto()

    missing_disulfide_partner    = enum.auto()

    unknown_hbond_acceptor = enum.auto()
    no_orient_atoms        = enum.auto()
    pseudobond_connection_change = enum.auto()
    rotlib_file            = enum.auto()
    duplicate_atom_name    = enum.auto()
    cutpoint_neighbour     = enum.auto()
    multiple_residue_bonds = enum.auto()
    abase2_nbrs            = enum.auto()
    icoord_depends_on_missing = enum.auto()
    no_usable_coords       = enum.auto()
    merge_with_next        = enum.auto()
    ld_chirality           = enum.auto()
    insufficient_mainchain = enum.auto()
    incompatible_w_polymer = enum.auto()
    aa_difference          = enum.auto()
    prepro_cyclic_pep      = enum.auto()
    missing_connection     = enum.auto()
    missing_bond           = enum.auto()
    not_being_packed       = enum.auto()
    no_hbond_deriv         = enum.auto()
    nu_conformer           = enum.auto()
    restype_configure      = enum.auto()

    bad_patch              = enum.auto()


def classify_pdb_diagnostic_log(log):
    ''' Classify given log PDB_diagnostic log and return one of PDB_Diagnostic_Codes
        refactored from rosetta/tools/chem_xrw/check_logs.py
    '''
    #log = log.lower()
    lines = log.split('\n')

    P = namedtuple('P', 'previous splited next line log')
    P.__new__.__defaults__ = (None, None, None, None, None)  # all arguments is optional

    def pattern_in(pattern, previous_line, line, next_line, splited, log):
        def in_(p, l):
            if p:
                p = p if type(p) == tuple else [p]
                for e in p:
                    if e not in l: return False
                return True
            return False

        return in_(pattern.previous, previous_line) or in_(pattern.line, line) or in_(pattern.next, next_line) or in_(pattern.splited, splited) or in_(pattern.log, log)


    error_map = OrderedDict( [ # Pattern → error_code, Note that order is slightly different from original code: it is now grouped by-pattern types for readability
        # Pattern Types:
        #   * line= -- The pattern is found in a line starting with ERROR:
        #   * previous= -- The pattern is found in a line right before one starting with ERROR:
        #   * next= -- The pattern is found in a line after one starting with ERROR:
        #   * log= -- The pattern is found anywhere in the log (which also contains a line starting with ERROR:)
        #              (Avoid when possible, as log is computationally more heavy.)

        #( P(previous='ace'), PDB_Diagnostic_Codes.ace ),

        ( P(line='unable to find desired variant residue:'),    PDB_Diagnostic_Codes.sugar_variant ),
        ( P(line=('atom name', 'not available in residue'), ),  PDB_Diagnostic_Codes.unknown_atom_name ),
        ( P(line=('Atom name not found in MutableResidueType.'), ),  PDB_Diagnostic_Codes.unknown_atom_name ),
        ( P(line='Unrecognized residue:'),                      PDB_Diagnostic_Codes.unrecognized_residue ),
        ( P(line='Unrecognized experimental technique string'), PDB_Diagnostic_Codes.unrecognized_experimental_technique ),
        ( P(line='Rosetta does not recognize the variant:'),    PDB_Diagnostic_Codes.unrecognized_variant ),
        ( P(line='Unrecognized compound token string'),         PDB_Diagnostic_Codes.unrecognized_compound ),

        ( P(line='Cannot parse CIF file. No atom block (chem_comp_atom) found for UNL#'),   PDB_Diagnostic_Codes.zero_atom_restype ),  # was: 'Cannot load in ResidueType for entry with no atoms.'   #
        ( P(line='SSBond records list multiple nonredundant disulfides for this residue!'), PDB_Diagnostic_Codes.multiple_disulfides ),


        ( P(next='with 3-letter code:'), PDB_Diagnostic_Codes.partly_recognized ),

        ( P(splited='fill_missing_atoms!'),                   PDB_Diagnostic_Codes.fill_missing_atoms ),
        ( P(splited='packed_rotno_conversion_data_current_'), PDB_Diagnostic_Codes.rotno ),

        ( P(log=('Cannot normalize xyzVector of length() zero', 'src/numeric/xyzVector.hh') ), PDB_Diagnostic_Codes.zero_length_xyzVector ),
        ( P(log=('Cannot create normalized xyzVector from vector of length() zero') ),         PDB_Diagnostic_Codes.zero_length_xyzVector ),

        ( P(line='Cannot reroot a disconnected ResidueType'),                        PDB_Diagnostic_Codes.reroot_disconnected ),
        ( P(line="Attempted to inappropriately reset ICOOR root atom"),              PDB_Diagnostic_Codes.reset_icoor_root ),
        ( P(line='atoms must be bonded to set atom base'),                           PDB_Diagnostic_Codes.atom_base_not_bonded ),
        ( P(next=('In chi','the base of the','rather than the', 'atom of the chi')), PDB_Diagnostic_Codes.base_of_chi ),
        ( P(line="set_atom_base: atom names don't exist"),                           PDB_Diagnostic_Codes.atom_base_not_found ),
        ( P(line='Assertion `first_base_atom != atom_vertex( "P" )` failed'),        PDB_Diagnostic_Codes.first_base_not_p ),

        ( P(previous='Residue outside res_map range'),                               PDB_Diagnostic_Codes.res_map_range),

        ( P(line="unrecognized stub atom id type"),                                  PDB_Diagnostic_Codes.unrecognized_stub_id ),

        ( P(line="Unable to add atom alias for non-existent atom"),                  PDB_Diagnostic_Codes.alias_missing_atom ),

        ( P(line="Can't find an atom to disulfide bond"),                            PDB_Diagnostic_Codes.missing_disulfide_partner ),
        ( P(line="Could not find disulfide partner for residue"),                    PDB_Diagnostic_Codes.missing_disulfide_partner ),

        ( P(line="unknown Hydrogen Bond acceptor type for"),                         PDB_Diagnostic_Codes.unknown_hbond_acceptor ),
        ( P(line="Could not find three atoms to use for orienting residue type"),    PDB_Diagnostic_Codes.no_orient_atoms ),
        ( P(line="Unable to handle change in the number of residue connections in the presence of pseudobonds"),  PDB_Diagnostic_Codes.pseudobond_connection_change ),

        ( P(line="Error!  Could not open rotamer library file"),                     PDB_Diagnostic_Codes.rotlib_file ),
        ( P(line=("Can't add atom named","as it already has one with that name.")),  PDB_Diagnostic_Codes.duplicate_atom_name ),

        ( P(line="Error in core::conformation::check_good_cutpoint_neighbour()"),    PDB_Diagnostic_Codes.cutpoint_neighbour ),
        ( P(line="dont add residue bonds more than once!"),                          PDB_Diagnostic_Codes.multiple_residue_bonds ),
        ( P(line="failed to set abase2 for acceptor atom, it has no nbrs"),          PDB_Diagnostic_Codes.abase2_nbrs ),
        ( P(line=("The internal coordinates of ResidueType","depend on a missing atom!")), PDB_Diagnostic_Codes.icoord_depends_on_missing ),
        ( P(line="No usable coordinates for mmCIF file"),                            PDB_Diagnostic_Codes.no_usable_coords ),
        ( P(line="has been indicated to merge with the next residue from the"),      PDB_Diagnostic_Codes.merge_with_next ),
        ( P(line="Needed atoms missing for detect_ld_chirality_from_polymer_residue"), PDB_Diagnostic_Codes.ld_chirality ),
        ( P(line="Insufficient mainchain atoms for residue"),                        PDB_Diagnostic_Codes.insufficient_mainchain ),
        ( P(previous=("Can't create a polymer bond","due to incompatible type")),    PDB_Diagnostic_Codes.incompatible_w_polymer ),
        ( P(line="Assertion `! found_aa_difference` failed"),                        PDB_Diagnostic_Codes.aa_difference ),
        ( P(line="The RamaPrePro term is incompatible with cyclic dipeptides"),      PDB_Diagnostic_Codes.prepro_cyclic_pep ),
        ( P(line="Residues which were assumed to be connected are not"),             PDB_Diagnostic_Codes.missing_connection ),
        ( P(line="Assertion `res2.is_bonded(res1)` failed"),                         PDB_Diagnostic_Codes.missing_bond ),
        ( P(line="Assertion `being_packed()` failed"),                               PDB_Diagnostic_Codes.not_being_packed ),
        ( P(line="Cannot compute derivative for hbond interaction"),                 PDB_Diagnostic_Codes.no_hbond_deriv ),
        ( P(line="No conformer with given nu angles found in this set"),             PDB_Diagnostic_Codes.nu_conformer ),
        ( P(line=("Residue type","is not appropriately configured.")),               PDB_Diagnostic_Codes.restype_configure ),

        ( P(log=("Patch","implies it can apply to residue type", "but actually applying it fails")), PDB_Diagnostic_Codes.bad_patch ),
    ] )


    for i, line in enumerate(lines):
        previous_line = lines[i-1] if i > 0 else ''
        next_line     = lines[i+1] if i < len(lines) - 1 else ''

        splited = line.replace(':', '').split()

        if 'ERROR' in splited:
            for pattern, error_code in error_map.items():
                if pattern_in(pattern, previous_line, line, next_line, splited, log):
                    return error_code

    #if '*****COMPLETED*****' in previous_line: return None
    # if 'assertion' in log  and  'failed.' in log: return PDB_Diagnostic_Codes.assert_segfault
    if 'exceeded the timeout and will be killed!' in log:
        return PDB_Diagnostic_Codes.exceed_timeout

    return PDB_Diagnostic_Codes.unknown


def main(args):
    ''' When invoked as standalone script run protein_data_bank_diagnostic for pdb's specified in given JSON file '''

    job_file_name = args[1]
    with open(job_file_name) as f: job = Job( **json.load(f) )

    error_logs = os.path.abspath(f'./error-logs')
    if not os.path.isdir(error_logs): os.makedirs(error_logs)

    results = dict(passed=[], failed=dict())
    for pdb, path in job.pdbs.items():
        print(f'Working on {pdb}...')
        res, output = execute('Scoring PDB...', job.command_line.format(input_file=path), return_='tuple')

        if res:
            code = classify_pdb_diagnostic_log(output)
            special_handled_errors = (PDB_Diagnostic_Codes.unknown, None)  # PDB_Diagnostic_Codes.misc_segfault removed from this list
            log_path = os.path.abspath(f'{pdb}.log' if code in special_handled_errors else f'{error_logs}/{pdb}.log')
            with open(log_path, 'w') as f: f.write( f'{pdb.upper()} process exit code: {code}\n\n{output}')

            results['failed'].setdefault(str(code), {})
            results['failed'][str(code)][pdb] = log_path

        else: results['passed'].append(pdb)


    with open(f'{job.name}.output.json', 'w') as f: json.dump(results, f, sort_keys=True, indent=2)


if __name__ == "__main__": main(sys.argv)
