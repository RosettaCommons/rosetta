#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   command.py
## @brief  scientific/command.py
## Python script for running multi-step scientific tests
## @author Sergey Lyskov

import os, json

# A bit of Python magic here, what we trying to say is this: from ../__init__ import *, but init is calculated from file location
import importlib.util, sys
importlib.util.spec_from_file_location(__name__, '/'.join(__file__.split('/')[:-2]) +  '/__init__.py').loader.exec_module(sys.modules[__name__])

_api_version_ = '1.1'

def symlink(source, dest):
    ''' Similar to os.symlink but if dest is alread exisist and if type of *source and *dest does not match or if link points to different location : remove dest first
    '''
    relative_source = os.path.relpath( os.path.abspath(source), os.path.dirname(dest) )

    if os.path.islink(dest):
        if os.readlink(dest) == relative_source: return
        else: os.unlink(dest)
    elif os.path.isdir(dest): shutil.rmtree(dest)
    elif os.path.isfile(dest): os.remove(dest)

    os.symlink(relative_source, dest)


def symlink_tree(source, dest):
    ''' Similar to symlink(...) above but recursivly recreate dir tree at dest and symlink all files
    '''
    source = os.path.abspath(source)
    for dir_name, dirs, files in os.walk(source):
        prefix = dir_name[len(source):] + '/'

        for d in dirs:
            dst = dest + prefix + d
            src = source + prefix + d

            #print(f'source:{source}, dir_name: {dir_name}, d: {d}')
            #print(f'src: {src}, dst: {dst}')

            if os.path.islink(src): symlink(src, dst)
            elif not os.path.isdir(dst): os.makedirs(dst)

        for f in files: symlink(source + prefix + f, dest + prefix + f)


def run_multi_step_test(test, rosetta_dir, working_dir, platform, config, hpc_driver, verbose, variants, python_packages):

    tests_dirs = 'tests sfxn_comparison_tests'.split()
    for tests_dir in tests_dirs:
        test_source_dir = f'{rosetta_dir}/tests/scientific/{tests_dir}/{test}'
        if os.path.isdir(test_source_dir): break
    else: raise BenchmarkError(f'Unknown scripts test: {test}! No dir {test!r} as found in {tests_dirs!r} tests dirs, aborting...')

    symlink_tree(test_source_dir, working_dir)

    # os.mkdir( f'{working_dir}/benchmark' )
    # symlink(rosetta_dir + '/tests/benchmark/__init__.py', working_dir + '/benchmark/__init__.py')
    # symlink(rosetta_dir + '/tests/benchmark/hpc_drivers', working_dir + '/benchmark/hpc_drivers')

    python_environment = local_python_install(platform, config)
    python_virtual_environment = setup_persistent_python_virtual_environment(python_environment, python_packages)

    multi_step_config = dict( config,
                              test = test,
                              rosetta_dir = rosetta_dir,
                              working_dir = working_dir,
                              platform = platform,
                              python_virtual_environment = python_virtual_environment._as_dict,
                              verbose = verbose,
                              variants = variants,
                              debug = 'debug' in variants,
    )

    #print('multi_step_config:', multi_step_config)
    with open(f'{working_dir}/{_multi_step_config_}', 'w') as f: json.dump(multi_step_config, f, sort_keys=True, indent=2)

    scripts = sorted( f for f in os.listdir(working_dir) if f[0].isdigit() and f.endswith('.py') )
    for script in scripts:
        #print(script)
        res, output = execute(f'Running {script}...', f'cd {working_dir} && {python_virtual_environment.activate} && python {script}', return_=tuple, add_message_and_command_line_to_output=True)  # source {ve.activate}

        if res:
            result = { _StateKey_ : _S_script_failed_,  _ResultsKey_ : {},
                       _LogKey_   : f'run_multi_step_test for {test} failed while running {script}...Aborting!\n\n{output}\n'
            }
            break

        if os.path.isfile(f'{working_dir}/{_multi_step_error_}'):
            with open(f'{working_dir}/{_multi_step_error_}') as f: result = json.load(f)
            break

    else:
        with open(f'{working_dir}/{_multi_step_result_}') as f: result = json.load(f)

    #if not config['emulation'] and os.path.isdir(python_virtual_environment_path): shutil.rmtree(python_virtual_environment_path)

    return result


def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    # map from test name to a space-separated-string-of-python-packages to be installed, for example: docking='numpy panda==0.23.4'
    # If package have not-yet-stable-api please make sure to SPECIFY THE EXACT VERSION of package to use so our testing-scripts
    # will not accidentally break when a new version of upstream package got released in the future
    tests = dict(
        _template_               = '',

        relax_cartesian  = 'numpy matplotlib==3.2',
        relax_fast       = 'numpy matplotlib==3.2',
        relax_fast_5iter = 'numpy matplotlib==3.2',

        # commented out because we're no longer running it continuously
        # has been superceded by glycan_dock and renamed to legacy_dock_glycans
#        dock_glycans     = 'numpy matplotlib==3.2',

        stepwise_rna_favorites = 'numpy matplotlib==3.2',
        rna_denovo_favorites   = 'numpy matplotlib==3.2',
        simple_cycpep_predict  = 'numpy matplotlib==3.2',
        peptide_pnear_vs_ic50  = 'numpy scipy matplotlib==3.2',

        design_fast            = 'numpy matplotlib==3.2',
        enzyme_design          = 'numpy matplotlib==3.2',
        fragments_picking      = 'numpy matplotlib==3.2',
        cofactor_binding_sites = 'numpy matplotlib==3.2',

        mp_f19_tilt_angle     = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',
        mp_f19_sequence_recovery    = 'numpy matplotlib==3.2',
        mp_f19_ddG_of_mutation      = 'numpy matplotlib==3.7', # Requires Python 3.8!
        mp_f19_decoy_discrimination = 'numpy matplotlib==3.2',

        mp_dock            = 'numpy matplotlib==3.2',
        mp_relax           = 'numpy matplotlib==3.2',
        mp_symdock         = 'numpy matplotlib==3.2',
        mp_lipid_acc       = 'numpy matplotlib==3.2',
        mp_domain_assembly = 'numpy matplotlib==3.2',

        sewing               = 'numpy matplotlib==3.2',

        antibody_grafting    = 'numpy matplotlib==3.2',
        antibody_h3_modeling = 'numpy matplotlib==3.2',
        antibody_snugdock    = 'numpy matplotlib==3.2 pandas==0.24.2',

        loop_modeling_ngk_12res           = 'numpy matplotlib==3.2',
        loop_modeling_kic_fragments_12res = 'numpy matplotlib==3.2',
        loop_modeling_kic_12res           = 'numpy matplotlib==3.2',
        loop_modeling_ccd_12res           = 'numpy matplotlib==3.2',

        make_fragments         = 'numpy matplotlib==3.2',
        ligand_docking	       = 'numpy matplotlib==3.2',
        ligand_scoring_ranking = 'numpy matplotlib==3.2 pandas==0.23.4 sklearn scipy ',
        mhc_epitope_energy     = 'numpy matplotlib==3.2',

        FlexPepDock            = 'numpy matplotlib==3.2',
        docking                = 'numpy matplotlib==3.2 pandas==0.24.2',
        docking_ensemble       = 'numpy matplotlib==3.2 pandas==0.24.2',
        ddg_ala_scan           = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',
        glycan_structure_prediction  = 'numpy matplotlib==3.2 pandas==0.24.2 seaborn',
        glycan_dock            = 'numpy matplotlib==3.2 pandas==0.24.2',
        RosettaCM              = 'numpy matplotlib==3.2',

        abinitio_RosettaNMR_rdc = 'numpy matplotlib==3.2',
        abinitio_RosettaNMR_pcs = 'numpy matplotlib==3.2',


	# LARGE-SCALE SCOREFUNCTION COMPARISONS============================

        sb_ref2015_relax_cartesian   = 'numpy matplotlib==3.2',
        sb_ref2015_relax_fast        = 'numpy matplotlib==3.2',
        sb_ref2015_relax_fast_5iter  = 'numpy matplotlib==3.2',

        sb_talaris14_relax_cartesian = 'numpy matplotlib==3.2',
        sb_talaris14_relax_fast      = 'numpy matplotlib==3.2',
        sb_talaris14_relax_fast_5iter= 'numpy matplotlib==3.2',

        sb_talaris13_relax_cartesian = 'numpy matplotlib==3.2',
        sb_talaris13_relax_fast      = 'numpy matplotlib==3.2',
        sb_talaris13_relax_fast_5iter= 'numpy matplotlib==3.2',

        sb_score12_relax_cartesian   = 'numpy matplotlib==3.2',
        sb_score12_relax_fast        = 'numpy matplotlib==3.2',
        sb_score12_relax_fast_5iter  = 'numpy matplotlib==3.2',

        sb_ref2015_loop_modeling_ngk_12res = 'numpy matplotlib==3.2',
        sb_ref2015_loop_modeling_kic_fragments_12res = 'numpy matplotlib==3.2',
        sb_ref2015_loop_modeling_kic_12res = 'numpy matplotlib==3.2',
        sb_ref2015_loop_modeling_ccd_12res = 'numpy matplotlib==3.2',

        sb_talaris14_loop_modeling_ngk_12res = 'numpy matplotlib==3.2',
        sb_talaris14_loop_modeling_kic_fragments_12res = 'numpy matplotlib==3.2',
        sb_talaris14_loop_modeling_kic_12res = 'numpy matplotlib==3.2',
        sb_talaris14_loop_modeling_ccd_12res = 'numpy matplotlib==3.2',

        sb_talaris13_loop_modeling_ngk_12res = 'numpy matplotlib==3.2',
        sb_talaris13_loop_modeling_kic_fragments_12res = 'numpy matplotlib==3.2',
        sb_talaris13_loop_modeling_kic_12res = 'numpy matplotlib==3.2',
        sb_talaris13_loop_modeling_ccd_12res = 'numpy matplotlib==3.2',

        sb_score12_loop_modeling_ngk_12res = 'numpy matplotlib==3.2',
        sb_score12_loop_modeling_kic_fragments_12res = 'numpy matplotlib==3.2',
        sb_score12_loop_modeling_kic_12res = 'numpy matplotlib==3.2',
        sb_score12_loop_modeling_ccd_12res = 'numpy matplotlib==3.2',

        sb_ligand_docking	= 'numpy matplotlib==3.2',

        sb_ref2015_docking      = 'numpy matplotlib==3.2 pandas==0.24.2',
        sb_talaris14_docking    = 'numpy matplotlib==3.2 pandas==0.24.2',
        sb_talaris13_docking    = 'numpy matplotlib==3.2 pandas==0.24.2',
        sb_score12_docking      = 'numpy matplotlib==3.2 pandas==0.24.2',

        sb_ref2015_fast_design   = 'numpy matplotlib==3.2',
        sb_talaris14_fast_design = 'numpy matplotlib==3.2',
        sb_talaris13_fast_design = 'numpy matplotlib==3.2',
        sb_score12_fast_design   = 'numpy matplotlib==3.2',

        sb_ref2015_mpddg         = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',
        sb_ref2015mem_mpddg      = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',
        sb_mpframework2012_mpddg = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',
        sb_franklin2019_mpddg    = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',

        sb_ref2015_mpddg_cartesian         = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',
        sb_ref2015mem_mpddg_cartesian      = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',
        sb_mpframework2012_mpddg_cartesian = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',
        sb_franklin2019_mpddg_cartesian    = 'numpy matplotlib==3.2 pandas==0.24.2 scipy==1.1.0',

        self   = '',
    )

    test, _, variants = test.partition('.')
    variants = list( set( variants.split('.') ) )

    if test in tests:
        return run_multi_step_test(
            test = test, python_packages = tests[test],
            rosetta_dir=repository_root, working_dir=working_dir, platform=platform,
            config=config, hpc_driver=hpc_driver, verbose=verbose, variants=variants,
        )

    else: raise BenchmarkError(f'Unknown scripts test: {test}!')
