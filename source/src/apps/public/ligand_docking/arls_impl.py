#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

import sys, os, re
from optparse import OptionParser, IndentedHelpFormatter

try: set
except: from sets import Set as set

# Better handle multiple paragraph descriptions.
class PreformattedDescFormatter (IndentedHelpFormatter):
    def format_description(self, description):
        return description.strip() + "\n" # Remove leading/trailing whitespace


def chain_from_iterables(iterables):
    for it in iterables:
        for el in it:
            yield el


def find_file(regex, paths, required=True):
    if isinstance(regex, str): regex = re.compile(regex)
    for path in paths:
        if not os.path.isdir(path): continue
        hits = [ f for f in os.listdir(path) if regex.match(f) is not None ]
        if len(hits) == 1:
            return os.path.join(path, hits[0])
        elif len(hits) > 1:
            if required: raise ValueError("Multiple hits for '%s' in '%s' (try using the --compile-tag option)" % (regex.pattern, path))
            else: return None
    if required: raise ValueError("'%s' not found" % regex.pattern)
    else: return None


class DockingCase(object):
    # Caching ensures that each file gets prepared only once,
    # and that cases that share files share DockingFile objects.
    # Classes of files are separated b/c e.g. ligands and cofactors are treated differently
    _protein_cache = {}
    _cofactor_cache = {}
    _ligand_cache = {}

    def __init__(self, paths, prot_cofs_lig):
        protein_exts = ['pdb']
        chemical_exts = ['mol', 'mol2', 'sdf']
        self.protein = self._pcache(DockingFile(paths, prot_cofs_lig[0], protein_exts))
        self.cofactors = [ self._ccache(DockingFile(paths, c, chemical_exts)) for c in prot_cofs_lig[1:-1] ]
        self.ligand = self._lcache(DockingFile(paths, prot_cofs_lig[-1], chemical_exts))

    def _pcache(self, docking_file):
        return self._protein_cache.setdefault(docking_file, docking_file)
    def _ccache(self, docking_file):
        return self._cofactor_cache.setdefault(docking_file, docking_file)
    def _lcache(self, docking_file):
        return self._ligand_cache.setdefault(docking_file, docking_file)


class DockingFile(object):
    def __init__(self, paths, base, exts):
        self.base = base
        ext_re = r'\.(%s)$' % '|'.join(exts) # == e.g. r'\.(pdb|mol|mol2|sdf)$'
        self.path = find_file(re.escape(base)+ext_re, paths)

    def __cmp__(self, other):
        return cmp(self.path, other.path)
    def __hash__(self):
        return hash(self.path)

    def __str__(self):
        return self.path
    def __repr__(self):
        return self.path


def find_programs(options):
    PATH_paths = os.environ['PATH'].split(os.pathsep)
    ros_paths = [os.path.join(options.mini,'bin')] + PATH_paths
    if options.compile_tag is not None:
        comp_tag = options.compile_tag.replace('.',r'\.')
    else:
        comp_tag = r'.+release'
    scr_paths = [os.path.join(options.mini,'src','apps','public','ligand_docking'),
        os.path.join(options.mini,'src','python','apps','public')] + PATH_paths
    oe_paths =  [os.path.join(options.openeye,'bin')] + PATH_paths
    if 'OE_LICENSE' in os.environ:
        oe_paths.append( os.path.join(os.path.dirname(os.environ['OE_LICENSE']),'bin') )
    progs = {
        'dock': find_file(r'ligand_dock\.'+comp_tag, ros_paths),
        'rpkmin': find_file(r'ligand_rpkmin\.'+comp_tag, ros_paths),
        'confidence': find_file(r'ligdock_confidence\.'+comp_tag, ros_paths, required=False),
        'extract': find_file(r'extract_atomtree_diffs\.'+comp_tag, ros_paths),
        'best_poses': find_file(r'select_best_unique_ligand_poses\.'+comp_tag, ros_paths, required=False),
        'assign_charges': find_file(r'assign_charges\.py', scr_paths, required=not options.skip_charges),
        'molfile_to_params': find_file(r'molfile_to_params\.py', scr_paths),
        'pdb_to_molfile': find_file(r'pdb_to_molfile\.py', scr_paths),
        'get_scores': find_file(r'get_scores\.py', scr_paths),
        'plot_funnels': find_file(r'plot_funnels\.R', scr_paths),
        'parallel': find_file(r'parallel\.py', scr_paths),
        'omega': find_file(r'omega2', oe_paths, required=not options.skip_omega),
    }
    if options.skip_omega: progs['omega'] = None
    if options.skip_charges: progs['assign_charges'] = None
    return progs


def read_input_list(infile):
    inputs = []
    paths = [ os.getcwd() ]
    for line in infile:
        line = line.rstrip()
        if line.startswith("#") or line == "": continue
        fields = line.split()
        if len(fields) < 2:
            raise ValueError("Not enough fields in list of input structures: "+line.rstrip())
        inputs.append( DockingCase(paths, fields) )
    return inputs


def get_procs_per_case(docking_cases, options):
    return max(1, int(options.njobs) // len(docking_cases))


def write_setup_script(outfile, options, progs, docking_cases): #{{{
    proteins = set([ d.protein for d in docking_cases ])
    cofactors = set( chain_from_iterables([d.cofactors for d in docking_cases]) )
    ligands = set([ d.ligand for d in docking_cases ])

    outfile.write('''#!/bin/bash
set -v

# Set up directory structure
mkdir -p ligand cofactor unbound input native

# Make ligands
pushd ligand
mkdir -p {fa,cen}/{conf1,confs,kins,withxtal}

''')
    if progs['omega'] is not None:
        outfile.write('omega="%s -includeInput -commentEnergy"\n' % progs['omega'])
    if options.num_procs > 1: outfile.write("%s -j%i <<HEREDOC\n" % (progs['parallel'], options.num_procs))
    for i, ligand in sampling(ligands):
        ligcode = "X%02X" % i # X01, X02, ..., XFF
        ligbase = ligand.base
        ligpath = ligand.path
        if progs['omega'] is not None:
            oldpath = ligpath
            ligpath = "%s.omega.mol2" % ligbase
            outfile.write("$omega -in %s -out %s -prefix _%s && " % (oldpath, ligpath, ligbase))
        if progs['assign_charges'] is not None:
            oldpath = ligpath
            ligpath = "%s.am1bcc.mol2" % ligbase
            outfile.write("%s %s < %s > %s && " % (options.python, progs['assign_charges'], oldpath, ligpath))
        outfile.write("%s %s -c -n%s -p%s -k%s.kin %s" % (options.python, progs['molfile_to_params'], ligcode, ligbase, ligbase, ligpath))
        for ext in ['fa', 'cen']:
            outfile.write(" && cat %s_????.%s.pdb | gzip -c > %s/withxtal/%s_confs.%s.pdb.gz" % (ligbase, ext, ext, ligbase, ext))
            # Force a _0002 conformer to exist, even if it's just a duplicate of _0001
            outfile.write(" && ( [ -f %s_0002.%s.pdb ] || cp %s_0001.%s.pdb %s_0002.%s.pdb )" % (ligbase, ext, ligbase, ext, ligbase, ext))
            outfile.write(" && mv %s_0001.%s.pdb %s/conf1/" % (ligbase, ext, ext))
            outfile.write(" && cat %s_????.%s.pdb | gzip -c > %s/%s_confs.%s.pdb.gz" % (ligbase, ext, ext, ligbase, ext))
            outfile.write(" && mv %s_????.%s.pdb %s/confs/" % (ligbase, ext, ext))
            outfile.write(" && echo 'PDB_ROTAMERS %s_confs.%s.pdb' >> %s.%s.params" % (ligbase, ext, ligbase, ext))
            outfile.write(" && cp %s.%s.params %s/withxtal/" % (ligbase, ext, ext))
            outfile.write(" && mv %s.%s.params %s/" % (ligbase, ext, ext))
            outfile.write(" && mv %s.%s.kin %s/kins/" % (ligbase, ext, ext))
        outfile.write("\n")
    if options.num_procs > 1: outfile.write("HEREDOC\n")
    outfile.write('''
popd

# Make cofactors
pushd cofactor
mkdir -p {fa,cen}/{conf1,confs,kins,withxtal}

''')
    if options.num_procs > 1: outfile.write("%s -j%i <<HEREDOC\n" % (progs['parallel'], options.num_procs))
    for i, ligand in sampling(cofactors):
        ligcode = "Q%02X" % i # Q01, Q02, ..., QFF
        ligbase = ligand.base
        ligpath = ligand.path
        # Generate conformers for cofactors, but by default only use xtal conf
        if progs['omega'] is not None:
            oldpath = ligpath
            ligpath = "%s.omega.mol2" % ligbase
            outfile.write("$omega -in %s -out %s -prefix _%s && " % (oldpath, ligpath, ligbase))
        if progs['assign_charges'] is not None:
            oldpath = ligpath
            ligpath = "%s.am1bcc.mol2" % ligbase
            outfile.write("%s %s < %s > %s && " % (options.python, progs['assign_charges'], oldpath, ligpath))
        outfile.write("%s %s -c -n%s -p%s -k%s.kin %s" % (options.python, progs['molfile_to_params'], ligcode, ligbase, ligbase, ligpath))
        for ext in ['fa', 'cen']:
            outfile.write(" && cat %s_????.%s.pdb | gzip -c > %s/withxtal/%s_confs.%s.pdb.gz" % (ligbase, ext, ext, ligbase, ext))
            # Force a _0002 conformer to exist, even if it's just a duplicate of _0001
            outfile.write(" && ( [ -f %s_0002.%s.pdb ] || cp %s_0001.%s.pdb %s_0002.%s.pdb )" % (ligbase, ext, ligbase, ext, ligbase, ext))
            outfile.write(" && mv %s_0001.%s.pdb %s/conf1/" % (ligbase, ext, ext))
            #outfile.write(" && cat %s_????.%s.pdb | gzip -c > %s/%s_confs.%s.pdb.gz" % (ligbase, ext, ext, ligbase, ext))
            outfile.write(" && mv %s_????.%s.pdb %s/confs/" % (ligbase, ext, ext))
            # Default is to generate the conformer libraries "just in case" but not enable them...
            #outfile.write(" && echo 'PDB_ROTAMERS %s_confs.%s.pdb' >> %s.%s.params" % (ligbase, ext, ligbase, ext))
            outfile.write(" && mv %s.%s.params %s/" % (ligbase, ext, ext))
            outfile.write(" && mv %s.%s.kin %s/kins/" % (ligbase, ext, ext))
        outfile.write("\n")
    if options.num_procs > 1: outfile.write("HEREDOC\n")
    outfile.write('''
popd

# Make "unbound" proteins (never prepacked, never a ligand or cofactor)
''')
    for protein in proteins:
        outfile.write("cp %s unbound/%s.pdb\n" % (protein.path, protein.base))
    outfile.write('''
# Make native and input PDBs, with ligands and cofactors
# May later be replaced by repacked PDB files
''')
    for docking_case in docking_cases:
        cofs = " ".join(["cofactor/fa/conf1/%s_0001.fa.pdb" % c.base for c in docking_case.cofactors])
        nat_files = " ".join([docking_case.protein.path, cofs, "ligand/fa/conf1/%s_0001.fa.pdb" % docking_case.ligand.base])
        inp_files = " ".join([docking_case.protein.path, cofs, "ligand/fa/confs/%s_0002.fa.pdb" % docking_case.ligand.base])
        outfile.write("cat %s > native/%s_%s.pdb\n" % (nat_files, docking_case.protein.base, docking_case.ligand.base))
        outfile.write("cat %s > input/%s_%s.pdb\n" % (inp_files, docking_case.protein.base, docking_case.ligand.base))
#}}}


def write_common_flags(outfile, options, docking_cases): #{{{
    database = options.database
    cofactors = set( chain_from_iterables([d.cofactors for d in docking_cases]) )
    ligands = set([ d.ligand for d in docking_cases ])
    extra_res_cen = " ".join(
        ["ligand/cen/%s.cen.params" % l.base for l in ligands]
        + ["cofactor/cen/%s.cen.params" % c.base for c in cofactors])
    extra_res_fa = " ".join(
        ["ligand/fa/%s.fa.params" % l.base for l in ligands]
        + ["cofactor/fa/%s.fa.params" % c.base for c in cofactors])
    outfile.write('''
-in
 -path
  ## On a large cluster, the database should be on a local scratch disk
  ## to avoid over-taxing NFS.
  #-database /scratch/ROSETTA/minirosetta_database
  ## "Fallback" database locations can also be specified,
  ## in case the primary database is missing on some nodes:
  #-database /work/davis/rosetta/rosetta_database
  ## Location provided by user:
  -database %(database)s
 -file
  ## You must supply .params files for any residue types (ligands)
  ## that are not present in the standard Rosetta database.
  ## Centroid residues are not needed for docking,
  ## but may be needed for other Rosetta protocols.
  #-extra_res_cen %(extra_res_cen)s
  -extra_res_fa %(extra_res_fa)s
-out
 ## These channels generate a LOT of output,
 ## especially with -improve_orientation.
 -mute protocols.geometry.RB_geometry core.scoring.dunbrack.SingleLigandRotamerLibrary core.scoring.rms_util core.pack.dunbrack.SingleLigandRotamerLibrary
 -file
  ## I prefer output structures with Rosetta numbering, from 1 to N residues.
  ## To keep the original PDB numbering, omit this flag:
  -renumber_pdb
-run
 ## Recording the SVN revision of the code in your output files
 ## makes it easier to figure out what you did later.
 -version
 ## MT is now the default random number generator, actually.
 #-rng mt19937
 ## Rosetta's default behavior is to pretend that atoms with zero occupancy
 ## don't exist at all, which can lead to whole residues disappearing...
 -ignore_zero_occupancy false
-packing
 ## Includes the input sidechain conformations as rotamers when repacking,
 ## but they can be "lost" if any other rotamer is chosen at that position.
 #-use_input_sc
 ## Instead, use -unboundrot to permanently add the rotamers from one or more
 ## PDB files.  Their rotamer ("Dunbrack") energies are also adjusted to be
 ## equal to the best rotamer at that position, thereby favoring the native.
 ## Since most sidechains do not change position much upon ligand binding,
 ## including this knowledge generally improves docking results:
 #-unboundrot ...
 ## Controls the number of (protein) rotamers used.
 -ex1
 -ex1aro
 -ex2
 ## Ensures that extra rotamers are used for surface residues too.
 -extrachi_cutoff 1
-docking
 -ligand
  ## Use soft-repulsive scoring during search (but not final minimization).
  ## This slightly improves search.  Hard-rep is used for final scoring,
  ## however, because it gives better discrimination.
  -soft_rep
  ## Like Rosetta++, only evaluate the Coloumb term between protein and ligand.
  -old_estat
''' % locals())
#}}}


def write_rpkmin_flags(outfile, options, docking_cases): #{{{
    outfile.write('''
-in
 -file
  ## If you weren't supplying the input file(s) on the command line,
  ## this is where you would put them:
  #-s ...
-out
 ## Number of trajectories to run (per input structure given with -s)
 -nstruct 10
 -path
  ## Where to write the PDB files.
  -pdb prepack/all_traj
-run
 ## Do nothing for 0 - N seconds on startup.  If many processes are started
 ## at once, this helps avoid too much I/O as they all load the database.
 #-random_delay 30
-packing
 ## If your PDB file does not have hydrogens in the right places, use this:
 -no_optH false
 ## You may also want to fix Asn/Gln/His sidechains:
 -flip_HNQ
''')
#}}}


def write_dock_flags(outfile, options, docking_cases): #{{{
    nstruct = 5000 // get_procs_per_case(docking_cases, options)
    outfile.write('''
-in
 -file
  ## If you weren't supplying the input file(s) on the command line,
  ## this is where you would put them:
  #-s ...
  ## To get meaningful RMS values in the output, either the input PDB must have
  ## the ligand in the correct place, or you must supply another file that does:
  #-native ...
-out
 ## Number of trajectories to run (per input structure given with -s)
 -nstruct %(nstruct)i
 ## With multiple processes, ensures a unique name for every output structure.
 ## Each process should get a different string here, so you can't really
 ## put it in the FLAGS.txt file, it has to go in the Condor script:
 #-suffix 3
 -path
  ## Where to write the silent.out file.  Specified in my Condor script.
  #-pdb work/jnk_pp_1_001/3
-run
 ## Do nothing for 0 - N seconds on startup.  If many processes are started
 ## at once, this helps avoid too much I/O as they all load the database.
 -random_delay 240
-packing
 ## If your PDB file does not have hydrogens in the right places, use this:
 #-no_optH false
 ## You may also want to fix Asn/Gln/His sidechains:
 #-flip_HNQ
-docking
 ## Flags to control the initial perturbation of the ligand,
 ## and thus how much of the binding pocket to explore.
 ##
 ## Random perturbation of up to N Angstroms in X, Y, and Z,
 ## drawn from a uniform distribution.
 ## This gives a cube, but points outside the sphere of radius N are discarded,
 ## resulting in uniform positional sampling within the sphere.
 -uniform_trans 5   # angstroms
 ## Randomize the starting orientation and conformer (respectively) of the
 ## ligand.  Unnecessary if using -improve_orientation.
 #-randomize2
 #-random_conformer
 ## An alternative to -uniform_trans and -randomize2.
 ## Random perturbations in X, Y, and Z and the 3 Euler angles,
 ## drawn from a Gaussian distribution.
 ## I *believe* that the three Gaussians in X, Y, and Z actually give a
 ## spherically isotropic distribution of positions, but large angles
 ## clearly give unreasonable results (use -randomize2 instead).
 ## This can be used for positional perturbation only (instead of
 ## -uniform_trans) by setting the angle component to zero.
 #-dock_pert 30 5  # rot degrees, trans angstroms
 -ligand
  ## Instead of choosing an orientation and conformer at random,  try all
  ## available conformers in N random orientations, and try to maximize
  ## shape complementarity to the protein (keeping in mind that sidechains
  ## may move later, etc etc.  See the paper or code for details.)
  ## We find this *significantly* improves search, and costs just a few seconds.
  -improve_orientation 1000
  ## Let ligand rotatable bonds minimize, with harmonic restraints
  ## where one standard deviation is 15 degrees:
  -minimize_ligand
  -harmonic_torsions 15
  ## Allows rotamers to "recombine" instead of always minimizing towards
  ## the last rotamer selected by the packer.
  -use_ambig_constraints
  ## Introduces random torsion tweaks in the MCM cycle; may help large ligands.
  ## Requires -use_ambig_constraints.
  -shear_moves 5
  ## In the final minimization, let the backbone minimize with harmonic
  ## restraints on the Calphas (stddev = 0.3 A).
  ## Only stretches of residues near the ligand are minimized,
  ## typically 20 - 40 in total.  (40 - 80 residues are repacked.)
  -minimize_backbone
  -harmonic_Calphas 0.3
  ## The 6-cycle protocol with repacks in cycles 1 and 4.
  -protocol abbrev2
  ## In each trajectory, one of the following points will be chosen at random
  ## as the starting point for the ligand centroid (followed by other
  ## perturbations like -uniform_trans, etc).
  #-start_from 20.0 28.4 26.0
  #-start_from 23.3 18.9 26.0
  #-start_from 28.7 14.0 19.6
  #-start_from 16.9 16.0 32.7
  #-start_from 21.3 10.8 31.1
  #-start_from 25.6  4.4 32.5
  ## By supplying multiple .params files that have the same 3-letter residue
  ## name, the same heavy atoms, but different proton configurations, you can
  ## sample different protomers / tautomers within the packing steps.
  ## In many cases this gives a negligible increase in run time (nice).
  ## However, I'm unaware of any cases (so far) where it actually helped
  ## to improve docking results, so use at your own risk.
  #-mutate_same_name3
''' % locals())
#}}}


def write_minnat_flags(outfile, options, docking_cases): #{{{
    outfile.write('''
## Flags for generating "minimized natives", if applicable.
## In docking benchmarks, these can be used to diagnose search problems vs scoring problems.
-out
 -nstruct 10
-run
 -random_delay 240
-packing
 -use_input_sc   # causes -use_ambig_constraints to include native rotamer!
-docking
 -ligand
  -minimize_ligand
  -harmonic_torsions 15
  -use_ambig_constraints
  -shear_moves 5
  -minimize_backbone
  -harmonic_Calphas 0.3
  -protocol abbrev2
''' % locals())
#}}}


def write_rpkmin_script(outfile, options, progs, docking_cases): #{{{
    # All unique protein/cofactor combinations
    proteins_and_cofactors = set([ tuple([d.protein]+d.cofactors) for d in docking_cases ])
    outfile.write('''#!/bin/bash
set -v

# Set up directory structure
mkdir -p prepack/{input,all_traj}

# Prepack and minimize protein receptors

''')
    if options.num_procs > 1: outfile.write("%s -j%i <<'HEREDOC'\n" % (progs['parallel'], options.num_procs))
    for p_and_c in proteins_and_cofactors:
        protein = p_and_c[0]
        cofactors = p_and_c[1:]
        base = "_".join([pc.base for pc in p_and_c])
        cof_files = " ".join(["cofactor/fa/conf1/%s_0001.fa.pdb" % c.base for c in cofactors])
        outfile.write("cat %s %s > prepack/input/%s.pdb" % (protein.path, cof_files, base))
        outfile.write(" && %s @COMMON.flags @RPKMIN.flags -s prepack/input/%s.pdb" % (progs['rpkmin'], base))
        outfile.write(" && egrep '^(ATOM|HETATM)' $(awk '/^pose / {print FILENAME, $NF}' prepack/all_traj/%s_????.pdb | sort -nk2 | head -1 | cut -d' ' -f1) > prepack/%s.pdb" % (base, base))
        outfile.write("\n")
    if options.num_procs > 1: outfile.write("HEREDOC\n")
    outfile.write('''
# Replace by old input PDB files with repacked ones
''')
    for d in docking_cases:
        base = "_".join([pc.base for pc in [d.protein]+d.cofactors])
        outfile.write("cat prepack/%s.pdb ligand/fa/confs/%s_0002.fa.pdb > input/%s_%s.pdb\n" % (base, d.ligand.base, d.protein.base, d.ligand.base))
#}}}


def write_tarball_pre_script(outfile, options, progs, docking_cases): #{{{
    if not options.local :
        outfile.write('''#!/bin/bash
#set -v  # this messes up the nice echo output at the end
shopt -s nullglob  # avoid error messages from tar about non-existant files

# Make tarball of files needed for docking
tar cvzf $(basename $(pwd)).tgz {cofactor,ligand}/{fa,cen}/{withxtal/,}*.params {cofactor,ligand}/{fa,cen}/{withxtal/,}*_confs.*.pdb* input/ native/ unbound/ *.flags [45]_*.sh

# Suggest copying
echo
echo "  Now copy this file to the cluster where you'll do docking."
echo
echo "    rsync -avxP $(basename $(pwd)).tgz syd:runs/$(basename $(pwd))/"
echo
echo "  *** If necessary, adjust path to database in COMMON.flags"
echo "  *** If necessary, adjust path to executable in 4_dock_condor.sh"
echo
''')
    else: # options.local = True
        outfile.write('''#!/bin/bash
echo
echo "  *** If necessary, adjust path to database in COMMON.flags"
echo "  *** If necessary, adjust path to executable in 4_dock_condor.sh"
echo
''')
#}}}


def write_dock_condor_script(outfile, options, progs, docking_cases): #{{{
    ligand_dock = progs['dock']
    #target_dir = os.path.join(os.getcwd(), options.target_dir)
    outfile.write('''#!/bin/bash

mkdir -p log

# Create the header for the CONDOR files
minnatfile="CONDOR.minnat"
dockfile="CONDOR.dock"

for f in "$minnatfile" "$dockfile"; do
    cat > "$f" <<'HEREDOC'
Executable = %(ligand_dock)s
Universe = vanilla
Requirements = Memory > 256
Notification = Never
copy_to_spool = false

Initialdir = .
HEREDOC
done

cat >> "$minnatfile" <<'HEREDOC'
Error =  log/err.minnat.$(Process).log
Log = log/condor.minnat.$(Process).log
Output = log/out.minnat.$(Process).log


HEREDOC

cat >> "$dockfile" <<'HEREDOC'
Error =  log/err.$(Process).log
Log = log/condor.$(Process).log
Output = log/out.$(Process).log


HEREDOC

# Function to set up each docking case
# Takes protein_name ligand_name num_procs
setup_docking_case()
{
    p="$1"
    l="$2"
    pl="${1}_${2}"
    mkdir -p "work/minnat/$pl"
    echo "Arguments = @COMMON.flags @MINNAT.flags -in:file:s native/$pl.pdb -in:file:native native/$pl.pdb -packing:unboundrot unbound/$p.pdb -out:path:pdb work/minnat/$pl"  >> "$minnatfile"
    echo "Queue" >> "$minnatfile"
    for((i=0;i<$3;i++)); do
        mkdir -p "work/$pl/$i"
        echo "Arguments = @COMMON.flags @DOCK.flags -in:file:s input/$pl.pdb -in:file:native native/$pl.pdb -packing:unboundrot unbound/$p.pdb -out:path:pdb work/$pl/$i -run:seed_offset $i -out:suffix _$i" >> "$dockfile"
        echo "Queue" >> "$dockfile"
    done
    echo >> "$dockfile"
}

# List the individual docking cases
''' % locals())
    procs_per_case = get_procs_per_case(docking_cases, options)
    outfile.write('procs_per_case=%i\n' % procs_per_case)
    for d in docking_cases:
        prot_lig = "%s_%s" % (d.protein.base, d.ligand.base)
        outfile.write('setup_docking_case "%s" "%s" $procs_per_case\n' % (d.protein.base, d.ligand.base))
    outfile.write('''
echo
echo "  condor_submit $minnatfile    # if you want minimized natives"
echo "  condor_submit $dockfile      # the actual docking calculations"
echo
''')
#}}}

def write_dock_bash_script(outfile, options, progs, docking_cases): #{{{
    ligand_dock = progs['dock']
    #target_dir = os.path.join(os.getcwd(), options.target_dir)
    outfile.write('''#!/bin/bash

Executable="%(ligand_dock)s"

mkdir -p log

# Create the header for the BASH files
minnatfile="BASH.minnat.sh"
dockfile="BASH.dock.sh"

# Make sure that neither of them already exist.
rm -f "$minnatfile"
rm -f "$dockfile"

# Function to set up each docking case
# Takes protein_name ligand_name num_procs
setup_docking_case()
{
    p="$1"
    l="$2"
    pl="${1}_${2}"
    mkdir -p "work/minnat/$pl"
    echo "${Executable} @COMMON.flags @MINNAT.flags -in:file:s native/$pl.pdb -in:file:native native/$pl.pdb -packing:unboundrot unbound/$p.pdb -out:path:pdb work/minnat/$pl > log/err.minnat.${pl}.log 2> log/out.minnat.${pl}.log"  >> "$minnatfile"
    echo >> "$minnatfile"
    for((i=0;i<$3;i++)); do
        mkdir -p "work/$pl/$i"
        echo "${Executable} @COMMON.flags @DOCK.flags -in:file:s input/$pl.pdb -in:file:native native/$pl.pdb -packing:unboundrot unbound/$p.pdb -out:path:pdb work/$pl/$i -run:seed_offset $i -out:suffix _$i > log/out.${pl}.${i}.log 2> log/err.${pl}.${i}.log" >> "$dockfile"
    done
    echo >> "$dockfile"
}

# List the individual docking cases
''' % locals())
    procs_per_case = get_procs_per_case(docking_cases, options)
    outfile.write('procs_per_case=%i\n' % procs_per_case)
    for d in docking_cases:
        prot_lig = "%s_%s" % (d.protein.base, d.ligand.base)
        outfile.write('setup_docking_case "%s" "%s" $procs_per_case\n' % (d.protein.base, d.ligand.base))
    outfile.write('''

# Make scripts executable
chmod u+x $minnatfile $dockfile

echo
echo "  nohup nice $minnatfile    # if you want minimized natives"
echo "  nohup nice $dockfile      # the actual docking calculations"
echo
echo " Or for multiple processors something like:"
echo
echo "  %s -j %s < $minnatfile  # if you want minimized natives"
echo "  %s -j %s < $dockfile    # the actual docking calculations"
''' % (progs['parallel'], options.num_procs, progs['parallel'], options.num_procs) )
#}}}


def write_concat_script(outfile, options, progs, docking_cases): #{{{
    outfile.write('''#!/bin/bash
mkdir -p out/minnat
for dir in work/minnat/*; do
    d=$(basename $dir)
    echo $d
    f="out/minnat/${d}_silent.out"
    mv "$dir"/silent*.out $f
    # Remove new file if it is contains no data (0 size)
    [ -s $f ] || rm $f
done

for dir in work/*; do
    d=$(basename $dir)
    [ "$d" = minnat ] && continue
    echo $d
    f="out/${d}_silent.out"
    # Remove old file first
    [ -f $f ] && rm $f
    find $dir -name 'silent*.out' -print0 | xargs -0 cat >> $f
    # Remove new file if it is contains no data (0 size)
    [ -s $f ] || rm $f
done

for f in out/*_silent.out; do
    echo -n "$f "
    a=$(grep -c '^POSE_TAG ' $f)
    b=$(grep -c '^END_POSE_TAG ' $f)
    if ((a == b)); then
        echo "ok $a $b"
    else
        echo "***ERROR*** $a $b"
    fi
done
''')
    if not options.local:
        outfile.write('''tar cvzf $(basename $(pwd))_out.tgz out/

# Suggest copying
echo
echo "  Now copy this file back to your computer for analysis."
echo
echo "    rsync -avxP $(basename $(pwd))_out.tgz dig:"
echo "    ssh dig"
echo "    cd path/to/project/"
echo "    tar xvzf ~/$(basename $(pwd))_out.tgz && rm ~/$(basename $(pwd))_out.tgz"
echo
''')
#}}}


def write_analysis_script(outfile, options, progs, docking_cases): #{{{
    outfile.write('''#!/bin/bash

# Check for result files
if [ ! -d out ]; then
    echo
    echo "  Can't find out/ directory.  Have you untarred the results?"
    echo
    exit
fi

# Make score files
mkdir -p {out,scores}/minnat
pushd out
pushd minnat
for f in *_silent.out; do
    echo $f
    %(get_scores)s < $f > ../../scores/minnat/${f/_silent.out/.tab}
done
popd
for f in *_silent.out; do
    echo $f
    %(get_scores)s < $f > ../scores/${f/_silent.out/.tab}
done
popd

# Make PDF plots if R is on the PATH
pushd scores
if which R > /dev/null; then
    R --vanilla < %(plot_funnels)s
else
    echo
    echo "  Can't find R on the PATH.  Skipping funnel plots..."
    echo
fi
popd
''' % progs)
    if progs['confidence'] is not None:
        outfile.write('''
# Estimate confidence in docking results
%(confidence)s @COMMON.flags -in:file:silent out/*.out
''' % progs)
    if progs['best_poses'] is not None:
        outfile.write('''
# Extract best-scoring docked poses, skipping duplicates
# Wow, the select_best_unique_ligand_poses program is awkward.  I should fix that.
mkdir -p pdbs
for f in out/*.out; do
    %(best_poses)s @COMMON.flags -docking:ligand:max_poses 10 -docking:ligand:min_rms 1.0 -out:path:pdb pdbs -in:file:silent $f -out:file:silent default.out
    tags=( $(head -1 cluster_rms.tab) )   # Bash array of selected tags
    %(extract)s @COMMON.flags -out:path:pdb pdbs -in:file:s $f -in:file:tags "${tags[@]}"
    # Prefix each extracted structure filename with its rank: 01, 02, 03, ...
    for ((i=0; i<${#tags[@]}; ++i)); do
        mv pdbs/${tags[$i]}.pdb pdbs/$(printf "%%02i_" $((i+1)))${tags[$i]}.pdb
    done
    # select_best_unique_ligand_poses seems to mangle the silentfile sometimes, so better to use the original
    #%(extract)s @COMMON.flags -out:path:pdb pdbs -in:file:s default.out
    rm default.out cluster_rms.tab
done

# Convert extracted PDBs into mol2 files, for people who prefer that
mkdir -p mol2_out
pdb_to_molfile="%(pdb_to_molfile)s --comment -oMOL2"
''' % progs)
        for d in docking_cases:
            prot_lig = "%s_%s" % (d.protein.base, d.ligand.base)
            if progs['assign_charges'] is not None: lig_path = "ligand/%s.am1bcc.mol2" % d.ligand.base
            elif progs['omega'] is not None: lig_path = "ligand/%s.omega.mol2" % d.ligand.base
            else: lig_path = d.ligand.path
            outfile.write('$pdb_to_molfile %s pdbs/%s*.pdb > mol2_out/%s.mol2\n' % (lig_path, prot_lig, prot_lig))
#}}}


def write_if_not_exists(filename, options, function, *func_args, **func_kwargs):
    filename = os.path.join(options.target_dir, filename)
    if os.path.exists(filename) and not options.clobber:
        print "File %s exists, use --clobber to overwrite" % filename
    else:
        try:
            outfile = open(filename, 'w')
            function(outfile, *func_args, **func_kwargs)
            outfile.close()
            if filename.endswith('.sh'): os.chmod(filename, 0750)
            print "Wrote %s" % filename
        except Exception, ex:
            print "Exception while writing %s: %s" % (filename, ex)
            # Remove partially complete junk file
            os.remove(filename)
            # Re-raise the exception
            raise # "raise ex" clobbers the original stack trace


def main(argv):
    '''
    Automatic RosettaLigand Setup (ARLS)

    Given a list of protein-ligand pairings, generates all the required input
    files for RosettaLigand docking, along with suggested flags and submit scripts.
    The scripts are numbered 1 - N and should be run in order.
    See the Doxygen documentation for important information on using RosettaLigand.

    The script takes one argument, a list of protein-ligand pairs to dock.
    Each line of the list should consist of the protein name, zero or more cofactor
    names, and the ligand name (whitespace separated).
    In the current directory should be files with the same name and an appropriate
    extension: .pdb for proteins, and one of (.mol, .sdf, .mol2) for ligands and
    cofactors. For example, the line "1abc cofactor mylig" could use files named
    "1abc.pdb", "cofactor.mol", and "mylig.mol2". No unnecessary work will be done
    even if a protein/cofactor/ligand appears in multiple lines.

    The protein file should only contain standard amino acid residues.
    The cofactor and ligand files should contain a single conformation
       of a single compound (unless using --skip-omega).
    The MiniRosetta database should be in ~/rosetta/rosetta_database;
        otherwise use --database.
    OpenEye's Omega should be in ~/openeye or on your PATH;
        otherwise use --openeye or --skip-omega.
    The OpenEye QUACPAC toolkit should be on your PYTHONPATH
        for charges to be assigned (else use --skip-charges).
    You can specify which clustering software to set up the docking runs for with
        --cluster. Specifying BASH will result in a standard shell script-type submission.
        for those without clustering software (It is recommended to also set --njobs.) You
        also can use the BASH output file with MOAB/TORQUE/PBS style clustering software.
    '''  # Preformatted

    def up_dir(path, num_steps=1):
        for i in xrange(num_steps):
            path = os.path.dirname(path)
        return path
    MINI_HOME = up_dir(os.path.abspath(sys.path[0]), 4)

    parser = OptionParser(usage="usage: %prog [<] LIST_OF_PAIRS.txt", formatter=PreformattedDescFormatter())
    parser.set_description(main.__doc__)
    # parser.add_option("-short", ["--long"],
    #   action="store|store_true|store_false",
    #   default=True|False|...
    #   type="string|int|float",
    #   dest="opt_name",
    #   help="store value in PLACE",
    #   metavar="PLACE",
    # )
    parser.add_option("-j", "--num-procs",
      default=1,
      type="int",
      help="Number of processors to use on local machine for non-CONDOR tasks (default: 1)",
    )
    parser.add_option("-m", "--mini",
      default=MINI_HOME,
      help="Directory where rosetta_source is found (default: %s)" % MINI_HOME,
    )
    parser.add_option("-d", "--database",
      default=os.path.join( os.path.expanduser("~"), "rosetta", "rosetta_database"),
      help="Directory where the Rosetta database is found (default: ~/rosetta/rosetta_database)",
    )
    parser.add_option("--openeye",
      default=os.path.join( os.path.expanduser("~"), "openeye"),
      help="Directory where OpenEye tools are found (default: ~/openeye)",
    )
    parser.add_option("-t", "--target-dir",
      help="Directory where output should be created (default: ./arls_work)",
      default="arls_work",
    )
    parser.add_option("--python",
      help="Python executable (2.5/2.6 required for Omega) ",
      default="/usr/bin/env python",
    )
    parser.add_option("--clobber",
      action="store_true",
      default=False,
      help="Overwrite pre-existing output files",
    )
    parser.add_option("--skip-omega",
      action="store_true",
      default=False,
      help="Do not use OpenEye's Omega to generate small molecule conformers.  Assumes small molecule file already contains all desired conformers.",
    )
    parser.add_option("--skip-charges",
      action="store_true",
      default=False,
      help="Do not assign partial charges with OpenEye's AM1BCC.  Assumes small molecule file is in .mol2 format and already has charges, or that you want default Rosetta charges.",
    )
    parser.add_option("--local",
      action="store_true",
      default=False,
      help="Run docking jobs locally (through current system), rather than on a remote cluster.",
    )
    parser.add_option("--cluster",
      default="CONDOR",
      help="Type of clustering software to use. Acceptable values: CONDOR, BASH",
    )
    parser.add_option("--njobs",
      default=250,
      help="Aproximate number of cluster jobs to split processing over (default: 250).",
    )
    parser.add_option("--compile-tag",
      default=None,
      help="""The compile tag (e.g. "static.linuxgccrelease") for Rosetta programs (default: use a very basic autodetect).""",
    )
    (options, args) = parser.parse_args(args=argv)

    if len(args) == 0:
        infile = sys.stdin
    elif len(args) == 1:
        infile = open(args[0])
    else:
        parser.print_help()
        print "Too many arguments!"
        return 1

    try:
        progs = find_programs(options)
        docking_cases = read_input_list(infile)

        if len(docking_cases) == 0:
            raise ValueError("No protein-ligand pairs were supplied!")
        if len(set([(d.protein, d.ligand) for d in docking_cases])) < len(docking_cases):
            raise ValueError("Some protein-ligand pairs are duplicated!")

        options.cluster = str(options.cluster).upper()
        if options.cluster not in ["CONDOR","BASH"]:
            raise ValueError('Invalid cluster type "%s"specified with --cluster'% options.cluster)
        if not os.path.isdir(options.mini):
            raise ValueError("Cannot find MiniRosetta source tree;  please specify with --mini")
        if not os.path.isdir(options.database):
            raise ValueError("Cannot find Rosetta database;  please specify with --database")
        if not options.skip_omega and not os.path.isdir(options.openeye):
            raise ValueError("Cannot find OpenEye tools;  please specify with --openeye (or use --skip-omega)")

        if not os.path.isdir(options.target_dir):
            os.mkdir(options.target_dir)

        write_if_not_exists("1_setup.sh", options,
            write_setup_script, options, progs, docking_cases)
        write_if_not_exists("COMMON.flags", options,
            write_common_flags, options, docking_cases)
        write_if_not_exists("RPKMIN.flags", options,
            write_rpkmin_flags, options, docking_cases)
        write_if_not_exists("DOCK.flags", options,
            write_dock_flags, options, docking_cases)
        write_if_not_exists("MINNAT.flags", options,
            write_minnat_flags, options, docking_cases)
        write_if_not_exists("2_prepack_minimize.sh", options,
            write_rpkmin_script, options, progs, docking_cases)
        write_if_not_exists("3_tarball_pre.sh", options,
            write_tarball_pre_script, options, progs, docking_cases)
        
        if options.cluster == "CONDOR":
            write_if_not_exists("4_dock_condor.sh", options,
                write_dock_condor_script, options, progs, docking_cases)
        elif options.cluster == "BASH":
            write_if_not_exists("4_dock_bash.sh", options,
                write_dock_bash_script, options, progs, docking_cases)
        else:
            raise ValueError("Shouldn't get here: invalid cluster type %s"%options.cluster)
        
        write_if_not_exists("5_concat.sh", options,
            write_concat_script, options, progs, docking_cases)
        write_if_not_exists("6_analyze_results.sh", options,
            write_analysis_script, options, progs, docking_cases)

    except ValueError, v:
        parser.print_help()
        print
        print v
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
    print "You should really be running the wrapper script, arls.py, instead!"
