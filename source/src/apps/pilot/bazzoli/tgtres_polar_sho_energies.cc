// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @brief prints the SHO energies of the polar atoms of selected residues.
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)
///
/// ARGUMENTS
///  ---> -s <POSE_PDB>
///  ---> -tgt_set <TGT_SET>, where <TGT_SET> is the path to a file describing the set of selected
///                           residues. The format of the file is specified in the comments to fuction
///                           load_set() (see below).
///
/// @details The output can be seen as a sequence of N blocks, where N is the number of residues in
///  <TGT_SET> that contain at least one polar heavy atom. The ith block contains atom
///  ids and sho energies of the ith such residue (i=1,...,N). Within each block, the jth line
///  contains the atom id and sho energy of the residue's jth polar heavy atom (j=0,...,H-1).
///

#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/types.hh>
#include <fstream>

static basic::Tracer TR( "apps.pilot.tgtres_polar_sho_energies" );

using core::Size;
using core::Real;
using basic::options::option;
using namespace basic::options::OptionKeys;

OPT_KEY( String, tgt_set )

///
/// @brief loads a set of residues from file. The file format is as follows:
///
///  C1 R1 I1\n
///  ...
///  CN RN IN\n
///
///   Here, Ci, Ri, and Ii indicate the chain identifier, residue index, and
///  insertion code (as specified in the pose's input PDB file) of the ith
///  residue in the set (i=1,...,N; N>=1).
///
/// @param[in] path to the input file.
/// @param[out] rset vector to hold the residues. The vector must be passed
///  empty.
/// @param[in] ps the pose.
///
/// @details: after this function has been called, rset[i] is the pose index
///  of the residue specified by the ith input line (i=1,...,N).
///
/// @details: blank chain identifiers and insertion codes must be specified
///  with the '_' character.
///
void load_set(std::string setf, utility::vector1<Size>& rset,
	core::pose::Pose &ps) {

	std::ifstream setfs(setf.c_str());
	if ( !setfs ) {
		TR << "can't open " << setf << std::endl;
		exit(0);
	}

	char cid;
	int idx;
	char ico;
	while ( setfs >> cid >> idx >> ico ) {

		if ( cid == '_' ) {
			cid = ' ';
		}

		if ( ico == '_' ) {
			ico = ' ';
		}

		rset.push_back(ps.pdb_info()->pdb2pose(cid, idx, ico));
	}
}


/// @brief prints an atom's identifier
///
/// @param[in] ato_idx index of the atom in its residue
/// @param[in] res_idx index of the residue in the pose
/// @param[in] ps the pose
///
void print_atom_info(
	Size ato_idx, Size res_idx, core::pose::Pose const& ps) {

	TR <<
		ps.pdb_info()->chain(res_idx) <<
		ps.pdb_info()->number(res_idx) <<
		"(" << res_idx << ')' <<
		ps.residue(res_idx).atom_name(ato_idx);
}


///
/// @brief accumulates the SHO energies of a residue's polar hydrogens that share a common base atom
///
/// @param[in] base_idx index of the base atom in the residue
/// @param[in] res_idx index of the residue in the pose
/// @param[in] ps the pose
/// @param[in] sho_meth energy method to compute SHO energies
///
/// @return the accumulated energy of such hydrogens
///
Real accum_hbdon_shoene(Size const base_idx, Size const res_idx, core::pose::Pose const& ps,
	core::scoring::geometric_solvation::ExactOccludedHbondSolEnergyCOP const& sho_meth) {

	core::conformation::Residue const& res = ps.residue(res_idx);

	core::Real tot = 0;
	for ( core::chemical::AtomIndices::const_iterator
			hnum = res.Hpos_polar().begin(), hnume = res.Hpos_polar().end();
			hnum != hnume; ++hnum ) {

		if ( res.atom_base( *hnum ) == base_idx ) {
			tot += sho_meth->compute_donor_atom_energy(res, res_idx, *hnum, ps);
		}
	}

	return tot;
}


/// MAIN

int main( int argc, char * argv [] )
{
	try {
		NEW_OPT( tgt_set, "set of residues the SHO energy of whose atoms is to be computed", "" );

		devel::init(argc, argv);

		// load pose from pdb file
		core::pose::Pose ps;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_pdb( ps, input_pdb_name );

		// score pose
		core::scoring::ScoreFunctionOP scorefxn(
			core::scoring::get_score_function());

		(*scorefxn)(ps);

		// build method to compute SHO energies
		core::scoring::geometric_solvation::ExactOccludedHbondSolEnergyCOP sho_meth(
			core::scoring::geometric_solvation::create_ExactSHOEnergy_from_cmdline());
		sho_meth->setup_for_scoring(ps, *scorefxn);

		// load target residue set
		utility::vector1<Size> targets;
		load_set(option[tgt_set], targets, ps);

		Size const N = targets.size();
		for ( Size i=1; i<=N; ++i ) {

			Size const residx = targets[i];
			core::conformation::Residue const& res = ps.residue(residx);

			Size const NHVY = res.nheavyatoms();
			for ( Size j = 1; j <= NHVY; ++j ) {

				bool atom_is_polar = false;
				core::Real atom_energy = 0;

				if ( res.heavyatom_is_an_acceptor(j) ) {
					atom_energy = sho_meth->compute_acceptor_atom_energy(res, residx, j, ps);
					atom_is_polar = true;
				}

				if ( res.heavyatom_has_polar_hydrogens(j) ) {
					atom_energy += accum_hbdon_shoene(j, residx, ps, sho_meth);
					atom_is_polar = true;
				}

				if ( atom_is_polar ) {
					print_atom_info(j, residx, ps);
					TR << " " << atom_energy << std::endl;
				}
			}
		}

		TR << "TASK COMPLETED" << std::endl;

	} // try
catch (utility::excn::Exception const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}
}


