// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Tells which atoms in a ligand have a sufficiently low SASA
///
/// @param[in] -s <PDB>, where <PDB> is the path to the PDB file containing
///   the ligand (and likely, the protein it is bound to).
///
/// @param[in] -extra_res_fa <PARAMS>, where <PARAMS> is the .params file
///  for the ligand.
///
/// @param[in] -lig_sasa_resfile <SASA>, where <SASA> is the path to a file
///  specifying the ligand atoms whose SASA is to be evaluated. The file
///  has the following format:
///
///   MAX_SASA\n
///   CID IDX\n
///   NAM_1\n
///   ...
///   NAM_N\n
///
///  There:
///  -MAX_SASA is a maximum SASA value that ligand atoms are checked against.
///  -CID is the chain identifier of the ligand in <PDB>.
///  -IDX is the index of the ligand in the chain, as specified in <PDB>.
///  -NAM_i is the name, as specified in <PDB>, of the ith ligand atom whose
///   SASA is to be evaluated (i=1,..,N).
///
/// @details
/// - If no atom names are present in <SASA>, the program evaluates SASA for
///  all heavy atoms of the ligand.
/// - The program outputs the names and SASA values of all atoms whose SASA
///   is to be evaluated, following the same order as in <SASA>, if any, or,
///   otherwise, Rosetta's heavy-atom order for the ligand. Atoms whose SASA
///   is greater than MAX_SASA have their names followed by the !HIGH! tag.
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/sasa.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
#include <core/pose/util.tmpl.hh>
#include <iomanip>
#include <string>
#include <fstream>
#include <stdlib.h>

using core::id::AtomID;
using core::Size;
using core::Real;
using core::pose::Pose;
using std::string;

static thread_local basic::Tracer TR( "apps.pilot.lig_low_sasa.main" );


/// @brief Returns the residue number of a residue in a pose.
///
/// @parm[in] pdbnum residue number of the residue in its PDB file.
/// @parm[in] pdbchn chain identifier of the residue in the PDB file.
/// @parm[in] ps pose that the residue has been loaded into.
///
Size get_pose_resnum(int const pdbnum, char const pdbchn, Pose& ps) {

	for ( Size j = 1; j <= ps.total_residue(); ++j ) {
		if ( ( ps.pdb_info()->chain(j) == pdbchn ) && (ps.pdb_info()->number(j) == pdbnum) ) {
			return j;
		}
	}

	// residue not found
	TR << "ERROR!! Could not find residue" << pdbnum << " and chain " << pdbchn << std::endl;
	exit(1);
}

OPT_KEY( String, lig_sasa_resfile )
using basic::options::option;
using basic::options::OptionKeys::lig_sasa_resfile;


int main( int argc, char * argv [] )
{

	try {

		NEW_OPT( lig_sasa_resfile, "ligand required SASA", "lig_sasa_resfile.txt" );

		devel::init(argc, argv);

		// create pose from pdb
		core::pose::Pose ps;
		string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_pdb( ps, input_pdb_name );


		// score pose
		core::scoring::ScoreFunctionOP scorefxn(core::scoring::get_score_function());
		(*scorefxn)(ps);

		// compute SASA for all atoms in the pose
		utility::vector1<Real> rsd_sasa( ps.n_residue(), 0.0 );

		core::id::AtomID_Map<Real> atom_sasa;
		core::pose::initialize_atomid_map( atom_sasa, ps, 0.0 );

		core::scoring::calc_per_atom_sasa( ps, atom_sasa, rsd_sasa, 1.4 );

		//// read ligand SASA requirements ////
		string const sasa_fnam = option[lig_sasa_resfile];
		std::ifstream sasa_fs(sasa_fnam.c_str());
		if ( !sasa_fs ) {
			TR << "can't open " << sasa_fnam << std::endl;
			return 0;
		}

		// maximum SASA
		Real max_sasa;
		sasa_fs >> max_sasa;
		Size const LINE_TERM_UB = 10;
		sasa_fs.ignore(LINE_TERM_UB, '\n');

		// chain id and number
		char lig_cid;
		int lig_idx;
		sasa_fs.get(lig_cid);
		sasa_fs >> lig_idx;
		sasa_fs.ignore(LINE_TERM_UB, '\n');

		// names of SASA relevant atoms. The ith name is stored at the ith
		// vector position (i=1,...,N, where N is the number of such atoms).
		utility::vector1<string> lig_sasa_atoms;
		string anam;
		while ( sasa_fs >> anam )
				lig_sasa_atoms.push_back(anam);

		Size const LIG_PS_IDX = get_pose_resnum(lig_idx, lig_cid, ps);
		core::conformation::Residue const& lig_res = ps.residue(LIG_PS_IDX);
		Size lig_nsasa = lig_sasa_atoms.size();

		// pick all heavy atoms if none are specified
		if ( lig_nsasa == 0 ) {
			lig_nsasa = lig_res.nheavyatoms();
			for ( Size i=1; i<=lig_nsasa; i++ ) {
				lig_sasa_atoms.push_back(lig_res.atom_name(i));
			}
		}

		// print atom names and their SASA values
		for ( Size i=1; i<=lig_nsasa; i++ ) {
			Size aidx = lig_res.atom_index(lig_sasa_atoms[i]);
			Real sasa = atom_sasa(AtomID(aidx, LIG_PS_IDX));
			TR << std::setw(4) << lig_sasa_atoms[i] << ": " << std::setw(10) << sasa;
			if ( sasa > max_sasa ) {
				TR << "  !HIGH!";
			}
			TR << std::endl;
		}

		// print average sasa of selected atoms
		Real tot = 0;
		for ( Size i=1; i<=lig_nsasa; i++ ) {
			Size aidx = lig_res.atom_index(lig_sasa_atoms[i]);
			tot += atom_sasa(AtomID(aidx, LIG_PS_IDX));
		}

		TR << std::endl;
		TR << "Average SASA of selected atoms: " << (tot/lig_nsasa) << std::endl;

		return 0;

	} // try
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}

} // main
