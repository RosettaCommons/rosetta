// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

///
/// @brief Computes the ligand RMSD between two poses that are different
/// 	instances of the same protein+ligand system. The RMSD is computed
// 		after superposition of the two instances.
///
/// @param[in] -s <REF>, where <REF> is the path to the PDB file containing the
/// 	pose that during the superposition remains fixed.
/// @param[in] -mov <MOV>, where <MOV> is the path to the PDB file containing
/// 	the pose that during the superposition is moved.
/// @param[in] -extra_res_fa <PARAM>, where <PARAM> is the path to the params
/// 	file describing the ligand residue type common to both poses.
///
/// @details The two poses are superposed by optimally fitting their CA traces
/// 	The resulting ligand RMSD is computed on heavy atoms only.
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)
///

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <utility/excn/Exceptions.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/file_data.hh>
#include <fstream>

using core::Size;
using core::pose::Pose;
using core::id::AtomID;
using core::conformation::Residue;


///
/// @brief superimposes a pose onto another by optimizing the fit of their CA traces
///
/// @param[in] mod_pose pose to be moved
/// @param[in] ref_pose pose that remains fixed
///
/// @details the two poses are assumed to be different instances of the same
/// 	protein+ligand complex.
///
void calpha_superimpose_pose(Pose& mod_pose, Pose const& ref_pose) {

	// map CA atoms in mod_pose to their images in ref_pose, when available
	Size const N = mod_pose.total_residue()-1;
	core::id::AtomID_Map< AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, core::id::BOGUS_ATOM_ID );
	for ( Size ii = 1; ii <= N; ++ii ) {
		Residue const& mres = mod_pose.residue(ii);
		if ( ! mres.has("CA") ) continue;
		for ( Size jj = 1; jj <= N; ++jj ) {
			Residue const& rres = ref_pose.residue(jj);
			if ( ! rres.has("CA") ) continue;
			if ( mod_pose.pdb_info()->chain(ii) != ref_pose.pdb_info()->chain(jj)) continue;
			if ( mod_pose.pdb_info()->number(ii) != ref_pose.pdb_info()->number(jj)) continue;
			AtomID const id1( mres.atom_index("CA"), ii );
			AtomID const id2( rres.atom_index("CA"), jj );
			atom_map.set( id1, id2 );
			break;
		}
	}

	core::scoring::superimpose_pose( mod_pose, ref_pose, atom_map );
}


///
/// @brief returns the ligand rmsd between two poses
///
/// @param[in] pose1 first pose
/// @param[in] pose2 second pose
///
/// @details no superposition of the poses is applied
///
/// @details it is assumed that in either pose the ligand is the last residue
///
core::Real ligand_rmsd (Pose const& pose1, Pose const& pose2) {

	using core::Real;
	using core::conformation::Atom;

	Size const N = pose1.total_residue();

 	Residue const & lig1 = pose1.residue(N);
	Residue const & lig2 = pose2.residue(N);

	Size const NLH1 = lig1.nheavyatoms();
	Size const NLH2 = lig2.nheavyatoms();

	assert(NLH1 == NLH2);

	Real rmsd = 0;
	for(Size i = 1; i <= NLH1; ++i) {

   	assert ( lig1.atom_name(i) == lig2.atom_name(i) );

		Atom const& ato1 = lig1.atom(i);
		Atom const& ato2 = lig2.atom(i);

		Real sqdx = ato1.xyz()(1) - ato2.xyz()(1);
		sqdx = pow(sqdx, 2);

		Real sqdy = ato1.xyz()(2) - ato2.xyz()(2);
		sqdy = pow(sqdy, 2);

		Real sqdz = ato1.xyz()(3) - ato2.xyz()(3);
		sqdz = pow(sqdz, 2);

		rmsd += sqdx + sqdy + sqdz;
	}
	return sqrt(rmsd/NLH1);
}


OPT_KEY( String, mov )
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                   MAIN                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char * argv [] )
{

try {

	NEW_OPT( mov, "PDB file containing the structure to be superposed", "");

	devel::init(argc, argv);

	// create wild-type pose
	Pose ref_ps;
	std::string const ref_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_pdb( ref_ps, ref_pdb_name );

	// create mutant pose
	Pose mov_ps;
	std::string const mov_pdb_name =
		basic::options::option[basic::options::OptionKeys::mov];
	core::import_pose::pose_from_pdb( mov_ps, mov_pdb_name );

	// superimpose
	calpha_superimpose_pose(mov_ps, ref_ps);

	std::cout << "ligand RMSD: " << ligand_rmsd(mov_ps, ref_ps) << std::endl;

}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
}

}

