// -*-
// mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t
// -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/** @page deepMove
	Read a PDB, move it, print it
	In this version, move using Jump methods
	Try:
		"readPDB.cc -in::file::s <pdb file> -in::path::database <DB root dir>"
*/

///
/// @file   apps/pilot/lemmon/deepMove.cc
///
/// @brief This is to illustrate packing residues in a PDB
/// @detail Run this script with the following arguments:
/// 1) in::file::s <list of one PDB with at least 2 residues to pack>
/// 2) in::path::database <list of one database root directory>
/// 3) in::file::extra_res_fa <list of extra params files, one per residue type>
/// @author Gordon Lemmon (glemmon@gmail.com)

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
//#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/Jump.hh>
#include <numeric/xyzVector.hh>
#include <core/types.hh> // for "Real" type


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    try {
	devel::init(argc, argv);

	utility::vector0<std::string> pdbs;
	{// process the options
		using namespace basic::options::OptionKeys;
		using basic::options::option;
		pdbs= option[in::file::s]();
	}
	core::pose::Pose pose; // starts NULL, coords *never* modified!
	{
		std::string pdb=pdbs[0];
		core::import_pose::pose_from_pdb(pose, pdb);
	}
	core::kinematics::Jump jump;
	int jump_num=1; // this assumes there is only 1 jump (2 molecules) in this PDB
	{
		jump= pose.jump(jump_num);
	}
	{
		int xAxis=1; // use 2-6 for y,z and xyz rotation
		{
			numeric::xyzVector< core::Real > translation  = jump.get_translation();
			double x= translation[xAxis];
			translation[xAxis]=  x+x*.2; // move the ligand an arbitrary amount along x axis
			jump.set_translation(translation);
		}
		// prefer the 3 lines below over the 4 lines above
		int const downstream=1;// anything else is upstream
		{
			core::Real const length= 100;// angstrom?
			jump.set_rb_delta(xAxis, downstream, length);
		}
		int alphaAxis=4;// 5&6 are beta and gamma (euler angles)
		{
			core::Real const radians=numeric::NumericTraits<Real>::pi();// not sure if this is degree or radians
			jump.set_rb_delta(alphaAxis, downstream, radians);
		}
		jump.fold_in_rb_deltas();// adjusts the center accordingly
	}
	pose.set_jump(jump_num, jump);
	// pose.set_jump(jumpNum, jump); // this form updates xyz coords of downstream atoms
	{
		const std::string output("output.pdb");
		pose.dump_pdb(output);
	}

    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}
