// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <devel/init.hh>


#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/kinematics/FoldTree.hh>
#include <protocols/rbsegment_moves/RMSVallData.hh>

#include <core/chemical/ChemicalManager.hh>



#include <core/conformation/ResidueFactory.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>


#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <protocols/evaluation/ChemicalShiftEvaluator.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <protocols/viewer/viewers.hh>
#include <basic/Tracer.hh>

#include <core/id/NamedStubID.hh>
#include <core/kinematics/Stub.hh>
#include <core/id/types.hh>
#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



int main(int argc, char **argv) {
	using namespace core::fragment;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;

	using core::Size;
	using std::string;

	std::cout << "USAGE:  \n";
	std::cout << "USAGE: " << argv[0] << "\n";
	std::cout << "               -in::file::native <CA-trace>\n";
	std::cout << "               -loops::vall_file <vall>\n";
	std::cout << "             [ -out::file::frag_prefix <prefix> ]\n";
	std::cout << "USAGE:  \n";
	devel::init( argc,argv );

	core::pose::Pose native_pose;
	core::import_pose::centroid_pose_from_pdb( native_pose, option[ in::file::native ]() );
	core::Size nres=native_pose.total_residue();
	std::string input_seq = native_pose.sequence();

	// load vall
	protocols::rbsegment_Moves::RMSVallData rms_vall( option[ OptionKeys::loops::vall_file ] );
	ConstantLengthFragSet frags3(3),frags9(9);

	// for each 3 mer get a frame
	for (int i =  1; i <= nres - 2; ++i) {
		FrameOP frame3_i = new core::fragment::Frame( i, 3 );
		utility::vector1< numeric::xyzVector< core::Real> > cas( 3 );
		for (int k=0; k<3; ++k)
			cas[k+1] = native_pose.residue(i+k).atom("CA").xyz();
		std::string frag_seq = input_seq.substr( i-1, 3 );
		rms_vall.get_frags( 200, cas, frag_seq, '-', frame3_i, 0.0 );
		frags3.add( frame3_i );
	}

	// for each 9 mer get a frame
	for (int i =  1; i <= nres - 8; ++i) {
		FrameOP frame9_i = new core::fragment::Frame( i, 9 );
		utility::vector1< numeric::xyzVector< core::Real> > cas( 9 );
		for (int k=0; k<9; ++k)
			cas[k+1] = native_pose.residue(i+k).atom("CA").xyz();
		std::string frag_seq = input_seq.substr( i-1, 9 );
		rms_vall.get_frags( 200, cas, frag_seq, '-', frame9_i, 0.0 );
		frags9.add( frame9_i );
	}

	// dump frame sets
	FragmentIO().write( option[ OptionKeys::out::file::frag_prefix ]()+"3", frags3 );
	FragmentIO().write( option[ OptionKeys::out::file::frag_prefix ]()+"9", frags9 );

	return 0;
}
