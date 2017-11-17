// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <devel/init.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/kinematics/FoldTree.hh>
#include <protocols/frags/RMSVallData.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/ResidueFactory.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <basic/database/open.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>


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


OPT_1GRP_KEY(Real, fpd, oversample)
OPT_1GRP_KEY(Boolean, fpd, skip3)
OPT_1GRP_KEY(Integer, fpd, nfrags)
OPT_1GRP_KEY(Integer, fpd, fraglen)


int main(int argc, char **argv) {
	try {
		using namespace core::fragment;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using basic::database::full_name;

		using core::Size;
		using std::string;

		NEW_OPT(fpd::oversample, "oversample rate; '-1=inf (slow!!!)", 10);
		NEW_OPT(fpd::fraglen, "fragment length", 9);
		NEW_OPT(fpd::nfrags, "fragments per position to generate", 25);

		std::cout << "USAGE:  \n";
		std::cout << "USAGE: " << argv[0] << "\n";
		std::cout << "               -in::file::native <CA-trace>\n";
		std::cout << "               -in::file::vall <vall>\n";
		std::cout << "               -fpd::oversample <vall>\n";
		std::cout << "               -fpd::nfrags <vall>\n";
		std::cout << "             [ -out::file::frag_prefix <prefix> ]\n";
		std::cout << "USAGE:  \n";
		devel::init( argc,argv );

		core::pose::Pose native_pose;
		core::import_pose::centroid_pose_from_pdb( native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);

		core::Size nres = native_pose.size();
		core::Size fraglength = option[ OptionKeys::fpd::fraglen ];
		int nfrags = option[ OptionKeys::fpd::nfrags ];
		std::string input_seq = native_pose.sequence();

		// load vall
		protocols::frags::RMSVallData rms_vall( full_name( option[ in::file::vall ][1] ) );
		ConstantLengthFragSet frags(fraglength);

		// for each 9 mer get a frame
		for (int i =  1; i <= nres - fraglength + 1; ++i) {
			FrameOP frame_i = new core::fragment::Frame( i, fraglength );
			utility::vector1< numeric::xyzVector< core::Real> > cas( fraglength );

			// check for cutpoints
			bool makefrags = true;
			for (int k=0; k<fraglength; ++k)
				if (native_pose.fold_tree().is_cutpoint( i+k ) || native_pose.residue(i+k).atom_index("CA") == 0)
					makefrags=false;

			if (makefrags) {
				for (int k=0; k<fraglength; ++k)
					cas[k+1] = native_pose.residue(i+k).atom("CA").xyz();
				std::string frag_seq = input_seq.substr( i-1, fraglength );
				rms_vall.get_frags( nfrags, cas, frag_seq, '-', frame_i, 0.0, option[ OptionKeys::fpd::oversample ]() );
				frags.add( frame_i );
			}
		}

		std::ostringstream oss;
		oss << option[ OptionKeys::out::file::frag_prefix ]() << fraglength;
		FragmentIO().write_data( oss.str() , frags );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
