// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief stupid test file for visual studio c++
/// @detailed

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/fragment/FragSet.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/ProteinSilentStruct.hh>


#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/io/ozstream.hh>

#include <string>

#include <numeric/random/random.hh>
#include <numeric/random/random.fwd.hh>


// option key includes

#include <devel/init.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char* argv [] ) {
	try {

	using std::min;
	using core::Size;
	using core::Real;
	using std::string;

	using namespace protocols::moves;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;

	using numeric::random::uniform;
	using numeric::random::random_range;

	devel::init( argc, argv );
	ScoreFunctionOP scorefxn( get_score_function() );

	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb(true);
	movemap->set_jump(true);

	core::pose::Pose fold_pose;
	MetaPoseInputStream input = streams_from_cmd_line();

	ResidueTypeSetCAP rsd_set =
		ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	input.fill_pose( fold_pose, *rsd_set );

	SilentFileData sfd;
	SilentStructOP ss = SilentStructFactory::get_instance()->get_silent_struct_out();
	string tag( "start" );
	ss->fill_struct( fold_pose, tag );
	sfd.write_silent_struct( *ss, option[ out::file::silent ]() );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
