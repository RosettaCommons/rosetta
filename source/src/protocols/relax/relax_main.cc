// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to relax
/// @detailed
/// @author Mike Tyka, Monica Berrondo
/// @author Roland A. Pache


#include <protocols/relax/RelaxProtocolBase.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh> //pba
// AUTO-REMOVED #include <core/pose/datacache/CacheableDataType.hh> //pba

// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/util.hh>
// AUTO-REMOVED #include <core/scoring/electron_density/util.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>


#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>

#include <protocols/moves/MoverContainer.hh>
// AUTO-REMOVED #include <protocols/electron_density/util.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/relax/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer tr( "protocols.relax" );

using namespace core;
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace relax {
////////////////////////////////////////////////////////////////////////////////////////////////////

int
Relax_main( bool ) {
	using namespace protocols::jobdist;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::moves::MoverOP protocol = generate_relax_from_cmd();
	protocols::jd2::set_native_in_mover( *protocol );

	// add constraints from cmd line
	if ( option[ OptionKeys::constraints::cst_fa_file ].user() || option[ OptionKeys::constraints::cst_file ].user()) {
			protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
			protocols::simple_moves::ConstraintSetMoverOP loadCsts( new protocols::simple_moves::ConstraintSetMover );
            if (option[ OptionKeys::constraints::cst_fa_file ].user()) {
                loadCsts->constraint_file( core::scoring::constraints::get_cst_fa_file_option() );
            } else {
                loadCsts->constraint_file( core::scoring::constraints::get_cst_file_option() );
            }
            seqmov->add_mover( loadCsts );
			seqmov->add_mover( protocol );
			protocol = seqmov;
	}

	// set pose for density scoring if a map was input
	// (potentially) dock map into density
	if ( option[ OptionKeys::edensity::mapfile ].user() ) {
		protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
		seqmov->add_mover( new protocols::electron_density::SetupForDensityScoringMover );
		seqmov->add_mover( protocol );
		protocol = seqmov;
	}

	// setup symmetry mover ... this should happen _before_ SetupForDensityScoringMover
	//   to avoid adding extra VRTs
	if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
	    seqmov->add_mover( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
	    seqmov->add_mover( protocol );
	    protocol = seqmov;
	}

	// superimpose input model to the native structure or to a supplied PDB file
	if ( option[ OptionKeys::relax::superimpose_to_file ].user() ||
	     option[ OptionKeys::relax::superimpose_to_native ].user()
	 ) {
			core::pose::Pose ref_pose;
			std::string ref_filename;
			if(  option[ OptionKeys::relax::superimpose_to_file ].user() ) ref_filename = option[ basic::options::OptionKeys::relax::superimpose_to_file ]();
			if(  option[ OptionKeys::relax::superimpose_to_native ].user() ) ref_filename =  option[ basic::options::OptionKeys::in::file::native ]();
			core::import_pose::pose_from_pdb( ref_pose, ref_filename );
			protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
			protocols::simple_moves::SuperimposeMoverOP sm  =  new protocols::simple_moves::SuperimposeMover;
			sm->set_reference_pose( ref_pose );
			seqmov->add_mover( sm );
			seqmov->add_mover( protocol );
			protocol = seqmov;
	}

	protocols::jd2::JobDistributor::get_instance()->go( protocol );

	return 0;
}



}
}
