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
/// @author James Thompson <tex@u.washington.edu>

// libRosetta headers


#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ReweightScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/rms_util.hh>

#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/PDBSilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <protocols/loop_build/LoopBuild.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.fwd.hh>
#include <protocols/comparative_modeling/ConstraintRemodelMover.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>
using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


/// @brief super-simple attempt to mimick the Modeller algorithm. Similar to FoldConstraints, except
/// moves are made based on minimization rather than based on fragment moves.

class ConstraintMinimizer : public protocols::moves::Mover {
	void apply( core::pose::Pose & pose ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring::constraints;
		using namespace protocols::abinitio;

		ConstraintSetOP orig_cstset
			= core::scoring::constraints::ConstraintIO::read_constraints(
			  core::scoring::constraints::get_cst_file_option(),
			  new core::scoring::constraints::ConstraintSet,pose
		);

		MaxSeqSepConstraintSetOP mss_cstset;
		mss_cstset = new MaxSeqSepConstraintSet( *orig_cstset, pose.fold_tree() );
		core::pose::Pose native_pose;
		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file( native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
		}

		protocols::comparative_modeling::ConstraintRemodelMover my_mover;
		for ( Size i = 1; i <= mss_cstset->largest_possible_sequence_sep( pose ); ++i ) {
			mss_cstset->set_max_seq_sep( i );
			std::cerr << "seqsep now set to " << i << std::endl;
			pose.constraint_set( mss_cstset );
			my_mover.apply( pose );

			if ( option[ james::debug ]() ) {
				core::io::silent::SilentFileData sfd;
				core::io::silent::SilentStructOP ss
					= core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
				ss->fill_struct( pose );
				ss->decoy_tag( "debug_" + string_of(i) );

				if ( option[ in::file::native ].user() ) {
					// calculate RMS
					core::Real rmsd  = core::scoring::CA_rmsd ( native_pose, pose );
					core::Real gdtmm = core::scoring::CA_gdtmm( native_pose, pose );

					ss->add_energy( "CA_rmsd", rmsd );
					ss->add_energy( "CA_gdtmm", gdtmm );
				}

				sfd.write_silent_struct( *ss, "james.debug" );
			} // option[ james::debug ]()
		} // for Size i
	} // apply
}; // ConstraintMinimizer

void*
my_main( void* ) {
	ConstraintMinimizer my_mover;
	protocols::jobdist::not_universal_main( my_mover );
	return 0;
}

int
main( int argc, char* argv [] )
{
	try {

	// options, random initialization
	devel::init( argc, argv );
	protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
