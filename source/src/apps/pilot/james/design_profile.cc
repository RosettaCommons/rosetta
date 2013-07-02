// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file design_profile.cc
/// @brief
/// @author James Thompson

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>

// option key includes

#include <basic/options/keys/james.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

class ComputeProfileMover : public protocols::moves::Mover {
	typedef utility::file::FileName FileName;

public:

	ComputeProfileMover() :
		protocols::moves::Mover("ComputeProfileMover"),
		n_designs_( 10 ),
		start_from_poly_ala_( true )
	{}

	ComputeProfileMover(
		core::Size n_designs,
		bool start_from_poly_ala = true
	) :
		protocols::moves::Mover("ComputeProfileMover"),
		n_designs_( n_designs ),
		start_from_poly_ala_( start_from_poly_ala )
	{}

	virtual ~ComputeProfileMover() {}

	Size n_designs() {
		return n_designs_;
	}

	void n_designs( core::Size n_designs ) {
		n_designs_ = n_designs;
	}

	/// @brief Uses fixed-backbone design to compute a profile for the given Pose.
	virtual void apply( core::pose::Pose & pose ) {
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
		protocols::simple_moves::PackRotamersMoverOP pack_mover
			= new protocols::simple_moves::PackRotamersMover;
		pack_mover->score_function( scorefxn );

		core::sequence::SequenceAlignment aln;
		std::string const nat_id( "native_seq" );
		core::sequence::Sequence native_seq( pose.sequence(), nat_id );
		aln.add_sequence( native_seq );

		//std::cout << "making " << n_designs() << " designs." << std::endl;
		for ( Size design_count = 1; design_count <= n_designs(); ++design_count ) {
			// make a copy of the Pose, convert all positions to alanine
			core::pose::Pose poly_ala_pose = pose;
			utility::vector1< core::Size > positions;
			for ( Size i = 1; i <= pose.total_residue(); ++i )
				positions.push_back( i );

			protocols::toolbox::pose_manipulation::construct_poly_ala_pose(
				poly_ala_pose, positions, false, false, false
			);

			//std::cout << "starting with sequence " << poly_ala_pose.sequence()
				//<< std::endl;

			// run fixbb design on the poly_ala_pose
			pack_mover->apply( poly_ala_pose );
			//core::Real design_score = (*scorefxn)(poly_ala_pose);

			// get back the sequence, add it to aln
			std::string const seq_id( "des_" + string_of( design_count) );
			core::sequence::Sequence design_seq( poly_ala_pose.sequence(), seq_id );
			aln.add_sequence( design_seq );
			//std::cout << "designed sequence " << design_seq
			//					<< " with score " << design_score << std::endl;
		} // n_designs

		//std::cout << "sequences are:" << std::endl << aln << std::endl;
	} // apply

	private:
		Size const n_designs_;
		bool const start_from_poly_ala_;
}; // ComputeProfileMover

int
main( int argc, char * argv [] ) {
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	devel::init( argc, argv );

	core::Size const n_designs( option[ james::n_designs ]() );
	ComputeProfileMover mover( n_designs );
	protocols::jobdist::not_universal_main( mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
	return 0;
}
