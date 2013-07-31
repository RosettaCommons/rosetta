// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief protocols for folding into density
/// @detailed
/// @author Frank DiMaio

#include <protocols/rbsegment_relax/OptimizeThreading.hh>
#include <protocols/rbsegment_relax/OptimizeThreadingCreator.hh>

#include <protocols/loops/loops_main.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Remarks.hh>

#include <protocols/rbsegment_relax/util.hh>
#include <protocols/rbsegment_relax/RBSegmentMover.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/loops/util.hh>
#include <protocols/moves/DsspMover.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace rbsegment_relax {

static basic::Tracer TR("protocols.rbsegment_relax.OptimizeThreading");

using namespace protocols;
using namespace core;


//////////////////
//////////////////
/// creator


std::string
OptimizeThreadingMoverCreator::keyname() const {
	return OptimizeThreadingMoverCreator::mover_name();
}

protocols::moves::MoverOP
OptimizeThreadingMoverCreator::create_mover() const {
	return new OptimizeThreadingMover;
}

std::string
OptimizeThreadingMoverCreator::mover_name() {
	return "OptimizeThreading";
}


//////////////////
//////////////////
/// mover

void OptimizeThreadingMover::apply( core::pose::Pose & pose ) {
	using namespace rbsegment_relax;

	core::Size nres = hybridization::get_num_residues_nonvirt( pose );
	if ( core::pose::symmetry::is_symmetric( pose ) )
		scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction( *scorefxn_ );

	// get sses & choose cutpoints
	core::scoring::dssp::Dssp dssp_obj( pose );
	dssp_obj.insert_ss_into_pose( pose );
	protocols::loops::Loops SSEs = protocols::loops::extract_secondary_structure_chunks( pose, "HE", 3, 6, 3, 4);
	utility::vector1< RBResidueRange > segments;
	core::Size prevcut = 1;
	for (int j=2; j<=(int)SSEs.num_loop(); ++j) {
		core::Size j0stop = SSEs[j-1].stop();
		core::Size j1start = SSEs[j].start();
		core::Size cutpoint = (j0stop+j1start)/2;  // randomize?
		segments.push_back( RBResidueRange( prevcut, cutpoint-1 ) );
		//TR << "Add segment: " << prevcut << " , " << cutpoint-1 << std::endl;
		prevcut = cutpoint;
	}
	segments.push_back( RBResidueRange( prevcut, nres ) );
	//TR << "Add segment: " << prevcut << " , " << nres << std::endl;

	// all continuous chains of chunks
	core::Size nsegments=segments.size();
	for (int i=1; i<=(int)nsegments; ++i) {
		for (int j=i+1; j<=(int)nsegments; ++j) {
			segments.push_back( RBResidueRange( segments[i].start(), segments[j].end() ) );
			//TR << "Add segment: " << segments[i].start() << " , " << segments[i].end() << std::endl;
		}
	}

	// mc loop
	core::Real best_score=99999, acc_score=99999;
	core::pose::Pose best_pose, acc_pose;
	utility::vector1< int > offsets( nres, 0 );
	for (int i=1; i<=(int)nsteps_; ++i) {
		// random segment
		core::Size seg_i = numeric::random::random_range(1, segments.size() );

		// apply movement & score pose
		SequenceShiftMover sshift( segments[seg_i], 4 ); // config max shift
		sshift.apply( pose );
		core::Real score_cst = (*scorefxn_)(pose);

		// apply shift to array and score
		core::Real score_aln = 0;
		for (int j=(int)segments[seg_i].start(); j<=(int)segments[seg_i].end(); ++j)
			offsets[j] += sshift.last_shift();
		for (int j=2; j<=(int)nres; ++j)
			score_aln += abs(offsets[j] - offsets[j-1]);

		// boltzmann
		core::Real temperature = 1.0;
		core::Real score = score_cst + weight_*score_aln;
		core::Real boltz_factor = ( acc_score - score ) / temperature;
		core::Real probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );
		if ( probability < 1 ) {
			if ( numeric::random::uniform() >= probability ) {
				// reject
				pose = acc_pose;
			} else {
				// accept, not new best
				acc_pose = pose;
				acc_score = score;
				TR << "Accept!    Step " << i << " score = " << score_cst << " + " << weight_ << " * " << score_aln << " = " << score << std::endl;
			}
		} else {
			TR << "New best!  Step " << i << " score = " << score_cst << " + " << weight_ << " * " << score_aln << " = " << score << std::endl;
			// new best
			acc_pose = pose;
			acc_score = score;
			best_pose = pose;
			best_score = score;
		}
	}

	pose = best_pose;
}

void OptimizeThreadingMover::parse_my_tag(
			utility::tag::TagPtr const tag,
			moves::DataMap & data,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ )
{
	if( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		scorefxn_ = (data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}

	nsteps_ = tag->getOption<core::Size>( "nsteps", 1000 );
	weight_ = tag->getOption<core::Real>( "weight", 0.1 );
}

}
}
