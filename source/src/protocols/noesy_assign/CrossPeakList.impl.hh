// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_CrossPeakList_impl_HH
#define INCLUDED_protocols_noesy_assign_CrossPeakList_impl_HH


// Unit Header
#include <protocols/noesy_assign/CrossPeakList.hh>

// Package Headers
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/DistanceScoreMover.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/jumping/JumpSample.hh>

#include <protocols/noesy_assign/StructureDependentPeakCalibrator.hh>
#include <protocols/noesy_assign/StructureIndependentPeakCalibrator.hh>

// Project Headers
#include <core/types.hh>
//#include <core/id/NamedAtomID.fwd.hh>
#include <core/chemical/AA.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/exit.hh>
// #include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
// #include <core/options/option.hh>
// #include <core/options/keys/abinitio.OptionKeys.gen.hh>
// #include <core/options/keys/run.OptionKeys.gen.hh>
//#include <core/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <list>
#include <map>
#include <cmath>


namespace protocols {
namespace noesy_assign {


template< class DecoyIterator >
void CrossPeakList::update_decoy_compatibility_score( DecoyIterator const& decoys_begin, DecoyIterator const& decoys_end ) {
	using namespace core;
	using namespace basic;
	basic::Tracer tr( "devel.noesy_assign.crosspeaks" );

	if ( decoys_begin == decoys_end ) return;

	Size ct_structures( 0 );
	for ( DecoyIterator iss = decoys_begin; iss != decoys_end; ++iss ) {
		++ct_structures;
	}
	tr.Info << "compute decoy compatibility from " << ct_structures << " structures" << std::endl;

	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );

	//compatibility with decoys.
	pose::Pose pose;
	decoys_begin->fill_pose( pose );
	protocols::jumping::JumpSample jumps( pose.fold_tree() );
	jumps.remove_chainbreaks( pose );
	DistanceScoreMoverOP cst_tool( new DistanceScoreMover( *this, pose, params.dcut_ ) );
	cst_tool->prepare_scoring( false /*not use for calibration */ );

	for ( DecoyIterator iss = decoys_begin; iss != decoys_end; ++iss ) {
		PROF_START( basic::NOESY_ASSIGN_DIST_MAKE_POSE );
		pose::Pose pose;
		iss->fill_pose( pose );
		protocols::jumping::JumpSample jumps( pose.fold_tree() );
		jumps.remove_chainbreaks( pose );
		PROF_STOP( basic::NOESY_ASSIGN_DIST_MAKE_POSE );
		tr.Debug << "score decoys " << iss->decoy_tag() << std::endl;
		cst_tool->apply( pose );
	}

	cst_tool->finalize_scoring();
	tr.Debug << "finished with decoy compatibility" << std::endl;
}


template < class DecoyIterator >
void CrossPeakList::calibrate( DecoyIterator const& begin, DecoyIterator const& end ) {
	using namespace core;
	using namespace basic;
	basic::Tracer tr( "protocols.noesy_assign.crosspeaks" );

	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	bool const structure_independent_calibration( params.calibration_target_ > 1.0 || begin == end );
	bool const elimination( params.calibration_eliminate_ && begin != end );

	if ( structure_independent_calibration ) {
		tr.Info << "structure independent calibration..."<<std::endl;
		PeakCalibratorMap calibrators( *this, PeakCalibratorOP( new StructureIndependentPeakCalibrator ) );
		calibrators.set_target_and_tolerance( params.calibration_target_, 0.1 );
		calibrators.do_calibration();
	};

	if ( !structure_independent_calibration || elimination ) {
		typedef utility::vector1< pose::PoseOP > Poses;
		Poses pose_cache;
		for ( DecoyIterator iss = begin; iss != end; ++iss ) {
			pose::PoseOP pose( new pose::Pose );
			iss->fill_pose( *pose ); //has to reread RDC file for each pose!
			protocols::jumping::JumpSample jumps( pose->fold_tree() );
			jumps.remove_chainbreaks( *pose );
			pose_cache.push_back( pose );
		}

		for ( Size cycles( params.calibration_cycles_ ); cycles >= 1; --cycles ) {
			PeakCalibratorMap calibrators( *this, PeakCalibratorOP( new StructureDependentPeakCalibrator( pose_cache, params.dcalibrate_ ) ) );
			if ( !structure_independent_calibration ) {
				tr.Info << "structure dependent calibration..."<<std::endl;
				calibrators.set_target_and_tolerance( params.calibration_target_, 0.005 );
				calibrators.do_calibration();
			}
			tr.Info << "structure dependent elimination and nudging ..."<< std::endl;
			calibrators.eliminate_violated_constraints(); //here also the nudging happens...
		}
	}
}

}
}

#endif
