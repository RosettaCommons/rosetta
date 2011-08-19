// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @detailed
///
///
/// @author Oliver Lange



// Unit Headers
#include <protocols/evaluation/RDC_Evaluator.hh>

// Package Headers
#include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/scoring/ScoreFunction.hh>

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>

// option key includes

// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <utility/options/keys/BooleanOptionKey.hh>

// C++ headers
#include <iterator>
#include <vector>

static basic::Tracer tr("protocols.evalution.RMSD");

namespace protocols {
namespace evaluation {
using namespace core;


SelectRDC_Evaluator::SelectRDC_Evaluator( std::list< Size > const& selection, std::string tag, std::string file )
  : SingleValuePoseEvaluator< Real >( "rdc"+tag ),
		selection_( selection ),
    tag_ ( tag ),
		rdc_file_( file )
{
	init_rdcs();
}

SelectRDC_Evaluator::SelectRDC_Evaluator( utility::vector1< Size> const& selection, std::string tag, std::string file )
  : SingleValuePoseEvaluator< Real >( "rdc"+tag ),
		tag_( tag ),
		rdc_file_ ( file )
{
	copy( selection.begin(), selection.end(), std::back_inserter( selection_ ) );
	init_rdcs();
}


SelectRDC_Evaluator::SelectRDC_Evaluator( core::pose::PoseCOP pose, std::string tag )
	: SingleValuePoseEvaluator< Real >( "rdc"+tag ),
		tag_( tag )
{
	find_existing_residues( pose, tag, selection_ );
	init_rdcs();
}

SelectRDC_Evaluator::SelectRDC_Evaluator( core::pose::Pose const& pose, std::string tag )
	: SingleValuePoseEvaluator< Real >( "rdc"+tag ),
		tag_( tag )
{
	find_existing_residues( new core::pose::Pose( pose ), tag, selection_ );
	init_rdcs();
}

Real
SelectRDC_Evaluator::apply( core::pose::Pose& pose ) const {
  core::Real rdc;
	scoring::ScoreFunction scorefxn;
  scorefxn.set_weight( scoring::rdc, 1 );
	pose::Pose test_pose( pose );
	scoring::store_RDC_in_pose( rdc_data_, test_pose );
	rdc = scorefxn( test_pose );
	return rdc;
}

void
SelectRDC_Evaluator::init_rdcs() {
	using namespace scoring;
	ResidualDipolarCoupling orig_rdcs( rdc_file_ );;
	ResidualDipolarCoupling::RDC_lines const& rdcs = orig_rdcs.get_RDC_data();
	if ( selection_.size() ) {
		ResidualDipolarCoupling::RDC_lines filtered;
		for ( ResidualDipolarCoupling::RDC_lines::const_iterator it = rdcs.begin(); it != rdcs.end(); ++it ) {
			core::scoring::ResidueSelection::const_iterator iter1 = find( selection_.begin(), selection_.end(), it->res1() );
			core::scoring::ResidueSelection::const_iterator iter2 = find( selection_.begin(), selection_.end(), it->res2() );
			if ( iter1 != selection_.end() && iter2 != selection_.end() ) {
				filtered.push_back( *it );
			}
		}
		rdc_data_ = new ResidualDipolarCoupling( filtered );
	} else {
		rdc_data_ = new ResidualDipolarCoupling( rdcs );
	}
}
}
}
