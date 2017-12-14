// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @details
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/simple_filters/RDC_Evaluator.hh>

// Package Headers
#include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/scoring/methods/ResidualDipolarCouplingEnergy.hh>


// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

// option key includes


#include <protocols/evaluation/util.hh>
#include <utility>
#include <utility/vector1.hh>

static basic::Tracer tr( "protocols.evalution.RMSD" );

namespace protocols {
namespace simple_filters {
using namespace core;

RDC_Evaluator::RDC_Evaluator( std::string const & tag ) :
	evaluation::SingleValuePoseEvaluator< Real >( tag )
{}

//RDC_Evaluator::RDC_Evaluator( utility::vector1< std::string > const& rdc_files, std::string tag) :
// evaluation::SingleValuePoseEvaluator< Real >( tag ),
// rdc_files_( rdc_files )
//{
// init_rdcs();
//}

static core::scoring::methods::ResidualDipolarCouplingEnergy energy_evaluator;
/// @brief evaluate pose
core::Real
RDC_Evaluator::apply( core::pose::Pose& pose ) const {
	return energy_evaluator.eval_dipolar( pose, rdc_data_ ); //const
}
/// @brief evaluate pose
SelectRDC_Evaluator::SelectRDC_Evaluator( std::list< Size > const & selection, std::string const & tag, std::string const & file )
: evaluation::SingleValuePoseEvaluator< Real >( "rdc"+tag ),
	selection_( selection ),
	tag_ ( tag ),
	rdc_file_( file )
{
	init_rdcs();
}

SelectRDC_Evaluator::SelectRDC_Evaluator( utility::vector1< Size> const& selection, std::string const & tag, std::string const & file )
: evaluation::SingleValuePoseEvaluator< Real >( "rdc"+tag ),
	tag_( tag ),
	rdc_file_ ( file )
{
	copy( selection.begin(), selection.end(), std::back_inserter( selection_ ) );
	init_rdcs();
}


SelectRDC_Evaluator::SelectRDC_Evaluator( core::pose::PoseCOP pose, std::string tag )
: evaluation::SingleValuePoseEvaluator< Real >( "rdc"+tag ),
	tag_( tag )
{
	evaluation::find_existing_residues( pose, tag, selection_ );
	init_rdcs();
}

SelectRDC_Evaluator::SelectRDC_Evaluator( core::pose::Pose const& pose, std::string tag )
: evaluation::SingleValuePoseEvaluator< Real >( "rdc"+tag ),
	tag_( tag )
{
	evaluation::find_existing_residues( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) ), tag, selection_ );
	init_rdcs();
}

Real
SelectRDC_Evaluator::apply( core::pose::Pose& pose ) const {
	return energy_evaluator.eval_dipolar( pose, *rdc_data_ ); //cons
}

void
SelectRDC_Evaluator::init_rdcs() {
	using namespace scoring;
	ResidualDipolarCoupling orig_rdcs( rdc_file_ );
	ResidualDipolarCoupling::RDC_lines const& rdcs = orig_rdcs.get_RDC_data();
	if ( selection_.size() ) {
		ResidualDipolarCoupling::RDC_lines filtered;
		for ( auto const & rdc : rdcs ) {
			core::scoring::ResidueSelection::const_iterator iter1 = find( selection_.begin(), selection_.end(), rdc.res1() );
			core::scoring::ResidueSelection::const_iterator iter2 = find( selection_.begin(), selection_.end(), rdc.res2() );
			if ( iter1 != selection_.end() && iter2 != selection_.end() ) {
				filtered.push_back( rdc );
			}
		}
		rdc_data_ = core::scoring::ResidualDipolarCouplingOP( new ResidualDipolarCoupling( filtered ) );
	} else {
		rdc_data_ = core::scoring::ResidualDipolarCouplingOP( new ResidualDipolarCoupling( rdcs ) );
	}
}
}
}
