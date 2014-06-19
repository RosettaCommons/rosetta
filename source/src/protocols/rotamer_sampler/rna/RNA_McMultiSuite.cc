// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_McMultiSuite.cc
/// @brief Markov chain sampler for mulitple RNA suite.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_McMultiSuite.hh>

// Package headers
#include <protocols/rotamer_sampler/rna/RNA_McSuite.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

using namespace core;

static basic::Tracer TR( "protocols.rotamer_sampler.rna.RNA_McMultiSuite" );

namespace protocols {
namespace rotamer_sampler {
namespace rna {

///////////////////////////////////////////////////////////////////////////
RNA_McMultiSuite::RNA_McMultiSuite():
	McComb()
{}
///////////////////////////////////////////////////////////////////////////
void  RNA_McMultiSuite::set_pucker_flip_rate( Real const setting ) {
	for ( Size i = 1; i <= suite_samplers_.size(); ++i )
			suite_samplers_[i]->set_pucker_flip_rate( setting );
}
///////////////////////////////////////////////////////////////////////////
void  RNA_McMultiSuite::set_pucker_flip_rate(
	Real const setting,
	Size const rotamer_id
) {
	runtime_assert( rotamer_id <= suite_samplers_.size() );
	suite_samplers_[rotamer_id]->set_pucker_flip_rate( setting );
}

///////////////////////////////////////////////////////////////////////////
void  RNA_McMultiSuite::set_gaussian_stdev( core::Real const setting ) {
	for ( Size i = 1; i <= suite_samplers_.size(); ++i )
			suite_samplers_[i]->set_gaussian_stdev( setting );
}
///////////////////////////////////////////////////////////////////////////
void  RNA_McMultiSuite::set_gaussian_stdev(
	Real const setting,
	Size const rotamer_id
) {
	runtime_assert( rotamer_id <= suite_samplers_.size() );
	suite_samplers_[rotamer_id]->set_gaussian_stdev( setting );
}
///////////////////////////////////////////////////////////////////////////
void  RNA_McMultiSuite::add_rotamer(
	McRotamerOP const & rotamer
) {
	RNA_McSuiteOP new_rotamer;
	new_rotamer =	dynamic_cast<RNA_McSuite *>( rotamer() );
	runtime_assert( new_rotamer );
	suite_samplers_.push_back( new_rotamer );
	McComb::add_rotamer( rotamer );
}
///////////////////////////////////////////////////////////////////////////
void  RNA_McMultiSuite::clear_rotamer() {
	suite_samplers_.clear();
	McComb::clear_rotamer();
}
///////////////////////////////////////////////////////////////////////////
void  RNA_McMultiSuite::set_init_from_pose( pose::Pose const & pose ) {
	for ( Size i = 1; i <= suite_samplers_.size(); ++i )
			suite_samplers_[i]->set_init_from_pose( pose );
}
///////////////////////////////////////////////////////////////////////////
} //rna
} //rotamer_sampler
} //protocols
