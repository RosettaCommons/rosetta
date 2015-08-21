// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rna/RNA_MC_MultiSuite.cc
/// @brief Markov chain sampler for mulitple RNA suite.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_MC_MultiSuite.hh>

// Package headers
#include <protocols/stepwise/sampler/rna/RNA_MC_Suite.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

using namespace core;

static thread_local basic::Tracer TR( "protocols.sampler.rna.RNA_MC_MultiSuite" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

///////////////////////////////////////////////////////////////////////////
RNA_MC_MultiSuite::RNA_MC_MultiSuite():
	MC_Comb()
{}
///////////////////////////////////////////////////////////////////////////
void  RNA_MC_MultiSuite::set_pucker_flip_rate( Real const setting ) {
	for ( Size i = 1; i <= suite_samplers_.size(); ++i ) {
		suite_samplers_[i]->set_pucker_flip_rate( setting );
	}
}
///////////////////////////////////////////////////////////////////////////
void  RNA_MC_MultiSuite::set_pucker_flip_rate(
	Real const setting,
	Size const rotamer_id
) {
	runtime_assert( rotamer_id <= suite_samplers_.size() );
	suite_samplers_[rotamer_id]->set_pucker_flip_rate( setting );
}

///////////////////////////////////////////////////////////////////////////
void  RNA_MC_MultiSuite::set_gaussian_stdev( core::Real const setting ) {
	for ( Size i = 1; i <= suite_samplers_.size(); ++i ) {
		suite_samplers_[i]->set_gaussian_stdev( setting );
	}
}
///////////////////////////////////////////////////////////////////////////
void  RNA_MC_MultiSuite::set_gaussian_stdev(
	Real const setting,
	Size const rotamer_id
) {
	runtime_assert( rotamer_id <= suite_samplers_.size() );
	suite_samplers_[rotamer_id]->set_gaussian_stdev( setting );
}
///////////////////////////////////////////////////////////////////////////
void  RNA_MC_MultiSuite::add_external_loop_rotamer(
	MC_StepWiseSamplerOP const & rotamer
) {
	RNA_MC_SuiteOP new_rotamer;
	new_rotamer = utility::pointer::dynamic_pointer_cast< RNA_MC_Suite >( rotamer );
	runtime_assert( new_rotamer != 0 );
	suite_samplers_.push_back( new_rotamer );
	MC_Comb::add_external_loop_rotamer( rotamer );
}
///////////////////////////////////////////////////////////////////////////
void  RNA_MC_MultiSuite::clear_rotamer() {
	suite_samplers_.clear();
	MC_Comb::clear_rotamer();
}
///////////////////////////////////////////////////////////////////////////
void  RNA_MC_MultiSuite::set_init_from_pose( pose::Pose const & pose ) {
	for ( Size i = 1; i <= suite_samplers_.size(); ++i ) {
		suite_samplers_[i]->set_init_from_pose( pose );
	}
}
///////////////////////////////////////////////////////////////////////////
} //rna
} //sampler
} //stepwise
} //protocols
