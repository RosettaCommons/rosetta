// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/rna/MC_RNA_MultiSuite.cc
/// @brief Markov chain sampler for mulitple RNA suite.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh>

// Package headers
#include <protocols/recces/sampler/rna/MC_RNA_Suite.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;

static THREAD_LOCAL basic::Tracer TR( "protocols.recces.sampler.rna.MC_RNA_MultiSuite" );

namespace protocols {
namespace recces {
namespace sampler {
namespace rna {

///////////////////////////////////////////////////////////////////////////
MC_RNA_MultiSuite::MC_RNA_MultiSuite():
	MC_Comb(),
	do_no_op_random_( false )
{
	set_name( "MC_RNA_MultiSuite" );
}
///////////////////////////////////////////////////////////////////////////
void MC_RNA_MultiSuite::operator++() {
	if ( do_no_op_random_ ) numeric::random::rg().uniform(); // increment
	MC_Comb::operator++();
}
///////////////////////////////////////////////////////////////////////////
void  MC_RNA_MultiSuite::set_pucker_flip_rate( Real const setting ) {
	for ( Size i = 1; i <= suite_samplers_.size(); ++i ) {
		suite_samplers_[i]->set_pucker_flip_rate( setting );
	}
}
///////////////////////////////////////////////////////////////////////////
void  MC_RNA_MultiSuite::set_pucker_flip_rate(
	Real const setting,
	Size const rotamer_id
) {
	runtime_assert( rotamer_id <= suite_samplers_.size() );
	suite_samplers_[rotamer_id]->set_pucker_flip_rate( setting );
}

///////////////////////////////////////////////////////////////////////////
void  MC_RNA_MultiSuite::set_gaussian_stdev( core::Real const setting ) {
	for ( Size i = 1; i <= suite_samplers_.size(); ++i ) {
		suite_samplers_[i]->set_gaussian_stdev( setting );
	}
}
///////////////////////////////////////////////////////////////////////////
void  MC_RNA_MultiSuite::set_gaussian_stdev(
	Real const setting,
	Size const rotamer_id
) {
	runtime_assert( rotamer_id <= suite_samplers_.size() );
	suite_samplers_[rotamer_id]->set_gaussian_stdev( setting );
}
///////////////////////////////////////////////////////////////////////////
void MC_RNA_MultiSuite::set_angle( pose::Pose const & pose ) {
	for ( Size i = 1; i<=suite_samplers_.size(); ++i ) {
		suite_samplers_[i]->set_angle( pose );
	}
}
///////////////////////////////////////////////////////////////////////////
void  MC_RNA_MultiSuite::add_rotamer(
	MC_RNA_SuiteOP rotamer
) {
	suite_samplers_.push_back( rotamer );
	MC_Comb::add_rotamer( rotamer );
}
///////////////////////////////////////////////////////////////////////////
void  MC_RNA_MultiSuite::clear_rotamer() {
	suite_samplers_.clear();
	MC_Comb::clear_rotamer();
}
///////////////////////////////////////////////////////////////////////////
void  MC_RNA_MultiSuite::set_init_from_pose( pose::Pose const & pose ) {
	for ( Size i = 1; i <= suite_samplers_.size(); ++i ) {
		suite_samplers_[i]->set_init_from_pose( pose );
	}
}
///////////////////////////////////////////////////////////////////////////
} //rna
} //sampler
} //recces
} //protocols
