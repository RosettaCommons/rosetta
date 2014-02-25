// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerComb.cc
/// @brief Aggregate of multiple rotamer samplers for sampling combinatorially.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/RotamerComb.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rotamer_sampler.RotamerComb" );

using namespace core;

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
RotamerComb::RotamerComb():
	RotamerBase(),
	is_empty_( false )
{}

RotamerComb::~RotamerComb(){}
///////////////////////////////////////////////////////////////////////////
void RotamerComb::init() {
	runtime_assert( !rotamer_list_.empty() );
	is_empty_ = false;
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->init();
		if ( !rotamer_list_[i]->not_end() ) {
			TR << "Got a null rotamer sampler!" << std::endl;
			is_empty_ = true;
		}
	}
	set_random( random() );
	set_init( true );
	reset();
}
///////////////////////////////////////////////////////////////////////////
void RotamerComb::reset() {
	runtime_assert( is_init() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->reset();
	}
}
///////////////////////////////////////////////////////////////////////////
void RotamerComb::operator++() {
	runtime_assert( not_end() );

	if ( random() ) {
		for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
			++( *rotamer_list_[i] );
		}
	} else {
		for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
			++( *rotamer_list_[i] );
			if ( rotamer_list_[i]->not_end() ) {
				break;
			} else {
				if ( i < rotamer_list_.size() ) rotamer_list_[i]->reset();
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////
bool RotamerComb::not_end() const {
	runtime_assert( is_init() );
	if ( is_empty_ ) return false;
	if ( random() ) return true;
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		if ( !rotamer_list_[i]->not_end() ) return false;
	}
	return true;
}
///////////////////////////////////////////////////////////////////////////
void RotamerComb::apply( Pose & pose ) {
	runtime_assert( is_init() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->apply( pose );
	}
}
///////////////////////////////////////////////////////////////////////////
void RotamerComb::set_random( bool const setting ) {
	RotamerBase::set_random( setting );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->set_random( setting );
		runtime_assert( rotamer_list_[i]->random() == setting );
	}
}
///////////////////////////////////////////////////////////////////////////
}
}
