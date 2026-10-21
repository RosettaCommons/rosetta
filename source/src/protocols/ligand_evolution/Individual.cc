// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/Individual.cc
/// @brief  Class definition for %Individual
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// unit headers
#include <protocols/ligand_evolution/Individual.hh>

// utility headers
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.ligand_evolution.Individual" );

namespace protocols {
namespace ligand_evolution {

Individual::Individual( LigandIdentifier const& identifier, utility::vector1< core::Size > const& parent_ids, std::string const& type_of_birth )
:
	identifier_( identifier ),
	parent_ids_( parent_ids ),
	type_of_birth_( type_of_birth )
{}

Individual::~Individual() = default;

core::Real Individual::score( std::string const& name ) const {
	if ( !is_scored_ ) { utility_exit_with_message( "Tried to access score but Individual is unscored." ); }
	return score_terms_.at( name );
}

void Individual::score( std::string const& name, core::Real score ) {
	score_terms_[ name ] = score;
}

core::Real Individual::score() const {
	if ( !is_scored_ ) { utility_exit_with_message( "Tried to access score but Individual is unscored." ); }
	return score_;
}

void Individual::score( core::Real score ) {
	is_scored_ = true;
	score_ = score;
}

bool Individual::is_scored() const {
	return is_scored_;
}

LigandIdentifier const& Individual::identifier() const {
	return identifier_;
}

bool Individual::id( core::Size id ) {
	if ( id_ == 0 ) {
		TR.Trace << "Set id " << id << std::endl;
		id_ = id;
		return true;
	} else {
		TR.Debug << "Prohibited overwrite of id " << id_ << " to id " << id << std::endl;
		return false;
	}
}

core::Size Individual::id() const {
	return id_;
}

std::map< std::string, core::Real > const& Individual::score_terms() const {
	return score_terms_;
}

utility::vector1< core::Size > const& Individual::parents() const {
	return parent_ids_;
}

std::string const& Individual::type_of_birth() const {
	return type_of_birth_;
}

}
}
