// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverOutputterContainer.cc
/// @brief a container class which holds a list of PeptideDeriverOutputter instances and delegates calls to those
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @author Orly Marcu (orly.marcu@mail.huji.ac.il)

// Unit header
#include <protocols/peptide_deriver/PeptideDeriverOutputterContainer.hh>
#include <protocols/peptide_deriver/PeptideDeriverFilter.hh>

// Project headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

static basic::Tracer TR( "protocols.peptide_deriver.PeptideDeriverOutputterContainer" );


namespace protocols {
namespace peptide_deriver {

PeptideDeriverOutputterContainer::PeptideDeriverOutputterContainer(): PeptideDeriverOutputter ()
{}

PeptideDeriverOutputterContainer::~PeptideDeriverOutputterContainer(){}

PeptideDeriverOutputterContainer::PeptideDeriverOutputterContainer( PeptideDeriverOutputterContainer const & src) {
	list_ = src.list_;
}

PeptideDeriverOutputterContainerOP
PeptideDeriverOutputterContainer::clone() const {
	return PeptideDeriverOutputterContainerOP( new PeptideDeriverOutputterContainer( *this ) );
}


void PeptideDeriverOutputterContainer::push_back(PeptideDeriverOutputterOP const item) {
	list_.push_back(item);
}

void PeptideDeriverOutputterContainer::clear() {
	list_.clear();
}

void PeptideDeriverOutputterContainer::begin_structure(core::pose::Pose const & pose, std::string const &name) {
	for ( PeptideDeriverOutputterOP const outputter : list_ ) {
		outputter->begin_structure(pose, name);
	}
}

void PeptideDeriverOutputterContainer::chain_pair_pose_prepared(core::pose::Pose const & pose) {
	for ( PeptideDeriverOutputterOP const outputter : list_ ) {
		outputter->chain_pair_pose_prepared(pose);
	}
}

void PeptideDeriverOutputterContainer::begin_receptor_partner_pair(char const receptor_chain_letter,
	char const partner_chain_letter, core::Real const total_isc,
	std::string const & options_string) {
	for ( PeptideDeriverOutputterOP const outputter : list_ ) {
		outputter->begin_receptor_partner_pair(receptor_chain_letter, partner_chain_letter, total_isc, options_string);
	}
}

void PeptideDeriverOutputterContainer::peptide_length(core::Size const pep_length) {
	for ( PeptideDeriverOutputterOP const outputter : list_ ) {
		outputter->peptide_length(pep_length);
	}
}

void PeptideDeriverOutputterContainer::peptide_entry(PeptideDeriverEntryType const entry_type, core::Real const total_isc,
	DerivedPeptideEntryCOP entry) {
	for ( PeptideDeriverOutputterOP const outputter : list_ ) {
		outputter->peptide_entry(entry_type, total_isc, entry);
	}
}

void PeptideDeriverOutputterContainer::end_receptor_partner_pair() {
	for ( PeptideDeriverOutputterOP const outputter : list_ ) {
		outputter->end_receptor_partner_pair();
	}
}

void PeptideDeriverOutputterContainer::end_structure() {
	for ( PeptideDeriverOutputterOP const outputter : list_ ) {
		outputter->end_structure();
	}
}

} //protocols
} //peptide_deriver






