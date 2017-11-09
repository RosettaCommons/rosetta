// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverBasicStreamOutputter.cc
/// @brief outputs a Peptiderive report to a stream in a basic, easily parsable format
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @author Orly Marcu (orly.marcu@mail.huji.ac.il)

// Unit headers
#include <protocols/peptide_deriver/PeptideDeriverBasicStreamOutputter.hh>
#include <protocols/peptide_deriver/PeptideDeriverFilter.hh>

// Project header
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility headers
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.hh>

// External headers
#include <boost/format.hpp>

// C++ headers
#include <string>

static THREAD_LOCAL basic::Tracer TR( "protocols.peptide_deriver.PeptideDeriverBasicStreamOutputter" );


namespace protocols {
namespace peptide_deriver {

PeptideDeriverBasicStreamOutputter::PeptideDeriverBasicStreamOutputter(utility::io::orstream & out, std::string prefix) {
	out_p_ = &out;
	prefix_ = prefix;

}

PeptideDeriverBasicStreamOutputter::~PeptideDeriverBasicStreamOutputter(){}

PeptideDeriverBasicStreamOutputter::PeptideDeriverBasicStreamOutputter( PeptideDeriverBasicStreamOutputter const & src ) {
	out_p_=src.out_p_;
	prefix_=src.prefix_;

}

PeptideDeriverBasicStreamOutputterOP
PeptideDeriverBasicStreamOutputter::clone() const {
	return PeptideDeriverBasicStreamOutputterOP( new PeptideDeriverBasicStreamOutputter( *this ) );
}

void PeptideDeriverBasicStreamOutputter::begin_structure(core::pose::Pose const &, std::string const &name) {
	(*out_p_) << prefix_ << "# structure: " << name << std::endl;
}

void PeptideDeriverBasicStreamOutputter::begin_receptor_partner_pair(char const receptor_chain_letter,
	char const partner_chain_letter, core::Real const total_isc,
	std::string const & options_string) {
	(*out_p_) << prefix_ << "# options: " << options_string << std::endl;
	(*out_p_) << prefix_ << "> chain_pair: receptor= " << receptor_chain_letter << " partner= " << partner_chain_letter << " total_isc= " << avoid_negative_zero(total_isc, 1e-7) << std::endl;
}

void PeptideDeriverBasicStreamOutputter::peptide_length(core::Size const pep_length) {
	(*out_p_) << prefix_ << ">> peptide_length: " << pep_length << std::endl;
}


void PeptideDeriverBasicStreamOutputter::peptide_entry(PeptideDeriverEntryType const entry_type,
	core::Real const total_isc, DerivedPeptideEntryCOP entry) {

	core::Real binding_contribution_fraction = entry->lin_isc / total_isc;

	std::string const binding_contribution_pct_str = binding_contribution_fraction >0 ? (boost::format("%1$.3f") % binding_contribution_fraction).str() : "*";

	(*out_p_) << prefix_
		<< "| "
		<< entry_type << " "
		<< entry->pep_start << " "
		<< avoid_negative_zero(entry->lin_isc, 1e-7) << " "
		<< binding_contribution_pct_str;

	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		CyclizedPeptideInfoOP cyc_info = entry->cyc_info_set[method];

		if ( cyc_info->is_cyclizable ) {  // NOTE : Yuval added this if; should be okay
			// NOTE : Yuval changed output so that it's unambiguous, in case only one of the cyclization methods is successful
			(*out_p_) << " " << CYCLIZATION_METHOD_NAMES[method]
				<< (cyc_info->cyc_comment.length()>0? " " : "") << cyc_info->cyc_comment
				<< (cyc_info->was_cyclic_model_created? (boost::format(" %1$.5f") % avoid_negative_zero(cyc_info->cyc_isc, 1e-3)).str() : "");
		}

	}
	(*out_p_) << std::endl;
}
void PeptideDeriverBasicStreamOutputter::end_receptor_partner_pair() {
	(*out_p_) << prefix_ << "# end chain pair" << std::endl;
}
void PeptideDeriverBasicStreamOutputter::end_structure() {
	(*out_p_) << prefix_ << "# end structure" << std::endl;
}

} //protocols
} //peptide_deriver






