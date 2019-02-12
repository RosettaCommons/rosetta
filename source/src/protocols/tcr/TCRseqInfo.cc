// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/TCRseqInfo.cc
/// @brief Class for parsing and numbering tcr a/b chain sequences
/// @author Ragul Gowthaman (ragul@umd.edu)

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

// Basic headers
#include <core/types.hh>
#include <basic/Tracer.hh>
// Protocol includes
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/util.hh>
// Option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/tcrmodel.OptionKeys.gen.hh>
// Utility Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

static basic::Tracer TR( "protocols.tcr.TCRseqInfo" );

using namespace core;
using namespace basic::options;

namespace protocols {
namespace tcr {

//TCRseqInfo::~TCRseqInfo() = default;

TCRseqInfo::TCRseqInfo( std::string const & aseq, std::string const & bseq ) {
	init();
	set_default();
	parse_tcrsegs(aseq,bseq);
}

TCRseqInfoOP
TCRseqInfo::clone() const{
	return utility::pointer::make_shared< TCRseqInfo >( * this);
}

TCRseqInfo::~TCRseqInfo() = default;

void TCRseqInfo::init() {
	cdr_posi_by_num_ = option[OptionKeys::tcrmodel::assign_cdr];
	if ( option[ OptionKeys::tcrmodel::assign_cdr1a ].user() ) {
		ausr_posi_.cdr1 = string_to_CDRbounds(option[OptionKeys::tcrmodel::assign_cdr1a]);
	}
	if ( option[ OptionKeys::tcrmodel::assign_cdr2a ].user() ) {
		ausr_posi_.cdr2hv4 = string_to_CDRbounds(option[OptionKeys::tcrmodel::assign_cdr2a]);
	}
	if ( option[ OptionKeys::tcrmodel::assign_cdr3a ].user() ) {
		ausr_posi_.cdr3 = string_to_CDRbounds(option[OptionKeys::tcrmodel::assign_cdr3a]);
	}
	if ( option[ OptionKeys::tcrmodel::assign_cdr1b ].user() ) {
		busr_posi_.cdr1 = string_to_CDRbounds(option[OptionKeys::tcrmodel::assign_cdr1b]);
	}
	if ( option[ OptionKeys::tcrmodel::assign_cdr2b ].user() ) {
		busr_posi_.cdr2hv4 = string_to_CDRbounds(option[OptionKeys::tcrmodel::assign_cdr2b]);
	}
	if ( option[ OptionKeys::tcrmodel::assign_cdr3b ].user() ) {
		busr_posi_.cdr3 = string_to_CDRbounds(option[OptionKeys::tcrmodel::assign_cdr3b]);
	}
	anarci_path_ = option[ OptionKeys::tcrmodel::anarci_path ];
	return;
}

void TCRseqInfo::set_default() {
	initialize_aho_numbers(aaho_posi_, baho_posi_);
	return;
}

void TCRseqInfo::parse_tcrsegs(std::string alpha_chain_sequence, std::string beta_chain_sequence){
	if ( cdr_posi_by_num_ ) {
		assign_CDRs_using_numbers(alpha_chain_sequence, ausr_posi_.cdr1, ausr_posi_.cdr2hv4, ausr_posi_.cdr3, aaho_posi_.cap, atcr_, aseq_posi_);
		assign_CDRs_using_numbers(beta_chain_sequence, busr_posi_.cdr1, busr_posi_.cdr2hv4, busr_posi_.cdr3, baho_posi_.cap, btcr_, bseq_posi_);
	} else {
		assign_achain_CDRs_using_REGEX(alpha_chain_sequence, atcr_, aseq_posi_);
		assign_bchain_CDRs_using_REGEX(beta_chain_sequence, btcr_, bseq_posi_);
	}
	if ( atcr_.truncdomain.empty() && !anarci_path_.empty() ) {
		assign_CDRs_using_anarci(alpha_chain_sequence, anarci_path_, aaho_posi_, aseq_posi_, atcr_);
	}
	if ( btcr_.truncdomain.empty() && !anarci_path_.empty() ) {
		assign_CDRs_using_anarci(beta_chain_sequence, anarci_path_, baho_posi_, bseq_posi_, btcr_);
	}
	adjust_position_for_chain(bseq_posi_, btcr_.truncdomain.length());
	return;
}

} // namespace tcr
} // namespace protocols


#endif // __ANTIBODY_GRAFTING__
