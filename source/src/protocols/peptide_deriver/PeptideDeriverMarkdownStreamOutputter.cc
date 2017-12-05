// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverMarkdownStreamOutputter.cc
/// @brief outputs a Markdown formatted Peptiderive report to a stream
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @author Orly Marcu (orlymarcu@gmail.com)

// Unit headers
#include <protocols/peptide_deriver/PeptideDeriverMarkdownStreamOutputter.hh>
#include <protocols/peptide_deriver/PeptideDeriverFilter.hh>

// Project headers
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility headers
#include <utility/version.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/ocstream.hh>

// External headers
#include <boost/format.hpp>

// C++ headers
#include <string>

static basic::Tracer TR( "protocols.peptide_deriver.PeptideDeriverMarkdownStreamOutputter" );


namespace protocols {
namespace peptide_deriver {

// BEGIN PeptideDeriverMarkdownStreamOutputter implementation

PeptideDeriverMarkdownStreamOutputter::PeptideDeriverMarkdownStreamOutputter(utility::io::orstream & out, std::string prefix) {
	out_p_ = &out;
	prefix_ = prefix;
	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		cyc_report_info_set_[method] = CyclizedReportInfoOP( new CyclizedReportInfo() );
	}
}

PeptideDeriverMarkdownStreamOutputter::PeptideDeriverMarkdownStreamOutputter( PeptideDeriverMarkdownStreamOutputter const & src) {
	out_p_=src.out_p_;
	prefix_=src.prefix_;
	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		cyc_report_info_set_[method] = src.cyc_report_info_set_[method];
	}
}

PeptideDeriverMarkdownStreamOutputter::~PeptideDeriverMarkdownStreamOutputter(){}



PeptideDeriverMarkdownStreamOutputterOP
PeptideDeriverMarkdownStreamOutputter::clone() const {
	return PeptideDeriverMarkdownStreamOutputterOP( new PeptideDeriverMarkdownStreamOutputter( *this ) );
}




void PeptideDeriverMarkdownStreamOutputter::clear_buffers() {
	header_.str("");
	best_linear_peptides_.str("");

	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		cyc_report_info_set_[method]->best_cyclic_peptides_.str("");
		cyc_report_info_set_[method]->cyclic_peptides_.str("");
	}
	all_peptides_.str("");
	footer_.str("");
}

void PeptideDeriverMarkdownStreamOutputter::begin_structure(core::pose::Pose const & pose, std::string const &name) {
	clear_buffers();

	header_ << prefix_ << "# Peptiderive report file " << name << std::endl
		<< prefix_ << std::endl
		<< prefix_ << "- Chains: " << pose.conformation().num_chains() << std::endl
		<< prefix_ << "- Rosetta version: " << utility::Version::package() << " " << utility::Version::version() << " commit " << utility::Version::commit() << " from " << utility::Version::url() << std::endl
		<< prefix_ << "- (*) in the 'Relative interface score (%)' column means positive values were calculated for these entries, indicating unfavorable interactions" << std::endl
		<< prefix_ << "- (**) in the 'Cyclized interface score' column means that a cyclized model was not constructed for this cyclizable peptide, since its energy contribution (in its linear form) was not significant" << std::endl
		<< prefix_ << "- Disulfide-cyclized peptide models have additional N- and C-terminal cysteine residues, not shown in the 'Sequence' column" << std::endl;



	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		cyc_report_info_set_[method]->best_cyclic_peptides_ << prefix_ << "## Best " << CYCLIZATION_METHOD_NAMES[method] << " cyclizable peptides for all chain pairs" << std::endl
			<< prefix_ << std::endl
			<< prefix_ << "| Receptor | Partner | Peptide length | Position | Interface score | Relative interface score (%) | Extra info     | Cyclized interface score | Sequence             |" << std::endl
			<< prefix_ << "|----------|---------|----------------|----------|-----------------|------------------------------|----------------|--------------------------|----------------------|" << std::endl;

		cyc_report_info_set_[method]->cyclic_peptides_ << prefix_ << "## All " << CYCLIZATION_METHOD_NAMES[method] << " cyclizable peptides" << std::endl
			<< prefix_ << std::endl;
	}

	best_linear_peptides_ << prefix_ << "## Best linear peptides for all chain pairs" << std::endl
		<< prefix_ << std::endl
		<< prefix_ << "| Receptor | Partner | Peptide length | Position | Interface score | Relative interface score (%) | Sequence             |" << std::endl
		<< prefix_ << "|----------|---------|----------------|----------|-----------------|------------------------------|----------------------|" << std::endl;

	all_peptides_ << prefix_ << "## All linear peptides" << std::endl
		<< prefix_ << std::endl;


}

void PeptideDeriverMarkdownStreamOutputter::end_structure() {
	footer_ << prefix_ << "*end of report*";

	( * out_p_ ) << header_.str() << std::endl << std::endl
		<< best_linear_peptides_.str() << std::endl << std::endl;
	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		( * out_p_ ) << cyc_report_info_set_[method]->best_cyclic_peptides_.str() << std::endl;
	}
	( * out_p_ ) << std::endl << all_peptides_.str() << std::endl;

	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		( * out_p_ ) << cyc_report_info_set_[method]->cyclic_peptides_.str() << std::endl;
	}
	( * out_p_ ) << std::endl << footer_.str() << std::endl;

	clear_buffers();
}

void PeptideDeriverMarkdownStreamOutputter::begin_receptor_partner_pair(char const receptor_chain_letter,
	char const partner_chain_letter, core::Real const total_isc,
	std::string const & /*options_string*/) {

	current_receptor_chain_letter_ = receptor_chain_letter;
	current_partner_chain_letter_ = partner_chain_letter;
	current_total_isc_ = total_isc;

	// all_peptides_ << "- Options: " << options_string << std::endl;
	// cyclic_peptides_ << std::endl << "- Options: " << options_string << std::endl;
}

void PeptideDeriverMarkdownStreamOutputter::peptide_length(core::Size const pep_length) {
	current_pep_length_ = pep_length;
	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		cyc_report_info_set_[method]->cyclic_peptide_encountered_for_current_pep_length = false;
	}
	all_peptides_ << prefix_ << ( boost::format("### Receptor= %1% Partner= %2% Peptide_length= %3%")
		% current_receptor_chain_letter_
		% current_partner_chain_letter_
		% current_pep_length_ ) << std::endl
		<< prefix_ << "- Total interface score: " <<  (boost::format("%1$-7.3f") % avoid_negative_zero(current_total_isc_, 1e-3) ) << std::endl
		<< prefix_ << std::endl
		<< prefix_ << "| Position | Interface score | Relative interface score (%) |"  << std::endl
		<< prefix_ << "|----------|-----------------|------------------------------|" << std::endl;

}


void PeptideDeriverMarkdownStreamOutputter::peptide_entry(PeptideDeriverEntryType entry_type, core::Real const total_isc, DerivedPeptideEntryCOP entry) {
	// TODO : maybe still consider calculating binding_contribution_fraction extenally from these function calls?
	core::Real binding_contribution_fraction = entry->lin_isc / total_isc;

	std::string const binding_contribution_pct_str = binding_contribution_fraction >0 ? (boost::format("%1$.2f") % (100*binding_contribution_fraction)).str() : "*";

	core::Real no_neg_zero_linear_isc = avoid_negative_zero(entry->lin_isc, 1e-3); // this should match precision of displayed value

	core::Size const PEPTIDE_CHAIN = 2;

	switch (entry_type) {
	case ET_BEST_LINEAR :
		best_linear_peptides_ << prefix_ << ( boost::format( "| %1$-8c | %2$-7c | %3$-14d | %4$-8d | %5$-15.3f | %6$-28s | %7$-20s |" )
			% current_receptor_chain_letter_
			% current_partner_chain_letter_
			% current_pep_length_
			% entry->pep_start
			% no_neg_zero_linear_isc
			% binding_contribution_pct_str
			% entry->lin_pose->chain_sequence(PEPTIDE_CHAIN) ) << std::endl;

		for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
			CyclizedPeptideInfoCOP cyc_info = entry->cyc_info_set[method];
			CyclizedReportInfoOP cyc_report_info = this->cyc_report_info_set_[method];
			std::string const cyclic_isc_string = (cyc_info->was_cyclic_model_created? (boost::format("%1$.3f") % avoid_negative_zero(cyc_info->cyc_isc, 1e-3)).str() : "**");

			// if the best linear is not the best cyclic but it can also be cyclized, we want to print this information in the best cyclic peptides table
			// TODO : assumes entry->lin_pose.chain_sequence(n) == entry->cyc_info_set[X]->cyc_pose.chain_sequence(n) for all n.
			if ( cyc_info->was_cyclic_model_created && entry->pep_start != cyc_report_info->best_cyclic_pep_start) {
				cyc_report_info->best_cyclic_peptides_ << prefix_ << ( boost::format( "| %1$-8c | %2$-7c | %3$-14d | %4$-8d | %5$-15.3f | %6$-24s | %7$-14s | %8$-24s | %9$-20s |" )
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_pep_length_
					% entry->pep_start
					% no_neg_zero_linear_isc
					% binding_contribution_pct_str
					% cyc_info->cyc_comment
					% cyclic_isc_string
					% entry->lin_pose->chain_sequence(PEPTIDE_CHAIN) ) << std::endl;
			}
		}// for method
		// the best linear peptide signals the end of a run for a certain chain pair in certain roles
		// add a separator at the end of the table for the 'general' peptides
		all_peptides_ << prefix_ << "---" << std::endl;
		break;

	case ET_BEST_CYCLIC :
		for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
			CyclizedPeptideInfoCOP cyc_info = entry->cyc_info_set[method];
			CyclizedReportInfoOP cyc_report_info = this->cyc_report_info_set_[method];
			if ( cyc_info->is_cyclizable ) {
				cyc_report_info->best_cyclic_pep_start = entry->pep_start;
				std::string const cyclic_isc_string = (cyc_info->was_cyclic_model_created? (boost::format("%1$.3f") % avoid_negative_zero(cyc_info->cyc_isc, 1e-3)).str() : "**");
				cyc_report_info->best_cyclic_peptides_ << prefix_ << ( boost::format( "| %1$-8c | %2$-7c | %3$-14d | %4$-8d | %5$-15.3f | %6$-28s | %7$-14s | %8$-24s | %9$-20s |" )
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_pep_length_
					% entry->pep_start
					% no_neg_zero_linear_isc
					% binding_contribution_pct_str
					% cyc_info->cyc_comment
					% cyclic_isc_string
					% cyc_info->cyc_pose->chain_sequence(PEPTIDE_CHAIN) ) << std::endl;

				cyc_report_info->best_cyclic_pep_start = entry->pep_start;
				cyc_report_info->cyclic_peptides_ << prefix_ << "----" << std::endl << std::endl;
			}
		}
		break;

	case ET_GENERAL :
		all_peptides_ << prefix_ << ( boost::format( "| %1$-8d | %2$-15.3f | %3$-28s |" )
			% entry->pep_start
			% no_neg_zero_linear_isc
			% binding_contribution_pct_str ) << std::endl;

		for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
			CyclizedPeptideInfoCOP cyc_info = entry->cyc_info_set[method];
			CyclizedReportInfoOP cyc_report_info = this->cyc_report_info_set_[method];

			std::string const cyclic_isc_string = (cyc_info->was_cyclic_model_created? (boost::format("%1$.3f") % avoid_negative_zero(cyc_info->cyc_isc, 1e-3)).str() : "**");

			if ( cyc_info->cyc_comment != "" ) {
				if ( !cyc_report_info->cyclic_peptide_encountered_for_current_pep_length ) {
					cyc_report_info->cyclic_peptides_ << prefix_ << ( boost::format("### Receptor= %1% Partner= %2% Peptide_length= %3%")
						% current_receptor_chain_letter_
						% current_partner_chain_letter_
						% current_pep_length_ ) << std::endl
						<< prefix_ << "- Total interface score: " << current_total_isc_ << std::endl << std::endl
						<< prefix_ << "| Position | Interface score | Relative interface score (%) | Cyclization info | Cyclized interface score |" << std::endl
						<< prefix_ << "|----------|-----------------|------------------------------|------------------|--------------------------|" << std::endl;
					cyc_report_info->cyclic_peptide_encountered_for_current_pep_length = true;
				}
				cyc_report_info->cyclic_peptides_ << prefix_ << ( boost::format( "| %1$-8d | %2$-15.3f | %3$-28s | %4$-16s | %5$-24s |" )
					% entry->pep_start
					% no_neg_zero_linear_isc
					% binding_contribution_pct_str
					% cyc_info->cyc_comment
					% cyclic_isc_string) << std::endl;
			}
		}

		break;

	default : // should never happen; means code wasn't written well
		utility_exit_with_message("Invalid enum value provided for entry_type");
		break;
	} // switch
}

void PeptideDeriverMarkdownStreamOutputter::end_receptor_partner_pair() {
	(*out_p_).flush();
}



} //protocols
} //peptide_deriver






