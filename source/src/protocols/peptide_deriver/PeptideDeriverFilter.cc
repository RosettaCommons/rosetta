// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  A filter that, for each dimer in a pose, outputs the peptide which contributes most to the interface.

/// @author Nir London
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @author Orly Marcu (orly.marcu@mail.huji.ac.il)
/// @date   Jul. 30, 2014

// TODO : consider whether JD2 has infrastructure to replace PeptideDeriverJobOutputter

// PeptideDeriverFilter


// Unit headers
#include <protocols/peptide_deriver/PeptideDeriverFilter.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/io/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/svn_version.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/util/disulfide_util.hh>
#include <numeric/xyzVector.hh>
#include <protocols/cyclic_peptide/PeptideCyclizeMover.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/jd2/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/DisulfideInsertionMover.hh>
#include <protocols/relax/AtomCoordinateCstMover.hh>



//Utility headers
#include <basic/MetricValue.hh>
#include <basic/options/keys/peptide_deriver.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility>
#include <utility/excn/Exceptions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/ocstream.hh>
#include <utility/tag/Tag.hh>

// RosettaScripts headers
#include <protocols/rosetta_scripts/util.hh>

// External headers
#include <boost/format.hpp>

// C++ headers
#include <cassert>
#include <string>
#include <limits>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
#include <protocols/peptide_deriver/PeptideDeriverFilterCreator.hh>

namespace protocols {
namespace peptide_deriver {

static THREAD_LOCAL basic::Tracer tracer( "protocols.peptide_deriver.PeptideDeriverFilter" );

// helper function: avoid +- signed on zero-values
core::Real avoid_negative_zero(core::Real const value, core::Real const threshold) {
	return ((value < 0) && (-value < threshold))? 0 : value;
}

// BEGIN DerivedPeptideEntry implementation
DerivedPeptideEntry::DerivedPeptideEntry(core::Real const lin_isc, core::Size const pep_start, core::pose::PoseCOP lin_pose)
{
	this->lin_isc = lin_isc;
	this->pep_start = pep_start;
	this->lin_pose = lin_pose;
	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		this->cyc_info_set[method] = CyclizedPeptideInfoOP(new CyclizedPeptideInfo(lin_pose));
	}
}


// END DerivedPeptideEntry implementation


// BEGIN CyclizedPeptideInfo implementation

CyclizedPeptideInfo::CyclizedPeptideInfo(core::pose::PoseCOP reference_pose) {
	// NOTE : this is a dummy value so that we don't dereference NULL; see below
	// avoid null pointers
	// TODO : solve the need for a dummy value for receptor_cyclic_peptide_pose
	// NOTE : receptor_peptide_pose is used here as dummy, since we can't pass NULL as reference

	this->pre_cyc_pose = reference_pose;
	this->cyc_pose = reference_pose;
	this->is_cyclizable = false;
	this->was_cyclic_model_created = false;
	this->cyc_comment = "";
}
// END CyclizedPeptideInfo implementation

// BEGIN PeptideDeriverOutputterContainer implementation

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

/*
void PeptideDeriverOutputterContainer::peptide_entry(core::pose::Pose const & pose, PeptideDeriverEntryType const entry_type, core::Size const pep_start,
core::Real const linear_isc, core::Real const binding_contribution_fraction, std::string const & disulfide_info, bool const was_cyclic_pep_modeled, core::pose::Pose const & cyclic_pose, core::Real const cyclic_isc) {
for ( PeptideDeriverOutputterOP const outputter : list_ ) {
outputter->peptide_entry(pose, entry_type, pep_start, linear_isc, binding_contribution_fraction, disulfide_info, was_cyclic_pep_modeled, cyclic_pose, cyclic_isc);
}
}
*/
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

// END PeptideDeriverOutputterContainer implementation


// BEGIN PeptideDeriverBasicStreamOutputter implementation

PeptideDeriverBasicStreamOutputter::PeptideDeriverBasicStreamOutputter(utility::io::orstream & out, std::string prefix) {
	out_p_ = &out;
	prefix_ = prefix;

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

// END PeptideDeriverBasicStreamOutputter implementation

// BEGIN PeptideDeriverMarkdownStreamOutputter implementation

PeptideDeriverMarkdownStreamOutputter::PeptideDeriverMarkdownStreamOutputter(utility::io::orstream & out, std::string prefix) {
	out_p_ = &out;
	prefix_ = prefix;
	for ( core::Size method = 0; method < NUM_CYCLIZATION_METHODS; ++method ) {
		cyc_report_info_set_[method] = CyclizedReportInfoOP( new CyclizedReportInfo() );
	}
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
		<< prefix_ << "- Rosetta version: " << core::minirosetta_svn_version() << " from " << core::minirosetta_svn_url() << std::endl
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
			if ( cyc_info->was_cyclic_model_created && entry->pep_start != cyc_report_info->best_cyclic_pep_start ) {
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

// END PeptideDeriverMarkdownStreamOutputter implementation


// BEGIN PeptideDeriverPoseOutputter implementation

void PeptideDeriverPoseOutputter::output_pose( core::pose::Pose & pose, std::string const & pose_name ) {
	(*scorefxn_)(pose);
	if ( protocols::jd2::jd2_used() ) {
		protocols::jd2::output_intermediate_pose(pose, pose_name );
	} else {
		std::string const file_name(pose_name + ".pdb");
		utility::io::ozstream out_stream( file_name );
		core::io::pdb::dump_pdb( pose, out_stream );
	}
}

PeptideDeriverPoseOutputter::PeptideDeriverPoseOutputter( bool const is_dump_best_peptide_pose,
	bool const is_dump_prepared_pose, bool const is_dump_cyclic_poses,
	core::scoring::ScoreFunctionCOP scorefxn) :
	is_dump_best_peptide_pose_(is_dump_best_peptide_pose),
	is_dump_prepared_pose_(is_dump_prepared_pose),
	is_dump_cyclic_poses_(is_dump_cyclic_poses),
	scorefxn_(std::move(scorefxn)) { }

void PeptideDeriverPoseOutputter::chain_pair_pose_prepared(core::pose::Pose const & pose) {
	// note: we use the term 'prepared' and 'minimized' interchangeably
	//       anticipating that this step might involve other things
	//       in the meantime, it makes sense to talk about the complex just before
	//       peptiderive goes over it

	// save a copy of the prepared pose aside, but we output only on begin_receptor_partner_pair()
	// since we want to include chain letter in output, and we don't get them from here
	current_chain_pair_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
	is_chain_pair_new_ = true;
}

void PeptideDeriverPoseOutputter::begin_receptor_partner_pair(char const receptor_chain_letter,
	char const partner_chain_letter, core::Real const,
	std::string const &) {
	current_receptor_chain_letter_ = receptor_chain_letter;
	current_partner_chain_letter_ = partner_chain_letter;

	// NOTE : this assumes the PeptideDeriver uses the same prepared
	//        structure when switching roles (receptor and partner)
	//        between partners in a chain pair, in case both roles are
	//        indeed evaluated.
	if ( is_dump_prepared_pose_ && is_chain_pair_new_ ) {
		std::string pose_name = ( boost::format("%1%%2%")
			% ( (current_receptor_chain_letter_ < current_partner_chain_letter_ )? current_receptor_chain_letter_ : current_partner_chain_letter_  )
			% ( (current_receptor_chain_letter_ < current_partner_chain_letter_ )? current_partner_chain_letter_  : current_receptor_chain_letter_ ) ).str();

		output_pose( *current_chain_pair_pose_, pose_name );
		is_chain_pair_new_ = false;
	}
}

void PeptideDeriverPoseOutputter::peptide_length(core::Size const pep_length) {
	current_peptide_length_ = pep_length;
}

void PeptideDeriverPoseOutputter::peptide_entry(PeptideDeriverEntryType const entry_type, core::Real const total_isc, DerivedPeptideEntryCOP entry) {

	tracer.Debug << "structure total_isc: " << total_isc << std::endl;

	// handles best linear peptide (first if statement) and all possible cyclic peptides generated for this peptide
	if ( is_dump_best_peptide_pose_ && entry_type == ET_BEST_LINEAR ) {
		std::string const lin_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_linear_linear_peptide")
			% current_receptor_chain_letter_
			% current_partner_chain_letter_
			% current_peptide_length_).str();
		core::pose::Pose temp_pose(*(entry->lin_pose));
		output_pose(temp_pose, lin_pose_name);

		for ( core::Size method = 0; method<NUM_CYCLIZATION_METHODS; ++method ) {
			if ( entry->cyc_info_set[method]->was_cyclic_model_created ) {
				std::string const cyc_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_linear_cyclic_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_).str();
				core::pose::Pose temp_cyclic_pose(*(entry->cyc_info_set[method]->cyc_pose));
				output_pose(temp_cyclic_pose, cyc_pose_name);
			}
		}
	} else {
		for ( core::Size method = 0; method<NUM_CYCLIZATION_METHODS; ++method ) {
			// handles best cyclic peptide (by any one of the cyclization methods) and its linear / other possible cyclization variants
			if ( is_dump_best_peptide_pose_ && entry_type == ET_BEST_CYCLIC && entry->cyc_info_set[method]->was_cyclic_model_created ) {
				std::string const lin_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_cyclic_linear_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_).str();
				core::pose::Pose temp_pose(*(entry->cyc_info_set[method]->pre_cyc_pose));
				output_pose(temp_pose, lin_pose_name);

				std::string const cyc_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_cyclic_cyclic_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_).str();
				core::pose::Pose temp_cyclic_pose(*(entry->cyc_info_set[method]->cyc_pose));
				output_pose(temp_cyclic_pose, cyc_pose_name);
			} else if ( is_dump_cyclic_poses_ && entry_type == ET_GENERAL && entry->cyc_info_set[method]->was_cyclic_model_created ) {
				// handles all non-best cyclic peptides that did meet the criteria to be generated
				std::string const pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_%4%_cyclic_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_
					% entry->cyc_info_set[method]->cyc_comment ).str();
				core::pose::Pose temp_pose(*(entry->cyc_info_set[method]->cyc_pose));
				output_pose(temp_pose, pose_name);
			}

		}
	}

}

// END PeptideDeriverPoseOutputter implementation


// -------------------------------------------------------------------------------- //
// ---                                                                          --- //
// ---                 PeptideDeriverFilter implementation                      --- //
// ---                                                                          --- //
// -------------------------------------------------------------------------------- //

// NOTE : based on LoopMover version
PeptideDeriverFilter & PeptideDeriverFilter::operator=( PeptideDeriverFilter const & rhs ) {
	if ( this == &rhs ) return *this;

	// NOTE : if Filter will have an operator=, call it:
	// Filter::operator=(rhs);
	PeptideDeriverFilter::assign(*this, rhs);
	return *this;
}

void PeptideDeriverFilter::parse_options() {
	// read user-defined options

	// the default value for pep_lengths is in options_rosetta.py
	// NOTE : the leading :: is to disambiguate peptide_deriver from the protocols:: namespace with the top-level one
	set_pep_lengths(basic::options::option[basic::options::OptionKeys::peptide_deriver::pep_lengths]());

	set_is_skip_zero_isc(basic::options::option[basic::options::OptionKeys::peptide_deriver::skip_zero_isc]());
	set_is_dump_peptide_pose(basic::options::option[basic::options::OptionKeys::peptide_deriver::dump_peptide_pose]());
	set_is_dump_report_file(basic::options::option[basic::options::OptionKeys::peptide_deriver::dump_report_file]());
	set_is_report_gzip(basic::options::option[basic::options::OptionKeys::peptide_deriver::report_gzip]());
	set_is_dump_prepared_pose(basic::options::option[basic::options::OptionKeys::peptide_deriver::dump_prepared_pose]());
	set_is_dump_cyclic_poses(basic::options::option[basic::options::OptionKeys::peptide_deriver::dump_cyclic_poses]());
	set_is_do_minimize(basic::options::option[basic::options::OptionKeys::peptide_deriver::do_minimize]());
	set_optimize_cyclic_threshold(basic::options::option[basic::options::OptionKeys::peptide_deriver::optimize_cyclic_threshold]());
	set_report_format( PeptideDeriverFilter::parse_report_format_string(basic::options::option[basic::options::OptionKeys::peptide_deriver::report_format]()) );

	utility::vector1<char> restrict_receptors_to_chains;

	for ( std::string const & chain_string : basic::options::option[basic::options::OptionKeys::peptide_deriver::restrict_receptors_to_chains]() ) {
		assert(chain_string.size() == 1);
		restrict_receptors_to_chains.push_back(chain_string[0]);

	}

	utility::vector1<char> restrict_partners_to_chains;
	for ( std::string const & chain_string : basic::options::option[basic::options::OptionKeys::peptide_deriver::restrict_partners_to_chains]() ) {
		assert(chain_string.size() == 1);
		restrict_partners_to_chains.push_back(chain_string[0]);
	}

	set_restrict_receptors_to_chains(restrict_receptors_to_chains);
	set_restrict_partners_to_chains(restrict_partners_to_chains);
}

PeptideDeriverFilter::PeptideDeriverFilter() :
	Filter("PeptideDeriverFilter") {
	// this is currently redundant, since we always read the options (parse_options())
	// and this is the default specified in options_rosetta.py, but I left it here
	// so that if options are moved out of the ctor, it would be clear that it's probably
	// a good idea to do this!
	pep_lengths_.push_back(10);

	parse_options();

	// the deriver uses the default score function
	scorefxn_deriver_ = core::scoring::get_score_function();

	// the minimizer uses the default score function, but with a full weight
	// on coordinate constraints (see PeptideDeriverFilter::minimize())
	scorefxn_minimizer_ = core::scoring::get_score_function();
	scorefxn_minimizer_->set_weight(core::scoring::coordinate_constraint, 1.0);
}

void
PeptideDeriverFilter::assign(PeptideDeriverFilter & lhs, PeptideDeriverFilter const & rhs) {
	lhs.scorefxn_deriver_          = rhs.scorefxn_deriver_         ;
	lhs.scorefxn_minimizer_        = rhs.scorefxn_minimizer_       ;
	lhs.pep_lengths_               = rhs.pep_lengths_              ;
	lhs.is_skip_zero_isc_          = rhs.is_skip_zero_isc_         ;
	lhs.is_dump_peptide_pose_      = rhs.is_dump_peptide_pose_     ;
	lhs.is_dump_report_file_       = rhs.is_dump_report_file_      ;
	lhs.is_report_gzip_            = rhs.is_report_gzip_           ;
	lhs.is_dump_prepared_pose_     = rhs.is_dump_prepared_pose_    ;
	lhs.is_dump_cyclic_poses_      = rhs.is_dump_cyclic_poses_     ;
	lhs.is_do_minimize_            = rhs.is_do_minimize_           ;
	lhs.optimize_cyclic_threshold_ = rhs.optimize_cyclic_threshold_;
	lhs.report_format_             = rhs.report_format_            ;
}

PeptideDeriverFilter::PeptideDeriverFilter(PeptideDeriverFilter const & rhs) :
	Filter(rhs) {
	PeptideDeriverFilter::assign(*this, rhs);
}

/// Minimize the pose, while using native constraints.
void
PeptideDeriverFilter::minimize(core::pose::Pose & pose) const {

	// set native constraints prior to minimization

	core::scoring::constraints::ConstraintSetOP cst_set(
		new core::scoring::constraints::ConstraintSet());
	core::scoring::func::HarmonicFuncOP spring(
		new core::scoring::func::HarmonicFunc(
		0, // mean
		1  // std. dev
		));
	core::conformation::Conformation const & conformation(pose.conformation());
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! pose.residue(i).is_protein() ) {
			continue;
		}
		core::id::AtomID CAi( pose.residue(i).atom_index("CA"), i );
		cst_set->add_constraint( core::scoring::constraints::ConstraintCOP(
			new core::scoring::constraints::CoordinateConstraint( CAi, CAi,
			conformation.xyz( CAi ), spring ) ) );
	}
	pose.constraint_set( cst_set );

	// minimize bb and sc
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	movemap->set_bb( true );
	movemap->set_chi( true );

	// TODO : consider using whatever is mentioned in -run:min_type and -run:min_tolerance
	//        and maybe change the defaults to this for Peptiderive app

	protocols::simple_moves::MinMover minimizer( movemap, scorefxn_minimizer_,
		"dfpmin_armijo_atol", 0.01, true /*nb_list*/ );

	// minimize pose
	minimizer.apply( pose );
	pose.remove_constraints();

}

// this filter doesn't do any filtering, only reporting
bool
PeptideDeriverFilter::apply(
	core::pose::Pose const & // pose - unused
) const {

	// we do nothing here
	// the interesting stuff happens in PeptideDeriverFilter::report()

	return true;
}

PeptideDeriverReportFormat
PeptideDeriverFilter::parse_report_format_string(
	std::string const & value) {
	if ( value == "markdown" ) {
		return PRF_MARKDOWN;
	}
	if ( value == "basic" ) {
		return PRF_BASIC;
	}
	throw utility::excn::EXCN_KeyError( "PeptideDeriverFilter::parse_report_format_string(std::string const value): value must be either 'markdown' or 'basic'." );
}

/// Verifies the pose has at least two chains. Minimizes the pose and detects
/// disulfides.
///
/// @param pose the pose to prepare.
void
PeptideDeriverFilter::prepare_pose(
	PeptideDeriverOutputter & ,
	core::pose::Pose & pose) const {
	// check number of chains
	if ( pose.conformation().num_chains() < 2 ) {
		throw utility::excn::EXCN_Msg_Exception(
			"Input PDB contains less than 2 chains.");
	}

	if ( is_do_minimize_ ) {
		// minimize pose
		tracer << "Start: minimization" << std::endl;
		minimize(pose);
		tracer << "End: minimization" << std::endl;
	}

	// detect SS
	tracer << "Start: detect disulfides" << std::endl;
	pose.conformation().detect_disulfides();
	tracer << "End: detect disulfides" << std::endl;

	// modify termini - add hydrogens where needed
	core::Size chain1_start = pose.conformation().chain_begin(1);
	core::Size chain1_end = pose.conformation().chain_end(1);
	core::Size chain2_start = pose.conformation().chain_begin(2);
	core::Size chain2_end = pose.conformation().chain_end(2);

	if ( pose.residue(chain1_start).has_variant_type(core::chemical::LOWERTERM_TRUNC_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::LOWERTERM_TRUNC_VARIANT, chain1_start );
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::LOWER_TERMINUS_VARIANT, chain1_start );
	}

	if ( pose.residue(chain1_end).has_variant_type(core::chemical::UPPERTERM_TRUNC_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::UPPERTERM_TRUNC_VARIANT, chain1_end );
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::UPPER_TERMINUS_VARIANT, chain1_end );
	}

	if ( pose.residue(chain2_start).has_variant_type(core::chemical::LOWERTERM_TRUNC_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::LOWERTERM_TRUNC_VARIANT, chain2_start );
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::LOWER_TERMINUS_VARIANT, chain2_start );
	}

	if ( pose.residue(chain2_end).has_variant_type(core::chemical::UPPERTERM_TRUNC_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::UPPERTERM_TRUNC_VARIANT, chain2_end );
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::UPPER_TERMINUS_VARIANT, chain2_end );
	}
}

// the report() method is where PeptideDeriver analyzes the pose
// and spits out information about derived peptides
void
PeptideDeriverFilter::report(std::ostream & out, core::pose::Pose const & orig_pose) const {

	tracer << "Start PeptideDeriver report" << std::endl;

	// determine what to output, where and how
	utility::io::ozstream file_out;
	utility::io::ocstream log_out( out );

	std::string const pose_name = ( protocols::jd2::jd2_used()? protocols::jd2::current_output_name() : "" );

	//std::string prefix;
	if ( is_dump_report_file_ ) {
		std::string path = ( basic::options::option[ basic::options::OptionKeys::out::path::all ].user()? basic::options::option[ basic::options::OptionKeys::out::path::all ]().path() : "" );
		std::string file_name = path + ( pose_name.length()? pose_name + "." : "" ) + "peptiderive.txt";
		if ( is_report_gzip_ ) {
			file_name += ".gz";
		}
		tracer << "Setting up report file output to " << file_name << std::endl;
		file_out.open( file_name );
	}

	PeptideDeriverOutputterContainer output;

	utility::io::orstream & output_stream = ( is_dump_report_file_?
		(utility::io::orstream &)file_out :
		(utility::io::orstream &)log_out );

	// this is simpler than creating a full-fledged factory. Sorry, design pattern people. --yuvals
	switch ( report_format_ ) {
	case PRF_MARKDOWN :
		output.push_back( PeptideDeriverOutputterOP( new PeptideDeriverMarkdownStreamOutputter( output_stream , /*
			prefix=*/"" ) ) );
		break;
	case PRF_BASIC :
		output.push_back( PeptideDeriverOutputterOP( new PeptideDeriverBasicStreamOutputter( output_stream , /*prefix=*/"" ) ) );
		break;
	}


	if ( is_dump_peptide_pose_ || is_dump_prepared_pose_ || is_dump_cyclic_poses_ ) {
		tracer << "Setting up pose output ("
			<< "peptides=" << is_dump_peptide_pose_ << ","
			<< "prepared=" << is_dump_prepared_pose_ << ","
			<< "cyclic=" << is_dump_cyclic_poses_ << ")"
			<< std::endl;
		output.push_back( PeptideDeriverOutputterOP( new PeptideDeriverPoseOutputter( is_dump_peptide_pose_,
			is_dump_prepared_pose_, is_dump_cyclic_poses_, scorefxn_deriver_ ) ) );
	}

	// prepare the PDB

	tracer << "Begin structure " << pose_name << std::endl;
	output.begin_structure(orig_pose, pose_name);

	// derive

	// if user didn't specify receptor/partner chains, we consider all chains as possible
	// receptor/partner chains (respectively).
	utility::vector1<core::Size> receptor_chain_indices = get_chain_indices(orig_pose, restrict_receptors_to_chains_);
	utility::vector1<core::Size> partner_chain_indices = get_chain_indices(orig_pose, restrict_partners_to_chains_);

	// since we have the 'both_ways' argument in PeptideDeriverFilter::derive_peptide()
	// we do some effort here to call the derivation function as little as possible.
	// as there is some set up done while preparing the two-monomer pose
	// so this figures out which poses can be run with both_ways = true, and which can't
	utility::vector1<core::Size> already_done_both_ways;

	// go through all requested chains
	for ( core::Size const receptor_chain : receptor_chain_indices ) {

		bool is_receptor_also_partner = partner_chain_indices.contains(receptor_chain);

		for ( core::Size const partner_chain : partner_chain_indices ) {
			if ( receptor_chain == partner_chain ) {
				continue;
			}

			// when both chains should be both receptors and partners, run only once with both_ways = true
			bool both_ways = (is_receptor_also_partner && receptor_chain_indices.contains(partner_chain));

			if ( both_ways ) {
				// check we've not been here already with the roles switched (see below, where already_done_both_ways is populated)
				if ( !already_done_both_ways.contains(partner_chain) ) {
					derive_peptide( output, orig_pose, receptor_chain, partner_chain, /* both_ways = */ true);
				}
			} else {
				derive_peptide( output, orig_pose, receptor_chain, partner_chain, /* both_ways = */ false );
			}
		}

		// avoid deriving from this chain when the time comes and other chains ask for it to be the partner
		if ( is_receptor_also_partner ) {
			already_done_both_ways.push_back(receptor_chain);
		}
	}


	tracer << "End structure " << pose_name << std::endl;
	output.end_structure();


}

utility::vector1<core::Size>
PeptideDeriverFilter::get_chain_indices(core::pose::Pose const & pose, utility::vector1<char> const & restrict_to_chains) {
	utility::vector1<core::Size> chain_indices;

	if ( restrict_to_chains.empty() ) {
		// no restriction mean return all chains
		core::Size chain_count = pose.conformation().num_chains();
		for ( core::Size i = 1; i <= chain_count; ++i ) {
			chain_indices.push_back(i);
		}
	} else {
		for ( char const chain_char : restrict_to_chains ) {
			chain_indices.push_back( core::pose::get_chain_id_from_chain( chain_char, pose ) );
		}
	}
	return chain_indices;
}


/// Pushes the residue indices belonging to the specified chain index into a given vector.
void
PeptideDeriverFilter::push_chain_residue_indices(core::pose::Pose const & many_chain_pose,
	core::Size const chain_index, utility::vector1<core::Size> & residue_indices ) {

	core::Size chain_begin = many_chain_pose.conformation().chain_begin( chain_index );
	core::Size chain_end = many_chain_pose.conformation().chain_end( chain_index );

	for ( core::Size i = chain_begin; i <= chain_end; ++i ) {
		residue_indices.push_back( i );
	}

}

/// cuts the specified chains from the original pose using core::io::pdb::pose_from_pose()
/// and calculates the complex (receptor and partner)
///
/// @param pose a many-chain pose
/// @param first_chain_index chain index of the receptor (and also partner, if both_ways is true)
/// @param second_chain_index chain index of the partner (and also receptor, if both_ways is true)
/// @param both_ways whether the two specified chains may play the role of both the partner and
///        the receptor. When false, the first chain is considered to be the receptor, and the
///        second one is considered to be the partner

void
PeptideDeriverFilter::derive_peptide(
	PeptideDeriverOutputter & output,
	core::pose::Pose const & pose,
	core::Size const first_chain_index,
	core::Size const second_chain_index,
	bool const both_ways ) const {

	// create a subpose with only the two chains
	utility::vector1<core::Size> residue_indices;

	PeptideDeriverFilter::push_chain_residue_indices( pose, first_chain_index, residue_indices );
	PeptideDeriverFilter::push_chain_residue_indices( pose, second_chain_index, residue_indices );

	core::pose::Pose chain_pair_pose;

	core::io::pose_from_pose( chain_pair_pose, pose, residue_indices );

	// prepare the chain pair pose
	prepare_pose(output, chain_pair_pose);

	// we shouldn't assume anything about the order of first_chain_index and second_chain_index
	// the chain given as first_chain_index might actually appear after the chain given as second_chain_index in the pose
	// in this case, first_chain_index will be chain 2 in the new pose
	//
	// two things to note:
	// 1. split_by_chain() will return the chains according to the ordering of their indices
	// 2. we assume a normal FoldTree, that is, chain-jump-chain. Constructing a new pose using pose_from_pose()
	//    aids this assumption, but if chains are messed up (e.g. a residue from one chain covalently connected
	//    to one from another chain), we will fail -- I honestly don't know if this is allowed in Rosetta.
	//    This assumption is employed in the call to get_jump_that_builds_residue(), and I think that's it.
	// - Yuval Sedan (yuval.sedan@mail.huji.ac.il), 2016-03-02
	core::Size const first_chain_new_index = (first_chain_index < second_chain_index? 1 : 2);
	core::Size const second_chain_new_index = (first_chain_index < second_chain_index? 2 : 1);

	// calculate the total energy of the chain pair pose
	// in any case, we expect the jump to come before chain #2 (though it might NOT be second_chain_new_index)
	core::Size const jump_id = chain_pair_pose.fold_tree().get_jump_that_builds_residue( chain_pair_pose.conformation().chain_begin(2) );
	core::Real const total_isc( calculate_interface_score(chain_pair_pose, jump_id) );


	tracer << "Chain pair prepared "
		<< chain_pair_pose.pdb_info()->chain(chain_pair_pose.conformation().chain_begin(first_chain_new_index))
		<< chain_pair_pose.pdb_info()->chain(chain_pair_pose.conformation().chain_begin(second_chain_new_index)) << std::endl;

	output.chain_pair_pose_prepared(chain_pair_pose);

	// derive the peptide
	utility::vector1<core::pose::PoseOP> chains = chain_pair_pose.split_by_chain();
	core::pose::PoseOP first_chain_pose = chains[first_chain_new_index];
	core::pose::PoseOP second_chain_pose = chains[second_chain_new_index];

	derive_peptide(output, *first_chain_pose, *second_chain_pose, total_isc);

	if ( both_ways ) {
		derive_peptide(output, *second_chain_pose, *first_chain_pose, total_isc);
	}
}

/// Goes through the residues in partner_pose with a sliding window (the size
/// of which is determined by PeptideDeriverFilter::pep_lengths_), and
/// calculates the energetic contribution of each such linear segment. This is
/// output for each position in the partner_pose, together with additional
/// data, such as the fraction of the energetic contribution to the total
/// interface score (@see PeptideDeriverFilter::calculate_interface_score()),
/// and putative SS bridge sites).
///
/// Finally, the best position is output again, with a special mark ('best')
/// to indicate that it is the best position.
///
/// Note: for lack of a better term, we call the chain from which the
///       peptide will be derived, the partner.
void
PeptideDeriverFilter::derive_peptide(
	PeptideDeriverOutputter & output,
	core::pose::Pose const & receptor_pose,
	core::pose::Pose const & partner_pose,
	core::Real const total_isc ) const {

	core::Size receptor_start = receptor_pose.conformation().chain_begin(1);
	//core::Size receptor_end = receptor_pose.conformation().chain_end(1);

	core::Size partner_start = partner_pose.conformation().chain_begin(1);
	core::Size partner_end = partner_pose.conformation().chain_end(1);

	char receptor_chain_letter = receptor_pose.pdb_info()->chain(receptor_start);
	char partner_chain_letter = partner_pose.pdb_info()->chain(partner_start);

	std::ostringstream options_string;

	// show options which have an effect on the way we're running
	options_string << "pep_lengths= "        << utility::join(pep_lengths_, ",") << " "
		<< "is_skip_zero_isc= "      << is_skip_zero_isc_                << " "
		<< "is_do_minimize= "        << is_do_minimize_                  ;

	tracer << "Begin receptor-partner pair " << receptor_chain_letter << partner_chain_letter << std::endl;
	output.begin_receptor_partner_pair(receptor_chain_letter, partner_chain_letter, total_isc, options_string.str());

	// this is used to skip discontinuities in the chain (and not derive discontinuous peptides)
	utility::vector1<core::Size> partner_cutpoints = partner_pose.fold_tree().cutpoints();
	if ( ! partner_cutpoints.empty() ) {
		tracer << "Partner cutpoints: " << utility::join(partner_cutpoints, ",") << std::endl;
	}

	// will be applied for each peptide with significant contribution to binding
	// where residues immediately before and after can be mutated to cysteins and form a disulfide bridge
	const core::Size PEPTIDE_CHAIN = 2;

	// TODO : perhaps we want the jump in the movemap that DisulfideInsertionMover uses to be user-defined (command-line option), to prevent the peptide from escaping the binding pocket
	protocols::simple_moves::DisulfideInsertionMoverOP disulfide_inserter( new protocols::simple_moves::DisulfideInsertionMover());
	disulfide_inserter->set_peptide_chain(PEPTIDE_CHAIN);

	core::scoring::ScoreFunctionOP scorefxn_N2C_minimize = core::scoring::get_score_function();
	scorefxn_N2C_minimize->set_weight(core::scoring::atom_pair_constraint, 0.1);
	scorefxn_N2C_minimize->set_weight(core::scoring::dihedral_constraint, 0.1);
	scorefxn_N2C_minimize->set_weight(core::scoring::angle_constraint, 0.1);
	//scorefxn_N2C_minimize->set_weight(core::scoring::coordinate_constraint, 0.001);

	for ( core::Size const pep_length : pep_lengths_ ) {

		if ( pep_length > (partner_end - partner_start + 1) ) {
			tracer << "Skip peptide length " << pep_length << " as it is larger than the chain" << std::endl;
			continue;
		}
		tracer << "Peptide length " << pep_length << std::endl;

		const core::Real CHECK_CYCLIZABLE_THRESHOLD_ISC = -0.01;
		const core::Real PRACTICALLY_ZERO_ISC = 1e-7;
		const core::Real UNLIKELY_ISC_VALUE = 2e6;


		DerivedPeptideEntryOP best_cyc_by_method[NUM_CYCLIZATION_METHODS];
		bool any_peptide_cyclizable_by_method[NUM_CYCLIZATION_METHODS];
		bool any_peptide_cyclic_model_created_by_method[NUM_CYCLIZATION_METHODS];


		// initialize arrays explicitly
		for ( core::Size method=0; method<NUM_CYCLIZATION_METHODS; ++method ) {
			best_cyc_by_method[method] = DerivedPeptideEntryOP( new DerivedPeptideEntry() );
			any_peptide_cyclizable_by_method[method] = false;
			any_peptide_cyclic_model_created_by_method[method] = false;
		}


		// DerivedPeptideEntry of the best linear initialized with high isc, 0 as pep_start position and a nullptr as the lin_pose
		// To be replaced by the first derived peptide, after isc comparison
		DerivedPeptideEntryOP best_lin(new DerivedPeptideEntry());


		output.peptide_length(pep_length);

		// keep track of the next chain break
		utility::vector1<core::Size>::const_iterator next_cutpoint_it = partner_cutpoints.begin();
		/*
		pseudo-code for better handling of non-protein residues :
		bool[partner_length] valid_start_points;
		foreach ( residue in partner_pose ) {
		if ( (! residue.is_protein()) || residue.is_cutpoint() ) {
		valid_start_points.set_range(residue.index-pep_length, pep_length, false); // make sure indexes work out
		}
		}
		*/

		for ( core::Size pep_start = partner_start; pep_start <= partner_end - pep_length + 1; ++pep_start ) {
			core::Size pep_end = pep_start + pep_length - 1;

			// TODO : skip windows with non-protein residues
			/*
			pseudo-code for better handling of non-protein residues :
			if ( ! valid_start_point[pep_start] ) continue;
			*/

			// skip sequences which have cutpoints
			// NOTE : for this to work, Rosetta must generate cutpoints where there is a chain discontinuity
			//        This is enforced by using the -in:missing_density_to_jump command line option, which is
			//        turned on by default in the Peptiderive app
			if ( (next_cutpoint_it != partner_cutpoints.end()) && (*next_cutpoint_it >= pep_start) && (*next_cutpoint_it < pep_end) ) {
				tracer << "Skip " << pep_start << "-" << *next_cutpoint_it << ", contains cutpoint" << std::endl;
				pep_start = *next_cutpoint_it; // it will be incremented further after 'continue'
				++next_cutpoint_it;
				continue;
			}


			// assemble a pose with the receptor and a cut-out peptide, starting at pep_start and ending at pep_end
			// here, the complex is the receptor and the derived peptide
			// NOTE : calculate_interface_score() needs to know the jump number, and for that we need to know the residue of the jump
			//        options were to use receptor_peptide_pose->fold_tree().get_residue_edge(receptor_pose.size() + 1), but we would have to
			//        know whether to look at the Edge's start() or stop() for the jump, and that depends on how build_receptor_peptide_pose() builds
			//        the pose. So it seemed preferable to just output it from build_receptor_peptide_pose().
			core::Size linear_jump_id;
			core::pose::PoseOP receptor_peptide_pose = build_receptor_peptide_pose(receptor_pose, partner_pose, pep_start, pep_end, linear_jump_id);

			// evaluate the newly created pose
			const core::Real linear_isc = calculate_interface_score(*receptor_peptide_pose, linear_jump_id);

			// skip zero linear_isc if requested
			if ( (is_skip_zero_isc_) && (fabs(linear_isc) < PRACTICALLY_ZERO_ISC) ) {
				tracer << "Skip " << pep_start << "-" << pep_end << ", zero (linear) isc" << std::endl;
				// we skipped a couple of residues, so make sure that the next_cutpoint_it is updated
				while ( (next_cutpoint_it != partner_cutpoints.end()) && (*next_cutpoint_it < pep_end) ) {
					++next_cutpoint_it;
				}
				pep_start = pep_end;
				continue;
			}

			DerivedPeptideEntryOP current_derived_peptide(new DerivedPeptideEntry(linear_isc, pep_start, receptor_peptide_pose));


			// This section checks and prints out closability of the peptide by head-to-tail (N-to-C termini) cyclization
			core::Length const MAX_DISTANCE = 5.0;
			core::Size pep_nter_idx = receptor_peptide_pose->conformation().chain_begin(PEPTIDE_CHAIN);
			core::Size pep_cter_idx = receptor_peptide_pose->conformation().chain_end(PEPTIDE_CHAIN);
			if ( receptor_peptide_pose->residue(pep_nter_idx).is_protein() && receptor_peptide_pose->residue(pep_cter_idx).is_protein() ) {
				core::Length N_to_C_dist = receptor_peptide_pose->residue( pep_nter_idx).xyz( "N" ).distance(receptor_peptide_pose->residue( pep_cter_idx ).xyz( "C" ));

				tracer.Debug << "N to C dist between positions " << pep_start << " and " << pep_end << " is: " << N_to_C_dist << std::endl;
				if ( (N_to_C_dist < MAX_DISTANCE) && ( (linear_isc / total_isc) >= optimize_cyclic_threshold_ ) ) {
					tracer << "calling PeptideCyclizeMover for partner positions: " << pep_start << " and " << pep_end << std::endl;
					CyclizedPeptideInfoOP N2C_pep_info = current_derived_peptide->cyc_info_set[CYCLIZATION_METHOD_END_TO_END];
					N2C_pep_info->is_cyclizable =true;
					any_peptide_cyclizable_by_method[CYCLIZATION_METHOD_END_TO_END] = true;
					core::pose::PoseOP N2C_cyclic_pose = generate_N2C_cyclic_peptide_protein_complex(receptor_peptide_pose, pep_nter_idx, pep_cter_idx, scorefxn_N2C_minimize);
					// accept or reject cyclic peptide - protein complex based on interface score
					core::Real N2C_cyclic_interface_score = calculate_interface_score(*N2C_cyclic_pose, linear_jump_id);

					core::scoring::EnergyMap cyclic_pose_emap (N2C_cyclic_pose->energies().total_energies());
					core::Real const cyc_pose_rama_energy = cyclic_pose_emap.get(core::scoring::score_type_from_name("rama_prepro"));
					core::Real const cyc_pose_omega_energy = cyclic_pose_emap.get(core::scoring::score_type_from_name("omega"));

					core::scoring::EnergyMap linear_pose_emap (receptor_peptide_pose->energies().total_energies());
					core::Real const linear_pose_rama_energy = linear_pose_emap.get(core::scoring::score_type_from_name("rama_prepro"));
					core::Real const linear_pose_omega_energy = linear_pose_emap.get(core::scoring::score_type_from_name("omega"));

					core::Real const delta_rama_score = cyc_pose_rama_energy - linear_pose_rama_energy;
					core::Real const delta_omega_score = cyc_pose_omega_energy - linear_pose_omega_energy;

					tracer << "changes in rama and omega due to cyclization: rama - " << delta_rama_score << "; omega - " << delta_omega_score << std::endl;
					tracer << "cyclic interface score of head-to-tail candidate is: " << N2C_cyclic_interface_score << std::endl;
					if ( N2C_cyclic_interface_score < 0 ) {
						N2C_pep_info->cyc_isc = N2C_cyclic_interface_score;
						N2C_pep_info->cyc_pose = N2C_cyclic_pose;
						N2C_pep_info->was_cyclic_model_created=true;
						std::stringstream N2C_info("");
						N2C_info << partner_chain_letter << "_" << pep_start << "-" << pep_end;
						N2C_pep_info->cyc_comment = N2C_info.str();
					}
				}
			}

			// This section checks and prints out closability of the peptide by a proposed putative wrapper disulfide bridge
			core::Size n_putative_cyd = pep_start - 1;
			core::Size c_putative_cyd = pep_end + 1;
			std::stringstream disulfide_info("");

			CyclizedPeptideInfoOP disulfide_pep_info = current_derived_peptide->cyc_info_set[CYCLIZATION_METHOD_DISULFIDE];

			if ( (linear_isc<CHECK_CYCLIZABLE_THRESHOLD_ISC) && (n_putative_cyd>=partner_start) && (c_putative_cyd<=partner_end) ) {
				protocols::simple_moves::DisulfideCyclizationViability cyclization_viable(disulfide_inserter->determine_cyclization_viability(partner_pose, n_putative_cyd, c_putative_cyd));

				disulfide_pep_info->is_cyclizable = ( (cyclization_viable == protocols::simple_moves::DCV_CYCLIZABLE) || (cyclization_viable == protocols::simple_moves::DCV_ALREADY_CYCLIZED) );


				if ( disulfide_pep_info->is_cyclizable ) {
					any_peptide_cyclizable_by_method[CYCLIZATION_METHOD_DISULFIDE] = true;
					disulfide_info << partner_chain_letter << "_" << n_putative_cyd << "-" << c_putative_cyd;

					// For peptides that can be closed and contribute more then a user defined fraction of the binding energy (default is 0.35)
					// mutate to cysteins and re-evaluate energy
					if ( (linear_isc / total_isc) >= optimize_cyclic_threshold_ ) {
						// NOTE : see note on linear_jump_id
						core::Size cyclic_jump_id;
						disulfide_pep_info->pre_cyc_pose = build_receptor_peptide_pose(receptor_pose, partner_pose, n_putative_cyd, c_putative_cyd, cyclic_jump_id);
						// create a copy of the linear pose which we are going to mutate and model
						// that is, for peptides for which we model the cyclic peptide, we also want to keep the expanded (+-1 residue) linear version
						core::pose::PoseOP receptor_cyclic_peptide_pose ( new core::pose::Pose( *(disulfide_pep_info->pre_cyc_pose)  ) );
						//define and create disulfide bond between the edges of the derived peptide using the DisulfideInsertionMover
						disulfide_inserter->apply(*receptor_cyclic_peptide_pose);
						// any_peptide_cyclic_model_created checks if any peptide in the pose has been cyclized
						// current_peptide_cyclic_model_created checks if the current peptide has been cyclized
						disulfide_pep_info->cyc_pose = receptor_cyclic_peptide_pose;
						if ( disulfide_inserter->get_last_move_status()==protocols::moves::MS_SUCCESS ) {
							any_peptide_cyclic_model_created_by_method[CYCLIZATION_METHOD_DISULFIDE] = true;
							disulfide_pep_info->was_cyclic_model_created = true;
							disulfide_pep_info->cyc_isc = calculate_interface_score(*receptor_cyclic_peptide_pose, cyclic_jump_id);
							disulfide_pep_info->cyc_comment =  disulfide_info.str();
						}
					}
				}
			}

			// output data for each cut-out peptide
			output.peptide_entry(ET_GENERAL, total_isc, current_derived_peptide);

			// save details for the linear peptide that scored best
			if ( linear_isc < best_lin->lin_isc ) {
				best_lin = current_derived_peptide;
			}

			// save details for the cyclizable peptide that scored best
			//
			// We only model cyclic peptides out of peptides that contribute enough to the interactions energy
			// (more than the fraction specified by optimize_cyclic_threshold_).
			//
			// There outputare three ways to define the scoring or the set from which the 'best' peptide may be taken:
			// 1. best scoring linear peptide pose
			// 2. best scoring linear peptide pose among the cyclization candidates (if any candidates exist)
			// 3. best scoring cyclic pose (if such were modeled among the candidates)
			// Right now we merge the 2nd and 3rd into
			// 2/3. the best cyclizable pose, which is the best scoring cyclic pose, if one is avaiable (i.e. was modeled), otherwise, the best linear pose among cyclization candidates
			//
			// Now, consider the following scenario:
			// Total interface score: -50
			// Position    Linear peptide interface score     Cyclic peptide interface score
			// 10          -10                                -4
			// 20          -20                                -3
			// 30          -30                                -3.5
			//
			// Let's assume all listed positions are cyclization candidates, and that they are the only cyclization candidates for this chain pair.
			// We might get different results here depending on the value of optimize_cyclic_threshold
			//
			// optimize_cyclic_threshold    Positions that produce cyclic peptides       Best cyclizable peptide       Based on
			// 0  -0.2                      10, 20, 30                                   10                            Cyclic peptide interface score
			// 0.2-0.4                      20, 30                                       30                            Cyclic peptide interface score
			// 0.4-0.6                      30                                           30                            Cyclic peptide interface score
			// 0.6-1                        N/A                                          30                            Linear peptide interface score
			//
			// Note that position 20 will never be considered best cyclizable, since if no cyclic model was produced, position 30 has a better linear peptide interface score.
			// If position 20 is compared using the cyclic peptide interface score, this means a cyclic model was produced for it, which means it would be produced for position 30 as well (has a better linear peptide interface score), and position 30 has a better cyclic interface score.
			// So, if a peptide B has both a better linear peptide interface score and a better cyclic peptide interface score than peptide A, peptide A will never be the best cyclizable peptide.

			// TODO : how safe is it to assume that the first CyclizationMethod enum value is 0?
			//        this is impotant since we map it between CYCLIZATION_METHOD_NAMES
			for ( core::Size method = 0; method<NUM_CYCLIZATION_METHODS; ++method ) {
				CyclizedPeptideInfoOP cyc_info = current_derived_peptide->cyc_info_set[method];
				CyclizedPeptideInfoOP best_cyc_info = best_cyc_by_method[method]->cyc_info_set[method];
				if ( cyc_info->is_cyclizable &&
						(
						(cyc_info->was_cyclic_model_created && (cyc_info->cyc_isc < best_cyc_info->cyc_isc)) ||
						((!any_peptide_cyclic_model_created_by_method[method]) && (current_derived_peptide->lin_isc < best_cyc_by_method[method]->lin_isc))
						) ) {
					best_cyc_by_method[method] = current_derived_peptide;
				}
			}
		} // for pep_start

		tracer << "Outputting best peptides" << std::endl;
		if ( best_lin->lin_isc < UNLIKELY_ISC_VALUE ) {
			output.peptide_entry(ET_BEST_LINEAR, total_isc, best_lin);
			for ( core::Size method = 0; method<NUM_CYCLIZATION_METHODS; ++method ) {
				if ( any_peptide_cyclizable_by_method[method] ) {
					// otherwise, cyc_pose_of_best_cyc is unset
					output.peptide_entry(ET_BEST_CYCLIC, total_isc, best_cyc_by_method[method]);
				}
			}
		}
	} // for each peptide_length

	tracer << "End receptor-partner pair" << std::endl;
	output.end_receptor_partner_pair();
}

/// Given a two-monomer pose and a score function - the interface score is the
/// score difference between the given complex and the unbound monomers,
/// 'unbound' meaning that they are spatially seperated.
///
/// @param pose A two-monomer pose.
/// TODO : this might have code in common with DdgScan. Consider using that.
core::Real
PeptideDeriverFilter::calculate_interface_score(core::pose::Pose const & pose, core::Size const jump_id) const {
	//make a new pose in which the monomers are far apart from one another
	core::pose::Pose unbound_pose(pose);

	core::Real const VERY_LARGE_TRANSLATION = 2e6;

	//create a Rigid body mover - the mover works on a rigid body jump
	protocols::rigid::RigidBodyTransMoverOP translate_away(
		new protocols::rigid::RigidBodyTransMover(unbound_pose, jump_id));

	//set the size of the move
	translate_away->step_size(VERY_LARGE_TRANSLATION);

	//calculate scores
	core::Real bound_energy = (*scorefxn_deriver_)(unbound_pose);  // before the move - as a complex.
	translate_away->apply(unbound_pose);
	core::Real unbound_energy = (*scorefxn_deriver_)(unbound_pose);  // after the move - unbound monomers.

	//answer
	return (bound_energy - unbound_energy);
}


/// Take a given linear peptide-protein complex and return a head-to-tail (N-ter to C-ter) cyclic peptide - receptor complex
/// @param receptor_peptide_pose linear peptide-protein complex
/// @param pep_nter_idx index of peptide start position in the given complex
/// @param pep_cter_idx index of peptide end position in the given complex
/// @param scorefxn_N2C_minimize scorefunction with constraints turned on to be used in complex minimization
core::pose::PoseOP
PeptideDeriverFilter::generate_N2C_cyclic_peptide_protein_complex(core::pose::PoseOP const receptor_peptide_pose, core::Size const pep_nter_idx, core::Size  const pep_cter_idx, core::scoring::ScoreFunctionOP const scorefxn_N2C_minimize) const {
	core::pose::PoseOP pose_for_N2C_cyclization( new core::pose::Pose(*receptor_peptide_pose));

	// setup for calling PeptideCyclizeMover to setup chemical bond and constraints on N-ter and C-ter residues
	char const cyclic_peptide_chain = pose_for_N2C_cyclization->pdb_info()->chain(pep_nter_idx);
	core::select::residue_selector::ChainSelectorCOP peptide_chain_selector (new core::select::residue_selector::ChainSelector( cyclic_peptide_chain));
	protocols::cyclic_peptide::PeptideCyclizeMoverOP cyclize_head_to_tail ( new protocols::cyclic_peptide::PeptideCyclizeMover () );
	cyclize_head_to_tail->set_selector(peptide_chain_selector);
	cyclize_head_to_tail->apply(*pose_for_N2C_cyclization);
	// minimization of the structure to actually close the bond
	core::kinematics::MoveMapOP N2C_minimize_movemap ( new core::kinematics::MoveMap );
	N2C_minimize_movemap->set_bb_true_range(pep_nter_idx, pep_cter_idx);
	protocols::simple_moves::MinMover N2C_minimizer( N2C_minimize_movemap, scorefxn_N2C_minimize, "dfpmin_armijo_atol", 0.01 /*tolerance*/, true /*nb_list*/ );
	N2C_minimizer.apply( *pose_for_N2C_cyclization );
	// additional bond declaration to correct the position of the amide proton and the carboxyl oxygen,
	// as the N-ter phi and C-ter psi are not sampled during minimization
	protocols::cyclic_peptide::DeclareBondOP declare_N2C_bond ( new protocols::cyclic_peptide::DeclareBond() );
	declare_N2C_bond->set(pep_nter_idx,"N",pep_cter_idx,"C",false /*add_termini*/,false /*run KIC*/, 0 /*KIC res1 */, 0 /*KIC res2*/, false /*rebuild_fold_tree*/);
	declare_N2C_bond->apply(*pose_for_N2C_cyclization);
	return pose_for_N2C_cyclization;
}



/// Given a two-monomer pose and a score function - the interface score (delta between energy of complex and of the seperated monomers) is calculated for each residue.
/// This signifies the energetic contribution of each residue to the complex binding.
/// @param pose a two-monomer pose.
void
PeptideDeriverFilter::calculate_per_residue_interface_score(
	core::pose::Pose & chain_pair_pose,
	std::set<core::Size> interface_residues,
	std::ostream & report_out,
	core::Size const jump_id) const {

	//seperate the given structure in space using a translate_away mover and score each of the states
	core::Real const VERY_LARGE_TRANSLATION = 2e6;
	core::pose::Pose unbound_pair(chain_pair_pose);
	protocols::rigid::RigidBodyTransMoverOP translate_away(new protocols::rigid::RigidBodyTransMover(unbound_pair, jump_id));
	translate_away->step_size(VERY_LARGE_TRANSLATION);
	translate_away->apply(unbound_pair);
	// make sure hbonds are included in the scoring function's per residue energies
	core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( scorefxn_deriver_->energy_method_options() ) );
	emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
	scorefxn_deriver_->set_energy_method_options( *emopts );
	(*scorefxn_deriver_)(chain_pair_pose);
	(*scorefxn_deriver_)(unbound_pair);

	//go over each interface residue and calculate its energetical contribution to binding
	for ( core::Size interface_residue : interface_residues ) {
		// delta between total energy in the bound and unbound states is in fact the interface score
		core::Real res_isc = chain_pair_pose.energies().residue_total_energy(interface_residue) - unbound_pair.energies().residue_total_energy(interface_residue);
		report_out << "residue " << interface_residue << " in chain " << chain_pair_pose.residue(interface_residue).chain() << " has interface score " << res_isc << std::endl;
	}

}

/// for a two-monomer pose, returns a set of interface residues
/// @param pose a two-monomer pose.
std::set<core::Size>
PeptideDeriverFilter::find_interface_residues(core::pose::Pose const & chain_pair_pose) {
	//find the interface residues using an InterfaceNeighborDefinition calculator and store them in the interface_residues set
	core::Size chain1_idx =1, chain2_idx =2;
	basic::MetricValue< std::set< core::Size > > mv_interface_set;
	core::pose::metrics::PoseMetricCalculatorOP if_residues_calculator(new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator(chain1_idx, chain2_idx));
	std::string InterfaceNeighborDefinition = "InterfaceNeighborDefinition_";
	core::pose::metrics::CalculatorFactory::Instance().register_calculator(InterfaceNeighborDefinition, if_residues_calculator);
	chain_pair_pose.metric( InterfaceNeighborDefinition, "interface_residues", mv_interface_set);
	std::set<core::Size> interface_residues = mv_interface_set.value();
	return interface_residues;
}

/// For two given full proteins poses, builds a new pose with the full receptor pose connected to the peptide found between start and end positions in the partner pose
/// @param receptor_pose the receptor chain to be maintained
/// @param partner_pose  the chain from which a peptide should be cut and connected to the receptor in a new pose
/// @param peptide_start starting position in the partner pose to cut from
/// @param peptide_end   last position in the partner pose to be cut out
/// @param jump_id       the number of the jump that was used to connect the two chains
core::pose::PoseOP PeptideDeriverFilter::build_receptor_peptide_pose(core::pose::Pose const & receptor_pose,
	core::pose::Pose const & partner_pose,
	core::Size peptide_start,
	core::Size peptide_end,
	core::Size & jump_id) const {
	// once created these will be the positions of the peptide start and end in the receptor_peptide_pose
	core::Size two_chain_peptide_start = receptor_pose.size() + 1;
	core::Size two_chain_peptide_end = two_chain_peptide_start + peptide_end - peptide_start;

	//copy receptor and start appending the peptide from its center of mass backward and forward
	core::pose::PoseOP receptor_peptide_pose( new core::pose::Pose(receptor_pose) );
	core::Size pep_cen_mass = core::pose::residue_center_of_mass(partner_pose, peptide_start, peptide_end);

	// from Flexpepdock's flags file; choose the residue in the receptor that is closest to the peptide center-of-mass
	// Atom(2) will be CA for canonical amino-acids and hopefully an atom that exists in any other residue type
	core::Size receptor_anchor_pos = core::pose::return_nearest_residue(receptor_pose, 1, receptor_pose.size(), partner_pose.residue(pep_cen_mass).atom(2).xyz());

	receptor_peptide_pose->append_residue_by_jump(partner_pose.residue(pep_cen_mass), receptor_anchor_pos, "", "", true);

	for ( core::Size i = pep_cen_mass ; i > peptide_start; --i ) {
		// tracer << "Prepending at " << two_chain_peptide_start << " " << partner_pose.residue(i-1) << std::endl;
		receptor_peptide_pose->conformation().safely_prepend_polymer_residue_before_seqpos(partner_pose.residue(i-1), two_chain_peptide_start, false);
	}

	for ( core::Size i = pep_cen_mass ; i < peptide_end; ++i ) {
		// tracer << "Appending at " << receptor_peptide_pose->size() << " " << partner_pose.residue(i+1) << std::endl;
		receptor_peptide_pose->conformation().safely_append_polymer_residue_after_seqpos(partner_pose.residue(i+1), receptor_peptide_pose->size(), false);
	}

	// see note for linear_jump_id
	jump_id = receptor_peptide_pose->fold_tree().get_jump_that_builds_residue(receptor_pose.size() + pep_cen_mass - peptide_start + 1);

	// NOTE : when enabling -in:missing_density_to_jump (which is enabled by default for the
	//        PeptideDeriver app), chain termini which are missing hydrogens will receive a
	//        UPPER/LOWERTERM_TRUNC_VARIANT. Without removing those before trying to add the
	//        UPPER/LOWER_TERMINUS_VARIANT the following error occurs:
	//        ERROR: unable to find desired variant residue: MET:NtermTruncation MET LOWER_TERMINUS_VARIANT
	if ( receptor_peptide_pose->residue(two_chain_peptide_start).has_variant_type(core::chemical::LOWERTERM_TRUNC_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue( *receptor_peptide_pose, core::chemical::LOWERTERM_TRUNC_VARIANT, two_chain_peptide_start );
	}
	if ( receptor_peptide_pose->residue(two_chain_peptide_end).has_variant_type(core::chemical::UPPERTERM_TRUNC_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue( *receptor_peptide_pose, core::chemical::UPPERTERM_TRUNC_VARIANT, two_chain_peptide_end );
	}
	core::pose::add_variant_type_to_pose_residue(*receptor_peptide_pose, core::chemical::LOWER_TERMINUS_VARIANT, two_chain_peptide_start);
	core::pose::add_variant_type_to_pose_residue(*receptor_peptide_pose, core::chemical::UPPER_TERMINUS_VARIANT, two_chain_peptide_end);

	// if a cysteine in the peptide was disulfide bonded to a residue that no longer exists in the protein-peptide pose, it is converted back to a regular cysteine
	for ( core::Size i = two_chain_peptide_start; i <= two_chain_peptide_end ; ++i ) {
		if ( receptor_peptide_pose->conformation().residue(i).has_variant_type( core::chemical::DISULFIDE ) ) {
			core::conformation::change_cys_state(i, "CYS", receptor_peptide_pose->conformation());
		}
	}
	receptor_peptide_pose->conformation().detect_disulfides();
	return receptor_peptide_pose;
}


void
PeptideDeriverFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & , // filters - unused
	protocols::moves::Movers_map const & , // movers - unused
	core::pose::Pose const & // pose - unused
) {

	// NOTE : we rely on the fact that all the options already have values
	//        from when parse_options() was called. This way command-line
	//        defaults are also effectively RosettaScripts defaults.

	if ( tag->hasOption("pep_lengths") ) {
		set_pep_lengths(utility::string_split( tag->getOption<std::string>("pep_lengths"), ',', core::Size() ));
	}

	if ( tag->hasOption("skip_zero_isc") ) {
		set_is_skip_zero_isc(tag->getOption<bool>("skip_zero_isc"));
	}

	if ( tag->hasOption("dump_peptide_pose") ) {
		set_is_dump_peptide_pose(tag->getOption<bool>("dump_peptide_pose"));
	}

	if ( tag->hasOption("dump_cyclic_poses") ) {
		set_is_dump_cyclic_poses(tag->getOption<bool>("dump_cyclic_poses"));
	}

	if ( tag->hasOption("dump_report_file") ) {
		set_is_dump_report_file(tag->getOption<bool>("dump_report_file"));
	}

	if ( tag->hasOption("report_gzip") ) {
		set_is_report_gzip(tag->getOption<bool>("report_gzip"));
	}

	if ( tag->hasOption("dump_prepared_pose") ) {
		set_is_dump_prepared_pose(tag->getOption<bool>("dump_prepared_pose"));
	}

	if ( tag->hasOption("do_minimize") ) {
		set_is_do_minimize(tag->getOption<bool>("do_minimize"));
	}

	if ( tag->hasOption("scorefxn_deriver") ) {
		set_scorefxn_deriver( protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_deriver", data, "" ) );
	}

	if ( tag->hasOption("optimize_cyclic_threshold") ) {
		set_optimize_cyclic_threshold( tag->getOption<core::Real>( "optimize_cyclic_threshold" ) );
	}

	if ( tag->hasOption("report_format") ) {
		set_report_format( PeptideDeriverFilter::parse_report_format_string(tag->getOption<std::string>( "report_format" )) );
	}

	utility::vector1<char> restrict_receptors_to_chains(
		utility::string_split( tag->getOption<std::string>("restrict_receptors_to_chains", ""), ',', char() ));

	utility::vector1<char> restrict_partners_to_chains(
		utility::string_split( tag->getOption<std::string>("restrict_partners_to_chains", ""), ',', char() ));

	set_restrict_receptors_to_chains(restrict_receptors_to_chains);
	set_restrict_partners_to_chains(restrict_partners_to_chains);
}

std::string PeptideDeriverFilter::name() const {
	return class_name();
}

std::string PeptideDeriverFilter::class_name() {
	return "PeptideDeriver";
}

void PeptideDeriverFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction format;
	format.name( "peptide_deriver_format" );
	format.base_type( xs_string );
	format.add_restriction( xsr_enumeration, "markdown" );
	format.add_restriction( xsr_enumeration, "basic" );
	xsd.add_top_level_element( format );
	AttributeList attlist;
	//Define report format string
	attlist
		+ XMLSchemaAttribute(
		"pep_lengths", xsct_nnegative_int_cslist,
		"Length(s) of derived peptides" )
		+ XMLSchemaAttribute(
		"skip_zero_isc", xsct_rosetta_bool,
		"Makes derivation go faster by skipping peptides with 0 interface score" )
		+ XMLSchemaAttribute(
		"dump_peptide_pose", xsct_rosetta_bool,
		"Output pose with peptide cut out (best one for each chain pair)" )
		+ XMLSchemaAttribute(
		"dump_cyclic_poses", xsct_rosetta_bool,
		"Output each cyclic peptide pose (those that are modeled; "
		"which is determined by -optimize_cyclic_threshold)" )
		+ XMLSchemaAttribute(
		"dump_report_file", xsct_rosetta_bool,
		"Send PeptideDeriver output to a file (.peptiderive.txt)" )
		+ XMLSchemaAttribute(
		"report_gzip", xsct_rosetta_bool,
		"Compress report" )
		+ XMLSchemaAttribute(
		"dump_prepared_pose", xsct_rosetta_bool,
		"Output each receptor-partner pose as PeptiDerive sees it, "
		"i.e. after preparation (minimization and disulfide detection)" )
		+ XMLSchemaAttribute(
		"do_minimize", xsct_rosetta_bool,
		"Perform minimization before everything." )
		+ XMLSchemaAttribute(
		"optimize_cyclic_threshold", xsct_rosetta_bool,
		"Value of peptide interface score percent of total isc "
		"from which to optimize cyclic peptide" )
		+ XMLSchemaAttribute(
		"report_format", "peptide_deriver_format",
		"The format of the report. Either basic (easily parsable format) "
		"or markdown (pretty, readable, but verbose format)" )
		+ XMLSchemaAttribute(
		"restrict_receptors_to_chains", xsct_chain_cslist,
		"Only use chains listed here as receptors. "
		"When empty, consider all chains." )
		+ XMLSchemaAttribute(
		"restrict_partners_to_chains", xsct_chain_cslist,
		"Only use chains listed here as partners. When empty, "
		"consider all chains. For each receptor-partner pair, a "
		"peptide is derived from the partner." );
	rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn_deriver" );


	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Filter for PeptiDerive, a simple application that derives from a given interface "
		"the linear stretch that contributes most of the binding energy (approximated as "
		"the score over the interface).",
		attlist );
}

std::string PeptideDeriverFilterCreator::keyname() const {
	return PeptideDeriverFilter::class_name();
}

protocols::filters::FilterOP
PeptideDeriverFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new PeptideDeriverFilter );
}

void PeptideDeriverFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PeptideDeriverFilter::provide_xml_schema( xsd );
}


} // namespace peptide_deriver
} // namespace protocols
