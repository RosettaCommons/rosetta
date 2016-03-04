// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief test suite for protocols::analysis::PeptideDeriverFilter
/// @author Yuval Sadan (yuval.sedan@mail.huji.ac.il)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/analysis/PeptideDeriverFilter.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>

// C++ headers
#include <fstream>

class PeptideDeriverUTracerOutputter : public protocols::analysis::PeptideDeriverOutputter {
private:
	utility::pointer::shared_ptr<test::UTracer> ut_;

public:
	PeptideDeriverUTracerOutputter() {
		ut_ = utility::pointer::shared_ptr<test::UTracer>( new test::UTracer("protocols/analysis/PeptideDeriverFilter.u") );
	}

	virtual void begin_structure(core::pose::Pose const & pose, std::string const &) {
		// note: current unit tests don't cover this code (called from report())
		(*ut_) << "begin_structure" << std::endl;
		(*ut_) << pose;
	}

	virtual void chain_pair_pose_prepared(core::pose::Pose const & /*pose*/) {
		(*ut_) << "chain_pair_pose_prepared" << std::endl;
	}

	virtual void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const /*total_isc*/,
		std::string const & options_string) {
		(*ut_) << "begin_receptor_partner_pair" <<
			" receptor=" << receptor_chain_letter <<
			" partner=" << partner_chain_letter <<
			" options=" << options_string <<
			std::endl;
	}

	virtual void peptide_length(core::Size const pep_length) {
		// we don't output here because we don't care in which order this method is called
		(*ut_) << "peptide_length " << pep_length << std::endl;
	}

	virtual void peptide_entry(core::pose::Pose const & /*pose*/,
		protocols::analysis::PeptideDeriverEntryType const entry_type, core::Size const pep_start,
		core::Real const /*linear_isc*/, core::Real const /*binding_contribution_fraction*/, std::string const & /*disulfide_info*/,
		bool const /*was_cyclic_pep_modeled*/, core::pose::Pose const & /*cyclic_pose*/,
		core::Real const /*cyclic_isc*/) {
		(*ut_) << "peptide_entry" <<
			" entry_type=" << entry_type <<
			" pep_start=" << pep_start <<
			std::endl;
		// we don't compare energies (was_cyclic_pep_modeled depends on energies, so we don't do that either)
	}

	virtual void end_receptor_partner_pair() {
		(*ut_) << "end_receptor_partner_pair" << std::endl;
	}

	virtual void end_structure() {
		// note: current unit tests don't cover this code (called from report())
		(*ut_) << "end_structure" << std::endl;
	}
};

class PeptideDeriverFilterTests : public CxxTest::TestSuite {
public:
	core::pose::PoseOP test_pose_;

	protocols::analysis::PeptideDeriverFilterOP peptiderive_;

	PeptideDeriverFilterTests() {}

	void setUp() {
		protocols_init_with_additional_options("-in:missing_density_to_jump 1");
		test_pose_ = core::pose::PoseOP( new core::pose::Pose() );

		// NOTE : this pose is PDB entry 2HLE modified such that
		// 1. only residues 144-157 from chain A and 114-134 from chain B are included
		// 2. minimized using the minimize app
		// 3. hydrogens and score terms removed from minimized file
		// 4. residue 117 from chain B removed
		// 5. residues 114-121 from chain B duplicated as 135-142 and shifted 100A in each axis
		// 6. residues 114-121 from chain B duplicated as chain C 1-12, shifted 100A in x, y and 121A in z
		//
		// Steps 1-2 are made so that processing time is minimal for the unit tests.
		//
		// Step 3 is made so that we handle missing hydrogens in input files properly (the
		// 'missing_density_to_jump' will tag termini as truncated, and there was a bug fix
		// for this, so this verifies we continue to handle this case correctly).
		//
		// Step 4 is made so that missing density is present, which should turn to jumps
		// and residues 114-116 should be skipped because no continuous peptides could be
		// derived from those positions.
		//
		// Step 5 is made so that a zero-isc peptide exists, and should be skipped when a 7-mer
		// should be extracted. This also has the advantage for different control flow for the
		// 10-mer, which should be skipped because the peptide is not long enough (and is an
		// edge case).
		//
		// Step 6 made so that we test the case for far-away chains C and A (where no peptide
		// will be derived at all), and so that we test more than two chains. The value of 121A
		// for the shift in the z-axis is 100A (the same as those residues shifted in chain B)
		// plus 21A, which makes sure chain C is shifted completely away from chain B (the width
		// in the z-axis of this segment is 18A), but still close enough so a B-C interaction
		// exists.
		core::import_pose::pose_from_file( *test_pose_, "protocols/analysis/2hle_remixed.pdb" , core::import_pose::PDB_file);

		peptiderive_ = protocols::analysis::PeptideDeriverFilterOP( new protocols::analysis::PeptideDeriverFilter() );

		peptiderive_->set_is_skip_zero_isc(false); // we want to count all residues
		peptiderive_->set_is_do_minimize(true);

		// the model above has exactly one cyclic peptide candidate; this make sure it gets modeled, so we cover the disulfide creation code as well
		peptiderive_->set_optimize_cyclic_threshold(-1);

		// the following only have an effect when report() is invoked.
		peptiderive_->set_is_dump_cyclic_poses(false);
		peptiderive_->set_is_dump_report_file(false);
		peptiderive_->set_is_report_gzip(false);
		peptiderive_->set_is_dump_peptide_pose(false);
		peptiderive_->set_is_dump_prepared_pose(false);
		peptiderive_->set_report_format(protocols::analysis::PRF_MARKDOWN);
	}

	void tearDown() {
		test_pose_ = NULL; // delete
	}

	void test_rosettascripts_options() {

		std::ifstream xml_stream( "protocols/analysis/test_peptiderive.xml" );

		utility::tag::TagCOP tag =  utility::tag::Tag::create( xml_stream );
		try {
			protocols::rosetta_scripts::RosettaScriptsParser parser;
			// NOTE : the following assumes that RosettaScriptParser::parse_protocol_tag() returns a pointer ParsedProtocol instance, even though the interface is to a MoverOP.
			protocols::rosetta_scripts::ParsedProtocolOP protocol(utility::pointer::dynamic_pointer_cast<protocols::rosetta_scripts::ParsedProtocol> (parser.parse_protocol_tag(tag)) );
			TS_ASSERT_EQUALS(protocol->size(), 1);
			protocols::rosetta_scripts::ParsedProtocol::MoverFilterPair pair = protocol->get_mover_filter_pair(1);
			TS_ASSERT_EQUALS(pair.filter().get_type(), "PeptideDeriverFilter");
			protocols::analysis::PeptideDeriverFilter const & filter( dynamic_cast<protocols::analysis::PeptideDeriverFilter const &> (pair.filter()) );

			TS_ASSERT_EQUALS(filter.get_is_skip_zero_isc(), false);
			TS_ASSERT_EQUALS(filter.get_pep_lengths().size(), 2);
			TS_ASSERT_EQUALS(filter.get_pep_lengths()[1], 7);
			TS_ASSERT_EQUALS(filter.get_pep_lengths()[2], 10);
			TS_ASSERT_DELTA(filter.get_optimize_cyclic_threshold(), 0.5, 1e-6);
			TS_ASSERT_EQUALS(filter.get_scorefxn_deriver()->get_name(), "soft_rep_design");
			TS_ASSERT_EQUALS(filter.get_is_do_minimize(), false);
			TS_ASSERT_EQUALS(filter.get_is_dump_prepared_pose(), true);
			TS_ASSERT_EQUALS(filter.get_is_dump_report_file(), true);
			TS_ASSERT_EQUALS(filter.get_is_report_gzip(), true);
			TS_ASSERT_EQUALS(filter.get_is_dump_cyclic_poses(), true);
			TS_ASSERT_EQUALS(filter.get_is_dump_peptide_pose(), true);
			TS_ASSERT_EQUALS(filter.get_restrict_receptors_to_chains().size(), 2);
			TS_ASSERT_EQUALS(filter.get_restrict_receptors_to_chains()[1], 'A');
			TS_ASSERT_EQUALS(filter.get_restrict_receptors_to_chains()[2], 'B');
			TS_ASSERT_EQUALS(filter.get_restrict_partners_to_chains().size(), 2);
			TS_ASSERT_EQUALS(filter.get_restrict_partners_to_chains()[1], 'B');
			TS_ASSERT_EQUALS(filter.get_restrict_partners_to_chains()[2], 'C');
		}
		catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_protocol() {
		PeptideDeriverUTracerOutputter ut_out;
		utility::vector1<core::Size> pep_lengths;
		pep_lengths.push_back(10);
		pep_lengths.push_back(7);
		peptiderive_->set_pep_lengths(pep_lengths);
		peptiderive_->derive_peptide( ut_out, *test_pose_, /*first_chain=*/1, /*second_chain=*/2, /*both_ways=*/true );

		// chain A and C are far apart. Here we look at what happens in this case (no peptide should be derived)
		peptiderive_->set_is_skip_zero_isc(true);
		peptiderive_->derive_peptide( ut_out, *test_pose_, /*first_chain=*/1, /*second_chain=*/3, /*both_ways=*/true );
		peptiderive_->set_is_skip_zero_isc(false);

	}

	void test_report_format_parsing() {
		TS_ASSERT_EQUALS(protocols::analysis::PeptideDeriverFilter::parse_report_format_string("basic"), protocols::analysis::PRF_BASIC);
		TS_ASSERT_EQUALS(protocols::analysis::PeptideDeriverFilter::parse_report_format_string("markdown"), protocols::analysis::PRF_MARKDOWN);
		TS_ASSERT_THROWS(protocols::analysis::PeptideDeriverFilter::parse_report_format_string(""), utility::excn::EXCN_KeyError);
	}

	void assert_peptiderive_filter_equal( protocols::analysis::PeptideDeriverFilter const & lhs,
		protocols::analysis::PeptideDeriverFilter const & rhs ) {
		TS_ASSERT_EQUALS(lhs.get_is_skip_zero_isc(), rhs.get_is_skip_zero_isc());
		TS_ASSERT_EQUALS(lhs.get_pep_lengths().size(), rhs.get_pep_lengths().size());
		for ( core::Size i = 1; i <= rhs.get_pep_lengths().size(); ++i ) {
			TS_ASSERT_EQUALS(lhs.get_pep_lengths()[i], rhs.get_pep_lengths()[i]);
		}
		TS_ASSERT_DELTA(lhs.get_optimize_cyclic_threshold(), rhs.get_optimize_cyclic_threshold(), 1e-6);
		TS_ASSERT_EQUALS(lhs.get_scorefxn_deriver()->get_name(), rhs.get_scorefxn_deriver()->get_name());
		TS_ASSERT_EQUALS(lhs.get_is_do_minimize(), rhs.get_is_do_minimize());
		TS_ASSERT_EQUALS(lhs.get_is_dump_prepared_pose(), rhs.get_is_dump_prepared_pose());
		TS_ASSERT_EQUALS(lhs.get_is_dump_report_file(), rhs.get_is_dump_report_file());
		TS_ASSERT_EQUALS(lhs.get_is_report_gzip(), rhs.get_is_report_gzip());
		TS_ASSERT_EQUALS(lhs.get_is_dump_cyclic_poses(), rhs.get_is_dump_cyclic_poses());
		TS_ASSERT_EQUALS(lhs.get_is_dump_peptide_pose(), rhs.get_is_dump_peptide_pose());
		TS_ASSERT_EQUALS(lhs.get_restrict_receptors_to_chains().size(), rhs.get_restrict_receptors_to_chains().size());
		for ( core::Size i = 1; i <= rhs.get_restrict_receptors_to_chains().size(); ++i ) {
			TS_ASSERT_EQUALS(lhs.get_restrict_receptors_to_chains()[i], rhs.get_restrict_receptors_to_chains()[i]);
		}
		TS_ASSERT_EQUALS(lhs.get_restrict_partners_to_chains().size(), rhs.get_restrict_partners_to_chains().size());
		for ( core::Size i = 1; i <= rhs.get_restrict_partners_to_chains().size(); ++i ) {
			TS_ASSERT_EQUALS(lhs.get_restrict_partners_to_chains()[i], rhs.get_restrict_partners_to_chains()[i]);
		}
		TS_ASSERT_EQUALS(lhs.get_report_format(), rhs.get_report_format());
	}

	void test_copy_ctor() {
		protocols::analysis::PeptideDeriverFilter another_peptiderive( *peptiderive_ );
		assert_peptiderive_filter_equal( another_peptiderive, *peptiderive_ );
	}

	void test_assignment() {
		protocols::analysis::PeptideDeriverFilter another_peptiderive;
		another_peptiderive = *peptiderive_;
		assert_peptiderive_filter_equal( another_peptiderive, *peptiderive_ );
	}

	// TODO : or not TODO ?
	// void test_commandline_options() { }
	// void test_report_null_input() {}
	// void test_report_single_chain_input() {}
	// void test_report_multiple_input() {}
	// void test_report_short_chains() {}
	// void test_report_has_expected_output() {}
	// void test_prepare_pose() {}
	// void test_apply_does_nothing() {}
};
