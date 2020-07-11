// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pdbinfo_manipulations/AlignPDBInfoToSequencesTests.cxxtest.hh
/// @brief  tests for AlignPDBInfoToSequences
/// @author Dan Farrell (danpf@uw.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1rpb.hh>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>

// Project Headers
#include <protocols/pdbinfo_manipulations/AlignPDBInfoToSequences.hh>
#include <protocols/pdbinfo_manipulations/AlignPDBInfoToSequencesUtil.hh>

// Core Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/json_utilities.hh>
#include <utility/io/izstream.hh>

static basic::Tracer TR("AlignPDBInfoToSequencesTests");


class AlignPDBInfoToSequencesTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options("-in:file:centroid -ignore_unrecognized_res -ignore_zero_occupancy false -missing_density_to_jump true -mute core.io core.conformation -packing::pack_missing_sidechains false");
	}

	void tearDown(){
	}

	utility::vector1<std::string>
	to_string_vec(std::string const & str_in) const {
		utility::vector1<std::string> str_vec_out;
		for ( char const & c : str_in ) str_vec_out.push_back(std::string(1, c));
		return str_vec_out;
	};

	void test_set_modes() {
		protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts;
		TS_ASSERT_THROWS_NOTHING(apts.set_mode_from_string("multiple"));
		TS_ASSERT_THROWS_NOTHING(apts.set_mode_from_string("single"));
		TS_ASSERT_THROWS(apts.set_mode_from_string("goober"), utility::excn::BadInput&);
	}

	void test_mode_enum() {
		TS_ASSERT_EQUALS(protocols::pdbinfo_manipulations::get_alignpdbinfotosequencesmode_from_string("single"), protocols::pdbinfo_manipulations::AlignPDBInfoToSequencesMode::single);
		TS_ASSERT_EQUALS(protocols::pdbinfo_manipulations::get_alignpdbinfotosequencesmode_from_string("multiple"), protocols::pdbinfo_manipulations::AlignPDBInfoToSequencesMode::multiple);
		TS_ASSERT_EQUALS(protocols::pdbinfo_manipulations::alignpdbinfotosequencesmode_to_string(protocols::pdbinfo_manipulations::AlignPDBInfoToSequencesMode::single), "single");
		TS_ASSERT_EQUALS(protocols::pdbinfo_manipulations::alignpdbinfotosequencesmode_to_string(protocols::pdbinfo_manipulations::AlignPDBInfoToSequencesMode::multiple), "multiple");
	}

	void test_real_example() {
		utility::vector1<protocols::pdbinfo_manipulations::SequenceSpecification> seq_specs;
		{
			utility::io::izstream instream("protocols/pdbinfo_manipulations/6q6h_seqs.json.gz");
			nlohmann::json json_result;
			instream >> json_result;
			utility::vector1<protocols::pdbinfo_manipulations::SequenceSpecification> const tmp_in(json_result.get<utility::vector1<protocols::pdbinfo_manipulations::SequenceSpecification>>());
			seq_specs.insert(seq_specs.end(),  tmp_in.begin(), tmp_in.end());
		}
		protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts;
		apts.set_mode_from_string("multiple");
		apts.set_target_sequences(seq_specs);

		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ) );

		core::pose::PoseOP const rcsb_pdb_pose(
			core::import_pose::pose_from_file(
			*residue_set,
			"protocols/pdbinfo_manipulations/6q6h.pdb.gz",
			false,
			core::import_pose::PDB_file
			)
		);

		core::pose::PoseOP const bad_aln_pdb_pose(
			core::import_pose::pose_from_file(
			*residue_set,
			"protocols/pdbinfo_manipulations/aln6q6h.pdb.gz",
			false,
			core::import_pose::PDB_file
			)
		);

		TS_ASSERT_DIFFERS(rcsb_pdb_pose->sequence(), bad_aln_pdb_pose->sequence());
		apts.apply(*bad_aln_pdb_pose);
		apts.apply(*rcsb_pdb_pose);
		std::set<std::tuple<char, int, char>> original, aligned;

		TS_ASSERT_EQUALS(rcsb_pdb_pose->size(), bad_aln_pdb_pose->size());
		for ( core::Size i=1; i<=rcsb_pdb_pose->size(); ++i ) {
			original.insert(std::make_tuple(rcsb_pdb_pose->pdb_info()->chain(i), rcsb_pdb_pose->pdb_info()->number(i), rcsb_pdb_pose->residue(i).name1()));
		}
		for ( core::Size i=1; i<=bad_aln_pdb_pose->size(); ++i ) {
			aligned.insert(std::make_tuple(bad_aln_pdb_pose->pdb_info()->chain(i), bad_aln_pdb_pose->pdb_info()->number(i), bad_aln_pdb_pose->residue(i).name1()));
		}
		for ( auto x : original ) {
			TS_ASSERT_EQUALS(aligned.count(x), 1);
		}
		for ( auto x : aligned ) {
			TS_ASSERT_EQUALS(original.count(x), 1);
		}
	}

	void test_SequenceSpecification() {
		std::string const base_sequence("CLGIGSCNDFAGCGYAVVCFW");
		std::string const   tobe_chains("AAAABBBBBBBBBCCCCCCDD"); // kept in line with split_chains_str
		std::string const    split_sequence("XXXXXXXXXXXXXXCLGIGSCNDFAGCGXXXXXXYAVVCFWXXXXX");
		std::string const  split_chains_str("XXXXXXXXXXXXXXAAAABBBBBBBBBCYYYYYYCCCCCDDZZZZZ");
		utility::vector1<std::string> const split_chains(to_string_vec(split_chains_str));
		utility::vector1<std::string> const split_idx(to_string_vec(split_sequence));

		{
			protocols::pdbinfo_manipulations::SequenceSpecification const ss(
				"AGGY",
				{"H", "H", "H", "H"},
				{"H", "H", "H", "H"},
				{"H", "H", "H", "H"}
			);
			TS_ASSERT_THROWS_NOTHING(ss.check_single_format());
		}
		{
			protocols::pdbinfo_manipulations::SequenceSpecification const ss(
				"AGGY",
				{"H", "H"},
				{"H", "H"},
				{"H", "H", "H", "H", "U"});
			TS_ASSERT_THROWS(ss.check_single_format(), utility::excn::BadInput&);
		}
		{
			protocols::pdbinfo_manipulations::SequenceSpecification const ss(
				"AGGY",
				{"H", "H", "H", "H", "U"});
			TS_ASSERT_THROWS(ss.check_single_format(), utility::excn::BadInput&);
		}
		{
			protocols::pdbinfo_manipulations::SequenceSpecification const ss(
				"AGGY",
				{"H", "H", "H", "H"},
				{"H", "H", "H", "H", "U"});
			TS_ASSERT_THROWS(ss.check_single_format(), utility::excn::BadInput&);
		}
		{
			protocols::pdbinfo_manipulations::SequenceSpecification const ss(
				"AGGY",
				{"H", "H", "H", "H"},
				{"H", "H", "H", "H"},
				{"H", "H", "H", "H", "U"});
			TS_ASSERT_THROWS(ss.check_single_format(), utility::excn::BadInput&);
		}
	}

	void test_throw_on_fail() {
		core::pose::Pose pose = pdb1rpb_pose();
		protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts;
		std::string const base_sequence("XLGIXSCNDFAGCGXAVXCFW");
		std::string const   tobe_chains("AAAABBBBBBBBBCCCCCCDD"); // kept in line with split_chains_str
		std::string const    split_sequence("XXXXXXXXXXXXXXCLGIGSCNDFAGCGXXXXXXYAVVCFWXXXXX");
		std::string const  split_chains_str("XXXXXXXXXXXXXXAAAABBBBBBBBBCYYYYYYCCCCCDDZZZZZ");
		utility::vector1<std::string> const split_chains(to_string_vec(split_chains_str));
		utility::vector1<std::string> const split_idx(to_string_vec(split_sequence));

		protocols::pdbinfo_manipulations::SequenceSpecification const ss(
			split_sequence,
			split_chains);

		apts.set_target_sequences(utility::vector1<protocols::pdbinfo_manipulations::SequenceSpecification>({ss}));
		apts.set_throw_on_fail(true);
		TS_ASSERT_THROWS(apts.apply(pose), utility::excn::BadInput&);
	}

	void test_single() {
		core::pose::Pose pose = pdb1rpb_pose();
		protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts;
		std::string const base_sequence("CLGIGSCNDFAGCGYAVVCFW");
		std::string const   tobe_chains("AAAABBBBBBBBBCCCCCCDD"); // kept in line with split_chains_str
		std::string const    split_sequence("XXXXXXXXXXXXXXCLGIGSCNDFAGCGXXXXXXYAVVCFWXXXXX");
		std::string const  split_chains_str("XXXXXXXXXXXXXXAAAABBBBBBBBBCYYYYYYCCCCCDDZZZZZ");
		utility::vector1<std::string> const split_chains(to_string_vec(split_chains_str));
		utility::vector1<std::string> const split_idx(to_string_vec(split_sequence));

		protocols::pdbinfo_manipulations::SequenceSpecification const ss(
			split_sequence,
			split_chains);

		apts.set_target_sequences(utility::vector1<protocols::pdbinfo_manipulations::SequenceSpecification>({ss}));
		TS_ASSERT_THROWS(apts.set_mode_from_string("failer"), utility::excn::BadInput&);
		apts.set_throw_on_fail(true);
		TS_ASSERT_THROWS(apts.apply(pose), utility::excn::BadInput&);
		apts.set_throw_on_fail(false);
		TS_ASSERT_THROWS_NOTHING(apts.set_mode_from_string("single"));
		TS_ASSERT_THROWS_NOTHING(apts.apply(pose));
		utility::vector1<core::Size> expected_numbers({
			15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 35, 36, 37, 38, 39, 40, 41});
		for ( core::Size i=1; i<=pose.size(); ++i ) {
			TS_ASSERT_EQUALS(pose.pdb_info()->number(i), expected_numbers[i]);
			TS_ASSERT_EQUALS(pose.pdb_info()->chain(i), tobe_chains.at(i-1));
		}
	}

	void test_multiple() {
		{ // test no doublers
			core::pose::PoseCOP const og_pdb_pose(
				utility::pointer::make_shared< core::pose::Pose >(
				*core::import_pose::pose_from_file(
				"protocols/pdbinfo_manipulations/5GRU_full.pdb",
				false,
				core::import_pose::PDB_file
				)
				)
			);
			core::pose::PoseOP const pose(
				core::import_pose::pose_from_file(
				"protocols/pdbinfo_manipulations/5GRU_full.pdb",
				false,
				core::import_pose::PDB_file
				)
			);

			std::string const chain_A(
				"KIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEIT"
				"PDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPL"
				"IAADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKV"
				"NYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELAKDPRIAAT"
				"MENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTLVPRGSAAALEHHHHHH");
			std::string const chain_H(
				"EVQLVESGGGLVQPGGSLRLSCAASGFNFSSSSIHWVRQAPGKGLEWVASISSSSGSTSYADSVKGRFTISADTSKNTAY"
				"LQMNSLTAEDTAVYYCARGYYYTGLWYPYAMYEFGMDYWGQGTLVTVSSGGGGSDIQLTQSTSSLPASLGDRVTISCRAG"
				"QDISNHLNWYQQKPDGTVKLLIYYTSRLHSGVPSRFSGSGSGTDYSLTISNLEQEDIATYFCQQGNTLPWTFGGGSKLEI"
				"KSRHHHHHH");
			std::string const chain_L(
				"QVQLKESGPGLVRPSQSLSLTCSVTGYSITSGYYWNWIRQFPGNKLEWMGYISYDGSNNYNPSLKGRISITRDTSKNQFF"
				"LKLNSVTTDDTATYYCARAYIGFAYWGQGTLVTVSSGGGGSDIQMTQSPSSLSASVGDRVTITCRASQSVSSAVAWYQQK"
				"PGKAPKLLIYSASSLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQSSSSLITFGQGTKVEIK");
			utility::vector1<std::string> const chain_L_chains({"Z"});

			utility::vector1<int> const L_residue_numbers = [&chain_L](){
				utility::vector1<int> ret;
				for ( int i = 0; i<int(chain_L.size()); ++i ) {
					if ( i < (int(chain_L.size())/2) ) { ret.push_back(i+1); }
					else { ret.push_back(i+100); }
				}
				return ret;
			}();

			protocols::pdbinfo_manipulations::SequenceSpecification const ss_A(
				chain_A,
				{"X"});
			std::cout << ss_A << std::endl;
			protocols::pdbinfo_manipulations::SequenceSpecification const ss_H(
				chain_H,
				{"Y"},
				utility::vector1<std::string>(),
				utility::vector1<std::string>(),
				{33});  // numbering should start at 33
			std::cout << ss_H << std::endl;
			protocols::pdbinfo_manipulations::SequenceSpecification const ss_L(
				chain_L,
				{"Z"},
				utility::vector1<std::string>(),
				utility::vector1<std::string>(),
				L_residue_numbers);
			std::cout << ss_L << std::endl;

			utility::vector1<protocols::pdbinfo_manipulations::SequenceSpecification> ss_inputs({
				ss_A, ss_H, ss_L});

			protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts;
			apts.set_mode_from_string("multiple");
			apts.set_target_sequences(ss_inputs);
			apts.apply(*pose);
			std::cout << apts.get_target_sequences() << std::endl;
			std::map<std::string, std::string> const old_to_new({{"A", "X"}, {"H", "Y"}, {"L", "Z"}});

			TS_ASSERT_EQUALS(pose->size(), og_pdb_pose->size());
			for ( core::Size i=1, X_count=0, Y_count=0, Z_count=0; i<=pose->size(); ++i ) {
				TS_ASSERT_EQUALS(
					std::string(1, pose->pdb_info()->chain(i)),
					old_to_new.at(std::string(1, og_pdb_pose->pdb_info()->chain(i)))
				);
				TS_ASSERT_DIFFERS(
					pose->pdb_info()->chain(i),
					og_pdb_pose->pdb_info()->chain(i)
				);
				if ( pose->pdb_info()->chain(i) == 'X' ) {
					++X_count;
					// this starts alignment at residue 7
					TS_ASSERT_EQUALS(X_count+6, pose->pdb_info()->number(i));
				} else if ( pose->pdb_info()->chain(i) == 'Y' ) {
					++Y_count;
					// this starts alignment at residue 2
					// so Y_count + 33(set start number)
					TS_ASSERT_EQUALS(Y_count+33, pose->pdb_info()->number(i));
				} else if ( pose->pdb_info()->chain(i) == 'Z' ) {
					++Z_count;
					TS_ASSERT_EQUALS(L_residue_numbers[Z_count], pose->pdb_info()->number(i));
				}
			}
		} // test no doublers
		{ // now try doubler
			core::pose::PoseOP const pose(
				core::import_pose::pose_from_file(
				"protocols/pdbinfo_manipulations/5GRU_full.pdb",
				false,
				core::import_pose::PDB_file
				)
			);
			core::Size const og_single_size(pose->size());

			pose->append_pose_by_jump(*pose->clone(), 1);

			core::pose::PoseCOP const og_pdb_pose(
				utility::pointer::make_shared< core::pose::Pose >(*pose)
			);

			std::string const chain_A(
				"KIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEIT"
				"PDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPL"
				"IAADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKV"
				"NYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELAKDPRIAAT"
				"MENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTLVPRGSAAALEHHHHHH");
			std::string const chain_H(
				"EVQLVESGGGLVQPGGSLRLSCAASGFNFSSSSIHWVRQAPGKGLEWVASISSSSGSTSYADSVKGRFTISADTSKNTAY"
				"LQMNSLTAEDTAVYYCARGYYYTGLWYPYAMYEFGMDYWGQGTLVTVSSGGGGSDIQLTQSTSSLPASLGDRVTISCRAG"
				"QDISNHLNWYQQKPDGTVKLLIYYTSRLHSGVPSRFSGSGSGTDYSLTISNLEQEDIATYFCQQGNTLPWTFGGGSKLEI"
				"KSRHHHHHH");
			std::string const chain_L(
				"QVQLKESGPGLVRPSQSLSLTCSVTGYSITSGYYWNWIRQFPGNKLEWMGYISYDGSNNYNPSLKGRISITRDTSKNQFF"
				"LKLNSVTTDDTATYYCARAYIGFAYWGQGTLVTVSSGGGGSDIQMTQSPSSLSASVGDRVTITCRASQSVSSAVAWYQQK"
				"PGKAPKLLIYSASSLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQSSSSLITFGQGTKVEIK");
			utility::vector1<std::string> const chain_L_chains({"Z"});

			protocols::pdbinfo_manipulations::SequenceSpecification const ss_A(
				chain_A,
				{"O", "P"},
				{"sqA1", "sqA2"});
			protocols::pdbinfo_manipulations::SequenceSpecification const ss_H(
				chain_H,
				{"Q", "R"},
				{"sqH1", "sqH2"});
			protocols::pdbinfo_manipulations::SequenceSpecification const ss_L(
				chain_L,
				{"T", "U"},
				{"sqL1", "sqL2"});

			utility::vector1<protocols::pdbinfo_manipulations::SequenceSpecification> ss_inputs({
				ss_A, ss_H, ss_L});

			protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts;
			apts.set_mode_from_string("multiple");
			apts.set_target_sequences(ss_inputs);
			apts.apply(*pose);
			std::map<std::string, utility::vector1<std::string>> const chain_old_to_new({{"A", {"O", "P"}}, {"H", {"Q", "R"}}, {"L", {"T", "U"}}});
			std::map<std::string, std::string> const chain_to_new_segid({{"O", "sqA1"}, {"P", "sqA2"}, {"Q", "sqH1"}, {"R", "sqH2"}, {"T", "sqL1"}, {"U", "sqL2"}});

			TS_ASSERT_EQUALS(pose->size(), og_pdb_pose->size());
			core::Size current_pos(1);
			for ( core::Size i=1; i<=pose->size(); ++i ) {
				TS_ASSERT_EQUALS(
					std::string(1, pose->pdb_info()->chain(i)),
					chain_old_to_new.at(std::string(1, og_pdb_pose->pdb_info()->chain(i)))[current_pos]
				);
				TS_ASSERT_DIFFERS(
					pose->pdb_info()->chain(i),
					og_pdb_pose->pdb_info()->chain(i)
				);

				TS_ASSERT_EQUALS(
					pose->pdb_info()->segmentID(i),
					chain_to_new_segid.at(std::string(1, pose->pdb_info()->chain(i)))
				);
				TS_ASSERT_DIFFERS(
					pose->pdb_info()->segmentID(i),
					og_pdb_pose->pdb_info()->segmentID(i)
				);

				if ( i == og_single_size ) {
					++current_pos;
				}
			}

		} // end doubler
	}

	void test_ctor_and_operators() {
		protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts1;
		TS_ASSERT_EQUALS(apts1.get_mode(), protocols::pdbinfo_manipulations::AlignPDBInfoToSequencesMode::unset);
		protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts2;
		apts2.set_mode_from_string("single");
		TS_ASSERT_EQUALS(apts2.get_mode(), protocols::pdbinfo_manipulations::AlignPDBInfoToSequencesMode::single);
	}

	void test_parse_target_tag() {
		protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts;
		apts.set_mode_from_string("single");
		{
			protocols::pdbinfo_manipulations::SequenceSpecification const t_ss(
				"MKVKIKCWNGVATWLWVANDENCGICRMAFNGCCPDCKVPGDDCPLVWGQCSHCFHMHCILKWLHAQQVQQHCPMCRQEWKFKE",
				{"C", "D"});
			utility::tag::TagCOP const tag = tagptr_from_string("<Target sequence=\"MKVKIKCWNGVATWLWVANDENCGICRMAFNGCCPDCKVPGDDCPLVWGQCSHCFHMHCILKWLHAQQVQQHCPMCRQEWKFKE\" chains=\"C,D\" />\n");
			protocols::pdbinfo_manipulations::SequenceSpecification const ss(apts.parse_target_tag(tag));
			TS_ASSERT_EQUALS(t_ss, ss);
		}
		{
			protocols::pdbinfo_manipulations::SequenceSpecification const t_ss(
				"MKVKIKCWNGVATWLWVANDENCGICRMAFNGCCPDCKVPGDDCPLVWGQCSHCFHMHCILKWLHAQQVQQHCPMCRQEWKFKE",
				{"C", "D"}, {"yyY", "XSd"}, {"A"}, {22,33,44});
			utility::tag::TagCOP const tag = tagptr_from_string("<Target sequence=\"MKVKIKCWNGVATWLWVANDENCGICRMAFNGCCPDCKVPGDDCPLVWGQCSHCFHMHCILKWLHAQQVQQHCPMCRQEWKFKE\" chains=\"C,D\" segmentIDs=\"yyY,XSd\" insCodes=\"A\" residue_numbers=\"22,33, 44\" />\n");
			protocols::pdbinfo_manipulations::SequenceSpecification const ss(apts.parse_target_tag(tag));
			TS_ASSERT_EQUALS(t_ss, ss);
		}
	}

	void test_dna_align() {
		core::pose::PoseOP const pose(
			core::import_pose::pose_from_file(
			"protocols/pdbinfo_manipulations/4gatA.pdb",
			false,
			core::import_pose::PDB_file
			)
		);

		// this should trigger a failure to align instance (a throw in source/src/core/sequence/Aligner.cc) which we should handle gracefully
		// I'm not sure if there's a way to prove that this is happening... but it is if you check the debug logs I swear!
		protocols::pdbinfo_manipulations::SequenceSpecification const ss(
			"MKVKIKCWNGVATWLWVANDENCGICRMAFNGCCPDCKVPGDDCPLVWGQCSHCFHMHCILKWLHAQQVQQHCPMCRQEWKFKE",
			{"C", "D"},
			utility::vector1<std::string>(),
			utility::vector1<std::string>(),
			{100});
		protocols::pdbinfo_manipulations::AlignPDBInfoToSequences apts;
		apts.set_mode_from_string("multiple");
		apts.set_target_sequences({ss});
		TS_ASSERT_THROWS_NOTHING(apts.apply(*pose));
	}
};
