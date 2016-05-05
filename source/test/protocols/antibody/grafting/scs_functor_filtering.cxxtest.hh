// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/grafting/scs_functor_filtering.cxxtest.hh
/// @brief  test suite for antibody scs_functor filtering code
/// @author Jeliazko Jeliazkov

#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <protocols/antibody/grafting/regex_based_cdr_detection.hh>
#include <protocols/antibody/grafting/scs_functor.hh>
#include <protocols/antibody/grafting/scs_blast.hh>
#include <protocols/antibody/grafting/util.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

using std::string;

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.Scs_functor_filtering_tests");
class Scs_functor_filtering_tests : public CxxTest::TestSuite
{
public:
	#ifdef __ANTIBODY_GRAFTING__
	// 2BRR, default values for testing
	// CDRs from SAbDab
	// H1: GYTFTNY
	// H2: NTYTGE
	// H3: DYYGSTYPYYAMDY
	// L1: RASSSVSSSYLH
	// L2: STSNLAS
	// L3: QQYSGYPYT
	string heavy_chain_sequence = "VQLEQSGPELKKPGETVKISCKASGYTFTNYGMNWVKQAPGKGLKWMGWINTYTGEPTYADDFKERFAFSLETSASAAYLQINNLKNEDTATYFCARDYYGSTYPYYAMDYWGQGTTVTVSS";
	string light_chain_sequence = "NVLTQSPAIMSASPGEKVTMTCRASSSVSSSYLHWYQQKSGASPKLWIYSTSNLASGVPARFSGSGSGTSYSLTISSVEAEDAATYYCQQYSGYPYTFGGGTKLEIK";
	#endif //__ANTIBODY_GRAFTING__
	
  void setUp() {
    // shared initialization goes here
		// need options for a specific filter...
		core_init_with_additional_options( "-antibody:grafting:exclude_homologs true -antibody:grafting:exclude_homologs_fr_cutoff 60 -out:level 500 -antibody:grafting:multi_template_graft true -antibody:grafting:n_multi_templates 3 -antibody:grafting:ocd_cutoff 5" );
		// sets sid on and cutoff to 80%
		// sets sid cutoff to 60% for FR
  } //setUp

  void tearDown() {
    // shared finalization goes here
  } //tearDown

  /// future filter tests should be generalized and only take filter type as input
  /// @brief specific filter test for SCS_BlastFilter_by_sequence_length
	void test_filter_by_sequence_length() {
	#ifdef __ANTIBODY_GRAFTING__
		using namespace protocols::antibody::grafting;
		
		// set up dummy AntibodySequence
		AntibodySequence as = AntibodySequence(heavy_chain_sequence, light_chain_sequence);
		RegEx_based_CDR_Detector().detect(as);

		// set up two BlastResults, a good one that should pass filtering and a dummy that should fail
		// in the sequence length filter, the dummy should have an length different from the good
		SCS_BlastResultOP dummy_br(new SCS_BlastResult);
		SCS_BlastResultOP good_br(new SCS_BlastResult);
		
		dummy_br->h1 = "GYTFTNYGH";
		dummy_br->h2 = "NTYTG";
		dummy_br->h3 = "DANHERTLK";
		dummy_br->frh = string(91,'X');
		dummy_br->l1 = "KDH";
		dummy_br->l2 = "STSNLASDHTPT";
		dummy_br->l3 = "QQYSGYPYTDDHN";
		dummy_br->frl = string(74,'X');
		dummy_br->pdb = "dummy";
		dummy_br->alignment_length = 92;
		dummy_br->bit_score = 0.3;
		dummy_br->bio_type = "human";
		dummy_br->light_type = "kappa";
		dummy_br->struct_source = "xray";

		
		good_br->h1 = "GYTFTNYGMN";
		good_br->h2 = "WINTYTGEPTYADDFKE";
		good_br->h3 = "DYYGSTYPYYAMDY";
		good_br->frh = string(63,'X'); // filter is either 61 if PDB ID is 2x7l or 63 otherwise
		good_br->l1 = "RASSSVSSSYLH";
		good_br->l2 = "STSNLAS";
		good_br->l3 = "QQYSGYPYT";
		good_br->frl = string(58,'X'); // either 60 for PDB ID 3h0t or 58
		good_br->pdb = "good";
		good_br->alignment_length = 96;
		good_br->bit_score = 0.001;
		good_br->bio_type = "human";
		good_br->light_type = "lambda";
		good_br->struct_source = "xray";
		
		// store BlastResults in ResultVector;
		SCS_ResultVector results;
		results.push_back(good_br);
		results.push_back(dummy_br);
		
		// store ResultVector in SCS_Results
		SCS_ResultsOP r(new SCS_Results);
		r->h1 = results;
		r->h2 = results;
		r->h3 = results;
		r->l1 = results;
		r->l2 = results;
		r->l3 = results;
		r->orientation = results; // what is stored here normally?
		r->frh = results;
		r->frl = results;

		// now apply single filter to ResultsOP and test
		SCS_FunctorCOP scs_filter( new SCS_BlastFilter_by_sequence_length );
		scs_filter->apply(as,r);
		
		// now test, the top result should be "good" whereas "dummy" should have been filtered out (i.e. len(r->x)=0)
		// orientation should remain unchanged.
		TS_ASSERT_EQUALS(r->h1.size(),1)
		TS_ASSERT_EQUALS(r->h2.size(),1)
		TS_ASSERT_EQUALS(r->h3.size(),1)
		TS_ASSERT_EQUALS(r->frh.size(),1)
		TS_ASSERT_EQUALS(r->l1.size(),1)
		TS_ASSERT_EQUALS(r->l2.size(),1)
		TS_ASSERT_EQUALS(r->l3.size(),1)
		TS_ASSERT_EQUALS(r->frl.size(),1)
		TS_ASSERT_EQUALS(r->orientation.size(),2)
		
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h1[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h2[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h3[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->frh[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l1[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l2[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l3[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->frl[0]->pdb),0);
		#endif //__ANTIBODY_GRAFTING__
	}
	
	/// @brief specific filter test for SCS_BlastFilter_by_alignment_length
	void test_filter_by_alignment_length() {
		#ifdef __ANTIBODY_GRAFTING__
		using namespace protocols::antibody::grafting;
		
		// set up dummy AntibodySequence
		AntibodySequence as = AntibodySequence(heavy_chain_sequence, light_chain_sequence);
		RegEx_based_CDR_Detector().detect(as);
		
		// set up two BlastResults, a good one that should pass filtering and a dummy that should fail
		// this filter only tests for alignment_length > const*region.size()
		// there is one alignment length per blast result, but six cdrs
		SCS_BlastResultOP dummy_br(new SCS_BlastResult);
		SCS_BlastResultOP good_br(new SCS_BlastResult);

		dummy_br->pdb = "dummy";
		good_br->pdb="good";
		
		// setting alignment length to 5 defines filter as
		// 5 < 0.7*len(h1)
		// 5 < 0.55*len(h2)
		// 5 < 0.1*len(h3)
		// 5 < 0.7*len(11)
		// 5 < 0.7*len(12)
		// 5 < 0.7*len(13)
		
		dummy_br->alignment_length=5;
		good_br->alignment_length=5;

		dummy_br->h1 = string(10,'X');
		dummy_br->h2 = string(11,'X');
		dummy_br->h3 = string(52,'X');
		
		dummy_br->l1 = string(12,'X');
		dummy_br->l2 = string(14,'X');
		dummy_br->l3 = string(15,'X');
		
		good_br->h1 = string(7,'X');
		good_br->h2 = string(8,'X');
		good_br->h3 = string(13,'X');
		
		good_br->l1 = string(5,'X');
		good_br->l2 = string(4,'X');
		good_br->l3 = string(1,'X');
		
		
		// store BlastResults in ResultVector;
		SCS_ResultVector results;
		results.push_back(dummy_br);
		results.push_back(good_br);
		
		// store ResultVector in SCS_Results
		SCS_ResultsOP r(new SCS_Results);
		r->h1 = results;
		r->h2 = results;
		r->h3 = results;
		r->l1 = results;
		r->l2 = results;
		r->l3 = results;
		r->orientation = results; // what is stored here normally?
		r->frh = results;
		r->frl = results;
		
		// now apply single filter to ResultsOP and test
		SCS_FunctorCOP scs_filter( new SCS_BlastFilter_by_alignment_length );
		scs_filter->apply(as,r);
		
		// now test, the top result should be "good" whereas "dummy" should have been filtered out (i.e. len(r->x)=0)
		// orientation should remain unchanged.
		TS_ASSERT_EQUALS(r->h1.size(),1)
		TS_ASSERT_EQUALS(r->h2.size(),1)
		TS_ASSERT_EQUALS(r->h3.size(),1)
		TS_ASSERT_EQUALS(r->frh.size(),2)
		TS_ASSERT_EQUALS(r->l1.size(),1)
		TS_ASSERT_EQUALS(r->l2.size(),1)
		TS_ASSERT_EQUALS(r->l3.size(),1)
		TS_ASSERT_EQUALS(r->frl.size(),2)
		TS_ASSERT_EQUALS(r->orientation.size(),2)
		
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h1[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h2[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h3[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l1[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l2[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l3[0]->pdb),0);
		#endif //__ANTIBODY_GRAFTING__
	}

	/// @brief specific filter test for SCS_BlastFilter_by_template_resolution
	void test_filter_by_template_resolution() {
		#ifdef __ANTIBODY_GRAFTING__
		using namespace protocols::antibody::grafting;
		
		// set up dummy AntibodySequence
		AntibodySequence as = AntibodySequence(heavy_chain_sequence, light_chain_sequence);
		RegEx_based_CDR_Detector().detect(as);
		
		// set up two BlastResults, a good one that should pass filtering and a dummy that should fail
		// this filter is easy to test as it elimnates templates with resolution worse than 2.8
		SCS_BlastResultOP dummy_br(new SCS_BlastResult);
		SCS_BlastResultOP good_br(new SCS_BlastResult);

		dummy_br->pdb = "dummy";
		good_br->pdb="good";
		
		dummy_br->resolution=2.9;
		good_br->resolution=2.1;
		
		// store BlastResults in ResultVector;
		SCS_ResultVector results;
		results.push_back(dummy_br);
		results.push_back(good_br);
		
		// store ResultVector in SCS_Results
		SCS_ResultsOP r(new SCS_Results);
		r->h1 = results;
		r->h2 = results;
		r->h3 = results;
		r->l1 = results;
		r->l2 = results;
		r->l3 = results;
		r->orientation = results; // what is stored here normally?
		r->frh = results;
		r->frl = results;
		
		// now apply single filter to ResultsOP and test
		// spoof flag/option since default sequence id cutoff is 0
		SCS_FunctorCOP scs_filter( new SCS_BlastFilter_by_template_resolution );
		scs_filter->apply(as,r);
		
		// now test, the top result should be "good" whereas "dummy" should have been filtered out (i.e. len(r->x)=0)
		// orientation should remain unchanged.
		TS_ASSERT_EQUALS(r->h1.size(),1)
		TS_ASSERT_EQUALS(r->h2.size(),1)
		TS_ASSERT_EQUALS(r->h3.size(),1)
		TS_ASSERT_EQUALS(r->frh.size(),1)
		TS_ASSERT_EQUALS(r->l1.size(),1)
		TS_ASSERT_EQUALS(r->l2.size(),1)
		TS_ASSERT_EQUALS(r->l3.size(),1)
		TS_ASSERT_EQUALS(r->frl.size(),1)
		TS_ASSERT_EQUALS(r->orientation.size(),1)
		
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h1[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h2[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h3[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l1[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l2[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l3[0]->pdb),0);
		#endif //__ANTIBODY_GRAFTING__
	}
	
	/// @brief specific filter test for SCS_BlastFilter_by_sequence_identity
	void test_filter_by_sequence_identity() {
		#ifdef __ANTIBODY_GRAFTING__
		using namespace protocols::antibody::grafting;
		
		// set up dummy AntibodySequence
		AntibodySequence as = AntibodySequence(heavy_chain_sequence, light_chain_sequence);
		RegEx_based_CDR_Detector().detect(as);
		
		// set up two BlastResults, a good one that should pass filtering and a dummy that should fail
		// this filter is annoying to test, it compares two sequences vs. a cutoff ( 80.0 ) so
		// to test we generate correct sequences which will be eliminated for being too identical
		SCS_BlastResultOP dummy_br(new SCS_BlastResult);
		SCS_BlastResultOP good_br(new SCS_BlastResult);
		
		dummy_br->h1 = "GYTFTNYVMN"; // 90% id
		dummy_br->h2 = "WINTYTGEPVYADDFKE";
		dummy_br->h3 = "DYYGSTYPFFAMDY";
		dummy_br->frh = "VQLEQSGPELKKPGETVKISCKASWVKQAPGKGLKWMGRFAFSLETSASAAYLQINNLKNEDTATYFCARWGQGTTVTVSS";
		dummy_br->l1 = "RASSSVSTSYLH";
		dummy_br->l2 = "STSNLAS";
		dummy_br->l3 = "QQYSGYPYT";
		dummy_br->frl = "NVLTQSPAIMSASPGEKVTMTCWYNNKSGASPKLWIYGVPARFSGSGSGTSYSLTISSVEAEDAATYYCFGGGTKLEIK";
		dummy_br->pdb = "dummy";
		dummy_br->alignment_length = 92;
		dummy_br->bit_score = 0.3;
		dummy_br->bio_type = "human";
		dummy_br->light_type = "kappa";
		dummy_br->struct_source = "xray";
		
		good_br->h1 = "VYTFTNYVMN"; // filtered on 80% id
		good_br->h2 = "HGTNVCSF"; // filtered on length (i.e. id =0)
		good_br->h3 = "DYYGSTYY"; // filtered on length
		good_br->frh = "ADPDNCTFVLRKTVETGLNDTRISYVKNVGRRPVNYWWWFAFSLETSDCPGIQILQMVPRDMNLDGDGVYTVVLERDAPHH"; // filtered on id--must 60% or less identical
		good_br->l1 = "KDHNVTTTR"; // filtered on length
		good_br->l2 = "LPTNLAS"; // filtered on id
		good_br->l3 = "WINRARIIT"; // filtered on id
		good_br->frl = string(58,'X'); //filtered on length
		good_br->pdb = "good";
		good_br->alignment_length = 96;
		good_br->bit_score = 0.001;
		good_br->bio_type = "human";
		good_br->light_type = "lambda";
		good_br->struct_source = "xray";
		
		// store BlastResults in ResultVector;
		SCS_ResultVector results;
		results.push_back(dummy_br);
		results.push_back(good_br);
		
		// store ResultVector in SCS_Results
		SCS_ResultsOP r(new SCS_Results);
		r->h1 = results;
		r->h2 = results;
		r->h3 = results;
		r->l1 = results;
		r->l2 = results;
		r->l3 = results;
		r->orientation = results; // what is stored here normally?
		r->frh = results;
		r->frl = results;
		
		// now apply single filter to ResultsOP and test
		// options for this filter are set in setUp/core_init
		SCS_FunctorCOP scs_filter( new SCS_BlastFilter_by_sequence_identity);
		scs_filter->apply(as,r);
		
		// now test, the top result should be "good" whereas "dummy" should have been filtered out (i.e. len(r->x)=0)
		// orientation should remain unchanged.
		TS_ASSERT_EQUALS(r->h1.size(),1)
		TS_ASSERT_EQUALS(r->h2.size(),1)
		TS_ASSERT_EQUALS(r->h3.size(),1)
		TS_ASSERT_EQUALS(r->frh.size(),1)
		TS_ASSERT_EQUALS(r->l1.size(),1)
		TS_ASSERT_EQUALS(r->l2.size(),1)
		TS_ASSERT_EQUALS(r->l3.size(),1)
		TS_ASSERT_EQUALS(r->frl.size(),1)
		TS_ASSERT_EQUALS(r->orientation.size(),2)
		
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h1[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h2[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->h3[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l1[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l2[0]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l3[0]->pdb),0);
		#endif //__ANTIBODY_GRAFTING__
	}
	
	/// @brief specific filter test for SCS_BlastFilter_by_outlier
	void test_filter_by_outlier() {
		#ifdef __ANTIBODY_GRAFTING__
		using namespace protocols::antibody::grafting;
		
		// set up dummy AntibodySequence
		AntibodySequence as = AntibodySequence(heavy_chain_sequence, light_chain_sequence);
		RegEx_based_CDR_Detector().detect(as);
		
		// set up two BlastResults, a good one that should pass filtering and a dummy that should fail
		// this filter is easy to test as it elimnates pdbs belonging to an outlier list in /info/outlier_list
		// to test get one outlier pdb per region and one good pdb so each region will have 7 passes and a fail?
		utility::vector0<SCS_BlastResultOP> dummy_brs;
		for (int i=0; i<7; ++i) {
			SCS_BlastResultOP br(new SCS_BlastResult);
			br->pdb = "dummy" + std::to_string(i);
			dummy_brs.push_back(br);
		}
		
		SCS_BlastResultOP good_br(new SCS_BlastResult);
		good_br->pdb="good";
		
		// set pdbs
		dummy_brs[0]->pdb="1hkl"; // h1 outlier
		dummy_brs[1]->pdb="1dl7"; // h2 outlier, also frh & frl
		dummy_brs[2]->pdb="1dl7"; // frh outlier, also frl & h2 outlier
		dummy_brs[3]->pdb="1it9"; // l1 outlier
		dummy_brs[4]->pdb="1gig"; // l2 outlier
		dummy_brs[5]->pdb="1nj9"; // l3 outlier
		dummy_brs[6]->pdb="1dl7"; // frl outlier, also frh & h2 outlier
		good_br->pdb="1gpo"; // no outliers
		
		// store BlastResults in ResultVector;
		SCS_ResultVector results;
		for (int i=0; i<7; ++i) {
			results.push_back(dummy_brs[i]);
		}
		results.push_back(good_br);
		
		// store ResultVector in SCS_Results
		SCS_ResultsOP r(new SCS_Results);
		r->h1 = results;
		r->h2 = results;
		r->l1 = results;
		r->l2 = results;
		r->l3 = results;
		r->frh = results;
		r->frl = results;
		
		// now apply single filter to ResultsOP and test
		// spoof flag/option since default sequence id cutoff is 0
		SCS_FunctorCOP scs_filter( new SCS_BlastFilter_by_outlier );
		scs_filter->apply(as,r);
		
		// now test, the top result should be "good" whereas "dummy" should have been filtered out (i.e. len(r->x)=0)
		TS_ASSERT_EQUALS(r->h1.size(),7)
		TS_ASSERT_EQUALS(r->h2.size(),5)
		TS_ASSERT_EQUALS(r->frh.size(),5)
		TS_ASSERT_EQUALS(r->l1.size(),7)
		TS_ASSERT_EQUALS(r->l2.size(),7)
		TS_ASSERT_EQUALS(r->l3.size(),7)
		TS_ASSERT_EQUALS(r->frl.size(),5)
		
		// some dummies are ok for other regions, let's see if they remain as well as the good
		TS_ASSERT_EQUALS(dummy_brs[2]->pdb.compare(r->h1[1]->pdb),0);
		TS_ASSERT_EQUALS(dummy_brs[4]->pdb.compare(r->h2[2]->pdb),0);
		TS_ASSERT_EQUALS(dummy_brs[5]->pdb.compare(r->frh[3]->pdb),0);
		TS_ASSERT_EQUALS(dummy_brs[5]->pdb.compare(r->l1[4]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l2[6]->pdb),0);
		TS_ASSERT_EQUALS(dummy_brs[1]->pdb.compare(r->l3[1]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->frl[4]->pdb),0);
		#endif //__ANTIBODY_GRAFTING__
	}
	
	/// @brief specific filter test for SCS_BlastFilter_by_template_bfactor
	void test_filter_by_template_bfactor() {
		#ifdef __ANTIBODY_GRAFTING__
		using namespace protocols::antibody::grafting;
		
		// set up dummy AntibodySequence
		AntibodySequence as = AntibodySequence(heavy_chain_sequence, light_chain_sequence);
		RegEx_based_CDR_Detector().detect(as);
		
		// set up two BlastResults, a good one that should pass filtering and a dummy that should fail
		// this filter is easy to test as it elimnates pdbs belonging to an outlier list in /info/list_bfactor50
		// to test get one outlier pdb per region and one good pdb so each region will have 6 passes and a fail?
		utility::vector0<SCS_BlastResultOP> dummy_brs;
		for (int i=0; i<6; ++i) {
			SCS_BlastResultOP br(new SCS_BlastResult);
			br->pdb = "dummy" + std::to_string(i);
			dummy_brs.push_back(br);
		}
		
		SCS_BlastResultOP good_br(new SCS_BlastResult);
		good_br->pdb="good";
		
		// set pdbs
		dummy_brs[0]->pdb="1a4k"; // h1 outlier
		dummy_brs[1]->pdb="3skj"; // h2 outlier
		dummy_brs[2]->pdb="15c8"; // h3 outlier
		dummy_brs[3]->pdb="1bbd"; // l1 outlier
		dummy_brs[4]->pdb="1ahw"; // l2 & h2 outlier
		dummy_brs[5]->pdb="1baf"; // l3 outlier
		good_br->pdb="12e8"; // no outliers
		
		// store BlastResults in ResultVector;
		SCS_ResultVector results;
		for (int i=0; i<6; ++i) {
			results.push_back(dummy_brs[i]);
		}
		results.push_back(good_br);
		
		// store ResultVector in SCS_Results
		SCS_ResultsOP r(new SCS_Results);
		r->h1 = results;
		r->h2 = results;
		r->h3 = results;
		r->l1 = results;
		r->l2 = results;
		r->l3 = results;
		
		// now apply single filter to ResultsOP and test
		SCS_FunctorCOP scs_filter( new SCS_BlastFilter_by_template_bfactor );
		scs_filter->apply(as,r);
		
		// now test, the top result should be "good" whereas "dummy" should have been filtered out
		TS_ASSERT_EQUALS(r->h1.size(),6)
		TS_ASSERT_EQUALS(r->h2.size(),5)
		TS_ASSERT_EQUALS(r->h3.size(),6)
		TS_ASSERT_EQUALS(r->l1.size(),6)
		TS_ASSERT_EQUALS(r->l2.size(),6)
		TS_ASSERT_EQUALS(r->l3.size(),6)
		
		// some dummies are ok for other regions, let's see if they remain as well as the good
		TS_ASSERT_EQUALS(dummy_brs[1]->pdb.compare(r->h1[0]->pdb),0);
		TS_ASSERT_EQUALS(dummy_brs[5]->pdb.compare(r->h2[3]->pdb),0);
		TS_ASSERT_EQUALS(dummy_brs[4]->pdb.compare(r->h3[3]->pdb),0);
		TS_ASSERT_EQUALS(dummy_brs[0]->pdb.compare(r->l1[0]->pdb),0);
		TS_ASSERT_EQUALS(dummy_brs[5]->pdb.compare(r->l2[4]->pdb),0);
		TS_ASSERT_EQUALS(good_br->pdb.compare(r->l3[5]->pdb),0);
		#endif //__ANTIBODY_GRAFTING__
	}
	
	/// @brief specific filter test for SCS_BlastFilter_by_OCD
	void test_filter_by_OCD() {
		#ifdef __ANTIBODY_GRAFTING__
		using namespace protocols::antibody::grafting;
		
		// set up dummy AntibodySequence
		AntibodySequence as = AntibodySequence(heavy_chain_sequence, light_chain_sequence);
		RegEx_based_CDR_Detector().detect(as);
		
		// set up two BlastResult vectors, a good one that should pass filtering and a dummy that should fail
		// intercalate good/bad (1st should be good, bad is always a below the distance cutoff
		utility::vector0<SCS_BlastResultOP> dummy_brs;
		utility::vector0<SCS_BlastResultOP> good_brs;
		
		for (int i=0; i<2; ++i) {
			SCS_BlastResultOP br(new SCS_BlastResult);
			br->pdb = "dummy" + std::to_string(i);
			dummy_brs.push_back(br);
		}
		
		for (int i=0; i<3; ++i) {
			SCS_BlastResultOP br(new SCS_BlastResult);
			br->pdb = "good" + std::to_string(i);
			good_brs.push_back(br);
		}
		
		// set pdbs good are OCD 5 away from each other, whereas "bad" will be filtered out
		good_brs[0]->pdb="1a2y"; // everything is compared to this
		dummy_brs[0]->pdb="15c8"; // OCD = 0.683
		good_brs[1]->pdb="12e8"; // OCD = 13.4
		dummy_brs[1]->pdb="1a14"; // OCD = 3.19
		good_brs[2]->pdb="1a3r"; // OCD = 8.74
		
		// store BlastResults in ResultVector;
		SCS_ResultVector results;
		results.push_back(good_brs[0]);
		results.push_back(dummy_brs[0]);
		results.push_back(good_brs[1]);
		results.push_back(dummy_brs[1]);
		results.push_back(good_brs[2]);
		
		// store ResultVector in SCS_Results
		SCS_ResultsOP r(new SCS_Results);
		r->orientation = results;
		
		// now apply single filter to ResultsOP and test
		// spoof flag/option since default sequence id cutoff is 0
		SCS_FunctorCOP scs_filter( new SCS_BlastFilter_by_OCD );
		scs_filter->apply(as,r);
		
		// now test, the top result should be "good" whereas "dummy" should have been filtered out (i.e. len(r->x)=0)
		TS_ASSERT_EQUALS(r->orientation.size(),3)
		
		// some dummies are ok for other regions, let's see if they remain as well as the good
		TS_ASSERT_EQUALS(good_brs[0]->pdb.compare(r->orientation[0]->pdb),0);
		TS_ASSERT_EQUALS(good_brs[1]->pdb.compare(r->orientation[1]->pdb),0);
		TS_ASSERT_EQUALS(good_brs[2]->pdb.compare(r->orientation[2]->pdb),0);
		#endif //__ANTIBODY_GRAFTING__
	}
};
