// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/grafting/regex_based_cdr_detection.cxxtest.hh
/// @brief  test suite for antibody regex_based_cdr_detection code and framework-trimming for Blast+ SCS
/// @author Sergey Lyskov

#include <protocols/antibody/grafting/regex_based_cdr_detection.hh>
#include <protocols/antibody/grafting/scs_blast.hh>
#include <protocols/antibody/grafting/util.hh>

#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/string_util.hh>

#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>

#include <fstream>

using std::string;

/* // Example data for tests debugging
string json_data = R"_qwe_(
{
"2brr" : {
"heavy_chain_sequence" : "AVQLEQSGPELKKPGETVKISCKASGYTFTNYGMNWVKQAPGKGLKWMGWINTYTGEPTYADDFKERFAFSLETSASAAYLQINNLKNEDTATYFCARDYYGSTYPYYAMDYWGQGTTVTVSSAKTTAPSVYPLAPVCGDTTGSSVTLGCLVKGYFPEPVTLTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVTSSTWPSQSITCNVAHPASSTKVDKKIEPRGG",
"light_chain_sequence" : "ENVLTQSPAIMSASPGEKVTMTCRASSSVSSSYLHWYQQKSGASPKLWIYSTSNLASGVPARFSGSGSGTSYSLTISSVEAEDAATYYCQQYSGYPYTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC",

"h1": "GYTFTNYGMN",
"h2": "WINTYTGEPTYADDFKE",
"h3": "DYYGSTYPYYAMDY",
"l1": "RASSSVSSSYLH",
"l2": "STSNLAS",
"l3": "QQYSGYPYT",
"frh": "ELKKPGETVKISCKASWVKQKWMGRFAFSLETSASAAYLQINNLKNEDTATYFCARWGQGTTV",
"frl": "IMSASPGEKVTMTCWYQQKLWIYGVPARFSGSGYSLTISSVEAEDAATYYCFGGGTKL",
"heavy": "AVQLEQSGPELKKPGETVKISCKASGYTFTNYGMNWVKQAPGKGLKWMGWINTYTGEPTYADDFKERFAFSLETSASAAYLQINNLKNEDTATYFCARDYYGSTYPYYAMDYWGQG",
"light": "ENVLTQSPAIMSASPGEKVTMTCRASSSVSSSYLHWYQQKSGASPKLWIYSTSNLASGVPARFSGSGSGTSYSLTISSVEAEDAATYYCQQYSGYPYTFGGG"
},
"2uyl" : {
"heavy_chain_sequence" : "QVQLEQPGAELVKPGASVKLSCKASGYTFTSNWINWVKQRPGQGLEWIGHISPGSSSTNYNEKFKSKATLTVDTSSSTAYMQLSSLTSDDSAVYYCGREETVRASFGNWGQGTLVTVSAAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWPSQTVTCNVAHPASSTKVDKKIVPRDC",
"light_chain_sequence" : "ELVMTQTPPSLPVSLGDQASISCRSSQSIVHSNGDTYLEWYLQKPGQSPKLLIYKVSNRFSGVPDRFSGSGSGTDFTLEISRVEAEDLGVYYCFQGSHVPRTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC",

"h1": "GYTFTSNWIN",
"h2": "HISPGSSSTNYNEKFKS",
"h3": "EETVRASFGN",
"l1": "RSSQSIVHSNGDTYLE",
"l2": "KVSNRFS",
"l3": "FQGSHVPRT",
"frh": "ELVKPGASVKLSCKASWVKQEWIGKATLTVDTSSSTAYMQLSSLTSDDSAVYYCGRWGQGTLV",
"frl": "SLPVSLGDQASISCWYLQKLLIYGVPDRFSGSGFTLEISRVEAEDLGVYYCFGGGTKL",
"heavy": "QVQLEQPGAELVKPGASVKLSCKASGYTFTSNWINWVKQRPGQGLEWIGHISPGSSSTNYNEKFKSKATLTVDTSSSTAYMQLSSLTSDDSAVYYCGREETVRASFGNWGQG",
"light": "ELVMTQTPPSLPVSLGDQASISCRSSQSIVHSNGDTYLEWYLQKPGQSPKLLIYKVSNRFSGVPDRFSGSGSGTDFTLEISRVEAEDLGVYYCFQGSHVPRTFGGG"
}
}
)_qwe_";
*/



class Regex_based_cdr_detection_tests : public CxxTest::TestSuite
{
public:
	Regex_based_cdr_detection_tests() {}

	// Shared initialization goes here.
	void setUp() { core_init(); }

	// Shared finalization goes here.
	void tearDown() {}

	// ------------------------------------------ //
	/// @brief test if CDR detection match knownw results from our DB
	void test_cdr_regions_detection() {
#ifdef __ANTIBODY_GRAFTING__

		using namespace protocols::antibody::grafting;
		using Value = utility::json_spirit::mValue;

		Value root;

		//if( utility::json_spirit::read(json_data, root) ) {

		std::ifstream data_stream("protocols/antibody/grafting/cdr-test-data.json");

		if( utility::json_spirit::read(data_stream, root) ) {

			for(auto &target : root.get_obj()) {
				std::cout << "Traget:" << target.first << std::endl;
				Value data { target.second };

				string heavy_chain_sequence { data.get_obj()["heavy_chain_sequence"].get_str() };
				string light_chain_sequence { data.get_obj()["light_chain_sequence"].get_str() };

				AntibodySequence as(heavy_chain_sequence, light_chain_sequence);

				RegEx_based_CDR_Detector().detect(as);

				TS_ASSERT_EQUALS( as.h1_sequence(), data.get_obj()["h1"].get_str() );
				TS_ASSERT_EQUALS( as.h2_sequence(), data.get_obj()["h2"].get_str() );
				TS_ASSERT_EQUALS( as.h3_sequence(), data.get_obj()["h3"].get_str() );

				TS_ASSERT_EQUALS( as.l1_sequence(), data.get_obj()["l1"].get_str() );
				TS_ASSERT_EQUALS( as.l2_sequence(), data.get_obj()["l2"].get_str() );
				TS_ASSERT_EQUALS( as.l3_sequence(), data.get_obj()["l3"].get_str() );

				FRH_FRL fr = calculate_frh_frl(as);

				string frh = fr.frh1 + fr.frh2 + fr.frh3 + fr.frh4;
				string frl = fr.frl1 + fr.frl2 + fr.frl3 + fr.frl4;

				TS_ASSERT_EQUALS( frh, data.get_obj()["frh"].get_str() );
				TS_ASSERT_EQUALS( frl, data.get_obj()["frl"].get_str() );


				// Chothia numberer tests
				AntibodyFramework heavy_fr = as.heavy_framework();
				AntibodyFramework light_fr = as.light_framework();

				trim_framework(as, heavy_fr, light_fr);

				AntibodyNumbering an( Chothia_Numberer().number(as, heavy_fr, light_fr) );

				string trimmed_heavy_sequence = heavy_fr.fr1 + as.h1_sequence() + heavy_fr.fr2 + as.h2_sequence() + heavy_fr.fr3 + as.h3_sequence() + heavy_fr.fr4;
				string trimmed_light_sequence = light_fr.fr1 + as.l1_sequence() + light_fr.fr2 + as.l2_sequence() + light_fr.fr3 + as.l3_sequence() + light_fr.fr4;

				string trimmed_heavy_numbering = utility::join( an.heavy.all(), " ");
				string trimmed_light_numbering = utility::join( an.light.all(), " ");

				TS_ASSERT_EQUALS( trimmed_heavy_sequence, data.get_obj()["trimmed_heavy_sequence"].get_str() );
				TS_ASSERT_EQUALS( trimmed_light_sequence, data.get_obj()["trimmed_light_sequence"].get_str() );

				TS_ASSERT_EQUALS( trimmed_heavy_numbering, data.get_obj()["trimmed_heavy_numbering"].get_str() );
				TS_ASSERT_EQUALS( trimmed_light_numbering, data.get_obj()["trimmed_light_numbering"].get_str() );
			}
		}
		else {
			TS_FAIL("Parsing of JSON data file finished with error!!!");
		}

#endif // __ANTIBODY_GRAFTING__
	}
};
