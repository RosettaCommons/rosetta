// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/grafting/json_based_cdr_detection.cxxtest.hh
/// @brief  test suite for antibody json_based_cdr_detection code and framework-trimming for Blast+ SCS
/// @author Sergey Lyskov


#include <utility/string_util.hh>
#include <basic/report.hh>

#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>

#include <fstream>

#include <utility/json_utilities.hh>
#include <protocols/antibody/grafting/json_based_cdr_detection.hh>
#include <protocols/antibody/grafting/util.hh>

static basic::Tracer TR("protocols.antibody.grafting.json_based_cdr_detection.cxxtest");

using std::string;


class Json_based_cdr_detection_tests : public CxxTest::TestSuite
{
public:
	Json_based_cdr_detection_tests() {}

	// Shared initialization goes here.
	void setUp() { core_init(); }

	// Shared finalization goes here.
	void tearDown() {}

	// ------------------------------------------ //
	/// @brief test if CDR defined by input json matches antibody sequence output
	void test_cdr_regions_detection() {
#ifdef __ANTIBODY_GRAFTING__
#ifdef _NLOHMANN_JSON_ENABLED_

		using namespace protocols::antibody::grafting;

		/// @brief heavy/light chain sequences for antibody PDBID:1EMT
		string heavy_chain_sequence( "QVHLQESGPELVRPGASVKISCKTSGYVFSSSWMNWVKQRPGQGLKWIGRIYPGNGNTNYNEKFKGKATLTADKSSNTAYMQLSSLTSVDSAVYFCATSSAYWGQGTLLTVSAAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSPRPSETVTCNVAHPASSTKVDKKIVPR" );
		string light_chain_sequence( "DIQMTQTTSSLSASLGDRVTFSCSASQDISNYLNWYQQKPDGTIKLLIYYTSSLRSGVPSRFSGSGSGTDYSLTINNLEPEDIATYFCQQYSRLPFTFGSGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC" );
		/// @brief input json defining CDR cutpoints
		string json_in_fname( "protocols/antibody/grafting/cdr_def_1emt.json" );

		AntibodySequence as(heavy_chain_sequence, light_chain_sequence);
		basic::ReportOP report = utility::pointer::make_shared<basic::Report>("cdr_def_1emt.report");
		Json_based_CDR_Detector( report, json_in_fname ).detect( as );

		TS_ASSERT_EQUALS( as.h1_sequence(), "GYVFSSSWMN" );
		TS_ASSERT_EQUALS( as.h2_sequence(), "RIYPGNGNTNYNEKFKG" );
		TS_ASSERT_EQUALS( as.h3_sequence(), "SSAY" );

		TS_ASSERT_EQUALS( as.l1_sequence(), "SASQDISNYLN" );
		TS_ASSERT_EQUALS( as.l2_sequence(), "YTSSLRS" );
		TS_ASSERT_EQUALS( as.l3_sequence(), "QQYSRLPFT" );

#endif // _NLOHMANN_JSON_ENABLED_
#endif // __ANTIBODY_GRAFTING__
		TS_ASSERT(true);
	}
};
