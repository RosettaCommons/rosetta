// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/json_based_cdr_detection.hh
/// @brief Json based detection of CRS's
/// @author Sergey Lyskov
/// @author Indigo King (indigo.c.king@gmail.com)
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)
/// @author Jeliazko Jeliazkov

/// @brief new class to set the CDR sequence position information in AntibodySequence.sequence.cdr.* based on input json numbers. the input json is parsed using the json parsing tools in utility/json_utitlties.hh. JSON format is checked every time and exit with informative message is called if input JSON is misformatted.
// note that this code is intentionally templated on the pattern in regex_based_cdr_detection class so that later development can call components of the regex predictor if only partial CDR information is supplied by the user. however, that addition will require a refactoring of the existing regex predictor code

//These provide the tokens for the ifdefs and must go first; transcluded from this file's header
//#include <protocols/antibody/grafting/util.hh>
//#include <utility/json_utilities.hh>

#include <protocols/antibody/grafting/json_based_cdr_detection.hh> //provides the ifdef tokens below
#ifdef __ANTIBODY_GRAFTING__
#ifdef _NLOHMANN_JSON_ENABLED_

#include <protocols/antibody/grafting/exception.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
#include <basic/report.hh>

#include <string>
#include <json.hpp>
#include <fstream>
#include <basic/options/option.hh>

namespace protocols {
namespace antibody {
namespace grafting {

using std::string;
using nlohmann::json;

static basic::Tracer TR("protocols.antibody.grafting");

Json_based_CDR_Detector::Json_based_CDR_Detector() : protocols::antibody::grafting::CDR_Detector() {}

// Delegate to the string constructor, assuming the options system has the string filename
Json_based_CDR_Detector::Json_based_CDR_Detector( basic::ReportOP rop ) : Json_based_CDR_Detector::Json_based_CDR_Detector( rop, basic::options::option[ basic::options::OptionKeys::antibody::json_cdr ].value()) {}

Json_based_CDR_Detector::Json_based_CDR_Detector( basic::ReportOP, std::string const json_filename ) : protocols::antibody::grafting::CDR_Detector()
{
	std::ifstream ifs( json_filename );
	json_in_ = nlohmann::json::parse( ifs );
}

Json_based_CDR_Detector::~Json_based_CDR_Detector() = default;

/// @brief get the "h3", "h1", and "h2" json input elements and set the associated AntibodySequence member variables
void Json_based_CDR_Detector::detect_heavy_chain(AntibodySequence &A)
{
	AntibodyChain & H(A.heavy);
	string const & heavy_chain_sequence(H.sequence);

	uint &H1_begin(H.cdr1.begin), &H1_end(H.cdr1.end), &H2_begin(H.cdr2.begin), &H2_end(H.cdr2.end), &H3_begin(H.cdr3.begin), &H3_end(H.cdr3.end);

	bool H1_detected(false), H2_detected(false), H3_detected(false);
	if ( json_in_.count( "h1" ) > 0 ) H1_detected = true;
	if ( json_in_.count( "h2" ) > 0 ) H2_detected = true;
	if ( json_in_.count( "h3" ) > 0 ) H3_detected = true;

	string heavy_first = ( heavy_chain_sequence.size() > 140) ? heavy_chain_sequence.substr(0, 70) : heavy_chain_sequence.substr(0, 60);

	// Identifying H1
	{
		if ( H1_detected ) {

			utility::verify_present_exactly_once_in_json( json_in_, "h1" );
			json H1js;
			utility::extract_nonempty_object_from_json( json_in_, "h1", H1js );
			utility::verify_present_exactly_once_in_json( H1js, "begin" );
			utility::verify_present_exactly_once_in_json( H1js, "end" );
			//might need to send in a platform::Real then static_cast to uint before assigning member indices
			platform::Real H1_begin_in( 0. );
			utility::extract_number_from_json( H1js, "begin", H1_begin_in );
			H1_begin = static_cast< uint >( H1_begin_in );

			platform::Real H1_end_in( 0. );
			utility::extract_number_from_json( H1js, "end", H1_end_in );
			H1_end = static_cast< uint >( H1_end_in );

			debug_assert( H1_end >= H1_begin );
			string H1( heavy_chain_sequence.substr( H1_begin, ( H1_end - H1_begin ) ) );

			TR << "H1 detected: " << H1 << " (" << H1.size() << " residues at positions: " << H1_begin << " to " << H1_end << ")" << std::endl;

		} else {
			//TR << "Trying to predict H1 NOT WORKING" << std::endl;
			//TODO: Call H1 RegEx predictor
			TR << "ERROR: H1 not detected in JSON" << std::endl;
			throw CREATE_EXCEPTION(_AE_cdr_detection_failed_, "All CDRs must be detected if using JSON detectr");
		}
	}

	// Identifying H3
	{
		if ( H3_detected ) {

			utility::verify_present_exactly_once_in_json( json_in_, "h3" );
			json H3js;
			utility::extract_nonempty_object_from_json( json_in_, "h3", H3js );
			utility::verify_present_exactly_once_in_json( H3js, "begin" );
			utility::verify_present_exactly_once_in_json( H3js, "end" );
			platform::Real H3_begin_in( 0. );
			utility::extract_number_from_json( H3js, "begin", H3_begin_in );
			H3_begin = static_cast< uint >( H3_begin_in );

			platform::Real H3_end_in( 0. );
			utility::extract_number_from_json( H3js, "end", H3_end_in );
			H3_end = static_cast< uint >( H3_end_in );

			debug_assert( H3_end >= H3_begin );
			string H3( heavy_chain_sequence.substr( H3_begin, ( H3_end - H3_begin ) ) );

			TR << "H3 detected: " << H3 << " (" << H3.size() << " residues at positions: " << H3_begin << " to " << H3_end << ")" << std::endl;

		} else {
			//TR << "Trying to predict H3 NOT WORKING" << std::endl;
			//TODO: Call H3 RegEx predictor
			TR << "ERROR: H3 not detected in JSON" << std::endl;
			throw CREATE_EXCEPTION(_AE_cdr_detection_failed_, "All CDRs must be detected if using JSON detectr");
		}
	}

	{
		if ( H2_detected ) {
			utility::verify_present_exactly_once_in_json( json_in_, "h2" );
			json H2js;
			utility::extract_nonempty_object_from_json( json_in_, "h2", H2js );
			utility::verify_present_exactly_once_in_json( H2js, "begin" );
			utility::verify_present_exactly_once_in_json( H2js, "end" );
			platform::Real H2_begin_in( 0. );
			utility::extract_number_from_json( H2js, "begin", H2_begin_in );
			H2_begin = static_cast< uint >( H2_begin_in );

			platform::Real H2_end_in( 0. );
			utility::extract_number_from_json( H2js, "end", H2_end_in );
			H2_end = static_cast< uint >( H2_end_in );

		} else {
			H2_begin = H1_end + 15 - 1;
			H2_end   = H3_begin - 33 + 1;
		}

		debug_assert( H2_end >= H2_begin );
		string H2( heavy_chain_sequence.substr(H2_begin, H2_end-H2_begin) );

		TR << "H2 detected: " << H2 << " (" << H2.size() << " residues at positions: " << H2_begin << " to " << H2_end << ")" << std::endl;
	}
}

/// @brief get the "l3", "l1", and "l2" json input elements and set the associated AntibodySequence member variables
void Json_based_CDR_Detector::detect_light_chain(AntibodySequence &A)
{
	AntibodyChain & L(A.light);
	string const & light_chain_sequence(L.sequence);

	uint &L1_begin(L.cdr1.begin), &L1_end(L.cdr1.end), &L2_begin(L.cdr2.begin), &L2_end(L.cdr2.end), &L3_begin(L.cdr3.begin), &L3_end(L.cdr3.end);

	bool L1_detected(false), L2_detected(false), L3_detected(false);
	if ( json_in_.count( "l1" ) > 0 ) L1_detected = true;
	if ( json_in_.count( "l2" ) > 0 ) L2_detected = true;
	if ( json_in_.count( "l3" ) > 0 ) L3_detected = true;
	std::string light_first = ( light_chain_sequence.size() > 130 )  ?  light_chain_sequence.substr( 0, 65 )  :  light_chain_sequence.substr( 0, 60 );

	// Identifying L1
	{
		if ( L1_detected ) {

			utility::verify_present_exactly_once_in_json( json_in_, "l1" );
			json L1js;
			utility::extract_nonempty_object_from_json( json_in_, "l1", L1js );
			utility::verify_present_exactly_once_in_json( L1js, "begin" );
			utility::verify_present_exactly_once_in_json( L1js, "end" );
			platform::Real L1_begin_in( 0. );
			utility::extract_number_from_json( L1js, "begin", L1_begin_in );
			L1_begin = static_cast< uint >( L1_begin_in );

			platform::Real L1_end_in( 0. );
			utility::extract_number_from_json( L1js, "end", L1_end_in );
			L1_end = static_cast< uint >( L1_end_in );

			debug_assert( L1_end >= L1_begin );
			string L1( light_chain_sequence.substr( L1_begin, ( L1_end - L1_begin ) ) );

			TR << "L1 detected: " << L1 << " (" << L1.size() << " residues at positions: " << L1_begin << " to " << L1_end << ")" << std::endl;


		} else {
			//TR << "Trying to predict L1 NOT WORKING" << std::endl;
			//TODO: Call L1 RegEx predictor
			TR << "ERROR: L1 not detected in JSON" << std::endl;
			throw CREATE_EXCEPTION(_AE_cdr_detection_failed_, "All CDRs must be detected if using JSON detectr");
		}
	}

	// Identifying L3
	{
		if ( L3_detected ) {

			utility::verify_present_exactly_once_in_json( json_in_, "l3" );
			json L3js;
			utility::extract_nonempty_object_from_json( json_in_, "l3", L3js );
			utility::verify_present_exactly_once_in_json( L3js, "begin" );
			utility::verify_present_exactly_once_in_json( L3js, "end" );
			platform::Real L3_begin_in( 0. );
			utility::extract_number_from_json( L3js, "begin", L3_begin_in );
			L3_begin = static_cast< uint >( L3_begin_in );

			platform::Real L3_end_in( 0. );
			utility::extract_number_from_json( L3js, "end", L3_end_in );
			L3_end = static_cast< uint >( L3_end_in );

			debug_assert( L3_end >= L3_begin );
			string L3( light_chain_sequence.substr( L3_begin, ( L3_end - L3_begin ) ) );

			TR << "L3 detected: " << L3 << " (" << L3.size() << " residues at positions: " << L3_begin << " to " << L3_end << ")" << std::endl;

		} else {
			//TR << "Trying to predict L3 NOT WORKING" << std::endl;
			//TODO: Call L1 RegEx predictor
			TR << "ERROR: L3 not detected in JSON" << std::endl;
			throw CREATE_EXCEPTION(_AE_cdr_detection_failed_, "All CDRs must be detected if using JSON detectr");
		}
	}

	{
		if ( L2_detected ) {

			utility::verify_present_exactly_once_in_json( json_in_, "l2" );
			json L2js;
			utility::extract_nonempty_object_from_json( json_in_, "l2", L2js );
			utility::verify_present_exactly_once_in_json( L2js, "begin" );
			utility::verify_present_exactly_once_in_json( L2js, "end" );
			platform::Real L2_begin_in( 0. );
			utility::extract_number_from_json( L2js, "begin", L2_begin_in );
			L2_begin = static_cast< uint >( L2_begin_in );

			platform::Real L2_end_in( 0. );
			utility::extract_number_from_json( L2js, "end", L2_end_in );
			L2_end = static_cast< uint >( L2_end_in );

		} else {
			L2_begin = L1_end + 16 - 1;
			L2_end   = L2_begin + 7;
		}

		debug_assert( L2_end >= L2_begin );
		string L2( light_chain_sequence.substr(L2_begin, L2_end-L2_begin) );

		TR << "L2 detected: " << L2 << " (" << L2.size() << " residues at positions: " << L2_begin << " to " << L2_end << ")" << std::endl;

	}
}

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // _NLOHMANN_JSON_ENABLED_
#endif // __ANTIBODY_GRAFTING__
