// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/grafting/json_based_cdr_detection.hh
/// @brief JSON based detection of CDRs
/// @author Indigo King (indigo.c.king@gmail.com)
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)
/// @author Jeliazko Jeliazkov

#ifndef INCLUDED_protocols_antibody_grafting_json_based_cdr_detection_HH
#define INCLUDED_protocols_antibody_grafting_json_based_cdr_detection_HH

//These provide the tokens for the ifdefs and must go first
#include <protocols/antibody/grafting/util.hh>
#include <utility/json_utilities.hh>

#ifdef __ANTIBODY_GRAFTING__
#ifdef _NLOHMANN_JSON_ENABLED_
#include <protocols/antibody/grafting/cdr_detection.hh>

#include <protocols/antibody/grafting/antibody_sequence.hh>

#include <basic/report.hh>
#include <json.hpp>

/// @brief new class to set the CDR sequence position information in AntibodySequence.sequence.cdr.* based on input json numbers. the input json is parsed using the json parsing tools in utility/json_utitlties.hh. JSON format is checked every time and exit with informative message is called if input JSON is misformatted.
// note that this code is intentionally templated on the pattern in regex_based_cdr_detection class so that later development can call components of the regex predictor if only partial CDR information is supplied by the user. however, that addition will require a refactoring of the existing regex predictor code

namespace protocols {
namespace antibody {
namespace grafting {


/// @brief Use input JSON (in same format as output report) to detect the CDR begin/end points
class Json_based_CDR_Detector : public protocols::antibody::grafting::CDR_Detector {
public:

	Json_based_CDR_Detector();
	Json_based_CDR_Detector( basic::ReportOP );
	Json_based_CDR_Detector( basic::ReportOP, std::string const json_filename );
	~Json_based_CDR_Detector() override;

	void detect_heavy_chain(AntibodySequence & ) override;
	void detect_light_chain(AntibodySequence & ) override;

private:
	nlohmann::json json_in_;
};

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // _NLOHMANN_JSON_ENABLED_
#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_json_based_cdr_detection_HH
