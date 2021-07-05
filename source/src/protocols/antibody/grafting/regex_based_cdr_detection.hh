// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/grafting/regex_based_cdr_detection.hh
/// @brief RegEx based detection of CRS's
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)
/// @author Jeliazko Jeliazkov

#ifndef INCLUDED_protocols_antibody_grafting_regex_based_cdr_detection_HH
#define INCLUDED_protocols_antibody_grafting_regex_based_cdr_detection_HH

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/antibody_sequence.fwd.hh>
#include <protocols/antibody/grafting/cdr_detection.hh>

#include <basic/report.fwd.hh>


namespace protocols {
namespace antibody {
namespace grafting {


/// @brief Use RegEx and antibody sequence information to detect CDR's
class RegEx_based_CDR_Detector : public protocols::antibody::grafting::CDR_Detector {
public:

	RegEx_based_CDR_Detector();

	RegEx_based_CDR_Detector( basic::ReportOP );
	~RegEx_based_CDR_Detector() override;

	void detect_heavy_chain(AntibodySequence &) override;
	void detect_light_chain(AntibodySequence &) override;
};

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_regex_based_cdr_detection_HH
