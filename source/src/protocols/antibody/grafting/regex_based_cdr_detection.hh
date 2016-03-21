// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/grafting/regex_based_cdr_detection.hh
/// @brief RegEx based detection of CRS's
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#ifndef INCLUDED_protocols_antibody_grafting_regex_based_cdr_detection_HH
#define INCLUDED_protocols_antibody_grafting_regex_based_cdr_detection_HH

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/antibody_sequence.hh>

#include <basic/report.hh>


namespace protocols {
namespace antibody {
namespace grafting {


/// @brief Base class for antibody CDR detector. Sub-class it to implement particular detection methods
class CDR_Detector : public basic::Reporter {
public:
	using Reporter::Reporter;

	virtual ~CDR_Detector() {}

	/// @brief Detect CDR's
	/// @throw _AE_cdr_detection_failed_ if for some of the loops detection failed
	virtual void detect(AntibodySequence &) = 0;
};


/// @brief Use RegEx and antibody sequence information to detect CDR's
class RegEx_based_CDR_Detector : public CDR_Detector {
public:
	using CDR_Detector::CDR_Detector;

	/// @brief Detect CDR's
	/// @throw _AE_cdr_detection_failed_ if for some of the loops detection failed
	virtual void detect(AntibodySequence &);

	void detect_heavy_chain(AntibodySequence &);
	void detect_light_chain(AntibodySequence &);
};



} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_regex_based_cdr_detection_HH
