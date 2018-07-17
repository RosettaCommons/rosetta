// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/grafting/cdr_detection.hh
/// @brief base class detection of CDRs
/// @author Indigo King (indigo.c.king@gmail.com)
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)
/// @author Jeliazko Jeliazkov

#ifndef INCLUDED_protocols_antibody_grafting_cdr_detection_HH
#define INCLUDED_protocols_antibody_grafting_cdr_detection_HH

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/antibody_sequence.hh>

#include <basic/report.hh>
#include <string> 
#include <utility/json_spirit/json_spirit.h>

/// @brief base class for CDR detection

namespace protocols {
namespace antibody {
namespace grafting {


/// @brief Base class for antibody CDR detector. Sub-class it to implement particular input file/cmdline or prediction
class CDR_Detector : public basic::Reporter {
public:
	using Reporter::Reporter;

	~CDR_Detector() override;

	CDR_Detector();

	CDR_Detector( basic::ReportOP );

	/// @brief Detect CDR's
	/// @throw _AE_cdr_detection_failed_ if for some of the loops detection is failed
	void detect(AntibodySequence & A);

	virtual void detect_heavy_chain( AntibodySequence & ) = 0;
	virtual void detect_light_chain( AntibodySequence & ) = 0;


};


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_cdr_detection_HH

