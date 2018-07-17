// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/cdr_detection.hh
/// @brief base class detection of CDRS's
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)
/// @author Jeliazko Jeliazkov
/// @author Indigo King (indigo.c.king@gmail.com)

/// @brief base class for CDR detection

#include <protocols/antibody/grafting/cdr_detection.hh>
#include <basic/report.hh>
#include <protocols/antibody/grafting/antibody_sequence.hh>
#include <string>
#include <utility/json_spirit/json_spirit.h>

#ifdef __ANTIBODY_GRAFTING__

namespace protocols {
namespace antibody {
namespace grafting {

CDR_Detector::CDR_Detector() : basic::Reporter() {}
CDR_Detector::CDR_Detector( basic::ReportOP ) : basic::Reporter() {}
CDR_Detector::~CDR_Detector()= default;

/// @brief Base class for antibody CDR detector. Sub-class it to implement particular input file/cmdline or prediction
/// @brief This function build a json spirit report while it calls the methods to detect the heavy+light chain CDR cutpoints, detect methods are virtual and defined in derived classes
void CDR_Detector::detect(AntibodySequence &A){
	*this << "CDR_Detector run with arguments:\n  heavy: " + A.heavy.sequence + "\n  light: " + A.light.sequence + "\n\n";
	set("heavy", A.heavy.sequence);  set("light", A.light.sequence);

	detect_heavy_chain(A);
	detect_light_chain(A);

	*this << "CDR_Detector results:\n";

	struct {
		std::string name;
		CDR_Bounds &cdr;
		std::string sequence;
	} J[] {
	{"h1", A.heavy.cdr1, A.h1_sequence()}, {"h2", A.heavy.cdr2, A.h2_sequence()}, {"h3", A.heavy.cdr3, A.h3_sequence()},
	{"l1", A.light.cdr1, A.l1_sequence()}, {"l2", A.light.cdr2, A.l2_sequence()}, {"l3", A.light.cdr3, A.l3_sequence()},
		};

	for ( auto &j : J ) {
		std::string line = j.name + ": " + j.sequence;
		line += line.size() < 8 ? "\t\t\t" : (line.size() < 16 ? "\t\t" : "\t");
		*this << "\t" + line << j.cdr << "\n";
		utility::json_spirit::Object cdr;

		cdr.push_back( utility::json_spirit::Pair("sequence", j.sequence) );
		cdr.push_back( utility::json_spirit::Pair("begin",    int(j.cdr.begin) ) );
		cdr.push_back( utility::json_spirit::Pair("end",      int(j.cdr.end)) );
		set(j.name, cdr);
	}
	*this << "\n";
}


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

