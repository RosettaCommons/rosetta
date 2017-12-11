// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/util.cc
/// @brief Utility functions for JD3.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd3/util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.jd3.util" );


namespace protocols {
namespace jd3 {


void
print_job_template(){



	TR << "The following is an empty (template) Job file:\n"
		<< "\n"
		<< "<JobDefinitionFile>\n"
		<< "\t<Job>\n"
		<< "\t\t<Input>\n"
		<< "\t\t</Input>\n"
		<< "\t\t<Output>\n"
		<< "\t\t</Output>\n"
		<< "\t\t<Options>\n"
		<< "\t\t</Options>\n"
		<< "\t</Job>\n"
		<< "</JobDefinitionFile>\n\n";

	TR << "Common Input Tags: \n "<<
		"<PDB filename=\"my_fname.pdb\"/>\n"<< std::endl;

	TR << "Common Ouptut Tags: \n" <<
		"<PDB filename_pattern=\"Some_job-specific_name_$\"/>\n" << std::endl;

	TR << std::endl;
	TR << "The rosetta_scripts application will now exit." << std::endl;
	TR.flush();

}



} //protocols
} //jd3


