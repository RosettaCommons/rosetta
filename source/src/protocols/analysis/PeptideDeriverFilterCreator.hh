// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  Creator class for PeptideDeriverFilterCreator

/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @date   Sep. 21, 2014


// PeptideDeriverFilterCreator


#ifndef INCLUDED_protocols_analysis_PeptideDeriverFilterCreator_hh
#define INCLUDED_protocols_analysis_PeptideDeriverFilterCreator_hh

// Package Headers
#include <protocols/filters/FilterCreator.hh>

namespace protocols {
namespace analysis {

class PeptideDeriverFilterCreator : public protocols::filters::FilterCreator {
public:
	// XRW TEMP  protocols::filters::FilterOP create_filter() const override;
	// XRW TEMP  std::string keyname() const override;
	protocols::filters::FilterOP create_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //namespace analysis
} //namespace protocols


#endif
// INCLUDED_protocols_analysis_PeptideDeriverFilterCreator_hh
