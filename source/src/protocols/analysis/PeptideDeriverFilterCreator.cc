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


#include <protocols/analysis/PeptideDeriverFilterCreator.hh>
#include <protocols/analysis/PeptideDeriverFilter.hh>

namespace protocols {
namespace analysis {

protocols::filters::FilterOP
PeptideDeriverFilterCreator::create_filter() const {
	protocols::filters::FilterOP filter( new protocols::analysis::PeptideDeriverFilter() );
	return filter;
}

std::string
PeptideDeriverFilterCreator::keyname() const {
	return "PeptideDeriver";
}

} //namespace analysis
} //namespace protocols
