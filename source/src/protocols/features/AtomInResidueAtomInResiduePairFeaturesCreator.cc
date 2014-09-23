// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/feature/AtomInResidueAtomInResiduePairFeaturesCreator.hh
/// @brief  Header for AtomInResidueAtomInResiduePairFeaturesCreator for the AtomInResidueAtomInResiduePairFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/AtomInResidueAtomInResiduePairFeaturesCreator.hh>

// Package Headers
#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/AtomInResidueAtomInResiduePairFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

AtomInResidueAtomInResiduePairFeaturesCreator::AtomInResidueAtomInResiduePairFeaturesCreator() {}
AtomInResidueAtomInResiduePairFeaturesCreator::~AtomInResidueAtomInResiduePairFeaturesCreator() {}
FeaturesReporterOP AtomInResidueAtomInResiduePairFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new AtomInResidueAtomInResiduePairFeatures );
}

std::string AtomInResidueAtomInResiduePairFeaturesCreator::type_name() const {
	return "AtomInResidueAtomInResiduePairFeatures";
}

} //namespace features
} //namespace protocols
