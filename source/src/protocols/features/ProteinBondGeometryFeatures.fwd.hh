// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProteinBondGeometryFeatures.fwd.hh
/// @brief  report protein bond geometry features
/// @author Patrick Conway

#ifndef INCLUDED_protocols_features_ProteinBondGeometryFeatures_fwd_hh
#define INCLUDED_protocols_features_ProteinBondGeometryFeatures_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace features {

class ProteinBondGeometryFeatures;
typedef utility::pointer::shared_ptr< ProteinBondGeometryFeatures > ProteinBondGeometryFeaturesOP;
typedef utility::pointer::shared_ptr< ProteinBondGeometryFeatures const > ProteinBondGeometryFeaturesCOP;

//class TorsionDatabase;
class BondAngleDatabase;
//class BondLengthDatabase;

//typedef  utility::pointer::shared_ptr< TorsionDatabase > TorsionDatabaseOP;
typedef  utility::pointer::shared_ptr< BondAngleDatabase > BondAngleDatabaseOP;
//typedef  utility::pointer::shared_ptr< BondLengthDatabase > BondLengthDatabaseOP;


} // namespace
} // namespace


#endif //include guard
