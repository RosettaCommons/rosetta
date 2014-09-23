// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinRMSDFeatures.fwd.hh
/// @brief  report RMSD similarity of a protein structure against supplied reference structure
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ProteinRMSDFeatures_fwd_hh
#define INCLUDED_protocols_features_ProteinRMSDFeatures_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace features{

class ProteinRMSDFeatures;
typedef utility::pointer::shared_ptr< ProteinRMSDFeatures > ProteinRMSDFeaturesOP;
typedef utility::pointer::shared_ptr< ProteinRMSDFeatures const > ProteinRMSDFeaturesCOP;

}//features
}//protocols

#endif // include guard
