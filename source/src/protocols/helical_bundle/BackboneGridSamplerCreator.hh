// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle/BackboneGridSamplerCreator.hh
/// @brief This class will create instances of Mover BackboneGridSampler for the MoverFactory.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_BackboneGridSamplerCreator_hh
#define INCLUDED_protocols_helical_bundle_BackboneGridSamplerCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace helical_bundle {

class BackboneGridSamplerCreator : public protocols::moves::MoverCreator {
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
};

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_BackboneGridSamplerCreator_hh

