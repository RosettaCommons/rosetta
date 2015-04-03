// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/HamiltonianExchangeCreator.hh
/// @brief This class will create instances of Mover HamiltonianExchange for the MoverFactory
/// @author Oliver Lange


#ifndef INCLUDED_protocols_canonical_sampling_HamiltonianExchangeCreator_hh
#define INCLUDED_protocols_canonical_sampling_HamiltonianExchangeCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace canonical_sampling {

/// @brief RosettaScripts factory for HamiltonianExchange.
class HamiltonianExchangeCreator : public protocols::moves::MoverCreator {
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
};

}
}

#endif

