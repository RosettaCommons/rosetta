// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/sampling/MonteCarloAssemblyMoverCreator.hh
///
/// @brief
/// @author Tim Jacobs



#ifndef INCLUDED_protocols_sewing_sampling_MonteCarloAssemblyMoverCreator_hh
#define INCLUDED_protocols_sewing_sampling_MonteCarloAssemblyMoverCreator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace sewing  {

class MonteCarloAssemblyMoverCreator : public protocols::moves::MoverCreator
{
public:
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static  std::string mover_name();
};

}
}

#endif

