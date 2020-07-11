// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constel/ResidueMask.hh
/// @brief Class implementing a boolean mask over the residues of a pose.
/// @author Andrea Bazzoli

#ifndef INCLUDED_ResidueMask_hh
#define INCLUDED_ResidueMask_hh

#include <core/pose/Pose.hh>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <core/types.hh>
#include <string>
#include <iosfwd>

namespace protocols {
namespace constel {

class ResidueMask : public utility::VirtualBase {

	utility::vector1<bool> mask;

public:
	ResidueMask(core::pose::Pose& ps, std::string const& fname);
	bool operator[](core::Size const i);
	void print(std::ostream& os) const;
};

typedef utility::pointer::shared_ptr< ResidueMask > ResidueMaskOP;

} // constel
} // protocols

#endif
