// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief Class implementing a boolean mask over the residues of a pose.
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#ifndef INCLUDED_ResidueMask_hh
#define INCLUDED_ResidueMask_hh

#include <core/pose/Pose.hh>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <string>
#include <iostream>

using utility::vector1;
using core::pose::Pose;
using core::Size;

namespace devel {
namespace constel {

class ResidueMask : public utility::pointer::ReferenceCount {

	vector1<bool> mask;

	public:
	ResidueMask(Pose& ps, std::string const& fname);
	bool operator[](Size const i);
	void print(std::ostream& os) const;
};

typedef utility::pointer::shared_ptr< ResidueMask > ResidueMaskOP;

} // constel
} // devel

#endif
