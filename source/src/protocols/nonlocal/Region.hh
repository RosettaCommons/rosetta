// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/Region.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_REGION_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_REGION_HH

// Unit header
#include <protocols/nonlocal/Region.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Project headers
#include <core/types.hh>

namespace protocols {
namespace nonlocal {

/// @class A continguous sequence of residues
class Region : public utility::pointer::ReferenceCount {
	typedef core::Size Size;

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Region();
	/// @brief Creates a new region with the specified start/stop residues
	Region(Size start_pos, Size stop_pos);

	/// @brief Returns the starting position of this region
	/// O(1)
	Size start() const;

	/// @brief Returns the stopping position of this region
	/// O(1)
	Size stop() const;

	/// @brief Returns the length of this region. Makes no assumption about
	/// directionality. That is, Region(3,5).length() == Region(5,3).length().
	/// O(1)
	Size length() const;

	/// @brief Returns true if start <= stop, false otherwise
	bool increasing() const;

	/// @brief Returns true if stop <= start, false otherwise
	bool decreasing() const;

private:
	/// @brief The starting position of the contiguous sequence
	Size start_;

	/// @brief The ending position of the contiguous sequence
	Size stop_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_REGION_HH_
