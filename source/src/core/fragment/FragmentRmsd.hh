// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/fragment/FragmentRmsd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_CORE_FRAGMENT_FRAGMENTRMSD_HH
#define INCLUDED_CORE_FRAGMENT_FRAGMENTRMSD_HH

// Unit header
#include <core/fragment/FragmentRmsd.fwd.hh>

// External headers
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace core {
namespace fragment {

class FragmentRmsd : public utility::pointer::ReferenceCount {
  typedef boost::unordered_map<core::Size, FrameCOP > FrameMap;

 public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FragmentRmsd();
  FragmentRmsd(FragSetCOP fragments);

  /// @brief Returns the kth fragment at the specified position
  /// in the fragment library.
  FragDataCOP fragment(core::Size position, core::Size k) const;

  /// @brief Returns the RMSD of the kth fragment at the specified position
  /// in the fragment library and pose.
  core::Real rmsd(core::Size position, core::Size k, const core::pose::Pose& reference) const;


 protected:
  /// @brief Returns the position'th frame in the fragment library
  FrameCOP frame(core::Size position) const;


 private:
  /// @brief Input fragment library
  FragSetCOP fragments_;

  /// @brief Position-indexable collection of fragments
  mutable FrameMap frames_;
};

}  // namespace fragment
}  // namespace core

#endif  // CORE_FRAGMENT_FRAGMENT_RMSD_HH_
