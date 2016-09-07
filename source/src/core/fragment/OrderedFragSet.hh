// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/ConstantLengthFragSet.hh
/// @brief  yields a simple implementation of a fragset
/// @author Oliver Lange ( olange@u.washington.edu)
/// @date   Wed Aug 22 12:08:31 2007


#ifndef INCLUDED_core_fragment_OrderedFragSet_HH
#define INCLUDED_core_fragment_OrderedFragSet_HH

// Unit Headers
#include <core/fragment/OrderedFragSet.fwd.hh>

// Package Headers
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.fwd.hh>
#ifdef WIN32
#include <core/fragment/Frame.hh> // WIN32 INCLUDE
#endif

// Project Headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// std Headers
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace fragment {

/* classic 9mer Frags would be in one of those */
class OrderedFragSet : public FragSet {
public:
	typedef FragSet Parent;
	typedef std::map< Size, FrameList > FrameMap;

	// ConstantLengthFragSet is a FragSet with only one frame per position!
public:

	OrderedFragSet();
	~OrderedFragSet() override;

	FragSetOP clone() const override;
	FragSetOP empty_clone() const override;

	Size region(
		kinematics::MoveMap const& mm,
		core::Size start,
		core::Size end, //not used
		core::Size , //min_overlap not used
		core::Size , //min_length not used
		FrameList &frames
	) const override;

	/// @brief Accessor for the Frame at the specified insertion position. Returns false if
	/// there is no frame at the specified position.
	Size frames( Size pos, FrameList &frames ) const override;

	ConstFrameIterator begin() const override;
	ConstFrameIterator end() const override;

	FrameIterator nonconst_begin() override;
	FrameIterator nonconst_end() override;

	bool empty() const override;

protected:

	void add_( FrameOP aframe ) override;


private:
	FrameMap frames_;

};

}
}

#endif
