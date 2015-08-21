// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragSet.hh
/// @brief  set of fragments
/// @author Oliver Lange ( olange@u.washington.edu)
/// @date   Wed Aug 22 12:08:31 2007


#ifndef INCLUDED_core_fragment_FragSetCollection_HH
#define INCLUDED_core_fragment_FragSetCollection_HH

// Unit Headers
#include <core/fragment/FragSetCollection.fwd.hh>

// Package Headers
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// Utility headers

#include <utility/vector1.hh>


// Package Headers

/* Just a mad thought: with fragments becoming ever more "Residue" like one might want to use the
packer to choose a combination of good fragments instead of makeing independent choices.
I guess, it is only a question of keeping the combinatorics in control...
maybe it makes sense to pack with only "unconfident" regions of the backbone flexible ..
*/

namespace core {
namespace fragment {

//* not too happy with this class... might be phased out.. */
// this interface extension yields an inexpensive way to pass several FragSets into a mover
class FragSetCollection : public FragSet {
public:
	typedef FragSet Parent;

public:
	FragSetCollection();
	~FragSetCollection();
	FragSetCollection( FragSetCollection const & );

	virtual FragSetOP clone() const;
	virtual FragSetOP empty_clone() const;

	void add_fragset( FragSetOP fragset );

	virtual Size region(
		kinematics::MoveMap const& move_map,
		Size start,
		Size end,
		Size min_overlap,
		Size min_length,
		FrameList &frames
	) const;

	virtual ConstFrameIterator begin() const;
	virtual ConstFrameIterator end() const;

	virtual FrameIterator nonconst_begin();
	virtual FrameIterator nonconst_end();

	virtual bool empty() const;


protected:
	virtual void add_( FrameOP );

private:
	typedef utility::vector1< FragSetOP > FragSetList;
	FragSetList fragset_list_;
};

}
}

#endif
