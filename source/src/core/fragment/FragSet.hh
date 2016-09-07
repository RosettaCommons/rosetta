// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/FragSet.hh
/// @brief  set of fragments
/// @author Oliver Lange ( olange@u.washington.edu)
/// @date   Wed Aug 22 12:08:31 2007
/// @author Roland A. Pache, PhD


#ifndef INCLUDED_core_fragment_FragSet_HH
#define INCLUDED_core_fragment_FragSet_HH

// Unit Headers
#include <core/fragment/FragSet.fwd.hh>

// Package Headers
#include <core/types.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/Frame.fwd.hh>

// Project headers
#include <core/kinematics/MoveMap.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <utility/vector1.hh>

#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif


// Package Headers

namespace core {
namespace fragment {

/// @detail The FragSet: (Interface Definition -- Virtual BaseClass )
/*
The FragSet is the basic interface to communicate with mover classes. Its main purpose is a fast access
to fragments by seqpos. More functionality might be coming.

To access fragments:

by region: --> region( .. )
all frames: FrameIterator it = fragset.begin() ... fragset.end()
all fragments: FragID_Iterator it = fragset.begin() ... fragset.end()

Two simple extensions are already provided:
ConstantLengthFragSet: fragments are stored by their start-sequence position, only one frame per position
OrderFragSet: fragments are stored by their start-sequence position, multiple frames per position
what might be useful: a map from position to all fragments that overlap with that position

*/


class FragSet : public utility::pointer::ReferenceCount {

	//how to iterate over fragments that fit certain search criteria ?
public:
	FragSet() : min_pos_ ( 100000 ), max_pos_( 0 ), max_frag_length_( 0 ), global_offset_ (0)
	{}

	~FragSet() override = default;

	/// @brief clone and copy the pointer of all Frames (Frames will not be copied)
	virtual FragSetOP clone() const = 0;

	/// @brief create an empty clone
	virtual FragSetOP empty_clone() const = 0;


	/// @brief iterate over contents
	virtual ConstFrameIterator begin() const = 0;
	virtual ConstFrameIterator end() const = 0;

	/// @brief iterate over contents
	virtual FrameIterator nonconst_begin() = 0;
	virtual FrameIterator nonconst_end() = 0;


	/// @brief appends frames at sequence position pos to frames, returns nr of frames added
	virtual Size frames( core::Size pos, FrameList &frames ) const {
		return region_simple( pos, pos, frames);
	}

	/// @brief returns fragments that exactly span seq_pos start...end
	virtual Size region_simple( core::Size start, core::Size end, FrameList &frames ) const {
		// return region( start, end, end-start+1, end-start+1, frames );
		return region_all( start, end, end-start+1, end-start+1, frames );
	}

	/// @brief return a list of frames that all sample the specified region, assume all motions are allowed
	virtual Size region_all(
		core::Size start,
		core::Size end,
		core::Size min_overlap,
		core::Size min_length,
		FrameList &frames
	) const;

	/// @brief the region thing has to be thought-over. How do we really want to sample fragments?
	/// for now, we ignore everything in this call and just return frags that have "start" as there specified start() entry.
	virtual Size region(
		kinematics::MoveMap const& move_map,
		core::Size start,
		core::Size end,
		core::Size min_overlap,
		core::Size min_length,
		FrameList &frames
	) const = 0;

	/// @brief returns the number and list of all fragment alignment frames that somehow overlap with the given region
	///(also allows those frames that start before the region and reach into it)
	virtual Size overlapping_with_region(
		kinematics::MoveMap const& mm,
		core::Size start,
		core::Size end,
		core::Size min_overlap,
		core::Size min_length,
		FrameList &frames
	) const;

	/// @brief InsertMap and InsertSize gives quick overview which residues can be affected by fragments.
	/// insert_map --- list of start-positions, insert_size corresponding list of longest fragment at position x
	virtual void generate_insert_map( kinematics::MoveMap const& mm, InsertMap &insert_map, InsertSize &insert_size) const;

	/// @brief returns the maximal sequence position that can be affected by fragments in this set
	Size max_pos() const
	{ return max_pos_; }

	/// @brief returns the first sequence position that can be affected by fragments in this set
	Size min_pos() const
	{ return min_pos_; }

	/// @brief returns the longest fragment stored in this FragSet.
	Size max_frag_length() const
	{ return max_frag_length_; }

	/// @brief add a single frame. if compatible frame is already in set the frames will be merged
	void add( FrameCOP aFrame );

	/// @brief add all Frames in list
	void add( FrameList const& frames );

	/// @brief add all Frames in FragSet
	void add( FragSet const& frames );

	/// @brief add single fragment
	void add( FragID const& );

	// put all fragments in FragID_list into this FragSet.
	// this function has the following effect:
	//      fragments that belong to the same frame are copied into a new frame
	//      the frame gets added. If all fragments of a frame are in the list, the frame is just added as is
	//
	void add( FragID_List & list ) { //need to work on const-correctness to make this a const-add
		insert_fragID_list( list );
	}

	/// @brief add all fragments in FragID_List
	void insert_fragID_list( FragID_List & list );

	/// @brief returns total size--> counts together all frags in each frame
	Size size() const;

	/// @brief counts number of frames ( slow! - it really counts )
	Size nr_frames() const;

	virtual bool empty() const = 0;

	friend std::ostream & operator<<(std::ostream & out, FragSet const & frags );

	/// @brief shift all frames in FragSet by offset
	virtual void shift_by( int offset );

	int global_offset() const{
		return global_offset_;
	}

	/// @brief resets global_offset of FragSet and shifts FragSet if necessary by calling shift_to
	void global_offset ( int );


	FragSetOP clone_shifted( int ) const;


protected:
	/// @brief setter for max_frag_length_
	void set_max_frag_length( Size setting ) {
		max_frag_length_ = setting;
	}

	void set_max_pos( Size pos ) { max_pos_=pos; }
	void set_min_pos( Size pos ) { min_pos_=pos; }

	/// @brief storage classes have to overload this one to add frames to their container
	virtual void add_( FrameOP aFrame ) = 0;


private:
	Size min_pos_;
	Size max_pos_;
	Size max_frag_length_;

	/// @brief global offset of current Fragmentset. Default = 0
	int global_offset_;

}; // class FragSet


} /* core */
} /* fragment */

#endif
