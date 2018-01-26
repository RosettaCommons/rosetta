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
/// @author Roland A. Pache, PhD


#ifndef INCLUDED_core_fragment_ConstantLengthFragSet_HH
#define INCLUDED_core_fragment_ConstantLengthFragSet_HH

// Unit Headers
#include <core/fragment/ConstantLengthFragSet.fwd.hh>

// Package Headers
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>

// Project Headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.fwd.hh>


// Package Headers

// std Headers

/* Just a mad thought: with fragments becoming ever more "Residue" like one might want to use the
packer to choose a combination of good fragments instead of makeing independent choices.
I guess, it is only a question of keeping the combinatorics in control...
maybe it makes sense to pack with only "unconfident" regions of the backbone flexible ..
*/

namespace core {
namespace fragment {

/* classic 9mer Frags would be in one of those */
/// @brief A set of fragments that contains fragments of a constant length read
/// in from a fragment file.
/// @note this object is a simple implementation of the FragSet
/// @note for custom fragments, check out https://robetta.bakerlab.org/fragmentsubmit.jsp
///
/// example:
///     ninemers = ConstantLengthFragSet(9)
///     ninemers.read_fragment_file("test9_fragments")
/// See also:
///     Pose
///     ClassicFragmentMover
class ConstantLengthFragSet : public FragSet {
	typedef FragSet Parent;
	// ConstantLengthFragSet is a FragSet with only one frame per position!
public:

	ConstantLengthFragSet( Size frag_length ) {
		set_max_frag_length( frag_length );
	};

	ConstantLengthFragSet( Size frag_length, std::string filename );

	ConstantLengthFragSet() {};

	~ConstantLengthFragSet() override;

	ConstantLengthFragSet( FragSet const& fragments );

	FragSetOP clone() const override
	{
		return FragSetOP( new ConstantLengthFragSet( *this ) );
	}

	FragSetOP empty_clone() const override
	{
		return FragSetOP( new ConstantLengthFragSet() );
	}

	/// @brief Loads the contents of  <filename>  into the ConstantLengthFragSet
	///
	/// example:
	///     ninemers.read_fragment_file("test9_fragments")
	/// See also:
	///     ConstantLengthFragSet
	///     Pose
	///     Pose.omega
	///     Pose.phi
	///     Pose.psi
	///     ClassicFragmentMover
	void read_fragment_file( std::string filename, Size top25 = 0, Size ncopies = 1, bool bAnnotation = false  );

	void read_fragment_stream( std::string const & filename, std::string const & first_line, std::istream & data, Size top25 = 0, Size ncopies = 1, bool bAnnotation = false  );

	// void print_fragments();
	/// @brief there is only one Frame per position, end / max_overlap are ignored
	Size region(
		kinematics::MoveMap const& mm,
		core::Size start,
		core::Size end,
		core::Size , //min_overlap not used
		core::Size , //min_length not used
		FrameList &frames
	) const override;

	/// @brief returns the number and list of all fragment alignment frames that somehow overlap with the given region
	///(also allows those frames that start before the region and reach into it)
	Size overlapping_with_region(
		kinematics::MoveMap const& mm,
		core::Size start,
		core::Size end,
		core::Size min_overlap,
		core::Size min_length,
		FrameList &frames
	) const override;


	//  /// @brief Accessor for the Frames at the specified insertion position. Returns 0=false if
	//  /// there is no frame at the specified position.
	//  virtual Size frames(
	//   Size pos,
	//   FrameList & frame
	//  );

	/// @brief iterate over contents
	ConstFrameIterator begin() const override;
	ConstFrameIterator end() const override;


	FrameIterator nonconst_begin() override;
	FrameIterator nonconst_end() override;

	bool empty() const override {
		return frames_.size()==0;
	}
protected:
	void add_( FrameOP aframe ) override;

private:
	FrameList frames_;
};

}
}

#endif
