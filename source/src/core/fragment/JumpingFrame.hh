// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/Frame.hh
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_JumpingFrame_HH
#define INCLUDED_core_fragment_JumpingFrame_HH

// Unit Headers
#include <core/fragment/JumpingFrame.fwd.hh>

// Package Headers
#include <core/fragment/Frame.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ STL Headers
#include <map>

#include <core/fragment/FragData.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>


namespace core {
namespace fragment {


/// @brief JumpingFrame is a discontinuous frame
/// i.e, the SRFDs stored in the FragData objects can be applied to residues anywhere
/// a 5 7 9 Frame of a FragData containing three BBTorsionSRFDs would change torsions of 5 7 9
/// a 5 32 2 Frame of a FragData containing two BBTorsionSRFD and on JumpSRFD would change torsions of 5 and 32 and the RT of jump_nr 2
/// note that in the latter case the 2 is not coding for a residue number!
///
/// what is meaning of start() and end() ? If some of the positions are actually jump_nr should we ignore them for the
/// start() / end() question ? or should we store which positions are jumps and translate into residue numbers from jump_nr ?
/// could a MixedFrame contain Frames MixedFrame -- > non-cont-Frame, JumpFrame
class NonContinuousFrame : public Frame {
  typedef utility::vector1< Size > PosList;
	typedef Frame Parent;
public:

  // c'stor
  NonContinuousFrame( Size start, Size end, Size length )
    : Frame( start, end, length ), pos_( length )
  {};

	/// @brief clone method, new frame with same alignment position, fragments are not copied!
 	virtual FrameOP clone() const;

//  -- cloneing with fragments is taken care of by base-class...

// 	/// @brief clone method, new frame with same alignment position, fragments are not copied!
// 	virtual NonContinuousFrameOP clone_with_frags() const;

// 	/// @brief clone method, new frame with same alignment position, one fragments is copied as template ( valid() == false )
// 	virtual NonContinuousFrameOP clone_with_template();


  /// @brief translate intra-frame position into sequence position. (trivial for base-class)
  virtual core::Size seqpos( core::Size intra_pos ) const { // BaseClass --> continuous frames
    runtime_assert( intra_pos <= length() );
    return pos_[ intra_pos ];
  }

	virtual bool moves_residue( core::Size pos ) const {
		PosList::const_iterator it = find( pos_.begin(), pos_.end(), pos );
		return it != pos_.end();
	}

  /// @brief true if frame is continuous
  virtual bool is_continuous() const
  { return false; };

  /// @brief assign sequence position or jump_nr to internal position pos
  void set_pos( Size intra_pos, Size setting ) {
    assert( intra_pos <= length() );
    pos_[ intra_pos ] = setting;
  }

	/// @brief shift to new start position ( change end accordingly )
	virtual void shift_to( core::Size setting );

	/// @brief shift frame by offset relative to current start position ( change end accordingly )
	virtual void shift_by( int offset );

	virtual void show( std::ostream& ) const;
	virtual void read( std::istream& );

	virtual bool align( core::id::SequenceMapping const& map );

protected:
	void show_pos( std::ostream &out ) const;
	PosList& pos() { return pos_ ; } ;
	PosList const& pos() const { return pos_; };
private:

  /// @brief stores the residue number's or jump_nr's associated with the SRFDs in FragData
  PosList pos_;
};


/// @brief JumpingFrame, so far there is nothing special about JumpingFrames.
/// but I think we might want to have additionally information like the start and end residues that belong to a certain jump_nr.!
/// okay: right now I require that the creator of a JumpingFrame sets start to the start-residue of the jump
class JumpingFrame : public NonContinuousFrame {
public:

	JumpingFrame()
		: NonContinuousFrame( 0.0, 0.0, 0.0 ) {};

  JumpingFrame( Size start, Size end, Size length )
		: NonContinuousFrame( start, end, length ) {};

	///@brief convience --- generate a matching FragData object dofs are unitialized!
	//FragDataOP generate_fragdata( SingleResidueFragDataOP frag_res_type, SingleResidueFragDataOP jump_frag_type  )

  /// @brief clone method, new frame with same alignment position, fragments are not copied!
  virtual FrameOP clone() const {
		JumpingFrameOP newFrame( new JumpingFrame( start(), end(), length() ) );
		newFrame->pos() = pos();
		return newFrame;
  }

	static std::string _static_type_name() {
		return "JUMPFRAME";
	}

	virtual std::string type() const {
		return _static_type_name();
	}
  /// fragment_as_pose ????




}; //JumpingFrame


}
}

#endif
