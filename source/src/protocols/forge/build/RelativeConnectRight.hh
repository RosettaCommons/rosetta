// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/RelativeConnectRight.hh
/// @brief version of ConnectRight instruction that depends upon results from
///  another BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_RelativeConnectRight_hh
#define INCLUDED_protocols_forge_build_RelativeConnectRight_hh

// unit headers
#include <protocols/forge/build/RelativeConnectRight.fwd.hh>

// package headers
#include <protocols/forge/build/ConnectRight.hh>
#include <protocols/forge/build/RelativeSequencePosition.fwd.hh>

#include <protocols/forge/build/BuildInstruction.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief version of ConnectRight instruction that depends upon results from
///  another BuildInstruction
class RelativeConnectRight : public ConnectRight {


private: // typedefs


	typedef ConnectRight Super;


public: // typedefs


	typedef Super::Size Size;

	typedef Super::ResidueTypeSetCAP ResidueTypeSetCAP;
	typedef Super::LengthEvent LengthEvent;
	typedef Super::MoveMap MoveMap;
	typedef Super::Pose Pose;

	typedef Super::Positions Positions;
	typedef Super::String String;

	typedef core::kinematics::Jump Jump;


public: // construct/destruct


	/// @brief default constructor
	RelativeConnectRight();


	/// @brief RelativeSequencePosition + position on-right jump constructor
	/// @param[in] rp RelativeSequencePosition defining the type of computation to perform.
	///  (will be cloned)
	/// @param[in] right_position connect at this position on 'pose_right'
	/// @param[in] pose_right connect this pose to the right of pose_left when
	///  modify( pose_left ) is called
	RelativeConnectRight(
		RelativeSequencePositionOP const & rp,
		Size const right_position,
		Pose const & pose_right
	);


	/// @brief copy constructor
	RelativeConnectRight( RelativeConnectRight const & rval );


	/// @brief default destructor
	virtual
	~RelativeConnectRight();


public: // assignment


	/// @brief copy assignment
	RelativeConnectRight & operator =( RelativeConnectRight const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	BuildInstructionOP clone() const;


public: // virtual Pose modification methods


	/// @brief do the actual work of modifying the Pose
	virtual
	void modify_impl( Pose & pose_left );


public: // instruction comparison


	/// @brief return set of any fixed positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @return always empty set, no fixed positions
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.  There is currently no way to
	///  represent the dependent fixed position, so we're forced to return an empty
	///  set.
	virtual
	Positions original_fixed_positions() const;


private: // data


	/// @brief function object used to compute the 'left_position' in
	///  ConnectRight
	RelativeSequencePositionOP rp_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_RelativeConnectRight_HH */
