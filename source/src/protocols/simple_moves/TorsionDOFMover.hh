// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/simple_moves/TorsionDOFMover.hh
/// @brief TorsionDOFMover header
/// @author Steven Lewis

#ifndef INCLUDED_protocols_simple_moves_TorsionDOFMover_hh
#define INCLUDED_protocols_simple_moves_TorsionDOFMover_hh

// Unit Headers
#include <protocols/simple_moves/TorsionDOFMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <core/id/AtomID.hh>

//#include <core/scoring/methods/MMTorsionEnergy.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

///TODO: De-duplicate shared code from RotateJumpAxisMover (angle picking code)

/// @details This mover rotates a specific AtomTree Torsion degree of freedom (any valid 4-body torsion).  It can rotate by a fixed amount, within a range, or randomly.  Optionally, the mover will attempt to internally score the move with MMTorsionEnergy (similar to check_rama in Small/ShearMover).  The mover will print a warning message at apply time if the specified DOF is bad.  For now this mover only allows one DOF; if you want to have it consider multiple DOF's that might be a good idea.  The DOF is determined by a set of 4 atoms; this allows the mover to check the validity of the DOF.  I found it conceptually simpler to think about the 4 atoms involved in the torsion than try to trace DOF_IDs.
class TorsionDOFMover : public protocols::moves::Mover {

public:
	/// @brief default ctor
	TorsionDOFMover();

	/// @brief constructor for random distribution (just needs torsion)
	TorsionDOFMover( core::id::AtomID const & atom1, core::id::AtomID const & atom2, core::id::AtomID const & atom3, core::id::AtomID const & atom4 );

	/// @brief constructor for range - these angles are in degrees, not radians!
	TorsionDOFMover( core::id::AtomID const & atom1, core::id::AtomID const & atom2, core::id::AtomID const & atom3, core::id::AtomID const & atom4, core::Angle const upper, core::Angle const lower );

	/// @brief constructor for single value - these angles are in degrees, not radians!
	TorsionDOFMover( core::id::AtomID const & atom1, core::id::AtomID const & atom2, core::id::AtomID const & atom3, core::id::AtomID const & atom4, core::Angle const angle );

	virtual ~TorsionDOFMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	//////////////////////////////////getters, setters////////////////////////////////////////
	/// @brief set range of desired change - on [180, -180) degrees
	void set_angle_range( core::Angle const upper, core::Angle const lower )
	{ upper_angle_ = upper; lower_angle_ = lower;}

	/// @brief return range of allowed angles
	void get_angle_range( core::Angle & upper, core::Angle & lower ) const
	{ upper = upper_angle_; lower = lower_angle_; }

	/// @brief change the torsion DOF under consideration
	void set_DOF( core::id::AtomID const & atom1, core::id::AtomID const & atom2, core::id::AtomID const & atom3, core::id::AtomID const & atom4 )
	{
		atom1_ = atom1;
		atom2_ = atom2;
		atom3_ = atom3;
		atom4_ = atom4;
	}

	/// @brief return DOF
	void get_DOF( core::id::AtomID & atom1, core::id::AtomID & atom2, core::id::AtomID & atom3, core::id::AtomID & atom4 ) const
	{
		atom1 = atom1_;
		atom2 = atom2_;
		atom3 = atom3_;
		atom4 = atom4_;
	}

	/// @brief (de)activate scoring check
	void check_mmt(bool const setting) { check_MMT_ = setting; }

	/// @brief getter for scoring check
	bool check_mmt() const { return check_MMT_; }

	/// @brief set temperature for scoring check
	void temp(core::Real const setting) { temp_ = setting; }

	/// @brief getter for temperature for scoring check
	core::Real temp() const { return temp_; }

	/// @brief set number of tries
	void tries(core::Size const setting) { tries_ = setting; }

	/// @brief getter for number of tries
	core::Size tries() const { return tries_; }

private:
	/// @brief calculate angle for perturbation - call to RNG
	core::Angle calc_angle();

	/// @brief calculate mmt score for the moving bond
	core::Energy score_torsion(core::pose::Pose & pose);

	/// @brief boltzmann calculation - is the new score acceptable?
	bool boltzmann( core::Energy const pre_score, core::Energy const post_score );

	//data

	/// @brief these atoms define the torsion
	core::id::AtomID atom1_, atom2_, atom3_, atom4_;

	/// @brief these angles define the range of transformations
	core::Angle upper_angle_, lower_angle_; //these angles are in degrees, not radians!

	/// @brief boolean - should we restrict moves based on MMTorsionEnergy?
	bool check_MMT_;

	/// @brief MMTorsionEnergy scorefunction
	core::scoring::ScoreFunctionOP mmt_;

	/// @brief temperature for accepting moves
	core::Real temp_;

	/// @brief number of attempts at finding a valid move
	core::Size tries_;

};//end TorsionDOFMover

}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_TorsionDOFMover_HH
