// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/kinematic_closure/KinematicPerturber.hh
/// @brief  Header file for KinematicPerturbers used by the kineamtic mover
/// @author Florian Richter, floric@u.washington.edu, march 2009

#ifndef INCLUDED_protocols_moves_kinematic_closure_KinematicPerturber_hh
#define INCLUDED_protocols_moves_kinematic_closure_KinematicPerturber_hh

#include <protocols/moves/KinematicMover.fwd.hh>
#include <protocols/moves/kinematic_closure/KinematicPerturber.fwd.hh>

// Rosetta Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/Ramachandran.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

// Utility Headers
#include <utility/LexicographicalIterator.hh> //not OPable; needs full header
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>

namespace protocols {
namespace moves {
namespace kinematic_closure {

/// @brief pure virtual base class for KinematicPerturber.  KinematicPerturbers determine HOW loops should be perturbed.  The base class contains a provision for determining WHERE they should be perturbed: MoveMap sensitivity.
class KinematicPerturber : public utility::pointer::ReferenceCount {

public:

	KinematicPerturber();
	~KinematicPerturber();

	virtual
	std::string perturber_type() const = 0;

	void
	set_kinmover( KinematicMoverCAP kinmover ){ kinmover_ = kinmover; }

	void
	set_movemap( core::kinematics::MoveMapCOP mm );

	core::kinematics::MoveMapCOP
	get_movemap() const ;

	/// @brief function that perturbs the chain, i.e. sets new values
	/// @brief for the torsions, bond angles and bond lengths
	/// @brief note: the torsions/angles/lengths that are input to this
	/// @brief function will be identical to the ones in the pose,
	/// @brief i.e. only the dof values that are explicitly set by
	/// @brief this function plus the pivots will have changed after loop closure
	virtual
	void
	perturb_chain(
		core::pose::Pose const & pose,
		utility::vector1< core::Real> & torsions,
		utility::vector1< core::Real> & bond_ang,
		utility::vector1< core::Real> & bond_len
  ) const = 0;


	/// @brief after the kinmover has closed the loop, the perturber needs
	/// @brief to put the solutions into the pose
	/// @brief note: the base class version of this function sets the torsions,
	/// @brief so any KinematicPerturber that only changes the torsions will not have
	/// @brief to implement this function
	virtual
	void
	set_pose_after_closure(
		core::pose::Pose & pose,
		utility::vector1< core::Real> const & torsions,
		utility::vector1< core::Real> const & , //bond_ang,
		utility::vector1< core::Real> const & ,//bond_len,
		bool //closure_successful
  ) const;

	virtual
	bool perturber_exhausted() const {
		return false; }

	void
	set_max_sample_iterations( core::Size sample_its ){
		max_sample_iterations_ = sample_its; }

	core::Size
	max_sample_iterations() const {
		return max_sample_iterations_; }


protected:

	//for access to the kinematic mover that owns this perturber
	KinematicMoverCAP kinmover_;

	core::Size max_sample_iterations_;

private:
	core::kinematics::MoveMapCOP movemap_;

};



/// @brief torsion sampling kinematic perturber
class TorsionSamplingKinematicPerturber : public KinematicPerturber {

public:

	typedef KinematicPerturber parent;

	TorsionSamplingKinematicPerturber( KinematicMoverCAP kinmover_in );

	~TorsionSamplingKinematicPerturber();

	std::string perturber_type() const {
		return "TorsionSampleKinematicPerturber"; }

	///@brief varies torsions always and bond angles sometimes.  Currently torsion varying will respect a movemap if present; angles do NOT look for a movemap.
	void
	perturb_chain(
		core::pose::Pose const & pose,
		utility::vector1< core::Real> & torsions,
		utility::vector1< core::Real> & bond_ang,
		utility::vector1< core::Real> & //bond_len
  ) const;


	void
	set_pose_after_closure(
		core::pose::Pose & pose,
		utility::vector1< core::Real> const & torsions,
		utility::vector1< core::Real> const & bond_ang,
		utility::vector1< core::Real> const & bond_len,
		bool closure_successful
  ) const;

	void
	set_vary_ca_bond_angles( bool vary_ca_bond_angles ) {
		vary_ca_bond_angles_ = vary_ca_bond_angles; }

	void
	set_sample_vicinity( bool sample_vicinity ) {
		sample_vicinity_ = sample_vicinity; }

	void
	set_degree_vicinity( core::Real degree_vicinity ) {
		degree_vicinity_ = degree_vicinity; }

private:

	bool vary_ca_bond_angles_;
	bool sample_vicinity_;
	core::Real degree_vicinity_;
	bool sample_omega_for_pre_prolines_;

	core::scoring::Ramachandran const & rama_;

};


/// @brief WARNING WARNING UNTESTED!!!! torsion sweeping kinematic perturber
/// @brief WARNING WARNING UNTESTED!!!! used to work in other implementation
/// @brief probably works now, but to make sure you shoud doublecheck
//  @brief if you plan to use
///@details Here's some commentary on why the TorsionSweepingKinematicPerturber is Fun and Awesome to use!  Basically, you take the loop as given and search for solutions nearby.  Its purpose is not to produce large changes, but instead to produce small perturbations to find the bottom of the current energy well.  It's not meant to be used for the original KIC protocol, in which the pivots are changed regularly: you can't sweep through nonpivot torsion space if the pivots don't stay the same.
class TorsionSweepingKinematicPerturber : public KinematicPerturber {

public:

	TorsionSweepingKinematicPerturber();
	~TorsionSweepingKinematicPerturber();

	std::string perturber_type() const {
		return "TorsionSweepingKinematicPerturber"; }

	///@brief movemap control NOT IMPLEMENTED in TorsionSweepingKP.  It is also NOT NEEDED because you can use set_nonpivot_res_to_sweep instead.
	void
	perturb_chain(
		core::pose::Pose const & pose,
		utility::vector1< core::Real > & torsions,
		utility::vector1< core::Real > & bond_ang,
		utility::vector1< core::Real > & //bond_len
  ) const;

	bool perturber_exhausted() const {
		return sweep_iterator_.at_end(); }

	void set_nonpivot_res_to_sweep( utility::vector1< Size > const & resids );
	void set_nonpivot_bb_torsion_id( utility::vector1< Size > const & bbtorids );
	void set_sweep_start_angle( utility::vector1< core::Real > const & angles_in_degrees );
	void set_sweep_step_size( utility::vector1< core::Real > const & angle_steps_in_degrees );
	void set_sweep_nsteps( utility::vector1< Size > const & nsteps );


private:

	utility::vector1< core::Size > nonpivot_res_to_sweep_;
	utility::vector1< core::Size > sweep_torsion_ids_;
	utility::vector1< core::Real > sweep_nonpivot_torsion_starts_;
	utility::vector1< core::Real > sweep_step_sizes_;
	mutable utility::LexicographicalIterator sweep_iterator_;

};



} // end namespace kinematic_closure
} // end namespace moves
} // end namespace protocols


#endif
