// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh
/// @brief  Header file for KinematicPerturbers used by the kineamtic mover
/// @author Florian Richter, floric@u.washington.edu, march 2009
/// @author Amelie Stein, amelie.stein@ucsf.edu, October 2012 -- refactoring vicinity sampling & new perturbers

#ifndef INCLUDED_protocols_loops_loop_closure_kinematic_closure_KinematicPerturber_hh
#define INCLUDED_protocols_loops_loop_closure_kinematic_closure_KinematicPerturber_hh

#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.fwd.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.fwd.hh>

// Rosetta Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/ppo_torsion_bin.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/Ramachandran.fwd.hh>
#include <core/scoring/Ramachandran2B.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

// Utility Headers
#include <utility/LexicographicalIterator.hh> //not OPable; needs full header
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>
#include <map>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace kinematic_closure {

/// @brief pure virtual base class for KinematicPerturber.  KinematicPerturbers determine HOW loops should be perturbed.  The base class contains a provision for determining WHERE they should be perturbed: MoveMap sensitivity.
class KinematicPerturber : public utility::pointer::ReferenceCount {

public:

	KinematicPerturber();
	virtual ~KinematicPerturber();

	virtual
	std::string perturber_type() const = 0;

	void
	set_kinmover( KinematicMoverCAP kinmover ){ kinmover_ = kinmover; }

	void
	set_movemap( core::kinematics::MoveMapCOP mm );

	core::kinematics::MoveMapCOP
	get_movemap() const ;

	/// @brief function that perturbs the chain, i.e. sets new values
	/// for the torsions, bond angles and bond lengths
	/// note: the torsions/angles/lengths that are input to this
	/// function will be identical to the ones in the pose,
	/// i.e. only the dof values that are explicitly set by
	/// this function plus the pivots will have changed after loop closure
	virtual
	void
	perturb_chain(
		core::pose::Pose const & pose,
		utility::vector1< core::Real> & torsions,
		utility::vector1< core::Real> & bond_ang,
		utility::vector1< core::Real> & bond_len
	) = 0; // needs to be non-const for Taboo Sampling


	/// @brief after the kinmover has closed the loop, the perturber needs
	/// to put the solutions into the pose
	/// note: the base class version of this function sets the torsions,
	/// so any KinematicPerturber that only changes the torsions will not have
	/// to implement this function
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

	void clear_torsion_string_stack() { } // only actually used for TabooSampling

	KinematicMoverCAP kinmover() const;

private:

	//for access to the kinematic mover that owns this perturber
	KinematicMoverCAP kinmover_;
	core::Size max_sample_iterations_;
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

	///@brief Varies torsions of a beta-amino acid residue based on minimia in the beta-amino acid Ramachandran cube.  This randomly picks a minimum, then chooses phi/theta/psi values randomly in a Gaussian centered in that minimum.
	void perturb_beta_residue (core::Real &phi, core::Real &theta, core::Real &psi, const core::Size beta_residue_type) const;

	///@brief Initialize positions of minima in the beta-amino acid Ramachandran cube.
	void initialize_betaresidue_minima (
		utility::vector1 < core::Real > &philist, //outputs -- will be cleared by this function
		utility::vector1 < core::Real > &thetalist, //outputs -- will be cleared by this function
		utility::vector1 < core::Real > &psilist, //outputs -- will be cleared by this function
		const core::Size mode //mode 1 initializes for beta-3-amino acids.
	) const;

	///@brief varies torsions always and bond angles sometimes.  Currently torsion varying will respect a movemap if present; angles do NOT look for a movemap.
	void
	perturb_chain(
		core::pose::Pose const & pose,
		utility::vector1< core::Real> & torsions,
		utility::vector1< core::Real> & bond_ang,
		utility::vector1< core::Real> & //bond_len
	);


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

private:

	bool vary_ca_bond_angles_;
	bool sample_omega_for_pre_prolines_;

	core::scoring::Ramachandran const & rama_;

};


/// @brief vicinity sampling kinematic perturber
/// @author Amelie Stein (just the refactoring)
class VicinitySamplingKinematicPerturber : public KinematicPerturber {

public:

	typedef KinematicPerturber parent;

	VicinitySamplingKinematicPerturber( KinematicMoverCAP kinmover_in );

	~VicinitySamplingKinematicPerturber();

	std::string perturber_type() const {
		return "VicinitySampleKinematicPerturber"; }

	///@brief varies torsions always and bond angles sometimes.  Currently torsion varying will respect a movemap if present; angles do NOT look for a movemap. -- note that the analytic closure and pivot selection currently do not respect movemaps though
	void
	perturb_chain(
		core::pose::Pose const & pose,
		utility::vector1< core::Real> & torsions,
		utility::vector1< core::Real> & bond_ang,
		utility::vector1< core::Real> & //bond_len
	);


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
		vary_ca_bond_angles_ = vary_ca_bond_angles;
	}


	void
	set_degree_vicinity( core::Real degree_vicinity ) {
		degree_vicinity_ = degree_vicinity;
	}

private:

	bool vary_ca_bond_angles_;
	core::Real degree_vicinity_;
	bool sample_omega_for_pre_prolines_;

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
	);

	bool perturber_exhausted() const { return sweep_iterator_.at_end(); }

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

}; // TorsionSweepingKinematicPerturber



/// @author Amelie Stein
/// @brief neighbor-dependent torsion sampling kinematic perturber -- uses rama2b for phi/psi lookup
class NeighborDependentTorsionSamplingKinematicPerturber : public KinematicPerturber {

public:

	typedef KinematicPerturber parent;

	NeighborDependentTorsionSamplingKinematicPerturber( KinematicMoverCAP kinmover_in );

	~NeighborDependentTorsionSamplingKinematicPerturber();

	std::string perturber_type() const {
		return "NeighborDependentTorsionSamplingKinematicPerturber"; }

	void
	perturb_chain(
		core::pose::Pose const & pose,
		utility::vector1< core::Real> & torsions,
		utility::vector1< core::Real> & bond_ang,
		utility::vector1< core::Real> & //bond_len
	);


	void
	set_pose_after_closure(
	  core::pose::Pose & pose,
	  utility::vector1< core::Real> const & torsions,
	  utility::vector1< core::Real> const & bond_ang,
	  utility::vector1< core::Real> const & bond_len,
	  bool closure_successful
	) const;

	void
	set_vary_ca_bond_angles( bool vary_ca_bond_angles ) { vary_ca_bond_angles_ = vary_ca_bond_angles; }


private:

	bool vary_ca_bond_angles_;
	bool sample_omega_for_pre_prolines_;

	core::scoring::Ramachandran2B const & rama_;

};


/// @brief torsion-restricted kinematic perturber (still samples randomly, but only within a given torsion bin)
/// @author Amelie Stein
class TorsionRestrictedKinematicPerturber : public KinematicPerturber {
public:

	typedef KinematicPerturber parent;

public:
	TorsionRestrictedKinematicPerturber(
		KinematicMoverCAP kinmover_in,
		core::conformation::torsion_bin_string const & torsion_bins
	);

	~TorsionRestrictedKinematicPerturber();

	std::string perturber_type() const {
		return "TorsionRestrictedKinematicPerturber"; }

	void
	perturb_chain(
				  core::pose::Pose const & pose,
				  utility::vector1< core::Real> & torsions,
				  utility::vector1< core::Real> & bond_ang,
				  utility::vector1< core::Real> & //bond_len
				  ) ;


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


private:

	bool vary_ca_bond_angles_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool sample_omega_for_pre_prolines_;
	core::conformation::torsion_bin_string predefined_torsions_;

	core::scoring::Ramachandran const & rama_;

};

/// @brief Taboo-sampling perturber base class, the two variants of which, share much
/// code in common, but interface with separate Ramachandran potentials.
class BaseTabooPerturber : public KinematicPerturber {
public:
	typedef KinematicPerturber parent;

public:
	// check if both c'tors are used -- if not, remove the unused one
	BaseTabooPerturber( KinematicMoverCAP kinmover_in );

	virtual ~BaseTabooPerturber();

	void
	perturb_chain(
		core::pose::Pose const & pose,
		utility::vector1< core::Real> & torsions,
		utility::vector1< core::Real> & bond_ang,
		utility::vector1< core::Real> & //bond_len
	);

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


	void clear_torsion_string_stack() {
		random_torsion_strings_.clear();
	}

	core::Size num_strings() const { return num_strings_; }

	core::conformation::torsion_bin_string
	next_torsion_string();

private:

	virtual
	void
	get_random_phi_psi_for_residue(
		core::pose::Pose const & pose,
		core::Size resid,
		core::conformation::ppo_torsion_bin torbin,
		core::Real & phi,
		core::Real & psi
	) const = 0;

	void
	refill_torsion_string_vector();

	/// @brief This function is implemented by the derived classes and is used to interface
	/// with the particular Ramachandran library that the derived classs relies upon
	/// it is invoked inside of refill_torsion_string_vector.
	virtual
	std::map< core::conformation::ppo_torsion_bin, core::Size >
	get_entries_per_torsion_bin(
		utility::vector1< core::chemical::AA > loop_seq,
		core::Size resid
	) const = 0;

private:
	bool vary_ca_bond_angles_;
	bool sample_omega_for_pre_prolines_;

	core::Size num_strings_;

	// holds a list of random torsion bin strings to be sampled, to ensure diversity -- accessed
	// (and filled, if necessary) by next_torsion_string( pose )
	utility::vector1< core::conformation::torsion_bin_string > random_torsion_strings_;

};

/// @brief Taboo-sampling kinematic perturber (still samples randomly, but only within a specific torsion bin, and
/// the Taboo sampler ensures that this torsion bin is varied in each iteration)
/// @author Amelie Stein
class TabooSamplingKinematicPerturber : public BaseTabooPerturber {
public:

	typedef BaseTabooPerturber parent;

public:

	TabooSamplingKinematicPerturber( KinematicMoverCAP kinmover_in );

	virtual ~TabooSamplingKinematicPerturber();

	std::string perturber_type() const {
		return "TabooSamplingKinematicPerturber"; }

private:

	virtual
	void
	get_random_phi_psi_for_residue(
		core::pose::Pose const & pose,
		core::Size resid,
		core::conformation::ppo_torsion_bin torbin,
		core::Real & phi,
		core::Real & psi
	) const;

	virtual
	std::map< core::conformation::ppo_torsion_bin, core::Size >
	get_entries_per_torsion_bin(
		utility::vector1< core::chemical::AA > loop_seq,
		core::Size resid
	) const;

private:

	core::scoring::Ramachandran const & rama_;

};



/// @brief Neighbor-dependent Taboo-sampling kinematic perturber (still samples randomly, but only
/// within a given torsion bin; the Taboo sampler ensures that this torsion bin is varied in each
/// iteration) that uses neighbor-dependent Ramachandran distributions (rama2b)
/// @author Amelie Stein
/// @date Mon May 21 11:39:26 PDT 2012
class NeighborDependentTabooSamplingKinematicPerturber : public BaseTabooPerturber {
public:

	typedef KinematicPerturber parent;

public:

	NeighborDependentTabooSamplingKinematicPerturber( KinematicMoverCAP kinmover_in );

	~NeighborDependentTabooSamplingKinematicPerturber();

	std::string perturber_type() const {
		return "NeighborDependentTabooSamplingKinematicPerturber"; }

private:

	virtual
	void
	get_random_phi_psi_for_residue(
		core::pose::Pose const & pose,
		core::Size resid,
		core::conformation::ppo_torsion_bin torbin,
		core::Real & phi,
		core::Real & psi
	) const;

	virtual
	std::map< core::conformation::ppo_torsion_bin, core::Size >
	get_entries_per_torsion_bin(
		utility::vector1< core::chemical::AA > loop_seq,
		core::Size resid
	) const;

private:

	core::scoring::Ramachandran2B const & rama_;

}; // NeighborDependentTabooSamplingKinematicPerturber

} // namespace kinematic_closure
} // namespace loop_closure
} // namespace moves
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_kinematic_closure_KinematicPerturber_hh
