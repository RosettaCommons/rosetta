// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (C) 199x-2008 Hebrew University, Jerusalem
//
/// @file   FlexPepDockingAbInitio.hh
///
/// @brief low-resolution agressive ab-initio folding and docking of peptide
/// @date August 5, 2008
/// @author Barak Raveh

#ifndef INCLUDED_protocols_flexpep_docking_FlexPepDockingAbInitio_hh
#define INCLUDED_protocols_flexpep_docking_FlexPepDockingAbInitio_hh


#include <protocols/flexpep_docking/FlexPepDockingFlags.hh>
#include <protocols/flexpep_docking/FlexPepDockingAbInitio.fwd.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.hh>
#include <string>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace flexpep_docking {

class FlexPepDockingAbInitio : public moves::Mover
{
public:

	/// @brief
	/// constructor for low resolution flexpible peptide docking
	//
	// @param[in] scorefxn_in
	//            The scoring function used for optimization
	// @param[in] rb_jump
	//            The FoldTree rigid body jump over
	//            which rigid-body pertrubations are made
	FlexPepDockingAbInitio(
		FlexPepDockingFlags flags_in,
		core::scoring::ScoreFunctionOP scorefxn_in,
		core::kinematics::MoveMapOP movemap_in,
		Size const rb_jump_in = 1
	);

	// empty destructor - for good inclusion of OP clasesses
	~FlexPepDockingAbInitio();

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	private:

	/// @brief initial setup for apply
	void setup_for_apply(core::pose::Pose& pose);

	// switch pose to centroid mode, if in full-atom
	void to_centroid ( core::pose::Pose & pose ) const;

	// switch pose to full-atom mode, if referncePose is full-atom
	// (using side-chains of referencePose)
	void to_allatom ( core::pose::Pose & pose, core::pose::Pose& referencePose ) const;

	// pose - the pose to perturb
	// cycles - # of MC cycles
	// acceptance_rate [Output] - the acceptance rate in the MC cycles
	void torsions_monte_carlo(
		core::pose::Pose & pose,
		const int cycles,
		double& acceptance_rate
	);

	// pose - the pose to loop close, containing the peptide
	void loopclosure_monte_carlo(
		core::pose::Pose & pose
	);

	// pose - the pose to perturb
	// cycles - # of MC cycles
	// trans_magnitude = random translation magnitude
	// rot_magnitude =  random rotation magnitude
	// acceptance_rate [Output] - the acceptance rate in the MC cycles
	void rigidbody_monte_carlo(
		core::pose::Pose & pose,
		const int cycles,
		const float trans_magnitude,
		const float rot_magnitude,
		double& acceptance_rate
	);

	/////////////////////////////////////////////////////////////////////////////
	// @brief convert all pSer residues to Asp, and save list of pSers to pSer_positions
	//
	// @param pose pose for conversion
	// @param pSer_positions[out] - positions that were converted from pSer to Asp
	/////////////////////////////////////////////////////////////////////////////
	void convertPSERtoASP(core::pose::Pose& pose, std::set<int>& pSer_positions);

	/////////////////////////////////////////////////////////////////////////////
	// @brief convert all pSer residues to Glu, and save list of pSers to pSer_positions
	//
	// @param pose pose for conversion
	// @param pSer_positions[out] - positions that were converted from pSer to Asp
	/////////////////////////////////////////////////////////////////////////////
        void convertPSERtoGLU(core::pose::Pose& pose, std::set<int>& pSer_positions);


	/////////////////////////////////////////////////////////////////////////////
	// @brief restore pSer residues if they were changed to Asp/Glu (based on pSer_positions)
	//
	// @param pose pose for conversion
	// @param pSer_positions - positions for pSer'ilation
	/////////////////////////////////////////////////////////////////////////////
	void restorePSER(core::pose::Pose& pose, std::set<int> const& pSer_positions);



private:
	FlexPepDockingFlags flags_; // all flags loaded from cmd-line options

	// energy scoring function for centroid optimization
	core::scoring::ScoreFunctionOP scorefxn_;

	// the flexpepdock protocol movemap // may change throughout the run
	core::kinematics::MoveMapOP movemap_;

	/// the jump number across which to do rigid_body transformations
	Size rb_jump_;

	/// whether or not to initialize the viewer (for opengl)
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool view_;

	moves::MonteCarloOP mc_;

	protocols::simple_moves::MinMoverOP minimizer_;

	// loop mover for modeling loop closure // TODO: this is a wrapper, use the loop modeller directly
	protocols::comparative_modeling::LoopRelaxMoverOP loop_relax_mover_;

        // set of fragments from a fragment library
	core::fragment::ConstantLengthFragSetOP fragset3mer_;
	core::fragment::ConstantLengthFragSetOP fragset9mer_;
	core::fragment::ConstantLengthFragSetOP fragset5mer_;

};

} // flexPepDocking
} // protocols

#endif

