// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file   FlexPepDockingProtocol.hh
///
/// @brief protocol for docking flexible peptides onto globular proteins
/// @date August 5, 2008
/// @author Barak Raveh

#ifndef INCLUDED_protocols_flexpep_docking_FlexPepDockingProtocol_hh
#define INCLUDED_protocols_flexpep_docking_FlexPepDockingProtocol_hh

#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/flexpep_docking/FlexPepDockingFlags.fwd.hh>
#include <protocols/flexpep_docking/FlexPepDockingPoseMetrics.hh>
#include <string>

#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace flexpep_docking {

class FlexPepDockingProtocol;
typedef utility::pointer::shared_ptr< FlexPepDockingProtocol > FlexPepDockingProtocolOP;
typedef utility::pointer::shared_ptr< FlexPepDockingProtocol const > FlexPepDockingProtocolCOP;

class FlexPepDockingProtocol : public moves::Mover
{
public:

	/// constructor
	FlexPepDockingProtocol(Size const rb_jump_in = 1);

	FlexPepDockingProtocol(Size const rb_jump_in, bool const fullatom, bool const view=false );

	// empty destructor - needed for proper inclusion of OP clasesses
	~FlexPepDockingProtocol();

	/// @brief setup that is called from constructor
	void set_default();

	/// @brief setter
	void view( bool view_in ) { view_=view_in; }

	virtual protocols::moves::MoverOP clone() const;

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);
private:

	void setup_foldtree( core::pose::Pose & pose );

	void minimize_only(
		core::pose::Pose & pose,
		const std::string & min_type,
		const float min_func_tol
	);

	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	///
	/// @brief prepack protocol for flexpepdock
	//         (0)  minimize all side-chains + peptide b.b., just to get a reference to best possible energy
	//         (i)   translate apart (by 1000A)
	//         (ii)  repack + minimize side-chains of protein
	//         (iii)  translate back (by 1000A)
	///
	//  @details
	//  Prepacking a docked structure
	//
	//  @param
	//  pose - the pose to prepack
	//  ppk_receptor - whether to prepack the receptor protein
	//  ppk_peptide - whether to prepack the lignad peptide
	////////////////////////////////////////////////////////////////////////////
	void prepack_only( core::pose::Pose & pose ,
		bool ppk_receptor, bool ppk_peptide );

	////////////////////////////////////////
	// set constraints to prevent receptor from changing too much from starting structure
	////////////////////////////////////////
	void set_receptor_constraints( core::pose::Pose & pose );

	void set_allowed_moves( );

	void random_peptide_phi_psi_perturbation( core::pose::Pose & pose );

	void extend_peptide( core::pose::Pose & pose );

	void place_peptide_on_binding_site ( core::pose::Pose & pose );

	void flip_in_pcs( core::pose::Pose & pose );

	void SlideIntoContact( core::pose::Pose & pose, core::Vector translate_axis );

	void random_rb_pert( core::pose::Pose & pose );

	void small_moves( core::pose::Pose & pose );

	void shear_moves( core::pose::Pose & pose );

	void backrub_move( core::pose::Pose & pose );

	void polyAla( core::pose::Pose & pose );

	bool check_filters( core::pose::Pose & pose);

	void calcMetrics( core::pose::Pose & pose );

	// randomly change peptide's jump residue every cycle
	// TODO: add a flag for this? (currently not used)
	void randomlySlidePeptideJump(core::pose::Pose & pose);

	void torsions_monte_carlo_minimize(
		core::pose::Pose & pose,
		const int cycles,
		const std::string & min_type,
		const float minimization_threshold,
		const float min_func_tol
	);

	void rigidbody_monte_carlo_minimize(
		core::pose::Pose & pose,
		const int cycles,
		const std::string & min_type,
		const float trans_magnitude,
		const float rot_magnitude,
		const float minimization_threshold,
		const float func_tol
	);

	void peptide_random_loop_model(
		core::pose::Pose & pose
	);


	/////////////////////////////////////////////////////////////////////////////
	// @brief mark the peptide residues in the native structure interface
	//
	// @param superpos_partner[out]
	//        An array of positions - true for peptide residues
	// @param native_interface_residues[out]
	//        An array of positions - true for interface peptide residues
	/////////////////////////////////////////////////////////////////////////////
	void markNativeInterface
	( ObjexxFCL::FArray1D_bool & superpos_partner,
		ObjexxFCL::FArray1D_bool & native_interface_residues) const;

	/////////////////////////////////////////////////////////////////////////////
	// @brief
	// add low resolution statistics to score_map (see make_statistics)
	//
	// @param[in] start_pose- the initial pose sent to this->apply()
	// @param[in] pose_after_lowres - the pose optimized at low-res by this->apply()
	/////////////////////////////////////////////////////////////////////////////
	void addLowResStatistics(
		core::pose::Pose const& start_pose,
		core::pose::Pose& pose_after_lowres ) const;


	/////////////////////////////////////////////////////////////////////////////
	// @brief
	// make statistics comparing the start pose, the final pose and the native
	// and updated final_pose with statistics
	//
	// @param[in] start_pose - the initial pose sent to this->apply()
	// @param[in] pose_after_lowres - the pose after low-res optimization (if applied),
	//                                and right before the hi-res optimization
	// @param[in] final_pose - the final pose optimized by this->apply()
	/////////////////////////////////////////////////////////////////////////////
	void storeJobStatistics(
		core::pose::Pose const& start_pose,
		core::pose::Pose& pose_after_lowres,
		core::pose::Pose& final_pose );


	/////////////////////////////////////////////////////////////////////////////
	// @brief
	// High-resolution protocol for docking a flexible peptide onto a globular
	// protein
	//
	// @param pose[in,out] an input pose conformation to be optimized
	void hires_fpdock_protocol(core::pose::Pose& pose);


private:
	/// information about the mode
	bool fullatom_;
	/// the jump number across which to do rigid_body transformations
	Size rb_jump_;
	/// whether or not to initialize the viewer (for opengl)
	bool view_;

private:

	bool is_fail_; // flag for failures in class operation

	FlexPepDockingFlags flags_; // all flags loaded from cmd-line options

	core::Size nres_receptor_;
	core::Size nres_peptide_;

	moves::MonteCarloOP mc_;

	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP scorefxn_lowres_; // for low-res optimization

	// metrics calculator object
	FlexPepDockingPoseMetrics fpdock_metrics_;
	std::map < std::string, core::Real > if_metrics_; // interface metrics

	// taskfactory for all-protein repacking operations
	core::pack::task::TaskFactoryOP allprotein_tf_;

	// task factory that restricts packing to the interface
	core::pack::task::TaskFactoryOP interface_tf_;

	// for packing the docking interface
	protocols::simple_moves::PackRotamersMoverOP interface_packer_;

	// for designing the peptide;
	//core::pack::task::operation::RestrictResidueToRepackingOP receptor_protector_oper_; // operation to prevent the receptor from being redesigned
	//core::pack::task::PackerTaskOP design_task_;
	//protocols::simple_moves::PackRotamersMoverOP design_mover_;

	// the flexpepdock protocol movemap // may change throughout the run
	core::kinematics::MoveMapOP movemap_;

	// movemap for the minimizer
	core::kinematics::MoveMapOP movemap_minimizer_;

	// loop mover for modeling loop closure
	protocols::comparative_modeling::LoopRelaxMoverOP loop_relax_mover_;

	// int num_jumps_;

	// flag to indicate when the flexible loops are stripped from the pose
	// this is useful for various tasks when we want to work only on the core
	// of the rigid body (like R.B. perturbation)
	bool is_trimmed_loops_;

};
} // flexPepDocking
} // protocols

#endif
