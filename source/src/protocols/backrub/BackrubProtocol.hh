// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file .../BackrubProtocol.hh
/// @brief
/// @author


#ifndef INCLUDED_protocols_backrub_BackrubProtocol_hh
#define INCLUDED_protocols_backrub_BackrubProtocol_hh

#include <protocols/backrub/BackrubProtocol.fwd.hh>

//Protocols
#include <protocols/backrub/BackrubMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.fwd.hh>

//Core
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
//Basic
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh> // AUTO IWYU For vector1


// Forward
namespace protocols {
namespace backrub {

class BackrubProtocol : public protocols::moves::Mover {

public:
	BackrubProtocol();
	BackrubProtocol(BackrubProtocol const & bp);
	~BackrubProtocol() override;

	void
	apply( core::pose::Pose& pose ) override;



	/// @brief Set the pivot residues for backrub.
	/// Will use contiguous residues as segments.
	/// Within these segments, backrub will occur.
	void
	set_pivot_residues( utility::vector1<core::Size> pivot_residues);

	/// @brief Set the type of atoms to use for backrub.  These are atom names.
	/// Default is to only use CA as pivots.
	void
	set_pivot_atoms( utility::vector1<std::string> pivot_atoms);


	/// @brief Set the movemap for the whole protocol.
	/// Uses contiguous backbone regions for backrub to set pivot residues
	/// Within these segments, backrub will occur.
	void
	set_movemap( core::kinematics::MoveMapCOP movemap);

	/// @brief Set the MoveMap that will only be used for the SmallMover.
	void
	set_movemap_smallmover(core::kinematics::MoveMapCOP movemap);

	/// @brief Set the taskfactory used for sidechain moves.
	/// NOTE: Clones the TF.  Sets C-Beta off as they are not compatable with branch angle optimization.
	void
	set_taskfactory( core::pack::task::TaskFactoryCOP tf);

	/// @brief Control whether the mover dumps the _last and _low intermediate structures.
	void
	set_dump_poses(bool);

	////////////////////////////////////////////////////////////////////////
	/// ScoreFunction
	///

	/// @brief Set the Scorefunction.
	/// NOTE: Clones the scorefxn.
	void
	set_scorefunction(core::scoring::ScoreFunctionCOP scorefxn);


	/// @brief Set a pre-configured Backrub Mover.
	void
	set_backrub_mover(protocols::backrub::BackrubMoverOP backrub_mover);

	void
	write_database();

	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	utility::vector1<core::Size>
	get_pivot_residues() const;

	void
	set_pivots_from_residue_subset( core::select::residue_selector::ResidueSubset residue_subset );

private:

	void
	read_cmd_line_options();

	///@brief Sets all member variables that are options that control behavior.
	/// These options can be set via parse_my_tag or command line options.
	void
	set_options(
		utility::vector1<core::Size> pivot_residues,
		utility::vector1<std::string> pivot_atoms,
		core::kinematics::MoveMapCOP minimize_movemap,
		core::kinematics::MoveMapCOP movemap_smallmover,
		core::pack::task::operation::TaskOperationCOP packing_operation,
		core::Size min_atoms,
		core::Size max_atoms,
		bool initial_pack,
		core::Real mm_bend_weight,
		core::Real sm_prob,
		core::Real sc_prob,
		core::Real sc_prob_uniform,
		core::Real sc_prob_withinrot,
		core::Real mc_kt,
		core::Size ntrials,
		bool trajectory,
		bool trajectory_gz,
		core::Size trajectory_stride
	);

	void
	finalize_setup(core::pose::Pose & pose);

private:

	core::scoring::ScoreFunctionOP scorefxn_;
	core::pack::task::TaskFactoryOP main_task_factory_;
	protocols::backrub::BackrubMoverOP backrubmover_;

	//Originally NOT Owning pointers!
	protocols::simple_moves::SmallMoverOP smallmover_;
	protocols::simple_moves::sidechain_moves::SidechainMoverOP sidechainmover_;
	protocols::minimization_packing::PackRotamersMoverOP packrotamersmover_;

private:
	//New
	core::Size ntrials_;
	core::kinematics::MoveMapCOP movemap_smallmover_;
	core::kinematics::MoveMapCOP minimize_movemap_;

	core::pack::task::operation::TaskOperationCOP packing_operation_;

	utility::vector1<core::Size> pivot_residues_;
	utility::vector1<std::string> pivot_atoms_;

	core::Size min_atoms_;
	core::Size max_atoms_;

	core::Real mm_bend_weight_;
	core::Real sm_prob_; //SmallMover probability.  Default 0.
	core::Real sc_prob_; //SideChain Refinement probability.

	core::Real sc_prob_uniform_; //probability of uniformly sampling chi angles
	core::Real sc_prob_withinrot_; //probability of sampling within the current rotamer
	core::Real mc_kt_;

	bool initial_pack_;
	bool trajectory_;
	bool trajectory_gz_;
	bool recover_low_;
	bool dump_poses_ = true; // Historical behavior
	core::Size trajectory_stride_;

	protocols::moves::MoverOP trajectory_apply_mover_;

	///@brief Residue selector to select pivot residues
	core::select::residue_selector::ResidueSelectorCOP pivots_residue_selector_;

};


} //backrub
} //protocols







#endif
