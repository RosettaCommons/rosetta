// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file .../BackrubProtocol.hh
/// @brief
/// @author


#ifndef INCLUDED_protocols_backrub_BackrubProtocol_hh
#define INCLUDED_protocols_backrub_BackrubProtocol_hh

#include <protocols/backrub/BackrubProtocol.fwd.hh>

//Protocols
#include <protocols/backrub/BackrubMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.fwd.hh>

//Core
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
//Basic
#include <basic/datacache/DataMap.fwd.hh>


// Forward
namespace protocols {
namespace backrub {

class BackrubProtocol : public protocols::moves::Mover {

public:
	BackrubProtocol();
	BackrubProtocol(BackrubProtocol const & bp);
	virtual ~BackrubProtocol();

	virtual void
	apply( core::pose::Pose& pose );



	///@brief Set the pivot residues for backrub.
	/// Will use contiguous residues as segments.
	/// Within these segments, backrub will occur.
	void
	set_pivot_residues( utility::vector1<core::Size> pivot_residues);

	///@brief Set the type of atoms to use for backrub.  These are atom names.
	/// Default is to only use CA as pivots.
	void
	set_pivot_atoms( utility::vector1<std::string> pivot_atoms);


	///@brief Set the movemap for the whole protocol.
	/// Uses contiguous backbone regions for backrub to set pivot residues
	/// Within these segments, backrub will occur.
	void
	set_movemap( core::kinematics::MoveMapCOP movemap);

	///@brief Set the MoveMap that will only be used for the SmallMover.
	void
	set_movemap_smallmover(core::kinematics::MoveMapCOP movemap);

	///@brief Set the taskfactory used for sidechain moves.
	/// NOTE: Clones the TF.  Sets C-Beta off as they are not compatable with branch angle optimization.
	void
	set_taskfactory( core::pack::task::TaskFactoryCOP tf);

	////////////////////////////////////////////////////////////////////////
	/// ScoreFunction
	///

	///@brief Set the Scorefunction.
	/// NOTE: Clones the scorefxn.
	void
	set_scorefunction(core::scoring::ScoreFunctionCOP scorefxn);


	///@brief Set a pre-configured Backrub Mover.
	void
	set_backrub_mover(protocols::backrub::BackrubMoverOP backrub_mover);

	void
	write_database();

	std::string
	get_name() const;

	virtual protocols::moves::MoverOP
	clone() const;

	virtual protocols::moves::MoverOP
	fresh_instance() const;

	//virtual void
	//parse_my_tag(
	// TagCOP tag,
	// basic::datacache::DataMap & data,
	// Filters_map const & filters,
	// moves::Movers_map const & movers,
	// core::pose::Pose const & pose );

private:

	void
	read_cmd_line_options();

	void
	finalize_setup(core::pose::Pose & pose);

private:

	core::scoring::ScoreFunctionOP scorefxn_;
	core::pack::task::TaskFactoryOP main_task_factory_;
	protocols::backrub::BackrubMoverOP backrubmover_;

	//Originally NOT Owning pointers!
	protocols::simple_moves::SmallMoverOP smallmover_;
	protocols::simple_moves::sidechain_moves::SidechainMoverOP sidechainmover_;
	protocols::simple_moves::PackRotamersMoverOP packrotamersmover_;

private:
	//New
	core::Size ntrials_;
	core::kinematics::MoveMapCOP movemap_smallmover_;
	core::kinematics::MoveMapCOP movemap_preminimization_;

	utility::vector1<core::Size> pivot_residues_;
	utility::vector1<std::string> pivot_atoms_;

	core::Size min_atoms_;
	core::Size max_atoms_;

	core::Real sm_prob_; //SmallMover probability.  Default 0.
	core::Real sc_prob_; //SideChain Refinement probability.

	core::Real sc_prob_uniform_; //probability of uniformly sampling chi angles
	core::Real sc_prob_withinrot_; //probability of sampling within the current rotamer
	core::Real mc_kt_;

	bool initial_pack_;
};


} //backrub
} //protocols







#endif







