// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RotamerTrialsMinMover.hh
/// @brief Wrapper mover for Rotamer-Trials with Minimization (based on RotamerTrialsMover)
/// @author Barak Raveh

#ifndef INCLUDED_protocols_simple_moves_RotamerTrialsMinMover_hh
#define INCLUDED_protocols_simple_moves_RotamerTrialsMinMover_hh

// Unit headers
#include <protocols/simple_moves/RotamerTrialsMinMover.fwd.hh>

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class RotamerTrialsMinMover : public protocols::moves::Mover {
public:

	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;

public:

	// default constructor
	RotamerTrialsMinMover();

	/// @brief constructor with PackerTask. use a PackerTask ONLY for fixed-sequence work.
	/// WARNING TO ANY DESIGNER WHO PASSES IN A TASK: YOUR DESIGN STEPS WILL BE UNDONE
	/// AS THIS TASK CONCEIVES OF THE INPUT SEQUENCE THAT CORRESPONDS TO THE ORIGINAL SEQUENCE
	///
	/// @param scorefxn_in The score function used for packing and minimization (which may be modified externally)
	/// @param task_in The task that will pass to the rotamers packer
	RotamerTrialsMinMover(
		ScoreFunctionCOP scorefxn_in,
		PackerTask & task_in
	);

	/// @brief constructor with TaskFactory for producing the packer task
	///
	/// @param scorefxn_in The score function used for packing and minimization (which may be modified externally)
	/// @param factory_in The task that will pass to the rotamers packer
	RotamerTrialsMinMover(
		ScoreFunctionCOP scorefxn_in,
		TaskFactoryCOP factory_in
	);

	void
	init();

	virtual ~RotamerTrialsMinMover();

	/// @brief Apply Rotamer-Trials with minimization to pose, using the score function
	///       and tasks provided by the constructor
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual void show(std::ostream & output=std::cout) const;

	//PackerTaskMover/RotamerTrialsMinMover needs to have a parent class that implements this?
	//bool task_is_valid( core::pose::Pose const & pose ) const;

	// setters
	void score_function( core::scoring::ScoreFunctionCOP sf );

	void task_factory( core::pack::task::TaskFactoryCOP tf );

	/// @brief Parse XML for RosettaScripts
	virtual void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	/// @brief Return a new mover instance (for RosettaScripts)
	virtual protocols::moves::MoverOP fresh_instance() const;
	/// @brief Return a copy of this mover instance (for RosettaScripts)
	virtual protocols::moves::MoverOP clone() const;

protected:

	/// @brief read access for derived classes
	ScoreFunctionCOP
	scorefxn() const;

	/// @brief read access for derived classes, pose needed to run TaskFactory
	PackerTaskOP
	task( core::pose::Pose const & pose ) const;

private:

	/// @brief RotamerTrailsMinMover does not have its own score function, rather, it shares
	/// one with other classes -- CAUTION: the score function is externally modifiable,
	/// but this is probably the way users expect this to behave when they use this mover).
	ScoreFunctionCOP scorefxn_;

	/// @brief use a PackerTask ONLY for fixed-sequence work.
	/// WARNING TO ANY DESIGNER WHO PASSES IN A TASK: YOUR DESIGN STEPS WILL BE UNDONE
	/// AS THIS TASK CONCEIVES OF THE INPUT SEQUENCE THAT CORRESPONDS TO THE ORIGINAL SEQUENCE
	///If a factory is present it overwrites this task with each call to apply()
	PackerTaskOP task_;

	/// @brief TaskFactory allows for nonconstant sequences to be used with RotamerTrialsMover
	///CAUTION: the factory is externally modifiable.
	TaskFactoryCOP factory_;
	bool nonideal_;
	bool cartesian_;
};  // class RotamerTrialsMinMover

std::ostream &operator<< (std::ostream &os, RotamerTrialsMinMover const &mover);


/// @brief Wrapper for Rotamer Trials with Minimization, which modifies only
///        rotamers whose energy changed by a given constant
class EnergyCutRotamerTrialsMinMover : public protocols::simple_moves::RotamerTrialsMinMover {
public:

	// default constructor
	EnergyCutRotamerTrialsMinMover();

	// constructor with arguments
	EnergyCutRotamerTrialsMinMover(
		ScoreFunctionCOP scorefxn_in,
		PackerTask & task_in,
		protocols::moves::MonteCarloOP mc_in,
		core::Real energycut_in
	);

	// constructor with arguments
	EnergyCutRotamerTrialsMinMover(
		ScoreFunctionCOP scorefxn_in,
		TaskFactoryCOP factory_in,
		protocols::moves::MonteCarloOP mc_in,
		core::Real energycut_in
	);

	virtual
	~EnergyCutRotamerTrialsMinMover();

public:

	/// @brief apply this mover to a pose
	virtual
	void
	apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

private:

	/// @brief selects a subset of residues to repack based on the per
	/// residue energies of the last accepted pose in the MC object.
	void
	setup_energycut_task(
		core::pose::Pose const & pose,
		protocols::moves::MonteCarlo const & mc,
		core::pack::task::PackerTask & task_in
	) const;

private:

	// data
	protocols::moves::MonteCarloOP mc_;
	core::Real energycut_;
};  // class EnergyCutRotamerTrialsMinMover

} // moves
} // protocols

#endif
