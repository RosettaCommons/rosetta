// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author Modified by Steven Lewis

#ifndef INCLUDED_protocols_simple_moves_RotamerTrialsMover_hh
#define INCLUDED_protocols_simple_moves_RotamerTrialsMover_hh

// Unit headers
#include <protocols/simple_moves/RotamerTrialsMover.fwd.hh>

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class RotamerTrialsMover : public protocols::moves::Mover {
public:

	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;
	typedef protocols::moves::MoverOP MoverOP;

public:

	// default constructor
	RotamerTrialsMover();

	/// @brief constructor with PackerTask. use a PackerTask ONLY for fixed-sequence work.
	/// WARNING TO ANY DESIGNER WHO PASSES IN A TASK: YOUR DESIGN STEPS WILL BE UNDONE
	/// AS THIS TASK CONCEIVES OF THE INPUT SEQUENCE THAT CORRESPONDS TO THE ORIGINAL SEQUENCE
	RotamerTrialsMover(
		ScoreFunctionCOP scorefxn_in,
		PackerTask & task_in
	);

	/// @brief constructor with TaskFactory
	RotamerTrialsMover(
		ScoreFunctionCOP scorefxn_in,
		TaskFactoryCOP factory_in
	);

	/// @brief copy constructor
	RotamerTrialsMover( RotamerTrialsMover const & rval );

	/// @brief destructor
	virtual ~RotamerTrialsMover();

	/// @brief clone this object
	virtual protocols::moves::MoverOP clone() const;

	/// @brief create this type of object
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual void show(std::ostream & output=std::cout) const;

	//PackerTaskMover/RotamerTrialsMover needs to have a parent class that implements this?
	//bool task_is_valid( core::pose::Pose const & pose ) const;

	// setters
	void score_function( core::scoring::ScoreFunctionCOP sf );
	void task_factory( core::pack::task::TaskFactoryCOP tf );

public:

	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

protected:

	/// @brief read access for derived classes
	ScoreFunctionCOP
	scorefxn() const;

	/// @brief read access for derived classes, pose needed to run TaskFactory
	PackerTaskCOP
	task( core::pose::Pose const & pose ) const;

private:

	/// @brief RotamerTrailsMover does not have its own score function, rather, it shares
	/// one with other classes -- CAUTION: the score function is externally modifiable.
	ScoreFunctionCOP scorefxn_;

	/// @brief use a PackerTask ONLY for fixed-sequence work.
	/// WARNING TO ANY DESIGNER WHO PASSES IN A TASK: YOUR DESIGN STEPS WILL BE UNDONE
	/// AS THIS TASK CONCEIVES OF THE INPUT SEQUENCE THAT CORRESPONDS TO THE ORIGINAL SEQUENCE
	///If a factory is present it overwrites this task with each call to apply()
	PackerTaskOP task_;

	/// @brief TaskFactory allows for nonconstant sequences to be used with RotamerTrialsMover
	///CAUTION: the factory is externally modifiable.
	TaskFactoryCOP factory_;

	/// @brief showing contents of PackerTask ( default false )
	bool show_packer_task_;
};  // class RotamerTrialsMover

std::ostream &operator<< (std::ostream &os, RotamerTrialsMover const &mover);


class EnergyCutRotamerTrialsMover : public protocols::simple_moves::RotamerTrialsMover {
public:

	// default constructor
	EnergyCutRotamerTrialsMover();

	// constructor with arguments
	EnergyCutRotamerTrialsMover(
		ScoreFunctionCOP scorefxn_in,
		PackerTask & task_in,
		protocols::moves::MonteCarloOP mc_in,
		core::Real energycut_in
	);

	// constructor with arguments
	EnergyCutRotamerTrialsMover(
		ScoreFunctionCOP scorefxn_in,
		TaskFactoryCOP factory_in,
		protocols::moves::MonteCarloOP mc_in,
		core::Real energycut_in
	);

	virtual
	~EnergyCutRotamerTrialsMover();

public:

	/// @brief apply this mover to a pose
	virtual
	void
	apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

protected:

	/// @brief selects a subset of residues to repack based on the per
	/// residue energies of the last accepted pose in the MC object.
	void
	setup_energycut_task(
		core::pose::Pose const & pose,
		protocols::moves::MonteCarlo const & mc,
		core::pack::task::PackerTask & task_in
	) const;

	protocols::moves::MonteCarloOP
	mc();

private:

	// data
	protocols::moves::MonteCarloOP mc_;
	core::Real energycut_;
};  // class EnergyCutRotamerTrialsMover

} // moves
} // protocols

#endif
