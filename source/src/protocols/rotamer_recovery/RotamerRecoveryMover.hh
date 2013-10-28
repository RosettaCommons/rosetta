// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RotamerRecoveryMover.hh
/// @brief A wrapper that measures how similar the rotamers are between before and after running the child mover
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author P. Douglas Renfrew

#ifndef INCLUDED_protocols_rotamer_recovery_RotamerRecoveryMover_hh
#define INCLUDED_protocols_rotamer_recovery_RotamerRecoveryMover_hh

// Unit Headers
#include <protocols/rotamer_recovery/RotamerRecovery.fwd.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>






namespace protocols {
namespace rotamer_recovery {

class RotamerRecoveryMover : public moves::Mover {
public: // type definitions

public: // constructors destructors

	/// @brief default constructor
	RotamerRecoveryMover();

	RotamerRecoveryMover(
		rotamer_recovery::RotamerRecoveryOP rotamer_recovery,
		core::scoring::ScoreFunctionOP scfxn,
		core::pack::task::TaskFactoryOP task_factory);

	~RotamerRecoveryMover();

	RotamerRecoveryMover( RotamerRecoveryMover const & src);

public: // functional interface

	virtual
	void
	register_options() const ;

	virtual
	void
	apply( core::pose::Pose & pose );

	virtual
	std::string
	get_name() const;

	/// @brief make a copy
	virtual moves::MoverOP clone() const;

	/// @brief make a copy but use default constructor
	virtual moves::MoverOP fresh_instance() const;

	virtual
	void
	parse_my_tag(
		TagCOP const /*tag*/,
		basic::datacache::DataMap & /*data*/,
		Filters_map const & /*filters*/,
		moves::Movers_map const & movers,
		Pose const &);

	///@brief this function informs the job distributor (august 08
	///vintage) whether this object needs to be freshly regenerated on
	///each use.
	virtual
	bool
	reinitialize_for_each_job() const;


	///@brief this function informs the job distributor (august 08
	///vintage) whether this object needs to be regenerated when the
	///input pose is about to change (for example, if the mover has
	///special code on the first apply() that is only valid for that
	///one input pose).
	virtual
	bool
	reinitialize_for_new_input() const;

	core::scoring::ScoreFunctionOP
	score_function();

	void
	score_function(
		core::scoring::ScoreFunctionOP);

	virtual
	void
	show() const;

	virtual
	void
	show(
		std::ostream & out
	) const;

private: // data

	rotamer_recovery::RotamerRecoveryOP rotamer_recovery_;
	core::scoring::ScoreFunctionOP scfxn_;
	core::pack::task::TaskFactoryOP task_factory_;
};

} // namespace rotamer_recovery
} // namespace protocols

#endif
