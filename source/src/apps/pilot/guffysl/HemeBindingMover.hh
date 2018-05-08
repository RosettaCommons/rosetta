// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/guffysl/heme_binding.cc
/// @brief Application to redesign an enzyme to bind an alternative ligand
/// @author Sharon Guffy

//Headers
#include <protocols/moves/Mover.hh>
#include <protocols/backrub/BackrubMover.fwd.hh>
#include <protocols/simple_moves/MinPackMover.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>




class HemeBindingMover: public protocols::moves::Mover{
public:
	//Default constructor
	HemeBindingMover();
	//Constructor with options
	HemeBindingMover(core::scoring::ScoreFunctionOP const & sfxn, protocols::backrub::BackrubMoverOP & brm,
		protocols::simple_moves::MinPackMoverOP & mpm, core::Size n_it = 1000, core::Real temperature = 0.6);
	//Constructor with different options
	HemeBindingMover(core::scoring::ScoreFunctionOP const & sfxn, core::pack::task::TaskFactoryOP const & tf, core::Size n_it, core::Real temperature = 0.6);
	//Copy constructor
	HemeBindingMover (HemeBindingMover const & other);
	//Destructor
	virtual ~HemeBindingMover();
	//Set private data
	void score_function(core::scoring::ScoreFunctionOP other_score_function);
	void backrub_mover(protocols::backrub::BackrubMoverOP br_mover);
	void minpack_mover(protocols::simple_moves::MinPackMoverOP mp_mover);
	void num_iterations(core::Size n_it);
	void temperature(core::Real temp);
	void task_factory(core::pack::task::TaskFactoryOP tf);
	//Get private data
	core::scoring::ScoreFunctionCOP score_function() const;
	protocols::backrub::BackrubMoverCOP backrub_mover() const;
	protocols::simple_moves::MinPackMoverCOP minpack_mover() const;
	core::Size num_iterations() const;
	core::Real temperature() const;
	core::pack::task::TaskFactoryCOP task_factory() const;
	//Virtual methods
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual void apply(core::pose::Pose & pose);
	virtual std::string get_name() const;
	//For RosettaScripts
	//Any other methods
private:
	//private data
	core::scoring::ScoreFunctionOP score_function_;
	//Store OPs to the other movers so we can set their properties directly
	protocols::backrub::BackrubMoverOP backrub_mover_;
	protocols::simple_moves::MinPackMoverOP minpack_mover_;
	core::Size num_iterations_;
	core::Real temperature_;
	core::pack::task::TaskFactoryOP tf_;
	protocols::simple_moves::PackRotamersMoverOP design_mover_;
};
//Define owning pointers
typedef utility::pointer::shared_ptr< HemeBindingMover > HemeBindingMoverOP;
typedef utility::pointer::shared_ptr< HemeBindingMover const > HemeBindingMoverCOP;
