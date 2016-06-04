// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/DesignRepackMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_simple_moves_DesignRepackMover_hh
#define INCLUDED_protocols_simple_moves_DesignRepackMover_hh

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>

// Unit headers
#include <protocols/simple_moves/DesignRepackMover.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @brief a pure virtual base class for movers which redesign and repack the interface
class DesignRepackMover : public protocols::moves::Mover
{
public:
	DesignRepackMover();
	DesignRepackMover( std::string const & name );
	void setup_packer_and_movemap( core::pose::Pose const & pose );
	protocols::moves::MoverOP clone() const = 0; // this is a pure virtual class that cannot be instantiated
	protocols::moves::MoverOP fresh_instance() const = 0;
	virtual void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	// managing minimization options
	void min_rb( utility::vector1< bool > const & min_rb ) { min_rb_ = min_rb; min_rb_set_ = true;}
	/// @brief in most cases, there would only be one rb dof making it useful to have a non-vector accessor to min_rb_. However, if the pose has multiple jumps, setting min_rb_ in this way might cause trouble in other parts of the code.
	void min_rb( bool const min_rb ) { min_rb_.clear(); min_rb_.push_back( min_rb ); min_rb_set_ = true; }
	utility::vector1< bool > const & min_rb() const { runtime_assert( min_rb_set() ); return min_rb_; }
	bool min_rb_set() const { return min_rb_set_; }
	void min_sc( utility::vector1< bool > const & min_sc ) { min_sc_ = min_sc; min_sc_set_ = true; }
	utility::vector1< bool > const & min_sc() const { runtime_assert( min_sc_set() ); return min_sc_; }
	bool min_sc_set() const { return min_sc_set_; }
	void min_bb( utility::vector1< bool > const & min_bb ) { min_bb_ = min_bb; min_bb_set_ = true; }
	utility::vector1< bool > const & min_bb() const { runtime_assert( min_bb_set() ); return min_bb_; }
	bool min_bb_set() const { return min_bb_set_; }
	bool optimize_foldtree() const { return optimize_foldtree_; }
	void optimize_foldtree( bool const opt ) { optimize_foldtree_ = opt; }
	/// @brief a dummy apply so that instantiation of this baseclass would be possible.
	virtual void apply( core::pose::Pose & ) {}
	virtual std::string get_name() const;
	void prevent_repacking( utility::vector1< core::Size > const &  p ) { prevent_repacking_ = p; }
	utility::vector1< core::Size > const & prevent_repacking() const { return( prevent_repacking_ ); }
	void restrict_to_repacking( utility::vector1< core::Size > const & p ) { restrict_to_repacking_ = p; }
	utility::vector1< core::Size > const & restrict_to_repacking() const { return( restrict_to_repacking_ ); }
	void design( bool const des ) { design_partner1_ = des; design_partner2_ = des; }
	bool design() const { return( design_partner1_ || design_partner2_ ); }
	void set_scorefxn_repack( core::scoring::ScoreFunctionCOP scorefxn );
	void set_scorefxn_minimize( core::scoring::ScoreFunctionCOP scorefxn );
	core::scoring::ScoreFunctionOP scorefxn_repack() const;
	core::scoring::ScoreFunctionOP scorefxn_minimize() const;
	core::pack::task::PackerTaskCOP task() const;
	core::pack::task::PackerTaskOP & task();
	/// @brief after fiddling with a task from outside this mover, clear it, or else, on the next iteration through
	/// the mover the changes will be remembered
	void clear_task();
	void clear_task_factory();
	void use_preset_task( bool const bt ) { use_preset_task_ = bt; }
	bool use_preset_task() const { return use_preset_task_; }
	void task_factory( core::pack::task::TaskFactoryOP p );
	core::pack::task::TaskFactoryOP & task_factory();
	core::pack::task::TaskFactoryOP task_factory() const;
	virtual ~DesignRepackMover();

protected:
	core::scoring::ScoreFunctionOP scorefxn_repack_;
	core::scoring::ScoreFunctionOP scorefxn_minimize_;
	bool repack_partner1_, repack_partner2_;
	bool design_partner1_, design_partner2_; // design or only repack?
	// curr_ variables are reset with every apply, whereas the min_ variables do not
	utility::vector1< bool > min_sc_, curr_min_sc_, min_rb_;
	utility::vector1< bool > min_bb_, curr_min_bb_, curr_min_rb_;
	bool min_rb_set_, min_sc_set_, min_bb_set_;
	utility::vector1< core::Size > target_residues_;
	core::Real interface_distance_cutoff_;
	bool repack_non_ala_; // do we redesign nonalanine positions?
	bool optimize_foldtree_; // do we want to optimize or keep the input fold tree or optimize it?
	bool automatic_repacking_definition_;
	core::pack::task::PackerTaskOP task_;
	bool use_preset_task_;// use as a basis for computing the packer task a task that was precomputed.
	utility::vector1< bool > allowed_aas_; //the default allowed aas in designable positions
	// The following parameters override the automatically decided parameters.
	utility::vector1< core::Size > prevent_repacking_; //residues that should not be repacked
	utility::vector1< core::Size > restrict_to_repacking_; //residues that should not be designed
	core::pack::task::TaskFactoryOP task_factory_; // sequence positions and residue-level tasks to apply when setup_packer_task is done
	bool symmetry_;
};

} // simple_moves
} // protocols

#endif /*INCLUDED_protocols_protein_interface_design_movers_DesignRepackMover_HH*/
