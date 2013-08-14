// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @detailed
///	  Contains currently: Relax Baseclass
///
///
/// @author Mike Tyka


#ifndef INCLUDED_protocols_relax_RelaxProtocolBase_hh
#define INCLUDED_protocols_relax_RelaxProtocolBase_hh

// Unit headers
#include <protocols/relax/RelaxProtocolBase.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

//// C++ headers
#include <string>

#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace relax {

class RelaxProtocolBase : public moves::Mover {
public:
	typedef moves::Mover parent;

	RelaxProtocolBase( core::scoring::ScoreFunctionOP );
	RelaxProtocolBase( std::string const & movername = "RelaxProtocol" );
	RelaxProtocolBase( std::string const & movername, core::scoring::ScoreFunctionOP );
	~RelaxProtocolBase();

	virtual	protocols::moves::MoverOP	fresh_instance() const { return clone(); };

	static void register_options();

  void apply_disulfides( core::pose::Pose & pose );

	// Default options -------------------------------------
	void set_defaults();
	void set_default_minimization_settings();
	void set_default_coordinate_settings();
	void set_default_movemap();

	// Public accessors
	core::kinematics::MoveMapOP get_movemap() { return movemap_; }
	const core::scoring::ScoreFunctionCOP get_scorefxn() const { return scorefxn_; }
	core::pack::task::TaskFactoryOP const & get_task_factory() const { return task_factory_; }
	
	bool cartesian() const { return cartesian_; }
	std::string min_type() const { return min_type_; }
	Size max_iter() const { return max_iter_; }
	bool dry_run() const { return dry_run_; }

	bool constrain_relax_to_native_coords() const { return constrain_relax_to_native_coords_; }
	bool constrain_relax_to_start_coords() const { return constrain_relax_to_start_coords_; }
	bool constrain_coords() const { return constrain_coords_; }
	bool explicit_ramp_constraints() const { return explicit_ramp_constraints_; }
	bool ramp_down_constraints() const { return ramp_down_constraints_; }
	bool constrain_relax_segments() const { return constrain_relax_segments_; }
	
	// Public mutators
	void set_movemap( core::kinematics::MoveMapOP movemap ) { movemap_ = movemap; }
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn ) { scorefxn_ = scorefxn; }
	void set_task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }

	void cartesian( bool newval ) { cartesian_ = newval; }
	void min_type( std::string min_type ) { min_type_ = min_type; }
	void max_iter( Size max_iter ) { max_iter_ = max_iter; }
	void dry_run( bool setting ) { dry_run_ = setting; }

	void constrain_relax_to_native_coords( bool constrain_relax_to_native_coords ) { constrain_relax_to_native_coords_ = constrain_relax_to_native_coords; }
	void constrain_relax_to_start_coords(  bool constrain_relax_to_start_coords ) { constrain_relax_to_start_coords_ = constrain_relax_to_start_coords; }
	void constrain_coords( bool constrain_coords ) { constrain_coords_ = constrain_coords; }
	void constrain_relax_segments( bool constrain_relax_segments ) { constrain_relax_segments_ = constrain_relax_segments; }
	void ramp_down_constraints( bool ramp_down_constraints ) {
		explicit_ramp_constraints_ = true;
		ramp_down_constraints_ =  ramp_down_constraints;
	}


protected:
	core::scoring::ScoreFunctionOP get_scorefxn() { return scorefxn_; }

	// Accessors -------------------------------------
	bool fix_omega() const { return fix_omega_; }
	bool minimize_bond_lengths() const { return minimize_bond_lengths_; }
	bool minimize_bond_angles() const { return minimize_bond_angles_; }
	int minimize_bondlength_subset() const { return minimize_bondlength_subset_; }
	int minimize_bondangle_subset() const { return minimize_bondangle_subset_; }
	bool limit_aroma_chi2() const { return limit_aroma_chi2_; }


	// Mutators -------------------------------------
	void fix_omega( bool fix_omega ) { fix_omega_ = fix_omega; }
  void minimize_bond_lengths( bool minimize_bond_lengths ) { minimize_bond_lengths_ = minimize_bond_lengths; }
  void minimize_bond_angles( bool minimize_bond_angles ) { minimize_bond_angles_ = minimize_bond_angles; }
  void minimize_bondangle_subset( int minimize_bondangle_subset ) { minimize_bondangle_subset_ = minimize_bondangle_subset; }
  void minimize_bondlength_subset( int minimize_bondlength_subset ) { minimize_bondlength_subset_ = minimize_bondlength_subset; }
	

	void initialize_movemap( core::pose::Pose const & pose, core::kinematics::MoveMap & movemap );
	void set_up_constraints( core::pose::Pose &pose,  core::kinematics::MoveMap & local_movemap );
	void output_debug_structure( core::pose::Pose & pose, std::string prefix );

private:
  // Essentially MoveMap settings
	bool fix_omega_;
	bool minimize_bond_lengths_;
	bool minimize_bond_angles_;
	int minimize_bondangle_subset_;
	int minimize_bondlength_subset_;

  // Constraint settings
	bool constrain_relax_to_native_coords_;
	bool constrain_relax_to_start_coords_;
	bool constrain_coords_;
	bool coord_constrain_sidechains_;
	bool explicit_ramp_constraints_;
	bool ramp_down_constraints_;
	bool constrain_relax_segments_;

	// The minimizer algorithm
  std::string min_type_;

  /// Do cartesian-space minimization?
	bool cartesian_;

	/// maximum minimizer iterations
	core::Size max_iter_;

	bool limit_aroma_chi2_;

	// The movemap
	core::kinematics::MoveMapOP movemap_;

	// Fullatom scoring function used
	core::scoring::ScoreFunctionOP scorefxn_;

	/// task factory for packer ( used only FastRelax as of 11/05/2010 )
	core::pack::task::TaskFactoryOP task_factory_;

	/// @brief Is this a dry run (i.e. do no cycles)
	bool dry_run_;
};

}
} // protocols

#endif
