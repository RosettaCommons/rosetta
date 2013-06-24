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
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>

//// C++ headers
#include <string>

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace relax {

class RelaxProtocolBase : public moves::Mover {
public:
	typedef moves::Mover parent;

public:

	RelaxProtocolBase( core::scoring::ScoreFunctionOP );
	RelaxProtocolBase( std::string const & movername = "RelaxProtocol" );
	RelaxProtocolBase( std::string const & movername, core::scoring::ScoreFunctionOP );
	~RelaxProtocolBase();

	virtual	protocols::moves::MoverOP	fresh_instance() const { return clone(); };

	static void register_options();

	// Default options -------------------------------------
	void set_defaults();
	void set_default_minimization_settings();
	void set_default_coordinate_settings();
	void set_default_movemap();

protected:

	void initialize_movemap( core::pose::Pose const & pose, core::kinematics::MoveMap & movemap );

public:
	// Accesors -------------------------------------

	void fix_omega( bool setting );
	void minimize_bond_lengths( bool setting );
	void minimize_bond_angles( bool setting );
	void minimize_bondlength_subset( int setting );
	void minimize_bondangle_subset( int setting );

	bool fix_omega() const;
	bool minimize_bond_lengths() const;
	bool minimize_bond_angles() const;
	int minimize_bondlength_subset() const;
	int minimize_bondangle_subset() const;

	bool constrain_relax_to_native_coords() const {  return    constrain_relax_to_native_coords_;}
	bool constrain_relax_to_start_coords() const {  return     constrain_relax_to_start_coords_;}
	bool constrain_coords() const {  return                    constrain_coords_;}
	bool explicit_ramp_constraints() const {  return           explicit_ramp_constraints_;}
	bool ramp_down_constraints() const {  return               ramp_down_constraints_;}
	bool constrain_relax_segments() const {  return            constrain_relax_segments_;}

	bool limit_aroma_chi2() const { return limit_aroma_chi2_; }

	void constrain_relax_to_native_coords( bool constrain_relax_to_native_coords ) {
		constrain_relax_to_native_coords_ = constrain_relax_to_native_coords;
	}
	void constrain_relax_to_start_coords(  bool constrain_relax_to_start_coords ) {
		constrain_relax_to_start_coords_ = constrain_relax_to_start_coords;
	}
	void constrain_coords( bool constrain_coords ) {
		constrain_coords_ = constrain_coords;
	}
	void ramp_down_constraints( bool ramp_down_constraints ) {
		explicit_ramp_constraints_ = true;
		ramp_down_constraints_ =  ramp_down_constraints;
	}
	void constrain_relax_segments( bool constrain_relax_segments ) {
		constrain_relax_segments_ = constrain_relax_segments;
	}

	// get/set cartesian minimization option
	void cartesian( bool newval ) { cartesian_ = newval; }
	bool cartesian( ) const { return cartesian_; }

	core::kinematics::MoveMapOP get_movemap();
	void set_movemap( core::kinematics::MoveMapOP movemap );

	void set_min_type( std::string min_type );
	void set_max_iter( Size max_iter );

public:
	void set_scorefxn( core::scoring::ScoreFunctionOP score );
	const core::scoring::ScoreFunctionCOP get_scorefxn() const;

	void set_task_factory( core::pack::task::TaskFactoryOP taskf );
	core::pack::task::TaskFactoryOP const & get_task_factory() const;

	void set_dry_run( bool setting ) {
		dry_run_ = setting;
	}

	bool dry_run() const {
		return dry_run_;
	}

protected:
	core::scoring::ScoreFunctionOP get_scorefxn();



public:
  void apply_disulfides( core::pose::Pose & pose );

protected:
	void set_up_constraints( core::pose::Pose &pose,  core::kinematics::MoveMap & local_movemap );

	void output_debug_structure( core::pose::Pose & pose, std::string prefix );





protected:  // Essentially MoveMap settings

	bool fix_omega_;
	bool minimize_bond_lengths_;
	bool minimize_bond_angles_;
	int minimize_bondangle_subset_;
	int minimize_bondlength_subset_;

protected:  // Constraint settings

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

protected:

	bool limit_aroma_chi2_;

private: // other private variables

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
