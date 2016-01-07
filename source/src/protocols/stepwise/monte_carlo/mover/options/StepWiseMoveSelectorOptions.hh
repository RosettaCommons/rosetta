// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/mover/options/StepWiseMoveSelectorOptions.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_mover_options_StepWiseMoveSelectorOptions_HH
#define INCLUDED_protocols_stepwise_monte_carlo_mover_options_StepWiseMoveSelectorOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/options/StepWiseMoveSelectorOptions.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

#if defined(WIN32) || defined(PYROSETTA)
#include <utility/tag/Tag.hh>
#endif

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {
namespace options {

class StepWiseMoveSelectorOptions: public virtual basic::resource_manager::ResourceOptions {

public:

	//constructor
	StepWiseMoveSelectorOptions();

	//destructor
	~StepWiseMoveSelectorOptions();

	StepWiseMoveSelectorOptions( StepWiseMoveSelectorOptions const & src );

public:


	/// @brief Initialize from the recursive "tag" structure.
	virtual
	void
	parse_my_tag( utility::tag::TagCOP ){}

	/// @brief The class name (its type) for a particular ResourceOptions instance.
	/// This function allows for better error message delivery.
	virtual
	std::string
	type() const{ return "StepWiseMoveSelectorOptions";}

	StepWiseMoveSelectorOptionsOP clone() const;

	void
	initialize_from_command_line();

	void set_allow_internal_hinge_moves( bool const & setting ){ allow_internal_hinge_moves_ = setting; }
	bool allow_internal_hinge_moves() const{ return allow_internal_hinge_moves_; }

	void set_allow_internal_local_moves( bool const & setting ){ allow_internal_local_moves_ = setting; }
	bool allow_internal_local_moves() const{ return allow_internal_local_moves_; }

	void set_add_delete_frequency( core::Real const & setting ){ add_delete_frequency_ = setting; }
	core::Real add_delete_frequency() const{ return add_delete_frequency_; }

	void set_from_scratch_frequency( core::Real const & setting ){ from_scratch_frequency_ = setting; }
	core::Real from_scratch_frequency() const{ return from_scratch_frequency_; }

	void set_docking_frequency( core::Real const & setting ){ docking_frequency_ = setting; }
	core::Real docking_frequency() const{ return docking_frequency_; }

	void set_submotif_frequency( core::Real const & setting ){ submotif_frequency_ = setting; }
	core::Real submotif_frequency() const{ return submotif_frequency_; }

	void set_switch_focus_frequency( core::Real const & setting ){ switch_focus_frequency_ = setting; }
	core::Real switch_focus_frequency() const{ return switch_focus_frequency_; }

	void set_skip_bulge_frequency( core::Real const & setting ){ skip_bulge_frequency_ = setting; }
	core::Real skip_bulge_frequency() const{ return skip_bulge_frequency_; }

	void set_vary_loop_length_frequency( core::Real const & setting ){ vary_loop_length_frequency_ = setting; }
	core::Real vary_loop_length_frequency() const{ return vary_loop_length_frequency_; }

	void set_allow_submotif_split( bool const & setting ){ allow_submotif_split_ = setting; }
	bool allow_submotif_split() const{ return allow_submotif_split_; }

	void set_force_submotif_without_intervening_bulge( bool const & setting ){ force_submotif_without_intervening_bulge_ = setting; }
	bool force_submotif_without_intervening_bulge() const{ return force_submotif_without_intervening_bulge_; }

	void set_filter_complex_cycles( bool const & setting ){ filter_complex_cycles_ = setting; }
	bool filter_complex_cycles() const{ return filter_complex_cycles_; }

private:

	bool allow_internal_hinge_moves_;
	bool allow_internal_local_moves_;
	core::Real add_delete_frequency_;
	core::Real from_scratch_frequency_;
	core::Real docking_frequency_;
	core::Real submotif_frequency_;
	core::Real switch_focus_frequency_;
	core::Real skip_bulge_frequency_;
	core::Real vary_loop_length_frequency_;
	bool filter_complex_cycles_;
	bool allow_submotif_split_;
	bool force_submotif_without_intervening_bulge_;

};

} //options
} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
