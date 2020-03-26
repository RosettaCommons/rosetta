// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/MC_Comb.hh
/// @brief Ensemble of Markov chain samplers for modeler combinatorially.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_sampler_MC_Comb_HH
#define INCLUDED_protocols_sampler_MC_Comb_HH

// Unit headers
#include <protocols/recces/sampler/MC_Comb.fwd.hh>

// Package headers
#include <protocols/recces/sampler/MC_Sampler.hh>

namespace protocols {
namespace recces {
namespace sampler {

class MC_Comb : public MC_Sampler {
public:
	MC_Comb();

	~MC_Comb() override{}

	/// @brief Initialization
	void init() override;

	/// @brief Reset to the first (or random if random()) rotamer
	void reset() override;

	/// @brief Move to next rotamer
	void operator++() override;

	/// @brief Update the DOFs
	void update() override;

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & pose ) override;

	/// @brief Get uniform modeler state
	void set_uniform_modeler( bool const setting ) override;

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_rotamer( MC_SamplerOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	toolbox::SamplerPlusPlusType type() const override { return toolbox::MC_COMB; }

	/// @brief set (optional) update_pose
	void
	set_update_pose( core::pose::PoseCOP setting ) override;

	/// @brief output summary of class
	void show( std::ostream & out, core::Size const indent = 0 ) const override;

	/// @brief return OP to the subsampler that controls exactly this torsion_id (assume only one).
	MC_SamplerOP
	find( core::id::TorsionID const & torsion_id ) override;

	core::Size
	num_rotamers() const { return rotamer_list_.size(); }

	MC_SamplerOP
	operator[]( core::Size const id ) { return rotamer_list_[ id ]; }

protected:
	utility::vector1<MC_SamplerOP> rotamer_list_;
};

} //sampler
} //recces
} //protocols

#endif

