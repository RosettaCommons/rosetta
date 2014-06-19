// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_McSugar.hh
/// @brief Markov chain sampler for sugar pucker.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_McSugar_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_McSugar_HH

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_McSugar.fwd.hh>

// Package headers
#include <core/chemical/rna/util.hh>
#include <protocols/rotamer_sampler/McRotamer.hh>
#include <protocols/rotamer_sampler/rna/RNA_SugarRotamer.fwd.hh>

using namespace core::chemical::rna;

namespace protocols {
namespace rotamer_sampler {
namespace rna {

class RNA_McSugar : public McRotamer {
public:

	RNA_McSugar(
		core::Size const rsd_id,
		core::Real const flip_rate = 0.1,
		PuckerState const init_pucker = NORTH
	);

	/// @brief Initialization
	void init();

	/// @brief Reset to current angle
	void reset() { stored_pucker_state_ = active_pucker_state_; }

	/// @brief Generate new active DOFs
	void operator++();

	/// @brief Update the stored DOFs
	void update() { active_pucker_state_ = stored_pucker_state_; }

	/// @brief Apply DOFs to pose
	void apply( core::pose::Pose & pose );

	/// @brief Set the flip rate of pucker
	void set_flip_rate( core::Real const setting ) {
		flip_rate_ = setting;
	}

	/// @brief Set if the sampler will skip pucker applying when input pose has
	//  same pucker assginment as sampler.
	void set_skip_same_pucker( bool const setting );

	/// @brief Set if using RNA_IdealCoord to sample puckers
	void set_idealize_coord( bool const setting );

	/// @brief Name of the class
	std::string get_name() const { return "RNA_McSugar"; }

	/// @brief Type of class (see enum in RotamerTypes.hh)
	virtual RotamerType type() const { return RNA_MC_SUGAR; }

private:
	PuckerState stored_pucker_state_, active_pucker_state_;
	core::Real flip_rate_;
	RNA_SugarRotamerOP sugar_rotamer_;
};

} //rna
} //rotamer_sampler
} //protocols

#endif
