// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/rna/MC_RNA_Suite.hh
/// @brief Markov chain sampler for RNA suite.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_sampler_rna_MC_RNA_Suite_HH
#define INCLUDED_protocols_sampler_rna_MC_RNA_Suite_HH

// Unit headers
#include <protocols/recces/sampler/rna/MC_RNA_Suite.fwd.hh>

// Package headers
#include <protocols/recces/sampler/MC_Comb.hh>
#include <protocols/recces/sampler/MC_OneTorsion.fwd.hh>
#include <protocols/recces/sampler/rna/MC_RNA_Sugar.fwd.hh>

// Project headers
#include <core/chemical/rna/util.hh>
#include <core/id/TorsionID.hh>


#ifdef WIN32
#include <protocols/recces/sampler/MC_OneTorsion.hh>
#include <protocols/recces/sampler/rna/MC_RNA_Sugar.hh>
#endif

// using namespace core::chemical::rna;
// Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
// by our coding convention due to problems it create on modern compilers and because of the name clashing.
// For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using

namespace protocols {
namespace recces {
namespace sampler {
namespace rna {

class MC_RNA_Suite : public MC_Comb {
public:

	MC_RNA_Suite( core::Size const rsd_id );

	/// @brief Initialization
	void init();

	/// @brief Clear all rotamer samplers stored in this sampler
	void clear_rotamer();

	/// @brief Set the standard deviation of Gaussian sampler
	void set_gaussian_stdev( core::Real const setting );

	/// @brief Set the angle range from the initial torsions
	void set_angle_range_from_init( core::Real const setting );

	/// @brief Set the flip rate of pucker
	void set_pucker_flip_rate( core::Real const setting );

	/// @brief Set if the sampler will skip pucker applying when input pose has
	//  same pucker assginment as sampler.
	void set_skip_same_pucker( bool const setting ) {
		set_and_reinit( skip_same_pucker_, setting );
	}

	/// @brief Set if using RNA_IdealCoord to sample puckers
	void set_idealize_coord( bool const setting ) {
		set_and_reinit( idealize_coord_, setting );
	}

	/// @brief Set angle ranges for A form rotamers
	void set_a_form_range( core::Real const setting ) {
		set_and_reinit( a_form_range_, setting );
	}

	/// @brief Set if only sample near A form rotamers
	void set_sample_near_a_form( bool const setting ) {
		set_and_reinit( sample_near_a_form_, setting );
	}

	/// @brief Set if sample backbone
	void set_sample_bb( bool const setting ) {
		set_and_reinit( sample_bb_, setting );
	}

	/// @brief Set if sample nucleoside
	void set_sample_lower_nucleoside( bool const setting ) {
		set_and_reinit( sample_lower_nucleoside_, setting );
	}

	/// @brief Set if sample nucleoside
	void set_sample_upper_nucleoside( bool const setting ) {
		set_and_reinit( sample_upper_nucleoside_, setting );
	}

	/// @brief Set starting torsions and pucker states
	void set_init_states(
		utility::vector1<core::Real> const & init_torsions, core::chemical::rna::PuckerState const init_pucker
	) {
		init_torsions_ = init_torsions;
		init_pucker_ = init_pucker;
		set_init( false );
	}

	/// @brief Set starting torsions and pucker states from pose
	void set_init_from_pose( core::pose::Pose const & pose );

	/// @brief Set the stored angle from a pose
	void set_angle( core::pose::Pose const & pose );

	/// @brief Name of the class
	std::string get_name() const { return "MC_RNA_Suite"; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::MC_RNA_SUITE; }

private:
	core::Size const rsd_id_;
	bool skip_same_pucker_, idealize_coord_, sample_near_a_form_, sample_bb_,
		sample_lower_nucleoside_, sample_upper_nucleoside_;
	core::Real pucker_flip_rate_, gaussian_stdev_, angle_range_, a_form_range_;
	core::chemical::rna::PuckerState init_pucker_;
	utility::vector1<core::Real> a_form_torsions_, init_torsions_;
	utility::vector1<core::id::TorsionID> torsion_ids_;
	utility::vector1<MC_OneTorsionOP> bb_samplers_;
	utility::vector1<MC_OneTorsionOP> chi_samplers_;
	utility::vector1<MC_RNA_SugarOP> sugar_samplers_;
};

} //rna
} //sampler
} //recces
} //protocols

#endif
