// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_McSuite.hh
/// @brief Markov chain sampler for RNA suite.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_McSuite_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_McSuite_HH

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_McSuite.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/McComb.hh>
#include <protocols/rotamer_sampler/McOneTorsion.fwd.hh>
#include <protocols/rotamer_sampler/rna/RNA_McSugar.fwd.hh>

// Project headers
#include <core/id/TorsionID.hh>

namespace protocols {
namespace rotamer_sampler {
namespace rna {

class RNA_McSuite : public McComb {
public:

	RNA_McSuite( core::Size const rsd_id );

	/// @brief Initialization
	void init();

	/// @brief Clear all rotamer samplers stored in this sampler
	void clear_rotamer();

	/// @brief Set the standard deviation of Gaussian sampler
	void set_gaussian_stdev( core::Real const setting );

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
		utility::vector1<core::Real> const & init_torsions,
		core::Size const init_pucker
	) {
		init_torsions_ = init_torsions;
		init_pucker_ = init_pucker;
		set_init( false );
	}

	/// @brief Set starting torsions and pucker states from pose
	void set_init_from_pose( core::pose::Pose const & pose );

	/// @brief Name of the class
	std::string get_name() const { return "RNA_McSuite"; }

	/// @brief Type of class (see enum in RotamerTypes.hh)
	virtual RotamerType type() const { return RNA_MC_SUITE; }

private:
	core::Size const rsd_id_;
	bool skip_same_pucker_, idealize_coord_, sample_near_a_form_, sample_bb_,
			 sample_lower_nucleoside_, sample_upper_nucleoside_;
	core::Real pucker_flip_rate_, gaussian_stdev_;
	core::Real const a_form_range_;
	core::Size init_pucker_;
	utility::vector1<core::Real> a_form_torsions_, init_torsions_;
	utility::vector1<core::id::TorsionID> torsion_ids_;
	utility::vector1<McOneTorsionOP> bb_samplers_;
	utility::vector1<McOneTorsionOP> chi_samplers_;
	utility::vector1<RNA_McSugarOP> sugar_samplers_;
};

}
}
}

#endif
