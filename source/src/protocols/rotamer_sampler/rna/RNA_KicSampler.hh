// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_AnalyticLoopCloseSampler.hh
/// @brief Sample and torsions and close an RNA loop.
/// @author Fang-Chieh Chou

// NOTE: The sampler makes the assumption that the sugar pucker of
// moving_suite + 1 is moving (with pucker and chi rotamer being sampled).
// This corresponds to the "Append" case in Parin's SWA code.
// I think the "Prepend" case is not closable using the current
// RNA_KinematicClosure class. But if we stick with "Append" then it does not
// matter.

#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_KicSampler_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_KicSampler_HH

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_KicSampler.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerBase.hh>
#include <protocols/rotamer_sampler/RotamerSizedComb.fwd.hh>
#include <protocols/rotamer_sampler/rna/RNA_KinematicCloser.fwd.hh>
#include <protocols/rotamer_sampler/rna/RNA_ChiRotamer.fwd.hh>
#include <protocols/rotamer_sampler/screener/RNA_TorsionScreener.fwd.hh>

namespace protocols {
namespace rotamer_sampler {
namespace rna {

class RNA_KicSampler : public RotamerBase {
public:
	RNA_KicSampler(
		core::pose::PoseOP const & ref_pose,
		core::Size const moving_suite,
		core::Size const chainbreak_suite
	);

	~RNA_KicSampler();

	/// @brief Initialization
	void init();

	/// @brief Reset to the first (or random if is_random()) rotamer.
	void reset();

	/// @brief Move to next rotamer
	void operator++();

	/// @brief Check if reach the end of rotamer list
	bool not_end() const;

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & pose );

	/// @brief Apply the current rotamer to ref_pose_
	void apply();

	/// @brief Set the random sampling state
	void set_random( bool const setting );

	/// @brief If the chain is closable (random sampling only)
	bool closable() const { return random_chain_closed_; }

	/// @brief Name of the class
	std::string get_name() const { return "RNA_KicSampler"; }

	// Set functions
	void set_verbose( bool const setting ) {
		set_and_reinit( verbose_, setting );
	}

	void set_extra_epsilon( bool const setting ) {
		set_and_reinit( extra_epsilon_, setting );
	}

	void set_extra_chi( bool const setting ) {
		set_and_reinit( extra_chi_, setting );
	}

	void set_skip_same_pucker( bool const setting ) {
		set_and_reinit( skip_same_pucker_, setting );
	}

	void set_idealize_coord( bool const setting ) {
		set_and_reinit( idealize_coord_, setting );
	}

	void set_torsion_screen( bool const setting ) {
		set_and_reinit( torsion_screen_, setting );
	}

	void set_pucker_state( core::Size const setting ) {
		set_and_reinit( pucker_state_, setting );
	}

	void set_base_state( core::Size const setting ) {
		set_and_reinit( base_state_, setting );
	}

	void set_bin_size( core::Real const setting ) {
		set_and_reinit( bin_size_, setting );
	}

	void set_fast( core::Real const setting ) {
		if ( setting ){
			extra_chi_ = false;
			extra_epsilon_ = false;
			set_and_reinit( bin_size_, 40.0 /*setting*/ );
		}
	}

	/// @brief Max # of step for trying in random samping
	void set_max_tries( core::Size const setting ) {
		set_and_reinit( max_tries_, setting );
	}

private:
	void get_next_valid_bb();

	//Disable copy constructor and assignment
  RNA_KicSampler( const RNA_KicSampler & );
  void operator=( const RNA_KicSampler & );

	core::pose::PoseOP const ref_pose_;
	core::Size const moving_suite_, chainbreak_suite_;
	core::Size pucker_state_, base_state_;
	core::Real bin_size_;
	core::Size max_tries_;
	bool verbose_, extra_epsilon_, extra_chi_, skip_same_pucker_,
			idealize_coord_, torsion_screen_, random_chain_closed_;

	RotamerSizedCombOP bb_rotamer_;
	RNA_KinematicCloserOP loop_closer_;
	RNA_ChiRotamerOP chi_rotamer_;
	screener::RNA_TorsionScreenerOP screener_;
};

}
}
}

#endif
