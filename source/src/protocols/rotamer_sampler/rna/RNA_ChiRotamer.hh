// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_ChiRotamer.hh
/// @brief Generate glycosidic chi rotamers for RNA.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_ChiRotamer_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_ChiRotamer_HH

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_ChiRotamer.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerOneTorsion.hh>

// Project headers
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>

namespace protocols {
namespace rotamer_sampler {
namespace rna {

class RNA_ChiRotamer : public RotamerOneTorsion {
public:
	RNA_ChiRotamer(
		core::Size const rsd_id,
		core::Size const pucker_state,
		core::Size const base_state
	);

	/// @brief Initialization
	void init();

	/// @brief Set the residue id
	void set_rsd_id( core::Size const setting ) {
		set_and_reinit( rsd_id_, setting );
	}

	/// @brief Get the pucker state (NORTH / SOUTH)
	core::Size pucker_state() { return pucker_state_; }

	/// @brief Set the pucker state (NORTH / SOUTH)
	void set_pucker_state( core::Size const setting ) {
		set_and_reinit( pucker_state_, setting );
	}

	/// @brief Get the base state (WHATEVER / ANTI / SYN)
	core::Size base_state() { return base_state_; }

	/// @brief Set the base state (WHATEVER / ANTI / SYN)
	void set_base_state( core::Size const setting ) {
		set_and_reinit( base_state_, setting );
	}

	/// @brief Set the bin_size (default: 20)
	void set_bin_size( core::Real const setting ) {
		set_and_reinit( bin_size_, setting );
	}

	/// @brief Set the max_range of sampling (default: +-20)
	void set_max_range( core::Real const setting ) {
		set_and_reinit( max_range_, setting );
	}

	/// @brief Set use extra chi (+-60)
	void set_extra_chi( bool const setting ) {
		if ( setting ) {
			set_max_range( 60 );
		} else {
			set_max_range( 20 );
		}
	}

	/// @brief Name of the class
	std::string get_name() const { return "RNA_ChiRotamer"; }

	/// @brief Type of class (see enum in RotamerTypes.hh)
	virtual RotamerType type() const { return RNA_CHI; }

private:
	core::Size rsd_id_, base_state_, pucker_state_;
	core::Real bin_size_, max_range_;
	bool extra_chi_;

	core::chemical::rna::RNA_FittedTorsionInfo const torsion_info_;
};

}
}
}

#endif
