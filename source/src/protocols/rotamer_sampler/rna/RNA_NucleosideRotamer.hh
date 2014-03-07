// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_NucleosideRotamer.hh
/// @brief Generate rotamers for one RNA nucleoside (pucker + glycosidic chi).
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_NucleosideRotamer_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_NucleosideRotamer_HH

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_NucleosideRotamer.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerSizedAny.hh>

// Project headers
#include <core/chemical/rna/RNA_FittedTorsionInfo.fwd.hh>
#include <core/chemical/rna/RNA_Util.hh>

using namespace core::chemical::rna;

namespace protocols {
namespace rotamer_sampler {
namespace rna {

class RNA_NucleosideRotamer : public RotamerSizedAny {
public:
	RNA_NucleosideRotamer(
		core::Size const rsd_id,
		PuckerState const pucker_state, //ANY_PUCKER, NORTH, SOUTH, NO_PUCKER
		ChiState const base_state //ANY_CHI, ANTI, SYN, NO_CHI
	);

	/// @brief Initialization wrapper
	void init();

	/// Set functions
	void set_extra_chi( bool const setting ) {
		set_and_reinit( extra_chi_, setting );
	}

	void set_skip_same_pucker( bool const setting ) {
		set_and_reinit( skip_same_pucker_, setting );
	}

	void set_idealize_coord( bool const setting ) {
		set_and_reinit( idealize_coord_, setting );
	}

	/// @brief Just sample the center chi torsion.
	void set_fast( bool const setting ) {
		set_and_reinit( fast_, setting );
	}

	void set_bin_size( core::Real const setting ) {
		set_and_reinit( bin_size_, setting );
	}

	/// @brief Name of the class
	std::string get_name() const { return "RNA_NucleosideRotamer"; }

	/// @brief Type of class (see enum in RotamerTypes.hh)
	virtual RotamerType type() const { return RNA_NUCLEOSIDE; }

private:
	core::Size const rsd_id_;
	ChiState base_state_;

	bool extra_chi_, skip_same_pucker_, idealize_coord_, fast_;

	core::Real bin_size_;

	utility::vector1<PuckerState> pucker_states_;
};

}
}
}

#endif
