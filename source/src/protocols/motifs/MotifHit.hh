// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MotifHit.hh
/// @brief Class declaration for class that holds information about a motif in the context of the search
/// @author sthyme (sthyme@gmail.com)

#ifndef INCLUDED_protocols_motifs_MotifHit_hh
#define INCLUDED_protocols_motifs_MotifHit_hh

// Unit Headers
#include <protocols/motifs/MotifHit.fwd.hh>

// Package Headers
#include <protocols/motifs/Motif.fwd.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers

namespace protocols {
namespace motifs {

class MotifHit : public utility::pointer::ReferenceCount {

public:

	// Constructor that defaults to empty allowed_types (aka, all allowed)
	MotifHit(
		Motif const & motif,
		Size const & vbpos,
		bool const passed_automorphism
	);

	// Destructor, copy constructor, clone method
	virtual ~MotifHit();
	MotifHit( MotifHit const & src );
	MotifHitOP clone() const;

	// Set final_test_
	void final_test(
		core::Real const & final_test
	);

	// Set the best rotamer and conformer
	void build_rotamer(
		core::conformation::Residue const & build_rotamer
	);

	void target_conformer(
		core::conformation::Residue const & target_conformer
	);

	void dump_motif_hit();

	// Accessors
	MotifCOP const & motifcop() const { return motifcop_; }
	Size const & vbpos() const { return vbpos_; }
	bool const & passed_automorphism() const { return passed_automorphism_; }
	core::Real const & final_test() const { return final_test_; }
	core::conformation::ResidueOP const & build_rotamer() const { return build_rotamer_; }
	core::conformation::ResidueOP const & target_conformer() const { return target_conformer_; }

private:
	// The motif that this class is based around
	MotifCOP motifcop_;
	// The pose residue that is the one that the motif is targeting
	// ie, for DNA this would be the closest DNA base of the correct type
	Size vbpos_;
	// Whether the automorphic version of the jump is a better hit
	bool passed_automorphism_;
	core::Real final_test_;
	// The best rotamer and conformer
	core::conformation::ResidueOP build_rotamer_;
	core::conformation::ResidueOP target_conformer_;
};

} // namespace motifs
} // namespace protocols

#endif // INCLUDED_protocols_motifs_MotifHit
