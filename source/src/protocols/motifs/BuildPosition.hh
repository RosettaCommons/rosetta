// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BuildPosition.hh
/// @brief Class declaration for position where motifs are built in MotifSearch
/// @author sthyme (sthyme@gmail.com)

#ifndef INCLUDED_protocols_motifs_BuildPosition_hh
#define INCLUDED_protocols_motifs_BuildPosition_hh

// Unit Headers
#include <protocols/motifs/BuildPosition.fwd.hh>

// Package Headers
#include <protocols/motifs/Motif.fwd.hh>
#include <protocols/motifs/MotifHit.fwd.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <set>

//Auto Headers
#ifdef WIN32
#include <protocols/motifs/MotifHit.hh>
#endif


namespace protocols {
namespace motifs {

class BuildPosition : public utility::pointer::ReferenceCount {

public:

	// Constructor that defaults to empty allowed_types (aka, all allowed)
	BuildPosition(
		Size const seqpos,
		utility::vector1< Size > const & target_positions
	);

	// Constructor that requires input of allowed_types for limiting aas
	BuildPosition(
		Size const seqpos,
		utility::vector1< Size > const & target_positions,
		std::set< std::string > const & allowed_types
	);

	// Destructor, copy constructor, clone method
	virtual ~BuildPosition();
	BuildPosition( BuildPosition const & src );
	BuildPositionOP clone() const;

	// Save this rotamer
	void keep_rotamer(
		core::conformation::Residue const & res
	);

	// Save this motif
	void keep_motif(
		Motif const & motif
	);

	// Save this motif hit
	void keep_motifhit(
		MotifHit const & motifhit
	);

	// Clears the best_motifs_ and best_rotamers_
	void clear_data();

	// Clears the rots only
	void clear_rots();

	// Accessors
	Size const & seqpos() const { return seqpos_; }
	utility::vector1< Size > const & target_positions() const { return target_positions_; }
	std::set< std::string > const & allowed_types() const { return allowed_types_; }
	core::pack::rotamer_set::Rotamers const & best_rotamers() const { return best_rotamers_; }
	MotifCOPs const & best_motifs() const { return best_motifs_; }
	MotifHitCOPs const & best_motifhits() const { return best_motifhits_; }

private:

	// The particular position associated with motifs and rotamers
	Size seqpos_;
	// Close potential targets for motifs (DNA bases in protein-DNA case)
	utility::vector1< Size > target_positions_;
	// Residues types (name3 string) allowed at this position
	std::set< std::string > allowed_types_;
	// Rotamers that pass cutoffs in the MotifSearch
	core::pack::rotamer_set::Rotamers best_rotamers_;
	// Motifs that pass cutoffs in the MotifSearch
	MotifCOPs best_motifs_;
	MotifHitCOPs best_motifhits_;

};

} // namespace motifs
} // namespace protocols

#endif // INCLUDED_protocols_motifs_BuildPosition
