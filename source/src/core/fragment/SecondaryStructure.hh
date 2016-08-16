// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief secondary structure will hold statistics about secondary structure predictions
/// sources can be from
///      - fragments
///      - psipred files ? other stuff
///
/// @details
///  from converting conformation_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange

#ifndef INCLUDED_core_fragment_SecondaryStructure_hh
#define INCLUDED_core_fragment_SecondaryStructure_hh

// Unit Headers
#include <core/fragment/SecondaryStructure.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
//#include <protocols/loops/LoopClass.hh>
#include <core/fragment/FragSet.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

namespace core {
namespace fragment {

/// @brief tiny helper class that knows the relative fractions of secondary structure  L,H,E
/// @detail
/// so far these fractions can be computed from a FragSet
/// other input strategies are conceivable but not implemented, yet: eg. psipred files, a bunch of poses,

class SecondaryStructure: public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~SecondaryStructure();

	SecondaryStructure() :
		total_residue_(0) {
	}

	/// @brief c'ctor that creates a SecondaryStructure as an average of several SecondaryStructure objects
	SecondaryStructure( utility::vector1<SecondaryStructureOP> &, utility::vector1<Real> &);

	/// @brief c'stor compute fractions from fragments
	SecondaryStructure(core::fragment::FragSet const& frags, core::Size nres =
		0, bool bJustCenterResidue = false) :
		total_residue_(nres) {
		compute_fractions(frags, bJustCenterResidue);
	}

	/// @brief c'stor compute fractions from pose ( well it won't be "fractions" )
	SecondaryStructure(core::pose::Pose const&);

	/// @brief return loop fraction at position
	core::Real loop_fraction(core::Size pos) const {
		runtime_assert( pos <= total_residue_ );
		return loop_fraction_(pos);
	}

	/// @brief return strand fraction at position
	core::Real strand_fraction(core::Size pos) const {
		runtime_assert( pos <= total_residue_ );
		return strand_fraction_(pos);
	}

	/// @brief alias for strand-fraction ...
	core::Real sheet_fraction( core::Size pos ) const {
		return strand_fraction( pos );
	}

	/// @brief helix fraction at position
	core::Real helix_fraction(core::Size pos) const {
		runtime_assert( pos <= total_residue_ );
		return std::max(0.0, 1.0 - loop_fraction_(pos) - strand_fraction_(pos));
	}

	/// @brief confidence at position
	core::Real confidence(core::Size pos) const {
		runtime_assert( pos <= total_residue_ );
		return confidence_(pos);
	}

	/// @brief sets secondary structure probabilities at a given position
	void set_fractions(core::Size pos, core::Real helix_fraction,
		core::Real sheet_fraction, core::Real loop_fraction, core::Real confidence=0) {
		runtime_assert( pos <= total_residue_ );

		core::Real s = helix_fraction + loop_fraction + sheet_fraction;
		if ( s != 1.0 ) {
			sheet_fraction = sheet_fraction / s;
			loop_fraction = loop_fraction / s;
		}

		loop_fraction_(pos) = (float)loop_fraction;
		strand_fraction_(pos) = (float) sheet_fraction;
		confidence_(pos)= (float) confidence;
	}

	/// @brief return loop fraction - FArray
	ObjexxFCL::FArray1D_float const& loop_fraction() const {
		return loop_fraction_;
	}

	/// @brief return strand fraction - FArray
	ObjexxFCL::FArray1D_float const& strand_fraction() const {
		return strand_fraction_;
	}

	/// @brief returns regions (in loop-class format) that belong to contiguous pieces of ss-structure
	//loops::Loops compute_ss_regions(
	// core::Real max_loop_frac = 0.3, core::Size min_length = 2
	//) const;

	/// @brief number of residues for which information is available
	core::Size total_residue() const {
		return total_residue_;
	}

	/// @brief returns the most probably secstruct at that position
	char secstruct(core::Size pos) const;

	/// @brief extends with pure 'L' at end until requested size is reached.
	void extend(core::Size);

	/// @brief read from file
	void read_from_file(std::string fn);

	/// @brief
	void show(std::ostream&) const;

	/// @brief write psipred format
	void write_psipred_ss2(std::ostream& os, std::string const& sequence) const;

	/// @brief write psipred format
	void read_psipred_ss2(std::istream& os);

	/// @brief write psipred format
	void read_psipred_ss2(std::string filename);

	/// @brief read talos+ format
	void read_talos_ss(std::istream& os);

	/// @brief read talos+ format
	void read_talos_ss(std::string filename);

private:
	void compute_fractions(core::fragment::FragSet const&,
		bool bJustCenterResidue);

	/// @brief store loop/strand fractions
	ObjexxFCL::FArray1D_float loop_fraction_;
	ObjexxFCL::FArray1D_float strand_fraction_;

	/// @brief store confidence values: used by talos, but not by psipred
	ObjexxFCL::FArray1D_float confidence_;

	/// @brief length of FArrays
	core::Size total_residue_;

};

/// @brief output operator
inline std::ostream & operator <<(std::ostream & os,
	SecondaryStructure const & t) {
	t.show(os);
	return os;
}

} //core
} //fragment

#endif

