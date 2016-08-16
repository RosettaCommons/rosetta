// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief A class to filter out constellations based on SASA (Solvent
///  Accessible Surface Area).
/// @author Andrea Bazzoli

#ifndef INCLUDED_FilterBySASA_hh
#define INCLUDED_FilterBySASA_hh


#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>
#include <map>
#include <string>

using core::Real;
using core::Size;
using core::pose::Pose;


namespace devel {
namespace constel {

class FilterBySASA {

private:

	/// @brief A table listing, for each amino acid type, the atoms whose SASA
	///  value is relevant to filtering.
	static std::map< char, utility::vector1<std::string> > aa_sasa_atoms;

	/// @brief A table holding the SASA values of all atoms in the pose to which
	///  constellations belong.
	static core::id::AtomID_Map<Real> atom_sasa;

	/// @brief Maximum allowed SASA value for a constellation atom.
	static Real MAX_ATOM_SASA;

public:

	static void init( Real const smax, Pose const& ps );

	/// @brief Tells whether a constellation has a sufficiently low per-atom SASA
	static bool has_low_per_atom_sasa( Pose const& ps,
		utility::vector1<Size> const& cnl );
};

} // constel
} // devel

#endif
