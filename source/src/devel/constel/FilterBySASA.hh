// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief A class to filter out constellations based on SASA (Solvent
/// 	Accessible Surface Area).
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
	/// 	value	is relevant to filtering.
	static std::map< char, utility::vector1<std::string> > aa_sasa_atoms;

	/// @brief A table holding the SASA values of all atoms in the pose to which
	/// 	constellations belong.
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
