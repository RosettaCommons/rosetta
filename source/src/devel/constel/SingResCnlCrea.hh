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

/// @brief Definition of a class to create single-residue constellations
/// @author jk
/// @author Andrea Bazzoli (bazzoli@ku.edu)
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>
#include <map>
#include <string>

namespace devel {
namespace constel {

using core::pose::Pose;

class SingResCnlCrea {

	/// @brief True if constellations must be output to file in stripped form
	static bool stripped_;

	/// @brief mutation-specific lists of constellation atoms whose occupancy must
	/// 	be zeroed 	when the constellation is output to file
	typedef std::pair<char, char> MutTyp;
	static std::map<MutTyp, std::string> zeroed_occ_;

	/// @brief sets occupancy to 0 for a residue's backbone atoms and hydrogen
	/// 	atoms. Sets occupancy to 1 for the residue's remaining atoms.
	static void zero_occ_bb_h(Pose& ps, core::Size seqpos);

public:

	static bool stripped() {return stripped_;}

	/// @brief initializes the static members of this class
	static void init(bool stripped);

	/// @brief Lists the amino acid types that a given amino acid type can be
	/// 	reduced to.
	static utility::vector1<char> list_allowable_mutations(
		char const starting_aa );

	/// @brief sets occupancy to zero for a residue's non-constellation atoms.
	static void zero_occ_for_deleted_atoms(Pose & pose, core::Size seqpos,
		char const target_aa);

	/// @brief sets occupancy to zero for constellation atoms that are not
	/// 	to be printed on output.
	static void strip_atoms(Pose & pose, core::Size seqpos,
		char const target_aa);
};

} // namespace constel
} // namespace devel
