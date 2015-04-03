// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/AtomBasedConstraintsScore.hh
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_AtomBasedConstraintsScore_hh
#define INCLUDED_protocols_frag_picker_scores_AtomBasedConstraintsScore_hh

#include <protocols/frag_picker/scores/AtomBasedConstraintsScore.fwd.hh>

// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.fwd.hh>
// mini
#include <numeric/xyzVector.hh>
#include <utility/vector1_bool.hh>
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief A base class for all scoring methods that need atom coordinates
/// @details The base class provides an access to atom coordinates from the current chunk
/// It takes care of caching proper atoms from the pose for fast access.
class AtomBasedConstraintsScore: public CachingScoringMethod {
public:

	/// @brief Prepare an atom-based score that utilizes some user-defined atoms
	/// @details User may provide names of atoms that will be cached when a new
	/// chunk is considered (i.e. at every do_caching() call)
	AtomBasedConstraintsScore(Size, Real, bool, Size, utility::vector1<
			std::string>, std::string);

	/// @brief Prepare an atom-based score that utilizes the following predefined atoms: N, CA, C, O and CB
	/// @details These atoms that will be cached when a new
	/// chunk is considered (i.e. at every do_caching() call)
	AtomBasedConstraintsScore(Size, Real, bool, Size, std::string);

	/// @brief In this case caching means copying coordinates of relevant atoms from a chunk's pose
	void do_caching(VallChunkOP);

	/// @brief Erases the internal array of coordinates
	void clean_up();

	/// @brief Returns true if a given atom has been successfully cached and one can use its coordinates
	/// @param residue_id the residue order number. The first is 1,
	///    the last one depends on the size of a pose i.e. the size of the current chunk
	/// @param atom_id the residue order number, the first is 1
	inline bool has_atom(Size residue_id, Size atom_id) {
		assert(residue_id<=atom_flags_.size());
		assert(atom_id<=atom_flags_[residue_id].size());
		return atom_flags_[residue_id][atom_id];
	}

	/// @brief Returns coordinates for a given atom
	/// @param residue_id the residue order number. The first is 1,
	///    the last one depends on the size of a pose i.e. the size of the current chunk
	/// @param atom_id the residue order number, the first is 1
	inline numeric::xyzVector<Real> get_atom_coordinates(Size residue_id,
			Size atom_id) {
		assert(residue_id<=chunk_atoms_xyz_.size());
		assert(atom_id<=chunk_atoms_xyz_[residue_id].size());
		return chunk_atoms_xyz_[residue_id][atom_id];
	}

	/// @brief Returns a map that defines all constrained atoms used by this scoring method
	/// @details Returned map defines the order of atoms as they are stored inside this object.
	/// Indexes defined y this map can be used as arguments for has_atom() and get_atom_coordinates() methods.
	inline std::map<std::string, Size> get_constrainable_atoms_map() {
		return constrainable_atoms_;
	}

	/// @brief returns an internal ID assigned to a given atom name
	/// @details this ID remains the same for all residues
        inline Size get_constrained_atom_id(std::string atom_name) {
		return constrainable_atoms_.find(atom_name)->second;
	}

	/// @brief returns a name of a constrained atom when its internal ID is known
	/// @details this is the oposite to get_constrained_atom_id(std::string)
	std::string get_constrained_atom_name(Size atom_id);

	/// @brief provides an access to the size of the length of a query sequence
	inline Size get_query_size() {
		return query_size_;
	}

private:
	Size query_size_;
	utility::vector1<utility::vector1<numeric::xyzVector<Real> > >
			chunk_atoms_xyz_;
	utility::vector1<utility::vector1<bool> > atom_flags_;
	Size get_atom_type(std::string atom_name);
	std::map<std::string, Size> constrainable_atoms_;
};

} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_AtomBasedConstraintsScore_HH */

