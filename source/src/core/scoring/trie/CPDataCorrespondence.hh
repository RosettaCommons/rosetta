// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/trie/trie.functions.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_trie_CPDataCorrespondence_hh
#define INCLUDED_core_scoring_trie_CPDataCorrespondence_hh

// Unit Headers
#include <core/scoring/trie/CPDataCorrespondence.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace trie {

class CPDataCorrespondence
{
public:
	CPDataCorrespondence();

	void n_entries( Size nentries );
	void resid_for_entry( Size entry, Size resid );
	void n_connpoints_for_entry( Size entry, Size nconnections );
	void connid_for_entry_connpoint( Size entry, Size connpoint, Size residue_connid );

	Size n_entries() const;
	Size resid_for_entry( Size entry ) const;
	Size n_connpoints_for_entry( Size entry ) const;
	Size connid_for_entry_connpoint( Size entry, Size connpoint ) const;

	void note_has_pseudobonds() { has_pseudobonds_ = true; }
	bool has_pseudobonds() const { return has_pseudobonds_; }

	Size max_connpoints_for_residue() const { return max_connpoints_for_residue_; }

private:

	Size n_entries_;
	Size max_connpoints_for_residue_;
	bool has_pseudobonds_;
	utility::vector1< Size > entry_2_resid_;
	utility::vector1< Size > nconnections_for_entry_;
	utility::vector1< utility::vector1< Size > > residue_connid_for_entry_connid_;
};

CPDataCorrespondence
create_cpdata_correspondence_for_rotamerset(
	conformation::RotamerSetBase const & rotset
);

CPDataCorrespondence
create_cpdata_correspondence_for_rotamer(
	conformation::Residue const & example_rotamer
);

} // namespace trie
} // namespace scoring
} // namespace core

#endif
