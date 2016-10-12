// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/VallChunkFilter.hh
/// @brief  defines two basic chunk filters based on pdb id
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#ifndef INCLUDED_protocols_frag_picker_PdbIdChunkFilter_hh
#define INCLUDED_protocols_frag_picker_PdbIdChunkFilter_hh

// package headers
#include <protocols/frag_picker/VallChunkFilter.hh>
#include <protocols/frag_picker/PdbIdChunkFilter.fwd.hh>

#include <core/types.hh>

#include <map>

#include <utility/vector1_bool.hh>
#include <iostream>


namespace protocols {
namespace frag_picker {

class PdbIdChunkFilter: public VallChunkFilter {
public:
	/// @brief Adds  a pdb id to the filter.
	/// @details What that PDB id means it depends on the filter implementation
	void add_pdb_id(std::string pdb_ids) {

		pdb_hash_[pdb_ids] = true;
		if ( pdb_ids.size() > 4 ) {
			if ( (pdb_ids[4]=='A') || (pdb_ids[4]=='a') ) {
				pdb_hash_[ pdb_ids.replace(4,1,1,'_') ] = true;
			}
			if ( pdb_ids[4]=='_' ) {
				pdb_hash_[ pdb_ids.replace(4,1,1,'A') ] = true;
			}
		}
	}

	/// @brief Adds a bunch of pdb ids to the filter
	void add_pdb_ids(utility::vector1<std::string> list_of_pdb_ids) {

		for ( core::Size i = 1; i <= list_of_pdb_ids.size(); i++ ) {
			pdb_hash_[list_of_pdb_ids[i]] = true;
			std::string & pdb_ids = list_of_pdb_ids[i];
			if ( pdb_ids.size() > 4 ) {
				if ( (pdb_ids[4]=='A') || (pdb_ids[4]=='a') ) {
					pdb_hash_[ pdb_ids.replace(4,1,1,'_') ] = true;
				}
				if ( pdb_ids[4]=='_' ) {
					pdb_hash_[ pdb_ids.replace(4,1,1,'A') ] = true;
				}
			}
		}
	}

	void load_pdb_id_from_file(std::string);

	void show_pdb_ids(std::ostream& out);

protected:
	std::map<std::string, bool> pdb_hash_;
};

/// @brief Accepts a chunk based on the pdb id of the source protein
/// @details If a given chunk comes from a registered pdb file then it will PASS the test.
/// Otherwise it will be rejected
class AllowPdbIdFilter: public PdbIdChunkFilter {
public:

	/// @brief say if a given chunk looks promising.
	/// @details Simply if its pdb id is on a list then it will pass the test
	bool test_chunk(VallChunkOP a_chunk) override;
};

/// @brief Denies a chunk based on the pdb id of the source protein
/// @details If a given chunk comes from a registered pdb file then it will FAIL the test.
/// Otherwise it will be accepted
class DenyPdbIdFilter: public PdbIdChunkFilter {
public:

	/// @brief say if a given chunk looks promising.
	/// @details Simply if its pdb id is on a list then it will be rejected
	bool test_chunk(VallChunkOP a_chunk) override;
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_PdbIdChunkFilter_HH */
