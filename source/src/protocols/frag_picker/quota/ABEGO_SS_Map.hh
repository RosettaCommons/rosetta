// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/quota/ABEGO_SS_Map.hh
/// @brief class for a quota pool based on secondary structure prediction and ABEGO torsion bin classification
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_quota_ABEGO_SS_Map_hh
#define INCLUDED_protocols_frag_picker_quota_ABEGO_SS_Map_hh

#include <protocols/frag_picker/quota/ABEGO_SS_Map.fwd.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// utility headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <string>

//Auto Headers
namespace protocols {
namespace frag_picker {
namespace quota {


Size torsion2big_bin_id(core::Real const phi,  core::Real const psi,  core::Real const omega);
Size abego_index(char);
Size ss_index(char);

/// @brief represents a single pool used by quota selector
class ABEGO_SS_Map : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~ABEGO_SS_Map();

	ABEGO_SS_Map(utility::vector1< std::pair<Size,Size> >);

	inline bool check_status(char ss_type,char abego_type) const {

	    return ss_abego_types_[ss_index(ss_type)][abego_index(abego_type)];
	}

	inline bool check_status(char ss_type,Size abego_type) const {

	    return ss_abego_types_[ss_index(ss_type)][abego_type];
	}

	inline bool check_status(Size ss_type,char abego_type) const {

	    return ss_abego_types_[ss_type][abego_index(abego_type)];
	}

	inline bool check_status(std::pair<Size,Size> bin_index) const {
	    return ss_abego_types_[bin_index.first][bin_index.second];
	}

	inline void set_status(std::pair<Size,Size> bin_index,bool new_status) {
    	    ss_abego_types_[bin_index.first][bin_index.second] = new_status;
	}

	std::string show_valid();

	char abego_char(Size abego_bin) { return all_abego_[abego_bin]; }

private:
	static char all_abego_[];
	static char all_ss_[];
	utility::vector1<utility::vector1<bool> > ss_abego_types_;
};

} // quota
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_quota_ABEGO_SS_Map_HH */
