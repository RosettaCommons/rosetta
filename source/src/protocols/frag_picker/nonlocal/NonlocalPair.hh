// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/nonlocal/NonlocalPair.hh
/// @brief  a nonlocal fragment pair.
/// @author David E. Kim (dekim@u.washington.edu)


#ifndef INCLUDED_protocols_frag_picker_nonlocal_NonlocalPair_hh
#define INCLUDED_protocols_frag_picker_nonlocal_NonlocalPair_hh

// unit headers
#include <protocols/frag_picker/nonlocal/NonlocalPair.fwd.hh>

// package headers
#include <protocols/frag_picker/Contact.hh>
#include <protocols/frag_picker/ContactTypes.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// type headers
#include <core/types.hh>

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>
#include <map>

namespace protocols {
namespace frag_picker {
namespace nonlocal {

typedef std::pair<Size,Size> PosPair;
typedef std::pair<FragmentCandidateOP, scores::FragmentScoreMapOP> Candidate;

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace protocols::frag_picker;

/// @brief  represents a nonlocal fragment pair.
/// @details NonlocalPair contains query start positions for each fragment of pair i and j, and their corresponding fragment candidate
class NonlocalPair: public utility::pointer::ReferenceCount {
public:
	NonlocalPair( Size query_pos_i, Size query_pos_j, Candidate candidate_i, Candidate candidate_j,
		Size candidate_i_rank, Size candidate_j_rank, utility::vector1<ContactOP> & contacts ) {
		query_position_i_ = query_pos_i;
		query_position_j_ = query_pos_j;
		candidate_i_ = candidate_i;
		candidate_j_ = candidate_j;
		candidate_i_rank_ = candidate_i_rank;
		candidate_j_rank_ = candidate_j_rank;
		contacts_ = contacts;
	}

	~NonlocalPair(){};

	utility::vector1<ContactOP> & get_contacts() {
		return contacts_;
	}

	Candidate get_candidate_i() {
		return candidate_i_;
	}

	Candidate get_candidate_j() {
		return candidate_j_;
	}

	Size get_candidate_i_rank() {
		return candidate_i_rank_;
	}

	Size get_candidate_j_rank() {
		return candidate_j_rank_;
	}

	Size get_query_pos_i() {
		return query_position_i_;
	}

	Size get_query_pos_j() {
		return query_position_j_;
	}

	void print(std::ostream& out) {
		out << "pair: " << query_position_i_ << " " << query_position_j_  << " " <<
			candidate_i_.first->get_residue(1)->resi() << " " << candidate_j_.first->get_residue(1)->resi() <<
			" " << candidate_i_rank_ << " " << candidate_j_rank_;
		std::map<ContactType, Size> contact_type_cnt;
		std::map<ContactType, Size>::iterator iter;
		for ( Size i=1; i<=contacts_.size(); ++i ) contact_type_cnt[contacts_[i]->type()]++;
		for ( iter = contact_type_cnt.begin(); iter != contact_type_cnt.end(); iter++ ) {
			out << " " << contact_name(iter->first) << " " << iter->second;
		}
		out << std::endl;
		candidate_i_.first->print_fragment(out);
		out << std::endl;
		candidate_j_.first->print_fragment(out);
		out << std::endl;
	}

private:
	Size query_position_i_;
	Size query_position_j_;
	Candidate candidate_i_;
	Candidate candidate_j_;
	Size candidate_i_rank_;
	Size candidate_j_rank_;
	utility::vector1<ContactOP> contacts_;
};

} // nonlocal
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_nonlocal_NonlocalPair_HH */
