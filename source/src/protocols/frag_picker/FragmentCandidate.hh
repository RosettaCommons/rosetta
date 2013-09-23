// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/FragmentCandidate.hh
/// @brief Something that might become a fragment if its scores will be good enough
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_FragmentCandidate_hh
#define INCLUDED_protocols_frag_picker_FragmentCandidate_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <core/fragment/FragData.hh>
#include <core/io/silent/SilentFileData.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/fragment/BBTorsionSRFD.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>


namespace protocols {
namespace frag_picker {

using namespace core;
using namespace core::fragment;

using ObjexxFCL::format::F;
using ObjexxFCL::format::I;

utility::vector1<FragmentCandidateOP> read_fragment_candidates(std::string,
		VallProviderOP, Size max_nfrags_per_pos = 900000000);

/// @brief Vector candidate says which X-mer from vall fits to a query sequence
/// @detailed Scores for a given fragment are stored separately in a FragmentScoreMap object
/// Therefore fragment containers hold std::pair<FragmentCandidateOP,FragmentScoreMapOP>
///
class FragmentCandidate: public utility::pointer::ReferenceCount {
public:

	FragmentCandidate(Size queryPosition, Size inChunkPosition,
			VallChunkOP chunk, Size fragmentLength) :
		chunk_(chunk) {

		assert(queryPosition>0);
		assert(inChunkPosition>0);
		queryResidueIndex_ = queryPosition;
		vallResidueIndex_ = inChunkPosition;
		fragmentLength_ = fragmentLength;
		pool_name_ = NULL;
	}

	virtual ~FragmentCandidate() { delete pool_name_; }

	/// @brief returns a pointer to the original chunk from vall the fragment comes from
	inline VallChunkOP get_chunk() const {
		return chunk_;
	}

	/// @brief returns a given residue from this fragment
	/// @detailed the irgument is in the range [1,fragmentLength]
	inline VallResidueOP get_residue(Size whichOne) const {
		runtime_assert( chunk_ );
		runtime_assert( whichOne >= 1 && whichOne <= fragmentLength_ );
		return chunk_->at(whichOne + vallResidueIndex_ - 1);
	}


	/// @brief creates a new string object that contains a sequence of this fragment
	inline std::string sequence() { return chunk_->get_sequence().substr(vallResidueIndex_-1, fragmentLength_); }

	/// @brief returns an integer key identifying a fragment
	inline Size key() const {
		return chunk_->at(vallResidueIndex_)->key();
	}

	/// @brief returns a PDB id of a protein from which the fragment has been extracted
	inline std::string get_pdb_id() const {
		return chunk_->get_pdb_id();
	}

	/// @brief returns a chain id of a protein from which the fragment has been extracted
	inline char get_chain_id() const {
		return chunk_->get_chain_id();
	}

	/// @brief returns the index of a very first residue in a query sequence that is covered by this fragment
	inline Size get_first_index_in_query() const {
		return queryResidueIndex_;
	}

	/// @brief returns the index of a very first residue in a Vall chunk that is covered by this fragment
	inline Size get_first_index_in_vall() const {
		return vallResidueIndex_;
	}

	/// @brief returns a vall index of a middle residue in this fragment
	/// @returns a position of the middle residue of this fragment in the vall sequence
	inline Size get_vall_middle_res_id() {
		return ( fragmentLength_/2 + 1) + ( vallResidueIndex_ ) - 1;
	}

	/// @brief returns a query index of a middle residue in this fragment
	/// @returns a position of the middle residue of this fragment in the query sequence
	inline Size get_query_middle_res_id() {
		return ( fragmentLength_/2 + 1) + ( queryResidueIndex_ ) - 1;
	}

	/// @brief returns the middle residue of this fragment candidate
	inline VallResidueOP get_middle_residue() const {
		return get_residue( fragmentLength_/2 + 1);
	}
	/// @brief returns secondary structure assigned to the middle residue of this fragment candidate
	/// @returns secondary structure of the middle residue of this fragment, as extracted from Vall data
	inline char get_middle_ss() const {

		return get_residue( fragmentLength_/2 + 1 )->ss();
	}

	/// @brief returns the length of this fragment
	inline Size get_length() const {
		return fragmentLength_;
	}

	inline FragDataOP get_frag_data() {

		AnnotatedFragDataOP fragdata = new AnnotatedFragData(get_pdb_id(),
				queryResidueIndex_);

		for (Size i = 1; i <= fragmentLength_; ++i) {
			fragdata->add_residue(
					chunk_->at(i + vallResidueIndex_ - 1)->bbtorsion_srfd());
		}

		if (fragdata->size() > 0)
			fragdata->set_valid();

		return fragdata;
	}

	/// @brief Prints fragment data, the output can be directly loaded to minirosetta
	void print_fragment(std::ostream& out, scores::FragmentScoreMapOP sc = NULL, scores::FragmentScoreManagerOP ms = NULL);

	/// @brief Prints fragment sequence, used for generating structure based sequence profiles
	void print_fragment_seq(std::ostream& out);

	/// @brief Prints fragment to silent struct
	void output_silent(core::io::silent::SilentFileData & sfd, const std::string sequence, const std::string silent_file_name, const std::string tag, scores::FragmentScoreMapOP sc, scores::FragmentScoreManagerOP ms);

	inline void set_pool_name(std::string pool_name) {
	    if(pool_name_!=NULL) delete pool_name_;
	    pool_name_= new std::string(pool_name);
	}

	inline std::string get_pool_name() {
	    if(pool_name_==NULL) return unknown_pool_name_;
	    else return *pool_name_;
	}

	bool same_chain( FragmentCandidateCOP fr );

protected:
	VallChunkOP chunk_;
	Size vallResidueIndex_;
	Size queryResidueIndex_;
	Size fragmentLength_;
private:
	static const std::string unknown_pool_name_;
	std::string* pool_name_;
};

inline std::ostream& operator<<(std::ostream& out, FragmentCandidate const& fr) {
	out << fr.get_pdb_id() << " " << fr.get_first_index_in_vall() << " : "
			<< fr.get_first_index_in_query() << " : ";
	return out;
}

inline std::ostream& operator<<(std::ostream& out, std::pair<
		FragmentCandidateOP, scores::FragmentScoreMapOP> const& pair) {

	out << pair.first->get_pdb_id() << " "
			<< pair.first->get_first_index_in_vall() << " : "
			<< pair.first->get_first_index_in_query() << " :";
	utility::vector1<Real> c = pair.second->get_score_components();
	for (Size i = 1; i <= c.size(); i++)
		out << " " << c.at(i);
	return out;
}

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_FragmentCandidate_HH */
