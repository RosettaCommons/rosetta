// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/FragmentPicker.hh
/// @brief  Fragment picker - the core part of picking machinery
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_FragmentPicker_hh
#define INCLUDED_protocols_frag_picker_FragmentPicker_hh

// unit headers
#include <protocols/frag_picker/FragmentPicker.fwd.hh>
// package headers
#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallChunk.fwd.hh>
#include <protocols/frag_picker/VallChunkFilter.hh>
#include <protocols/frag_picker/CandidatesCollector.hh>
#include <protocols/frag_picker/quota/QuotaCollector.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <protocols/frag_picker/scores/PValuedFragmentScoreManager.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// utility headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/fragment/SecondaryStructure.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/fragment/ConstantLengthFragSet.hh>

// C++
#include <string>
#include <map>
#include <sstream>

namespace protocols {
namespace frag_picker {

class QuotaDebug : public std::ostringstream {

public:
    Size nFrags_;

    QuotaDebug(Size nFrags) { nFrags_ = nFrags; }
    utility::vector1<std::string> tags_;
    std::map<std::string,Size> tag_map_;
    void setup_summary(quota::QuotaCollector* collector_);
    void write_summary();
    void log(Size,Size,utility::vector1<Real>);
    Size max_pools();

private:
	QuotaDebug(QuotaDebug const &);
};


/// @brief The core of the fragment picking machinery
/// @detailed FragmentPicker class does:\n
///    - know about query data: sequence, sequence profile, secondary structure.
///      All other query data must be loaded directly to the relevant scoring methods
///    - provide slots for 'plugable' parts of the machinery, like chunk filters, scoring methods and so on.
///    - pick fragments
class FragmentPicker: public utility::pointer::ReferenceCount {
public:

	utility::vector1<Size> frag_sizes_;
	Size n_frags_;
	Size n_candidates_;
	std::string prefix_;
	FragmentSelectingRuleOP selector_;
	scores::FragmentScoreManagerOP scores_;

	FragmentPicker() {
		scores_ = new scores::FragmentScoreManager();
		max_frag_size_ = 0;
	}

	FragmentPicker(std::string fragment_score_manager_type) {

		if (fragment_score_manager_type.compare("PValuedFragmentScoreManager") == 0)
		    scores_ = (scores::FragmentScoreManager*) new scores::PValuedFragmentScoreManager();
		else
		    scores_ = new scores::FragmentScoreManager();
		max_frag_size_ = 0;
	}

	~FragmentPicker();

	// Command line processing and high-level stuff -----------------
	void parse_command_line();

	void parse_quota_command_line();

	void set_up_quota_nnmake_style();

	void set_up_ss_abego_quota();

	// Picking machinery -----------------
	/// @brief sets a collector where fragment candidates will be kept until final selection
	void set_candidates_collector(Size frag_size,
			CandidatesCollectorOP sink) {
		candidates_sink_.insert(std::pair<Size, CandidatesCollectorOP>(
				frag_size, sink));
	}

	/// @brief returns a pointer to the candidates collector currently used
	/// @details one may need this instance to access the candidates that have been found by the picker
	CandidatesCollectorOP get_candidates_collector(Size frag_size) {
		return candidates_sink_[frag_size];
	}

	/// @brief adds a chunk filter that will be used to screen chunks before they are cut into fragments
	void add_chunk_filter(VallChunkFilterOP filter) {
		filters_.push_back(filter);
	}

	/// @brief returns a pointer to a scoring manager
	scores::FragmentScoreManagerOP get_score_manager() {
		return scores_;
	}

	// vall stuff -----------------
	/// @brief reads a vall file
	void read_vall( std::string const & fn );
	void read_vall( utility::vector1< std::string > const & fns );

	/// @brief sets vall data
	void set_vall(VallProviderOP vall_chunks) {
		chunks_ = vall_chunks;
	}

	/// @returns a pointer to Vall provider
	inline VallProviderOP get_vall() {
		return chunks_;
	}

	// query sequence -----------------

	/// @brief sets the query sequence
	/// @detailed Well, it is a sequence profile, but the sequence can be extracted from it
	void set_query_seq(core::sequence::SequenceProfileOP query_sequence) {

		query_profile_ = query_sequence;
		query_seq_as_string_ = query_profile_->sequence();
		set_picked_positions(1,query_sequence->length());
	}

	/// @brief sets the query sequence
	void set_query_seq(std::string & query_sequence) {

		if (query_profile_ == 0)
			query_profile_ = new core::sequence::SequenceProfile();
		query_profile_->sequence(query_sequence);
		query_seq_as_string_ = query_sequence;
		set_picked_positions(1,query_profile_->length());
	}

	/// @brief Returns the sequence we are picking fragments for (as a string)
	 std::string& get_query_seq_string() {
		return query_seq_as_string_;
	}

	core::sequence::SequenceProfileOP get_query_seq() {
		return query_profile_;
	}

	/// @brief Sets the sequence object we are picking fragments for
	void set_query_profile(core::sequence::SequenceProfileOP profile) {
		query_profile_ = profile;
		query_seq_as_string_ = query_profile_->sequence();
		set_picked_positions(1,profile->length());
	}

	// query secondary structure -----------------
	/// @brief Returns the query secondary structure as a SecondaryStructure object
	core::fragment::SecondaryStructureOP get_query_ss(std::string prediction_name) {
		if (query_ss_profile_.find(prediction_name) != query_ss_profile_.end())
			return query_ss_profile_.find(prediction_name)->second;
		else {

			return 0;
		}
	}

	/// @brief Reads a bunch of ss predicitons from files and plugs them into the picker
	void read_ss_files(utility::vector1<std::string> sec_str_input);

	/// @brief Returns the query secondary structure as a string
	inline std::string & get_query_ss_string(std::string prediction_name) {
		return query_ss_as_string_.find(prediction_name)->second;
	}

	/// @brief Sets the query secondary structure
	void add_query_ss(std::string, std::string);

	/// @brief Sets the query secondary structure
	inline Size count_query_ss() {
		return query_ss_profile_.size();
	}

	/// @brief Identifies if the ss file is psipred or talos, calls appropriate reader
	void read_ss_file(std::string const &, std::string);

	// other stuff -----------------
	inline Size size_of_query() {
		assert(query_seq_as_string_.length() == query_profile_->length());
		return query_seq_as_string_.length();
	}

	/// @brief Asks the picker to pick fragments from a given range in a query sequence
	void set_picked_positions(Size,Size);

	/// @brief Asks the picker to pick fragments for given positions in a query sequence
	void set_picked_positions(utility::vector1<Size>);

	/// @brief picks fragment candidates.
	/// @detailed These basically become fragments if pass the final selection.
	/// Fragment candidates are stored in a container that a user must plug into the picker
	void pick_candidates();

	void pick_candidates(Size i_pos,Size frag_len);


	/// @brief Calculates total score for a given vector of small scores
	/// (FragmentScoreMap object pointer)
	/// @detailed FragmentScoreManager that is stored inside the picker is used
	/// for this calculation. It particular it provides weights
	double total_score(scores::FragmentScoreMapOP);

	/// @brief Picks fragments and saves them into proper files - independently for each query position.
	/// @detailed This protocol scores all vall data against a given position
	///	and keeps all the candidates unless they fail an energy filter
	///	When all candidates for a positions are scored, it selects fragments
	///	for the position and proceeds to the next position. Bounded queue is not used
	void keep_all_protocol();
	/// @brief Picks fragments and saves them into proper files - uses bounded queue.
	/// @detailed This protocol scores all vall data against all query positions
	///	and keeps a limited number of candidates per position using a bounded queue.
	///	When all candidates for all positions are scored, it selects final fragments.
	void bounded_protocol();

	void quota_protocol();

	void nonlocal_pairs_protocol();

	// save results
	void save_candidates();
	void save_fragments();

	/// @brief How long is the longest fragment?
	Size max_frag_len() { return *std::max_element( frag_sizes_.begin(), frag_sizes_.end() ); }

	// Delegators from FragmentScoreManager ---------------
	void show_scoring_methods(std::ostream & out) {
		scores_->show_scoring_methods(out);
	}

	/// @brief adds a new scoring method to the scoring scheme
	void add_scoring_method(
			scores::FragmentScoringMethodOP scoring_term,
			Real weight) {
		scores_->add_scoring_method(scoring_term, weight);
	}

	// Convert to a FragSet
	utility::vector1<core::fragment::ConstantLengthFragSetOP> getFragSet(int residueInPose_);

	static QuotaDebug log_25_;
	static QuotaDebug log_200_;
private:

        void setup_summary(quota::QuotaCollector* collector_);
	void write_summary();

	core::sequence::SequenceProfileOP query_profile_;
	std::string query_seq_as_string_;
	std::map<std::string, core::fragment::SecondaryStructureOP>
			query_ss_profile_;
	std::map<std::string, std::string> query_ss_as_string_;
	std::map<Size, CandidatesCollectorOP> candidates_sink_;
	VallProviderOP chunks_;
	utility::vector1<VallChunkFilterOP> filters_;
	Size max_frag_size_;

	utility::vector1<Size> query_positions_;

	/// @brief Reads query secondary structure prediction from a PsiPred file
	void read_psipred_ss2(std::string const &, std::string);

	/// @brief Reads query secondary structure prediction from a Talos+ file
	void read_talos_ss(std::string const &, std::string);
};



} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_FragmentPicker_HH */

