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
#include <protocols/frag_picker/ContactCounts.hh>
#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallChunk.fwd.hh>
#include <protocols/frag_picker/VallChunkFilter.hh>
#include <protocols/frag_picker/CandidatesCollector.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <protocols/frag_picker/scores/AtomPairConstraintsScore.hh>
#include <protocols/frag_picker/scores/PValuedFragmentScoreManager.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/SidechainContactDistCutoff.hh>
#include <protocols/frag_picker/nonlocal/NonlocalPair.hh>
#include <protocols/frag_picker/ContactTypes.hh>
#include <protocols/frag_picker/quota/QuotaCollector.fwd.hh>

// type headers
#include <core/types.hh>

// core headers
#include <core/fragment/SecondaryStructure.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++
#include <string>
#include <map>
#include <sstream>
#include <set>

namespace protocols {
namespace frag_picker {

typedef std::pair<Size,Size> PosPair;
typedef std::pair<FragmentCandidateOP, scores::FragmentScoreMapOP> Candidate;
typedef utility::vector1<Candidate> Candidates;
typedef std::map<Size, CandidatesCollectorOP> CandidatesSink;

class QuotaDebug : public std::ostringstream {

public:
	Size nFrags_;
	QuotaDebug(Size nFrags) { nFrags_ = nFrags; }
	utility::vector1<std::string> tags_;
	std::map<std::string,Size> tag_map_;
	void setup_summary(quota::QuotaCollector const & collector_);
	void write_summary();
	void log(Size,Size,utility::vector1<Real>);
	Size max_pools();

private:
	QuotaDebug(QuotaDebug const &);
};

/// @brief The core of the fragment picking machinery
/// @details FragmentPicker class does:\n
///    - know about query data: sequence, sequence profile, secondary structure.
///      All other query data must be loaded directly to the relevant scoring methods
///    - provide slots for 'plugable' parts of the machinery, like chunk filters, scoring methods and so on.
///    - pick fragments
class FragmentPicker: public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< FragmentPicker >
{
public:

	// constructors

	FragmentPicker() {
		scores_.push_back(utility::pointer::shared_ptr<class protocols::frag_picker::scores::FragmentScoreManager>( new scores::FragmentScoreManager() ));
		CandidatesSink storage;
		candidates_sinks_.push_back(storage);
		max_frag_size_ = 0;
		max_threads_ = 1;
		prefix_ = "frags"; // umm... let's not make hidden files
		contacts_min_seq_sep_ = 12;
		contacts_dist_cutoff_squared_ = 81;
		nonlocal_min_contacts_per_res_ =  1.0;
	}

	FragmentPicker(std::string fragment_score_manager_type) {
		if ( fragment_score_manager_type.compare("PValuedFragmentScoreManager") == 0 ) {
			scores_.push_back(utility::pointer::shared_ptr<class protocols::frag_picker::scores::FragmentScoreManager>( new scores::PValuedFragmentScoreManager() ));
		} else {
			scores_.push_back(utility::pointer::shared_ptr<class protocols::frag_picker::scores::FragmentScoreManager>( new scores::FragmentScoreManager() ));
		}
		CandidatesSink storage;
		candidates_sinks_.push_back(storage);
		max_frag_size_ = 0;
		max_threads_ = 1;
		prefix_ = "frags"; // umm... let's not make hidden files
		contacts_min_seq_sep_ = 12;
		contacts_dist_cutoff_squared_ = 81;
		nonlocal_min_contacts_per_res_ =  1.0;
	}

	// destructor

	virtual ~FragmentPicker();

public:
	inline FragmentPickerCOP get_self_ptr() const { return shared_from_this(); }
	inline FragmentPickerOP  get_self_ptr() { return shared_from_this(); }
	inline FragmentPickerCAP get_self_weak_ptr() const { return FragmentPickerCAP( shared_from_this() ); }
	inline FragmentPickerAP  get_self_weak_ptr() { return FragmentPickerAP( shared_from_this() ); }

public:

	// public methods

	// Command line processing and high-level stuff -----------------
	void parse_command_line();

	/// @brief Picks fragments and saves them into proper files - independently for each query position.
	/// @details This protocol scores all vall data against a given position
	/// and keeps all the candidates unless they fail an energy filter
	/// When all candidates for a positions are scored, it selects fragments
	/// for the position and proceeds to the next position. Bounded queue is not used
	void keep_all_protocol();

	/// @brief Picks fragments and saves them into proper files - uses bounded queue.
	/// @details This protocol scores all vall data against all query positions
	/// and keeps a limited number of candidates per position using a bounded queue.
	/// When all candidates for all positions are scored, it selects final fragments.
	void bounded_protocol();

	void quota_protocol();

	void fragment_contacts( Size const fragment_size, utility::vector1<Candidates> const & fragment_set );

	// multi-threaded task
	void nonlocal_pairs_at_positions( utility::vector1<Size> const & positions, Size const & fragment_size, utility::vector1<bool> const & skip,
		utility::vector1<Candidates> const & fragment_set, utility::vector1<nonlocal::NonlocalPairOP> & pairs );

	void nonlocal_pairs( Size const fragment_size, utility::vector1<Candidates> const & fragment_set );

	// these should be private methods but some classes directly access these


	/// @brief returns a pointer to a scoring manager
	scores::FragmentScoreManagerOP get_score_manager(Size index=1) {
		assert(index <= scores_.size());
		return scores_[index];
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
	/// @details Well, it is a sequence profile, but the sequence can be extracted from it
	void set_query_seq(core::sequence::SequenceProfileOP query_sequence);

	/// @brief sets the query sequence
	void set_query_seq(std::string & query_sequence);

	/// @brief Returns the sequence we are picking fragments for (as a string)
	std::string& get_query_seq_string() {
		return query_seq_as_string_;
	}

	core::sequence::SequenceProfileOP get_query_seq() {
		return query_profile_;
	}

	/// Returns the solvent accessibility prediction
	utility::vector1<Real> & get_query_sa_prediction() {
		return query_sa_prediction_;
	}

	/// Returns the phi prediction
	utility::vector1<Real> & get_query_phi_prediction() {
		return query_phi_prediction_;
	}

	/// Returns the psi prediction
	utility::vector1<Real> & get_query_psi_prediction() {
		return query_psi_prediction_;
	}

	/// Returns the phi prediction confidence
	utility::vector1<Real> & get_query_phi_prediction_conf() {
		return query_phi_prediction_conf_;
	}

	/// Returns the psi prediction confidence
	utility::vector1<Real> & get_query_psi_prediction_conf() {
		return query_psi_prediction_conf_;
	}

	/// Returns residue depth values
	utility::vector1<Real> & get_query_residue_depth() {
		return query_residue_depth_;
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
		if ( query_ss_profile_.find(prediction_name) != query_ss_profile_.end() ) {
			return query_ss_profile_.find(prediction_name)->second;
		} else {
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

	/// @brief Reads spine-x phi, psi, and solvent accessibility predictions
	void read_spine_x(std::string const & file_name);

	/// @brief Reads DEPTH residue depth values
	void read_depth(std::string const & file_name);

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
	/// @details These basically become fragments if pass the final selection.
	/// Fragment candidates are stored in a container that a user must plug into the picker

	/// multi-threaded task
	void pick_chunk_candidates(utility::vector1<VallChunkOP> const & chunks, Size const & index);

	void pick_candidates();

	void pick_candidates(Size i_pos,Size frag_len);


	/// @brief Calculates total score for a given vector of small scores
	/// (FragmentScoreMap object pointer)
	/// @details FragmentScoreManager that is stored inside the picker is used
	/// for this calculation. It particular it provides weights
	double total_score(scores::FragmentScoreMapOP f, Size index=1);

	// save results
	void save_candidates();
	void save_fragments();

	/// @brief How long is the longest fragment?
	Size max_frag_len() { return *std::max_element( frag_sizes_.begin(), frag_sizes_.end() ); }

	// Delegators from FragmentScoreManager ---------------
	void show_scoring_methods(std::ostream & out, Size index=1) {
		if ( index > scores_.size() ) return;
		scores_[index]->show_scoring_methods(out);
	}

	/// @brief adds a new scoring method to the scoring scheme
	void add_scoring_method(
		scores::FragmentScoringMethodOP scoring_term,
		Real weight, Size index=1) {
		if ( index > scores_.size() ) return;
		scores_[index]->add_scoring_method(scoring_term, weight);
	}

	// Convert to a FragSet
	utility::vector1<core::fragment::ConstantLengthFragSetOP> getFragSet(int residueInPose_);

	/// @brief Reads query secondary structure prediction from a PsiPred file
	void read_psipred_ss2(std::string const &, std::string);

	/// @brief Reads query secondary structure prediction from a Talos+ file
	void read_talos_ss(std::string const &, std::string);

	// should be private but some classes directly access these
	utility::vector1<Size> frag_sizes_;

	static QuotaDebug log_25_;
	static QuotaDebug log_200_;

	void parse_quota_command_line();

	void set_up_quota_nnmake_style();

	void set_up_ss_abego_quota();

	bool is_valid_chunk( VallChunkOP chunk );

	bool is_valid_chunk( Size const frag_len, VallChunkOP chunk );

	// Output fragments
	void output_fragments( Size const fragment_size, utility::vector1<Candidates> const & final_fragments );

	// Picking machinery -----------------
	/// @brief sets a collector where fragment candidates will be kept until final selection
	void set_candidates_collector(Size frag_size,
		CandidatesCollectorOP sink, Size index = 1) {
		if ( index == 0 ) {
			candidates_sink_.insert(std::pair<Size, CandidatesCollectorOP>(frag_size, sink));
		} else {
			candidates_sinks_[index].insert(std::pair<Size, CandidatesCollectorOP>(frag_size, sink));
		}
	}

	/// @brief returns a pointer to the candidates collector currently used
	/// @details one may need this instance to access the candidates that have been found by the picker
	CandidatesCollectorOP get_candidates_collector(Size frag_size, Size index=1) {
		if ( index == 0 ) {
			return candidates_sink_[frag_size];
		} else {
			return candidates_sinks_[index][frag_size];
		}
	}

	/// @brief adds a chunk filter that will be used to screen chunks before they are cut into fragments
	void add_chunk_filter(VallChunkFilterOP filter) {
		filters_.push_back(filter);
	}

	Size n_candidates_;
	Size n_frags_;
	std::string prefix_;
	FragmentSelectingRuleOP selector_;

private:
	Size max_threads_;

	Size max_frag_size_;

	core::sequence::SequenceProfileOP query_profile_;
	std::string query_seq_as_string_;
	std::map<std::string, core::fragment::SecondaryStructureOP> query_ss_profile_;
	std::map<std::string, std::string> query_ss_as_string_;
	utility::vector1<Size> query_positions_;

	VallProviderOP chunks_;
	utility::vector1<VallChunkFilterOP> filters_;

	CandidatesSink candidates_sink_;

	// for multi-threaded picking there should only be one candidate sink per thread
	utility::vector1<CandidatesSink>  candidates_sinks_;

	// for multi-threaded scoring there should only be one FragmentScoreManager per thread due to cached chunk scores
	utility::vector1<scores::FragmentScoreManagerOP> scores_;


	// for frag contacts
	Size contacts_min_seq_sep_;
	std::set<ContactType> contact_types_;
	Real contacts_dist_cutoff_squared_;
	utility::vector1<Real> contacts_dist_cutoffs_squared_;
	SidechainContactDistCutoffOP sidechain_contact_dist_cutoff_;

	// for nonlocal contacts
	Real nonlocal_min_contacts_per_res_;

	// phi,psi,sa predictions
	utility::vector1<Real> query_sa_prediction_;
	utility::vector1<Real> query_phi_prediction_;
	utility::vector1<Real> query_psi_prediction_;
	utility::vector1<Real> query_phi_prediction_conf_;
	utility::vector1<Real> query_psi_prediction_conf_;

	// residue depth
	utility::vector1<Real> query_residue_depth_;

	// atom pair constraint contact map
	utility::vector1<utility::vector1<Real> >  atom_pair_constraint_contact_map_;

	void
	output_pair_counts(
		Size const fragment_size,
		Size const neighbors,
		std::map<std::pair<Real,ContactType>, ContactCountsOP> contact_counts
	);
};


} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_FragmentPicker_HH */

