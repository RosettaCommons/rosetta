// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/legacy_sewing/hashing/Hasher.cc
///
/// @author Tim Jacobs
/// @author Doo Nam Kim (box_length part)

//Unit headers
#include <protocols/legacy_sewing/hashing/Hasher.hh>

//Package headers
#include <basic/Tracer.hh>

//Package headers
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

//Utility headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/HomogeneousTransform.hh>

#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

//External headers
#include <boost/unordered_map.hpp>
#include <boost/math/special_functions/round.hpp>

namespace protocols {
namespace legacy_sewing  {

static THREAD_LOCAL basic::Tracer TR("protocols.legacy_sewing.Hasher");

Hasher::Hasher(){
	hash_map_.clear();
}


HashMap const &
Hasher::hash_map() const {
	return hash_map_;
}


ScoreResult
Hasher::score_one(
	Model const & m1,
	SewResidue const & m1_basis,
	Model const & m2,
	SewResidue const & m2_basis,
	core::Size box_length // Doonam introduced this argument to allow user to used customized box_length
){
	Model transformed_m1 = transform_model(m1, m1_basis);
	Model transformed_m2 = transform_model(m2, m2_basis);

	//Get a map from the m2 to m1
	hash_model(transformed_m2, m2_basis);
	ScoreResults alignment_scores;

	if ( box_length == 3 ) {
		score_basis(alignment_scores, transformed_m1, m1_basis, true);
	} else if ( box_length == 5 ) { // Doonam introduced this box_length == 5 case
		score_basis_125(alignment_scores, transformed_m1, m1_basis, true);
	} else {
		TR << "box_length should be either 3 or 5!!" << std::endl;
		exit(1);
	}

	if ( TR.Debug.visible() ) {
		TR << "alignment_scores.size(): " << alignment_scores.size() << std::endl;
	}
	runtime_assert(alignment_scores.size() == 1);
	//runtime_assert(alignment_scores.begin()->second.segment_matches.size() == 2);
	if ( TR.Debug ) {
		std::map< SegmentPair, AtomMap > segment_matches = alignment_scores.begin()->second.segment_matches;
		TR.Debug << "After scoring once." << std::endl;
		TR.Debug << "Models: " << m1.model_id_ << " " << m2.model_id_ << std::endl;
		TR.Debug << "Basis residues: " << m1_basis.resnum_ << " " << m2_basis.resnum_ << std::endl;
		TR.Debug << "\tNumber of matched segments: " << segment_matches.size() << std::endl;
		std::map< SegmentPair, AtomMap >::const_iterator it = segment_matches.begin();
		std::map< SegmentPair, AtomMap >::const_iterator it_end = segment_matches.end();
		for ( ; it != it_end; ++it ) {
			TR.Debug << "\tSegments " << it->first.first << " and " << it->first.second << " have " << it->second.size() << " overlapping atoms." << std::endl;
		}
	}

	return *alignment_scores.begin();

}


///@details Insert the features  into the HashMap. This is done by
///iterating through all possible basis sets (coordinate frames defined by atom positions), and for each set
///transform all features to the corresponding local coordinates and then insert then into the
///hash table. The complexity of this operation should be O(m^4) where m is the number of features.
///The total size of the HashMap should be O(N * m^4) where N is the number of hashed objects.
void
Hasher::insert(
	Model const & model
) {
	utility::vector1<SewSegment> segments = model.segments_;
	for ( utility::vector1<SewSegment>::const_iterator seg_it = segments.begin(); seg_it != segments.end(); ++seg_it ) {

		//Don't score the segments that we aren't hashing
		if ( ! seg_it->hash_ ) continue;

		utility::vector1<SewResidue> residues = seg_it->residues_;
		for ( core::Size basis_i=1; basis_i<=residues.size(); ++basis_i ) {
			SewResidue basis_residue = residues[basis_i];
			if ( basis_residue.basis_atoms_.size() < 3 ) {
				utility_exit_with_message("Each SewResidue must have at least 3 atoms to be used as a coordinate frame");
			}
			Model transformed_model = transform_model(model, basis_residue);
			hash_model(transformed_model, basis_residue);
		}
	}
}

ScoreResults
Hasher::score(
	Model const & model, // it1->second (Model itself) from legacy_sewing_hasher
	core::Size num_segment_matches,
	core::Size min_segment_score,
	core::Size max_clash_score,
	bool store_atoms,
	core::Size box_length // Doonam introduced this argument to allow user to used customized box_length
) const {
	TR << "[First score function in Hasher.cc]" << std::endl;

	std::set<core::Size> all_segments;
	utility::vector1<SewSegment>::const_iterator it= model.segments_.begin();

	// For 5-ss-models, this 'for loop' iterates 5 times
	for ( ; it!=model.segments_.end(); ++it ) {
		all_segments.insert(it->segment_id_);
	}
	return score(model, num_segment_matches, min_segment_score, max_clash_score, all_segments, store_atoms, box_length);
}


///@details Tally the score of each model/basis_set in the HashMap against the input pose.
///This is very similar to the insert function, but instead of populating the hash with the transformed features,
///you tally the number of HashMap hits corresponding to each structure/basis_set pair.
ScoreResults
Hasher::score(
	Model const & model, // it1->second (Model itself) from legacy_sewing_hasher
	core::Size num_segment_matches,
	core::Size min_segment_score,
	core::Size max_clash_score,
	std::set<core::Size> const & score_segments, //all_segments with all segment_ids
	bool store_atoms,
	core::Size box_length // Doonam introduced this argument to allow user to used customized box_length
) const {
	TR << "[Second score function in Hasher.cc]" << std::endl;

	ScoreResults alignment_scores; //typedef std::map< BasisPair, HashResult > ScoreResults;
	//Basis is a struct with elements of 'model_id  and resnum'
	//HashResult is a struct with elements of 'segment_matches, segment_match_counts and clash_count'

	utility::vector1<SewSegment> segments = model.segments_;
	utility::vector1<SewResidue> all_residues;
	utility::vector1<SewSegment>::const_iterator it = model.segments_.begin();
	for ( ; it != model.segments_.end(); ++it ) {

		//Don't score the segments that we aren't hashing, if it is not hashed, it is just a linker segment!
		if ( ! it->hash_ ) continue;

		if ( score_segments.find(it->segment_id_) != score_segments.end() ) {
			all_residues.insert(all_residues.begin(), it->residues_.begin(), it->residues_.end());;
		}
	}
	TR << "A size of all residues in all hashed segments: " << all_residues.size() << std::endl;

	for ( core::Size basis_i=1; basis_i<=all_residues.size(); ++basis_i ) {

		SewResidue basis_residue = all_residues[basis_i];
		//SewResidue is a struct with elements of 'resnum_, residue_type_ and num_neighbors_'

		Model transformed_model = transform_model(model, basis_residue);
		//HomogenousTransform all features into the local coordinate frame

		runtime_assert_msg ( (box_length == 3) || (box_length == 5), "box_length should be either 3 or 5 as of 2015/December" );

		//  if ( TR.Debug.visible() ) {
		//   TR.Debug << box_length*box_length*box_length << " boxes for neighborhood lookup " << std::endl;
		//  }

		if ( box_length == 3 ) {
			score_basis(alignment_scores, transformed_model, basis_residue, store_atoms);
		} else { // (box_length == 5)
			score_basis_125(alignment_scores, transformed_model, basis_residue, store_atoms);
		}

		// trim the given ScoreResults based on
		// the number of segments that match between two models,
		// the number of atom matches for each of these segments (will be compared against to min_hash_score), and the
		// clash score (number of hits between atoms of different atom_types)
		trim_scores(alignment_scores, num_segment_matches, min_segment_score, max_clash_score);

	}//for each basis_residue (from query)

	//Keep only the segment matches between two models that has the most aligned atoms.
	alignment_scores = remove_duplicates(alignment_scores);

	TR << "Done scoring model id: " << model.model_id_ << std::endl;

	return alignment_scores;
}// the 2nd score fn


void
Hasher::score_basis(
	ScoreResults & alignment_scores, // here no element, just typedef std::map< BasisPair, HashResult > ScoreResults;
	Model const & transformed_model,
	SewResidue const & basis_residue,
	bool store_atoms
) const {
	// TR << "Hasher::score_basis" << std::endl;

	utility::fixedsizearray1<HashMap::const_iterator, 27> hit_its(hash_map_.end()); //put here for speed
	Basis reference_bp(transformed_model.model_id_, basis_residue.resnum_);
	ModelConstIterator<SewSegment> model_it = transformed_model.model_begin();

	for ( ; model_it != transformed_model.model_end(); ++model_it ) {
		SewAtom const & cur_atom = *model_it.atom();
		//SewAtom is a struct with elements of 'atomno_ and coords_'

		//iterate through hits in the hash map and tally the score for the model/residue alignment pairs
		HashKey key = generate_key(cur_atom);
		neighborhood_lookup(key, hit_its);
		for ( core::Size hit_it_i = 1; hit_it_i<= hit_its.size(); ++hit_it_i ) { // hit_its.size() is 27 for 3 box_length
			if ( hit_its[hit_it_i] == hash_map_.end() ) continue;
			for ( utility::vector1<HashValue>::const_iterator it = hit_its[hit_it_i]->second.begin(); it!=hit_its[hit_it_i]->second.end(); ++it ) {
				HashValue const & cur_hit = *it;
				if ( transformed_model.model_id_ != cur_hit.model_id ) { //Don't score the model against itself
					BasisPair basis_pair = std::make_pair(reference_bp, Basis(cur_hit.model_id, cur_hit.basis_resnum));
					if ( cur_atom.atomno_ == cur_hit.atomno ) {
						using namespace basic::options;
						if ( ! option[OptionKeys::legacy_sewing::score_between_opposite_terminal_segments].user() ) {
							option[OptionKeys::legacy_sewing::score_between_opposite_terminal_segments].value( 0 );
						}
						bool score_between_opposite_terminal_segments = option[OptionKeys::legacy_sewing::score_between_opposite_terminal_segments];
						// <begin> identifying the 1st and the last segment
						utility::vector1<SewSegment> segments = transformed_model.segments_;
						utility::vector1<SewSegment>::const_iterator it= transformed_model.segments_.begin();
						core::Size segment_id_1st = 9999;
						core::Size segment_id_last = 9999;
						for ( ; it != transformed_model.segments_.end(); ++it ) {
							//Don't score the segments that we aren't hashing, if it is not hashed, it is just a linker segment!
							if ( it->hash_ ) {
								if ( segment_id_1st == 9999 ) {
									segment_id_1st = it->segment_id_;
								} else {
									segment_id_last = it->segment_id_;
								}
							}
						}// <end> identifying the 1st and the last segment
						if (
								(!score_between_opposite_terminal_segments) ||
								(
								((segment_id_1st == model_it.segment()->segment_id_) && (segment_id_last == cur_hit.segment_id))
								||
								((segment_id_last == model_it.segment()->segment_id_) && (segment_id_1st == cur_hit.segment_id))
								)
								) {
							SegmentPair segment_pair = std::make_pair(model_it.segment()->segment_id_, cur_hit.segment_id);

							std::map<SegmentPair, core::Size>::iterator seg_pair_it = alignment_scores[basis_pair].segment_match_counts.find(segment_pair);
							if ( seg_pair_it == alignment_scores[basis_pair].segment_match_counts.end() ) {
								alignment_scores[basis_pair].segment_match_counts[segment_pair] = 0;
							}
							alignment_scores[basis_pair].segment_match_counts[segment_pair]++;
							if ( store_atoms ) {
								core::id::AtomID query_atom(cur_atom.atomno_, model_it.residue()->resnum_);
								core::id::AtomID hit_atom(cur_hit.atomno, cur_hit.resnum);
								alignment_scores[basis_pair].segment_matches[segment_pair].insert(std::make_pair(query_atom, hit_atom));
							}
						}
					} else { //if ( cur_atom.atomno_ == cur_hit.atomno ) {
						alignment_scores[basis_pair].clash_count++;
					}
				}
			}//foreach hit
		}//foreach hit_iterator
	}//model iterator
}//score_basis


//Doonam added this score_basis_125 fn, to allow 125 boxes instead of typical 25 boxes
void
Hasher::score_basis_125(
	ScoreResults & alignment_scores,
	Model const & transformed_model,
	SewResidue const & basis_residue,
	bool store_atoms
) const {

	utility::fixedsizearray1<HashMap::const_iterator, 125> hit_its(hash_map_.end()); //put here for speed

	Basis reference_bp(transformed_model.model_id_, basis_residue.resnum_);

	ModelConstIterator<SewSegment> model_it = transformed_model.model_begin();


	for ( ; model_it != transformed_model.model_end(); ++model_it ) {

		SewAtom const & cur_atom = *model_it.atom();

		//iterate through hits in the hash map and tally the score for the model/residue alignment pairs
		HashKey key = generate_key(cur_atom);

		neighborhood_lookup_125(key, hit_its);

		for ( core::Size hit_it_i = 1; hit_it_i<= hit_its.size(); ++hit_it_i ) {
			if ( hit_its[hit_it_i] == hash_map_.end() ) continue;

			for ( utility::vector1<HashValue>::const_iterator it = hit_its[hit_it_i]->second.begin(); it!=hit_its[hit_it_i]->second.end(); ++it ) {
				HashValue const & cur_hit = *it;

				if ( transformed_model.model_id_ != cur_hit.model_id ) { //Don't score the model against itself
					BasisPair basis_pair = std::make_pair(reference_bp, Basis(cur_hit.model_id, cur_hit.basis_resnum));

					if ( cur_atom.atomno_ == cur_hit.atomno ) {
						SegmentPair segment_pair = std::make_pair(model_it.segment()->segment_id_, cur_hit.segment_id);

						std::map<SegmentPair, core::Size>::iterator seg_pair_it = alignment_scores[basis_pair].segment_match_counts.find(segment_pair);
						if ( seg_pair_it == alignment_scores[basis_pair].segment_match_counts.end() ) {
							alignment_scores[basis_pair].segment_match_counts[segment_pair] = 0;
						}

						alignment_scores[basis_pair].segment_match_counts[segment_pair]++;
						if ( store_atoms ) {
							core::id::AtomID query_atom(cur_atom.atomno_, model_it.residue()->resnum_);
							core::id::AtomID hit_atom(cur_hit.atomno, cur_hit.resnum);
							alignment_scores[basis_pair].segment_matches[segment_pair].insert(std::make_pair(query_atom, hit_atom));
						}
					} else {
						alignment_scores[basis_pair].clash_count++;
					}
				} //if( transformed_model.model_id_ != cur_hit.model_id ) { //Don't score the model against itself
			}//foreach hit
		}//foreach hit_iterator
	}//model iterator

} //score_basis_125


///@details trim the given ScoreResults based on the number of segments that match
///between two models, the number of atom matches for each of these segments, and the
///clash score (number of hits between atoms of different atom_types)
void
Hasher::trim_scores(
	ScoreResults & scores,
	core::Size num_segment_matches,
	core::Size min_segment_score,
	core::Size max_clash_score
) const{
	TR << "Hasher::trim_scores" << std::endl;
	using namespace core;
	using namespace basic::options;

	// Doonam introduced this disregard_num_segment_matches option for development purpose
	if ( ! option[OptionKeys::legacy_sewing::disregard_num_segment_matches].user() ) {
		option[OptionKeys::legacy_sewing::disregard_num_segment_matches].value( 0 );
	}
	bool disregard_num_segment_matches = option[OptionKeys::legacy_sewing::disregard_num_segment_matches];
	//TR << "disregard_num_segment_matches: " << disregard_num_segment_matches << std::endl;
	//  TR << "[attention] [condition] max_clash_score: " << max_clash_score << std::endl;
	//  TR << "[attention] [condition] min_segment_score: " << min_segment_score << std::endl;
	//  TR << "[attention] [condition] num_segment_matches: " << num_segment_matches << std::endl;

	ScoreResults::iterator it = scores.begin();
	//typedef std::map< BasisPair, HashResult > ScoreResults;
	ScoreResults::iterator it_end = scores.end();

	while ( it!=it_end ) {
		// TR << "it!=it_end" << std::endl;

		//Ensure that there are no clashes
		bool erase=false;

		// it->second represents HashResult which is the struct with elements of 'segment_matches, segment_match_counts and clash_count'
		if ( it->second.clash_count > max_clash_score ) {
			TR << "[reason of being erased] it->second.clash_count > max_clash_score" << std::endl;
			TR << "[being erased] it->second.clash_count: " << it->second.clash_count << std::endl;
			erase=true;
		} else if ( it->second.segment_match_counts.size() != num_segment_matches ) {
			if ( (!disregard_num_segment_matches) ) {
				//Ensure hits are between the specified number segments (and only the specified number of segments, to prevent clashes)
				TR << "[reason of being erased] it->second.segment_match_counts.size() != num_segment_matches" << std::endl;
				TR << "[being erased] it->second.segment_match_counts.size(): " << it->second.segment_match_counts.size() << std::endl;
				erase=true;
			}
		} else {
			//Ensure each matched segment contains at least the minimum number of atoms in the same bin
			std::map<SegmentPair, core::Size> const & segment_matches = it->second.segment_match_counts;
			std::map<SegmentPair, core::Size>::const_iterator seg_it = segment_matches.begin();
			std::map<SegmentPair, core::Size>::const_iterator seg_it_end = segment_matches.end();
			std::set<core::Size> source_segments;
			std::set<core::Size> target_segments;
			for ( ; seg_it != seg_it_end; ++seg_it ) {
				source_segments.insert(seg_it->first.first);
				target_segments.insert(seg_it->first.second);
				if ( seg_it->second < min_segment_score ) {
					TR << "[attention] seg_it->second (superimposed # of atoms): " << seg_it->second << std::endl;
					TR << "[reason of being erased] seg_it->second < min_segment_score" << std::endl;
					erase=true;
					break;
				}
			}

			//Delete any matches where one segment from one model matches more than one segment of another model
			if ( source_segments.size() != target_segments.size() ) {
				TR << "[reason of being erased] source_segments.size() != target_segments.size()" << std::endl;
				erase=true;
			}
		}

		// just for devel purpose
		if ( TR.Debug.visible() ) {
			std::map<SegmentPair, core::Size> const & segment_matches = it->second.segment_match_counts;
			std::map<SegmentPair, core::Size>::const_iterator seg_it = segment_matches.begin();
			std::map<SegmentPair, core::Size>::const_iterator seg_it_end = segment_matches.end();
			std::set<core::Size> source_segments;
			std::set<core::Size> target_segments;
			for ( ; seg_it != seg_it_end; ++seg_it ) {
				source_segments.insert(seg_it->first.first);
				target_segments.insert(seg_it->first.second);
				if ( seg_it->second < min_segment_score ) {
					TR << "[attention] seg_it->second (superimposed # of atoms): " << seg_it->second << std::endl;
					TR << "[reason of being erased] seg_it->second < min_segment_score" << std::endl;
					erase=true;
					break;
				}
			}
		}

		if ( !erase ) {
			++it;
		} else {
			scores.erase(it++);
		}
	}
}//trim_scores

///@details Keep only the segment matches between two models that has the most aligned atoms.
ScoreResults
Hasher::remove_duplicates(
	ScoreResults const & scores
) const {
	TR << "Hasher::remove_duplicates" << std::endl;
	//there is almost certainly a better way to do this....
	typedef std::pair<int, std::set<core::Size> > score_node;

	std::map< std::pair< score_node, score_node >, std::pair<core::Size, ScoreResults::const_iterator> > best_hits;

	ScoreResults::const_iterator it = scores.begin();
	ScoreResults::const_iterator it_end = scores.end();
	for ( ; it != it_end; ++it ) {
		score_node model_1_node;
		model_1_node.first = it->first.first.model_id;

		score_node model_2_node;
		model_2_node.first = it->first.second.model_id;

		core::Size pair_score = 0;
		std::map<SegmentPair, core::Size> segment_matches = it->second.segment_match_counts;
		std::map<SegmentPair, core::Size>::const_iterator seg_it = segment_matches.begin();
		std::map<SegmentPair, core::Size>::const_iterator seg_it_end = segment_matches.end();
		for ( ; seg_it != seg_it_end; ++seg_it ) {
			pair_score += seg_it->second;
			model_1_node.second.insert(seg_it->first.first);
			model_2_node.second.insert(seg_it->first.second);
		}

		std::pair<score_node, score_node> cur_node_pair = std::make_pair(model_1_node, model_2_node);

		std::map< std::pair< score_node, score_node >, std::pair<core::Size, ScoreResults::const_iterator> >::const_iterator find_it =
			best_hits.find(cur_node_pair);
		if ( find_it == best_hits.end() || find_it->second.first < pair_score ) {
			best_hits[cur_node_pair] = std::make_pair(pair_score, it);
		}
	}

	ScoreResults trimmed_scores;
	std::map< std::pair< score_node, score_node >, std::pair<core::Size, ScoreResults::const_iterator> >::const_iterator best_it = best_hits.begin();
	std::map< std::pair< score_node, score_node >, std::pair<core::Size, ScoreResults::const_iterator> >::const_iterator best_end = best_hits.end();
	for ( ; best_it != best_end; ++best_it ) {
		trimmed_scores.insert( *best_it->second.second );
	}
	return trimmed_scores;

}//remove_duplicates


//remove edges between segments that both have 'next' or 'previous' segments
void
Hasher::remove_connection_inconsistencies(
	std::map< int, Model > const & models,
	ScoreResults & scores
) const {
	TR << "[remove_connection_inconsistencies]" << std::endl;
	ScoreResults::iterator scores_it = scores.begin();
	ScoreResults::iterator scores_it_end = scores.end();
	while ( scores_it != scores_it_end ) {
		BasisPair bp = scores_it->first;

		int model_id_1 = bp.first.model_id;
		int model_id_2 = bp.second.model_id;
		TR << "model_id_1: " << model_id_1 << std::endl;
		TR << "model_id_2: " << model_id_2 << std::endl;

		std::map<SegmentPair, core::Size> segment_matches = scores_it->second.segment_match_counts;
		std::map<SegmentPair, core::Size>::const_iterator seg_it = segment_matches.begin();
		std::map<SegmentPair, core::Size>::const_iterator seg_it_end = segment_matches.end();
		bool erase = false;
		for ( ; seg_it != seg_it_end; ++seg_it ) {
			bool has_next_1 = false;
			bool has_next_2 = false;
			bool has_previous_1 = false;
			bool has_previous_2 = false;

			Model model_1 = models.find(model_id_1)->second;
			Model model_2 = models.find(model_id_2)->second;

			for ( core::Size i=1; i<=model_1.segments_.size(); ++i ) {
				if ( model_1.segments_[i].segment_id_ == seg_it->first.first ) {
					if ( model_1.segments_.has_next(i) ) { has_next_1 = true; }
					if ( model_1.segments_.has_previous(i) ) { has_previous_1 = true; }
					break;
				}
			}
			for ( core::Size i=1; i<=model_2.segments_.size(); ++i ) {
				if ( model_2.segments_[i].segment_id_ == seg_it->first.second ) {
					if ( model_2.segments_.has_next(i) ) { has_next_2 = true; }
					if ( model_2.segments_.has_previous(i) ) { has_previous_2 = true; }
					break;
				}
			}

			if ( (has_next_1 && has_next_2) || (has_previous_1 && has_previous_2) ) {
				erase = true;
				break;
			}
		}
		if ( erase ) {
			scores.erase(scores_it++);
		} else {
			++scores_it;
		}
	}

}


///@details when doing a lookup, look in each of the quarter angstrom bins surrounding
///the query key. This should prevent issues of close matches being missed due to being
///across bin boundaries.
//utility::vector1<HashValue>
void
Hasher::neighborhood_lookup(
	HashKey const & key,
	utility::fixedsizearray1<HashMap::const_iterator, 27> & hit_its
) const {

	for ( int i=-1; i<=1; ++i ) {
		for ( int j=-1; j<=1; ++j ) {
			for ( int k=-1; k<=1; ++k ) {
				HashKey modified_key = key;
				modified_key[1]+=i;
				modified_key[2]+=j;
				modified_key[3]+=k;
				core::Size index = ((i+1) * 3 * 3) + ((j+1) * 3) + (k+1);
				index++; //1-indexed array
				hit_its[index] = hash_map_.find(modified_key);
			}
		}
	}
}


///@details when doing a lookup, look in each of the quarter angstrom bins surrounding
///the query key. This should prevent issues of close matches being missed due to being
///across bin boundaries.
//utility::vector1<HashValue>
// Doonam introduced this neighborhood_lookup_125 fn to deal with 125 box length case
void
Hasher::neighborhood_lookup_125(
	HashKey const & key,
	utility::fixedsizearray1<HashMap::const_iterator, 125> & hit_its
) const {
	int index = 0; // just initial value
	int box_number = 0;

	for ( int i=-2; i<=2; ++i ) {
		for ( int j=-2; j<=2; ++j ) {
			for ( int k=-2; k<=2; ++k ) {
				box_number++;

				HashKey modified_key = key;
				// The HashKey is the 3 indices that describe the location the feature in discretized space, and the atomno of the point being described.
				// The hasher is responsible for creating a (hopefully) uniform distribution of keys in the table.

				modified_key[1]+=i;
				modified_key[2]+=j;
				modified_key[3]+=k;

				index = ((i+2) * 5 * 5) + ((j+2) * 5) + (k+2);
				index++; //1-indexed array
				hit_its[index] = hash_map_.find(modified_key);
			}
		}
	}
}


///@details Construct HomogenousTransform using the 3 points in the BasisSet. Use this
///HomogenousTransform to transform all features into the local coordinate frame
Model
Hasher::transform_model(
	Model const & model_in,
	SewResidue const & basis_residue
) const {

	utility::vector1<SewAtom> basis_atoms = basis_residue.basis_atoms_;
	numeric::HomogeneousTransform<core::Real> ht(basis_atoms[1].coords_,
		basis_atoms[2].coords_, basis_atoms[3].coords_);

	Model model = model_in;
	ModelIterator<SewSegment> model_it = model.model_begin();
	ModelIterator<SewSegment> erase_it;
	for ( ; model_it != model.model_end(); ++model_it ) {
		if ( model_it.residue()->resnum_ != basis_residue.resnum_ && model_it.segment()->hash_ ) {
			model_it.atom()->coords_ = ht.to_local_coordinate(model_it.atom()->coords_);
		} else {
			erase_it = model_it;
		}
	}
	erase_it.segment()->residues_.erase(erase_it.residue());

	return model;
}

///@details For each of the transformed features, insert the appropriate
///key-value pair into the hash table. The key is the 3D-voxel (or bin)
///corresponding to the basis set. The value is the residues number for the
///residue that generated basis set that the transformed atom coordinates are in frame of
void
Hasher::hash_model(
	Model const & transformed_model,
	SewResidue const & basis_residue
) {
	//For each set of transformed query_residues, count up and save hits to other models
	ModelConstIterator<SewSegment> model_it = transformed_model.model_begin();
	for ( ; model_it != transformed_model.model_end(); ++model_it ) {
		HashValue value;
		value.model_id = transformed_model.model_id_;
		value.basis_resnum = basis_residue.resnum_;

		value.segment_id = model_it.segment()->segment_id_;
		value.resnum = model_it.residue()->resnum_;

		value.atomno = model_it.atom()->atomno_;

		HashKey const & key = generate_key(*model_it.atom());
		hash_map_[key].push_back(value);
	}
}


HashKey
Hasher::generate_key(
	SewAtom const & atom
) const {
	HashKey key;
	key[1] = boost::math::iround(atom.coords_.x()*4);
	key[2] = boost::math::iround(atom.coords_.y()*4);
	key[3] = boost::math::iround(atom.coords_.z()*4);
	return key;
}



///@details write the hash table to disk!
void
Hasher::write_to_disk(std::string filename) const {

	utility::io::ozstream file;
	file.open(filename);
	for ( HashMap::const_iterator it = hash_map_.begin(); it != hash_map_.end(); ++it ) {
		HashKey const & key = it->first;
		file << "KEY " << key[1] << " " << key[2] << " " << key[3] << std::endl;

		utility::vector1<HashValue> const & values = it->second;
		for ( core::Size i=1; i<=values.size(); ++i ) {
			file << "VALUE " << values[i].model_id << " " << values[i].basis_resnum << " " << values[i].segment_id << " " << values[i].resnum << " " << values[i].atomno << std::endl;
		}
	}
	file.close();
}



///@details read the hash table from disk. This function clears the contents of the hash map
void
Hasher::read_from_disk(std::string filename) {

	utility::io::izstream file(filename);
	if ( !file.good() ) {
		utility_exit_with_message("Could not find Hasher file with name: " + filename);
	}

	hash_map_.clear();

	HashKey key;
	HashValue value;
	std::string line;
	while ( getline( file, line ) ) {

		utility::vector1<std::string> tokens = utility::string_split(line);
		debug_assert(tokens.size() > 0);

		if ( tokens[1]=="KEY" ) {
			debug_assert(tokens.size() == 4);
			key[1]=utility::string2int(tokens[2]);
			key[2]=utility::string2int(tokens[3]);
			key[3]=utility::string2int(tokens[4]);
		} else if ( tokens[1]=="VALUE" ) {
			debug_assert(tokens.size() == 6);
			value.model_id = utility::string2int(tokens[2]);
			value.basis_resnum = utility::string2int(tokens[3]);
			value.segment_id = utility::string2int(tokens[4]);
			value.resnum = utility::string2int(tokens[5]);
			value.atomno = utility::string2int(tokens[6]);
			hash_map_[key].push_back(value);
		}
	}
	file.close();
}



} //legacy_sewing namespace
} //protocols namespace
