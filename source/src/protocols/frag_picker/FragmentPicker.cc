// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/FragmentPicker.cc
/// @brief  Fragment picker - the core part of picking machinery
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

// unit headers
#include <protocols/frag_picker/FragmentPicker.hh>

#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/VallChunkFilter.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/CandidatesCollector.hh>
#include <protocols/frag_picker/GrabAllCollector.hh>
#include <protocols/frag_picker/quota/QuotaSelector.hh>
#include <protocols/frag_picker/quota/QuotaConfig.hh>
#include <protocols/frag_picker/quota/QuotaCollector.hh>
#include <protocols/frag_picker/quota/ABEGO_SS_Config.hh>
#include <protocols/frag_picker/quota/ABEGO_SS_Pool.hh>
#include <protocols/frag_picker/quota/SecondaryStructurePool.hh>
#include <protocols/frag_picker/BoundedCollector.hh>
#include <protocols/frag_picker/PdbIdChunkFilter.hh>
#include <protocols/frag_picker/BestTotalScoreSelector.hh>
#include <protocols/frag_picker/CustomScoreSelector.hh>
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <protocols/frag_picker/scores/SecondarySimilarity.hh>
#include <protocols/frag_picker/scores/PartialSecondarySimilarity.hh>
#include <protocols/frag_picker/scores/ProfileScoreL1.hh>
#include <protocols/frag_picker/scores/RamaScore.hh>
#include <protocols/frag_picker/scores/CSScore.hh>
#include <protocols/frag_picker/scores/ABEGO_SS_Score.hh>
#include <protocols/frag_picker/scores/TorsionBinSimilarity.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/TorsionBinIO.hh>

#include <protocols/frag_picker/Contact.hh>
#include <protocols/frag_picker/ContactCounts.hh>
#include <protocols/frag_picker/nonlocal/NonlocalPair.hh>

#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/id/SequenceMapping.hh>

#include <core/fragment/FragID.hh> // required for windows build
#include <core/fragment/SecondaryStructure.hh>
#include <core/fragment/ConstantLengthFragSet.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/model_quality/rms.hh>
#include <core/fragment/util.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <utility>
#include <sstream>
#include <iostream>
#include <fstream>

#include <ObjexxFCL/format.hh>


#if defined MULTI_THREADED && defined CXX11
#include <thread>
#elif defined USE_BOOST_THREAD
// Boost headers
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#endif

namespace protocols {
namespace frag_picker {

using namespace core;
using namespace core::fragment;
using namespace protocols::frag_picker;
using namespace protocols::frag_picker::scores;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer tr( "protocols.frag_picker.FragmentPicker" );

FragmentPicker::~FragmentPicker() {}

void FragmentPicker::bounded_protocol() {
	tr.Info << "pick fragments using bounded protocol..." << std::endl;
	pick_candidates();
	save_fragments();
}

void FragmentPicker::quota_protocol() {
	using namespace ObjexxFCL;
	tr.Info << "pick fragments using quota protocol..." << std::endl;
	pick_candidates();

	const bool skip_merge = (candidates_sinks_.size() == 1) ? true : false;
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
		Size fragment_size = frag_sizes_[iFragSize];
		quota::QuotaCollectorOP c = (skip_merge) ? 
			utility::pointer::dynamic_pointer_cast<quota::QuotaCollector> (candidates_sinks_[1][fragment_size]) :
			utility::pointer::dynamic_pointer_cast<quota::QuotaCollector> (candidates_sink_[fragment_size]); // merged storage
		if (c == 0)
			utility_exit_with_message("Cant' cast candidates' collector to QuotaCollector. Is quota set up correctly?");
		log_25_.setup_summary(*c);
		log_200_.setup_summary(*c);
		Size maxqpos = size_of_query() - fragment_size + 1;

		utility::vector1<Candidates> final_fragments(maxqpos); // final fragments

		for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) {
			Size qPos = query_positions_[iqpos];
			if ( qPos > maxqpos) continue;
			Candidates dummy_input, final_out;
			if (!skip_merge) {  // merge candidates
				for (Size i=1;i<=candidates_sinks_.size();++i)
					candidates_sink_[fragment_size]->insert(qPos, candidates_sinks_[i][fragment_size]);
			}
			// select the final fragments
			quota::QuotaSelector selector(n_frags_, qPos, c );
			selector.select_fragments(dummy_input,final_out);
			final_fragments[qPos] = final_out;
		}
		log_25_.write_summary();
		log_200_.write_summary();

		output_fragments( fragment_size, final_fragments );
	}
}

void FragmentPicker::keep_all_protocol() {
	using namespace ObjexxFCL;
	tr.Info << "pick fragments using keep-all protocol..." << std::endl;
	// memory usage makes the following options impractical
	if (max_threads_ > 1)
		tr.Warning << "Ignoring -j option for keep_all_protocol" << std::endl;
	if (option[frags::nonlocal_pairs].user())
		tr.Warning << "Ignoring -nonlocal_pairs option for keep_all_protocol" << std::endl;

	CompareTotalScore comparator(get_score_manager());
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
		Size fragment_size = frag_sizes_[iFragSize];

		utility::io::ozstream out_file;
		if (option[frags::describe_fragments].user()) {
			std::string describe_name = option[frags::describe_fragments]()+"." + string_of(n_frags_) +"."+string_of(fragment_size)+"mers";
			out_file.open(describe_name.c_str());
		}

		std::string out_file_name = prefix_ + "." + string_of(n_frags_) + "." + string_of(fragment_size) + "mers";
		utility::io::ozstream output(out_file_name);
		CandidatesCollectorOP storage = get_candidates_collector(fragment_size);

		Size maxqpos = size_of_query() - fragment_size + 1;
		for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) {
			Size qPos = query_positions_[iqpos];
			if ( qPos > maxqpos ) continue;

			Candidates out;
			pick_candidates(qPos,fragment_size);
			Candidates candidates = storage->get_candidates(qPos);
			std::sort(candidates.begin(),candidates.end(),comparator);
			selector_->select_fragments(candidates, out);
			if(out.size() == 0) continue;
			output << "position: " << I(12, qPos) << " neighbors:   " << I(10,out.size()) << std::endl << std::endl;
			FragmentScoreManagerOP ms = get_score_manager();
			if( ms->if_late_scoring_for_zeros() )  {
				for (Size fi = 1; fi <= out.size(); ++fi)
					ms->score_zero_scores(out[fi].first,out[fi].second);
			}
			for (Size fi = 1; fi <= out.size(); ++fi) {
				out[fi].first->print_fragment(output);
				output << std::endl;
			}
			if (option[frags::describe_fragments].user()) {
				get_score_manager()->describe_fragments(out, out_file);
			}
			tr.Info << "Collected candidates of size "<<fragment_size<<" at pos"<<qPos<<std::endl;
			storage->clear();
			tr.Debug<< storage->count_candidates()<<" candidates left in a sink after flushing"<<std::endl;
		}
		output.close();
		out_file.close();
	}
	tr.Info<<std::endl;
}



void FragmentPicker::fragment_contacts( Size const fragment_size, utility::vector1<Candidates> const & fragment_set ) {
	using namespace ObjexxFCL;

	// how many neighboring residues shall we also find contacts for?
	Size neighbors = option[frags::contacts::neighbors]();

	// output all contacts?
	std::set<ContactType>::iterator it;
	utility::io::ozstream output_all_contacts;
	bool output_all = option[frags::contacts::output_all]();
	if (output_all) {
		std::string scale_factor = "0";
		for (it=contact_types_.begin(); it!=contact_types_.end(); it++)
			if (*it == CEN) scale_factor = string_of(sidechain_contact_dist_cutoff_->scale_factor());
		replace( scale_factor.begin(), scale_factor.end(), '.', '_' );
		const std::string out_file_name_all_contacts = prefix_ + "." + string_of(contacts_min_seq_sep_) + "." + string_of(sqrt(contacts_dist_cutoff_squared_)) + "." + scale_factor +
				"." + string_of(n_frags_) + "." + string_of(fragment_size) + "mers.contacts";
		output_all_contacts.open(out_file_name_all_contacts.c_str());
	}

	// initialize contact counts
	std::map<std::pair<Real,ContactType>, ContactCountsOP> contact_counts;
	for (it=contact_types_.begin(); it!=contact_types_.end(); it++) {
		if (*it == CEN) {
			std::pair<Real,ContactType> p(0,*it);
			contact_counts[p] = protocols::frag_picker::ContactCountsOP( new ContactCounts() );
		} else {
			for (Size i=1; i<=contacts_dist_cutoffs_squared_.size();++i) {
				std::pair<Real,ContactType> p(contacts_dist_cutoffs_squared_[i],*it);
				contact_counts[p] = protocols::frag_picker::ContactCountsOP( new ContactCounts() );
			}
		}
	}

	// output all header
	if (output_all)
		output_all_contacts << "# i j type dist cutoff frag_pos frag_rank" << std::endl;

	for (Size iqpos = 1; iqpos <= query_positions_.size()-fragment_size+1; ++iqpos) {
		Size qPosi = query_positions_[iqpos];
		Candidates const & outi = fragment_set[qPosi];
		for (Size fi = 1; fi <= outi.size(); ++fi) { // loop through selected fragments at qPosi

      // for neighboring contacts
      // chunks can be different chains
      VallChunkOP chunk = outi[fi].first->get_chunk();
      int cPos_offset = outi[fi].first->get_first_index_in_vall() - qPosi;

			for (Size i=1; i<=fragment_size;++i) {
				VallResidueOP ri = outi[fi].first->get_residue(i);
				Size q_pos_i = qPosi + i - 1;
				for (Size j=i+1; j<=fragment_size;++j) {
					Size q_pos_j = qPosi + j - 1;

					// skip local contacts relative to query
					if (std::abs(int(q_pos_i-q_pos_j)) < (int)contacts_min_seq_sep_) continue;
					// skip local contacts relative to fragments
					if (std::abs(int(ri->resi() - outi[fi].first->get_residue(j)->resi() )) < (int)contacts_min_seq_sep_) continue;

					for (it=contact_types_.begin(); it!=contact_types_.end(); it++) {

						// pair distance
						Real distance_squared = ri->distance_squared(outi[fi].first->get_residue(j), *it);

						// distance cutoff
						Real cutoff_dist_squared = (*it == CEN) ?
								sidechain_contact_dist_cutoff_->get_cutoff_squared( ri->aa(), outi[fi].first->get_residue(j)->aa() ) :
								contacts_dist_cutoff_squared_;

						if (distance_squared <= cutoff_dist_squared) {

							// output all row
							if (output_all)
								output_all_contacts << q_pos_i << " " << q_pos_j << " " << contact_name(*it) << " " <<
										format::F(5, 2, sqrt(distance_squared)) << " " << format::F(5, 2, sqrt(cutoff_dist_squared)) << " " << qPosi << " "<< fi << std::endl;

							// sorry for the copy and paste code below
							if (*it == CEN) {
								// iterate contact counts
								std::pair<Real,ContactType> p(0,*it);
								std::pair<Size,Size> querypair(q_pos_i, q_pos_j);
								contact_counts[p]->iterate(querypair);

								// iterate neighboring contact counts
								if (neighbors > 0) {
									int m_min_tmp = q_pos_i-neighbors;
									Size m_min = (m_min_tmp >= 1) ? m_min_tmp : 1;
									Size m_max = q_pos_i+neighbors;
									int n_min_tmp = q_pos_j-neighbors;
									Size n_min = (n_min_tmp >= 1) ? n_min_tmp : 1;
									Size n_max = q_pos_j+neighbors;
									// m == query position i
									for (Size m = m_min; m <= m_max; ++m) {
										if (m > size_of_query()) continue;
										// chunk_i = chunk position i
										Size chunk_i = cPos_offset + m;
										if (chunk_i < 1 || chunk_i > chunk->size()) continue;
										// n == query position j
										for (Size n = n_min; n <= n_max; ++n) {
											if (n > size_of_query()) continue;
											if (m == q_pos_i && n == q_pos_j) continue;
											// chunk_j = chunk position j
											int chunk_j = cPos_offset + n;
											if (chunk_j < 1 || chunk_j > (int)chunk->size()) continue;

											// skip local contacts relative to query
											if (std::abs(int(m-n)) < (int)contacts_min_seq_sep_) continue;
											// skip local contacts relative to fragments
											if (std::abs(int( chunk->at(chunk_i)->resi() - chunk->at(chunk_j)->resi() )) < (int)contacts_min_seq_sep_) continue;

											Real dist_squared = chunk->at(chunk_i)->distance_squared(chunk->at(chunk_j), *it);
											if (dist_squared <= sidechain_contact_dist_cutoff_->get_cutoff_squared( chunk->at(chunk_i)->aa(), chunk->at(chunk_j)->aa() )) {
												 std::pair<Size,Size> neighbor_pair(m, n);
												 contact_counts[p]->iterate_neighbor(querypair, neighbor_pair);
											}
										}
									}
								}

							} else {
								for (Size cdi=1; cdi<=contacts_dist_cutoffs_squared_.size();++cdi) {
									if (distance_squared < contacts_dist_cutoffs_squared_[cdi]) {
										// iterate contact counts
										std::pair<Real,ContactType> p(contacts_dist_cutoffs_squared_[cdi],*it);
										std::pair<Size,Size> querypair(q_pos_i, q_pos_j);
										contact_counts[p]->iterate(querypair);

										// iterate neighboring contact counts
										if (neighbors > 0) {
											int m_min_tmp = q_pos_i-neighbors;
											Size m_min = (m_min_tmp >= 1) ? m_min_tmp : 1;
											Size m_max = q_pos_i+neighbors;
											int n_min_tmp = q_pos_j-neighbors;
											Size n_min = (n_min_tmp >= 1) ? n_min_tmp : 1;
											Size n_max = q_pos_j+neighbors;
											// m == query position i
											for (Size m = m_min; m <= m_max; ++m) {
												if (m > size_of_query()) continue;
												// chunk_i = chunk position i
												int chunk_i = cPos_offset + m;
												if (chunk_i < 1 || chunk_i > (int)chunk->size()) continue;
												// n == query position j
												for (Size n = n_min; n <= n_max; ++n) {
													if (n > size_of_query()) continue;
													if (m == q_pos_i && n == q_pos_j) continue;
													// chunk_j = chunk position j
													int chunk_j = cPos_offset + n;
													if (chunk_j < 1 || chunk_j > (int)chunk->size()) continue;

													// skip local contacts relative to query
													Size m_n_sep = std::abs(int(m - n));  // sequence separation
													if (m_n_sep < contacts_min_seq_sep_) continue;
													// skip local contacts relative to fragments
													if (std::abs(int( chunk->at(chunk_i)->resi() - chunk->at(chunk_j)->resi() )) < (int)contacts_min_seq_sep_) continue;

													Real dist_squared = chunk->at(chunk_i)->distance_squared(chunk->at(chunk_j), *it);
													if (dist_squared <= contacts_dist_cutoffs_squared_[cdi]) {
														 std::pair<Size,Size> neighbor_pair(m, n);
														 contact_counts[p]->iterate_neighbor(querypair, neighbor_pair);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if (output_all) output_all_contacts.close();

	// now output pair counts - sorry for the copy and paste code here
	for (it=contact_types_.begin(); it!=contact_types_.end(); it++) {
		if (*it == CEN) {
			std::string scale_factor = string_of(sidechain_contact_dist_cutoff_->scale_factor());
			replace( scale_factor.begin(), scale_factor.end(), '.', '_' );
			const std::string out_file_name_contacts = prefix_ + "." + contact_name(*it) + "." + string_of(contacts_min_seq_sep_) + "." + scale_factor + "."  +
					string_of(n_frags_) + "." + string_of(fragment_size) + "mers.contacts";
			utility::io::ozstream output_contacts(out_file_name_contacts);
			output_contacts << "# i j count";
			if (neighbors > 0) output_contacts << " neighbors_" << neighbors << "_i_j_count";
			output_contacts << std::endl;
			std::pair<Real,ContactType> p(0,*it);
			std::map<std::pair<Size,Size>, Size> query_counts = contact_counts[p]->counts();
			std::map<std::pair<Size,Size>, Size>::iterator iter;
			for ( iter = query_counts.begin(); iter != query_counts.end(); iter++ ) {
				std::pair<Size,Size> query_pair = iter->first;
				output_contacts << query_pair.first << " " << query_pair.second << " " << iter->second;
				if (neighbors > 0 && contact_counts[p]->neighbor_counts_exist(query_pair)) {
					std::map<std::pair<Size,Size>, Size> neighbor_counts = contact_counts[p]->neighbor_counts(query_pair);
					std::map<std::pair<Size,Size>, Size>::iterator neigh_iter;
					for ( neigh_iter = neighbor_counts.begin(); neigh_iter != neighbor_counts.end(); neigh_iter++ ) {
						std::pair<Size,Size> neighbor_pair = neigh_iter->first;
						output_contacts << " " << neighbor_pair.first << " " << neighbor_pair.second << " " << neigh_iter->second;
					}
				}
				output_contacts << std::endl;
			}
			output_contacts.close();
		} else {
			for (Size i=1; i<=contacts_dist_cutoffs_squared_.size();++i) {
				const std::string out_file_name_contacts = prefix_ + "." + contact_name(*it) + "." + string_of(contacts_min_seq_sep_) + "." + string_of(sqrt(contacts_dist_cutoffs_squared_[i])) +
						"." + string_of(n_frags_) + "."  + string_of(fragment_size) + "mers.contacts";
				utility::io::ozstream output_contacts(out_file_name_contacts);
				output_contacts << "# i j count";
				if (neighbors > 0) output_contacts << " neighbors_" << neighbors << "_i_j_count";
				output_contacts << std::endl;
				std::pair<Real,ContactType> p(contacts_dist_cutoffs_squared_[i],*it);
				std::map<std::pair<Size,Size>, Size> query_counts = contact_counts[p]->counts();
				std::map<std::pair<Size,Size>, Size>::iterator iter;
				for ( iter = query_counts.begin(); iter != query_counts.end(); iter++ ) {
					std::pair<Size,Size> query_pair = iter->first;
					output_contacts << query_pair.first << " " << query_pair.second << " " << iter->second;
					if (neighbors > 0 && contact_counts[p]->neighbor_counts_exist(query_pair)) {
						std::map<std::pair<Size,Size>, Size> neighbor_counts = contact_counts[p]->neighbor_counts(query_pair);
						std::map<std::pair<Size,Size>, Size>::iterator neigh_iter;
						for ( neigh_iter = neighbor_counts.begin(); neigh_iter != neighbor_counts.end(); neigh_iter++ ) {
							std::pair<Size,Size> neighbor_pair = neigh_iter->first;
							output_contacts << " " << neighbor_pair.first << " " << neighbor_pair.second << " " << neigh_iter->second;
						}
					}
					output_contacts << std::endl;
				}
				output_contacts.close();
			}
		}
	}
}


// should be thread safe
void FragmentPicker::nonlocal_pairs_at_positions( utility::vector1<Size> const & positions, Size const & fragment_size, utility::vector1<bool> const & skip,
				utility::vector1<Candidates> const & fragment_set, utility::vector1<nonlocal::NonlocalPairOP> & pairs ) {

	Size const maxjqpos = size_of_query()-fragment_size+1;

	// loop through query positions, qPosi
	for (Size p = 1; p <= positions.size(); ++p) {
		Size const qPosi = positions[p];
		Candidates const & outi = fragment_set[qPosi];  // candidates at i
		Size const minjqpos = qPosi+fragment_size+contacts_min_seq_sep_-1;
		// loop through nonlocal query positions, qPosj
		for (Size jqpos = 1; jqpos <= query_positions_.size(); ++jqpos) {
			Size qPosj = query_positions_[jqpos];
			if (qPosj > maxjqpos || qPosj < minjqpos) continue;
			bool skip_it = true;
			for (Size i=0; i<fragment_size;++i) {
				if (!skip[qPosi+i] || !skip[qPosj+i]) {
					skip_it = false;
					break;
				}
			}
			if (skip_it) continue;
			Candidates const & outj = fragment_set[qPosj]; // candidates at j
			for (Size fi = 1; fi <= outi.size(); ++fi) { // loop through selected fragments at qPosi
				for (Size fj = 1; fj <= outj.size(); ++fj) { // loop through selected fragments at qPosj
					if (!outi[fi].first->same_chain( outj[fj].first )) continue; // skip if not from same pdb chain
					//if (outi[fi].first->get_residue(1)->resi() >= outj[fj].first->get_residue(1)->resi()) continue; // skip inverse pairs
					//if (std::abs(int(outi[fi].first->get_residue(1)->resi()-outj[fj].first->get_residue(1)->resi())) < min_pdb_seq_sep) continue; // skip if too local in PDB
					if (std::abs(int(outi[fi].first->get_residue(1)->resi()-outj[fj].first->get_residue(1)->resi())) < (int)fragment_size) continue; // skip overlapping fragments in PDB
					Size qpi = qPosi; // query position i in fragment
					utility::vector1<ContactOP> contacts;
					bool skip = false;
					bool has_good_constraint = false;
					bool has_constraints = (atom_pair_constraint_contact_map_.size() > 0) ? true : false;
					for (Size i=1; i<=fragment_size;++i) {
						VallResidueOP ri = outi[fi].first->get_residue(i);
						Size qpj = qPosj; // query position j in fragment
						for (Size j=1; j<=fragment_size;++j) {
							if (skip) continue;
							if (std::abs(int(qpi-qpj)) < (int)contacts_min_seq_sep_) continue;
							// skip local contacts relative to fragments
							if (std::abs(int( ri->resi()-outj[fj].first->get_residue(j)->resi() )) < (int)contacts_min_seq_sep_) continue;
							std::set<ContactType>::iterator it;
							for (it=contact_types_.begin(); it!=contact_types_.end(); it++) {
								// contact distance cutoff
								Real cutoff_dist_squared = (*it == CEN) ?
										sidechain_contact_dist_cutoff_->get_cutoff_squared( ri->aa(), outj[fj].first->get_residue(j)->aa() ) :
										contacts_dist_cutoff_squared_;
								// contact distance
								Real dist_squared = ri->distance_squared(outj[fj].first->get_residue(j), *it);
								if (has_constraints && atom_pair_constraint_contact_map_[qpi][qpj] > 0) {
									if (dist_squared > atom_pair_constraint_contact_map_[qpi][qpj]) {
										skip = true;
										continue;
									} else {
										has_good_constraint = true;
									}
								}
								if (dist_squared <= cutoff_dist_squared) contacts.push_back(ContactOP( new Contact( qpi, qpj, dist_squared, *it ) ));
							}
							qpj++;
						}
						qpi++;
					}
					if (!skip && contacts.size() > 0 && (!has_constraints || has_good_constraint)) {
						// save all fragment pairs with contacts
						nonlocal::NonlocalPairOP pair( new nonlocal::NonlocalPair( qPosi, qPosj, outi[fi], outj[fj], fi, fj, contacts ) );
						pairs.push_back(pair);
					} // contact
				} // fi
			} // fj
		} // jqpos
	}

}

void FragmentPicker::nonlocal_pairs( Size const fragment_size, utility::vector1<Candidates> const & fragment_set ) {
	using namespace ObjexxFCL;

	// always print ca coords
	bool orig_opt = option[frags::write_ca_coordinates]();
	option[frags::write_ca_coordinates].value(true);

	// how many neighboring residues shall we also find contacts for?
	Size neighbors = option[frags::contacts::neighbors]();

	// native
	bool has_native = false;
	core::pose::PoseOP nativePose;
	if (option[in::file::native].user()) {
		nativePose = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::native]());
		has_native = true;
	} else if (option[in::file::s].user()) {
		nativePose = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::s]()[1]);
		has_native = true;
	}

	// atom pair constraints?
  if (option[constraints::cst_file].user()) {
    tr.Info << "Reading constraints from: "
        << option[constraints::cst_file]()[1] << std::endl;
		// initialize
		atom_pair_constraint_contact_map_.resize(size_of_query());
		for ( Size qi = 1; qi <= size_of_query(); qi++) {
			for ( Size qj = 1; qj <= size_of_query(); qj++) {
				atom_pair_constraint_contact_map_[qi].push_back( 0 );
			}
		}

		utility::io::izstream data(option[constraints::cst_file]()[1].c_str());
		if (!data) {
			utility_exit_with_message("[ERROR] Unable to open constraints file: "
        + option[constraints::cst_file]()[1]);
		}
		std::string line;
		getline(data, line); // header line
		std::string tag;
		Size n_constr = 0;
		while (!data.fail()) {
			char c = data.peek();
			if (c == '#' || c == '\n') {
				getline(data, line); //comment
				continue;
			}
			data >> tag;
			if (data.fail()) {
				tr.Debug << option[constraints::cst_file]()[1]
					<< " end of file reached" << std::endl;
				break;
			}
			if (tag == "AtomPair") {
				std::string name1, name2, func_type;
				Size id1, id2;
				data >> name1 >> id1 >> name2 >> id2 >> func_type;
				tr.Debug << "read: " << name1 << " " << id1
					<< " " << name2 << " " << id2 << " func: " << func_type
					<< std::endl;
				if (id1 <= size_of_query() && id2 <= size_of_query()) {
					atom_pair_constraint_contact_map_[id1][id2] = 81.0; // hard code this for a casp10 hack
					atom_pair_constraint_contact_map_[id2][id1] = 81.0;
					n_constr++;
				}
			}
		}
		tr.Info << n_constr << " constraints loaded from a file" << std::endl;
  }


	// skip positions from an input alignment if one exists
	utility::vector1<bool> skip_position( size_of_query(), false );
	//bool has_positions_to_skip = false;
	if (option[ in::file::alignment ].user()) {
		utility::vector1<core::sequence::SequenceAlignment> alns =
			core::sequence::read_aln( option[ cm::aln_format ](), option[ in::file::alignment ]()[1] );
		tr.Info << "Input alignment used to skip aligned positions: " << std::endl;
		tr.Info << alns[1] << std::endl;
		Size const query_idx( 1 );
		Size const templ_idx( 2 );
		Size nres = size_of_query();
		core::id::SequenceMapping mapping_(
			alns[1].sequence_mapping( query_idx, templ_idx )
		);
		for ( Size resi = 1; resi <= nres; resi++ ) {
			Size t_resi = mapping_[ resi ];
			bool const gap_exists( t_resi == 0 ); // query residue maps to a gap
			if ( !gap_exists ) {
				skip_position[resi] = true;
				//has_positions_to_skip = true;  set but never used ~Labonte
			}
		}
	}

	// initialize contact counts
	std::map<std::pair<Real,ContactType>, ContactCountsOP> contact_counts;
	std::set<ContactType>::iterator it;
	for (it=contact_types_.begin(); it!=contact_types_.end(); it++) {
		if (*it == CEN) {
			std::pair<Real,ContactType> p(0,*it);
			contact_counts[p] = protocols::frag_picker::ContactCountsOP( new ContactCounts() );
		} else {
			for (Size i=1; i<=contacts_dist_cutoffs_squared_.size();++i) {
				std::pair<Real,ContactType> p(contacts_dist_cutoffs_squared_[i],*it);
				contact_counts[p] = protocols::frag_picker::ContactCountsOP( new ContactCounts() );
			}
		}
	}


	Real const min_contacts = (nonlocal_min_contacts_per_res_ < 1.0) ? 1.0 : nonlocal_min_contacts_per_res_*(Real)fragment_size;
	Size const maxiqpos = size_of_query()-(contacts_min_seq_sep_-1)-fragment_size-fragment_size+1;

	time_t time_start = time(NULL);

	utility::vector1<utility::vector1<Size> > qPosi_to_run( max_threads_ );
	Size positions_cnt = 0;
	for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) {
		Size qPosi = query_positions_[iqpos];
		if (qPosi > maxiqpos) continue;
		positions_cnt++;
	}
	const Size qPosi_per_thread = positions_cnt/max_threads_;
	Size thread = 1;
	for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) {
		Size qPosi = query_positions_[iqpos];
		if (qPosi > maxiqpos) continue;
		qPosi_to_run[thread].push_back( qPosi );
		if (qPosi_to_run[thread].size() >= qPosi_per_thread && thread < max_threads_) ++thread;
	}
	utility::vector1<utility::vector1<nonlocal::NonlocalPairOP> > thread_pairs(max_threads_);

#if (defined MULTI_THREADED && defined CXX11) || defined USE_BOOST_THREAD

#if defined MULTI_THREADED && defined CXX11
	utility::vector1<std::thread> threads;
#elif defined USE_BOOST_THREAD
	boost::thread_group threads;
#endif
	tr.super_mute(true); // lets suppress tracer output when running multi threads
	for (Size j = 1; j <= max_threads_; ++j) {
		if (qPosi_to_run[j].size() > 0) {
			std::cout << "thread: " << j << " - " << qPosi_to_run[j].size() << " positions -";
			for (Size pos = 1; pos <= qPosi_to_run[j].size(); ++pos) std::cout << " " << qPosi_to_run[j][pos];
			std::cout << std::endl;
#if defined MULTI_THREADED && defined CXX11
			threads.push_back( std::thread( boost::bind(
				&FragmentPicker::nonlocal_pairs_at_positions,
				this,
				boost::ref(qPosi_to_run[j]),
				fragment_size,
				boost::ref(skip_position),
				boost::ref(fragment_set),
				boost::ref(thread_pairs[j]))));
			//&FragmentPicker::nonlocal_pairs_at_positions, this, qPosi_to_run[j], fragment_size, skip_position, fragment_set, thread_pairs[j]));
#elif defined USE_BOOST_THREAD
			threads.create_thread(boost::bind(&FragmentPicker::nonlocal_pairs_at_positions, this, boost::ref(qPosi_to_run[j]), fragment_size, boost::ref(skip_position),
				boost::ref(fragment_set), boost::ref(thread_pairs[j])));
#endif
		}
	}
#if defined MULTI_THREADED && defined CXX11
	for (auto& th : threads) th.join();
#elif defined USE_BOOST_THREAD
	threads.join_all();
#endif
	tr.super_mute(false);

#else // (defined MULTI_THREADED && defined CXX11) || defined USE_BOOST_THREAD

	// single thread
	nonlocal_pairs_at_positions( qPosi_to_run[1], fragment_size, skip_position, fragment_set, thread_pairs[1] );

#endif

	// silent output
	std::string scale_factor = "0";
	for (it=contact_types_.begin(); it!=contact_types_.end(); it++)
	    if (*it == CEN) scale_factor = string_of(sidechain_contact_dist_cutoff_->scale_factor());
	replace( scale_factor.begin(), scale_factor.end(), '.', '_' );
	const std::string silent_out_file_name = prefix_ + "." + string_of(contacts_min_seq_sep_) + "." + string_of(sqrt(contacts_dist_cutoff_squared_)) + "." + scale_factor +
			"." + string_of(n_frags_) + "." + string_of(fragment_size) + "mers.nonlocal_pairs.out";
	core::io::silent::SilentFileData sfd;

	// contacts output
	utility::io::ozstream contacts_output_all;
	bool output_all = option[frags::contacts::output_all]();
	if (output_all) {
		const std::string contacts_out_file_name = prefix_ + "." + string_of(contacts_min_seq_sep_) + "." + string_of(sqrt(contacts_dist_cutoff_squared_)) + "." + scale_factor +
				"." + string_of(n_frags_) + "." + string_of(fragment_size) + "mers.nonlocal_pairs.contacts";
		contacts_output_all.open(contacts_out_file_name.c_str());
		// contacts output header
		contacts_output_all << "# i j type dist frag_i frag_j rank_i rank_j" << std::endl;
	}

	// get score manager
	FragmentScoreManagerOP ms = get_score_manager();

	// now output the thread pairs and merge the contacts maps
	for (Size j = 1; j <= thread_pairs.size(); ++j) {
		for (Size k = 1; k <= thread_pairs[j].size(); ++k) {

			Size qPosi = thread_pairs[j][k]->get_query_pos_i(); // fragment i pos
			Size qPosj = thread_pairs[j][k]->get_query_pos_j();	// fragment j pos

			// for neighboring contacts
			// chunks can be different chains
			VallChunkOP chunki = thread_pairs[j][k]->get_candidate_i().first->get_chunk();
			VallChunkOP chunkj = thread_pairs[j][k]->get_candidate_j().first->get_chunk();
			int cPosi_offset = thread_pairs[j][k]->get_candidate_i().first->get_first_index_in_vall() - qPosi;
			int cPosj_offset = thread_pairs[j][k]->get_candidate_j().first->get_first_index_in_vall() - qPosj;

			// get contacts
			utility::vector1<ContactOP> contacts = thread_pairs[j][k]->get_contacts();

			if (option[frags::nonlocal::output_silent]()) {
				// check if enough contacts exist to output fragment pair
				std::map<ContactType, Size> contact_type_cnt;
				std::map<ContactType, Size>::iterator iter;
				bool output_pair = false;
				for (it=contact_types_.begin(); it!=contact_types_.end(); it++)
					contact_type_cnt[*it] = 0;
				for (Size i=1; i<=contacts.size(); ++i) contact_type_cnt[contacts[i]->type()]++;
				for ( iter = contact_type_cnt.begin(); iter != contact_type_cnt.end(); iter++ ) {
					if ((Real)iter->second >= min_contacts) {
						output_pair = true;
						break;
					}
				}

				if (output_pair) {

					// make pose from frag
					utility::vector1<fragment::FragDataCOP> fragdatapair;
					fragdatapair.push_back(thread_pairs[j][k]->get_candidate_i().first->get_frag_data());
					fragdatapair.push_back(thread_pairs[j][k]->get_candidate_j().first->get_frag_data());
					std::string const & sequence = get_query_seq_string().substr(qPosi-1,fragment_size) +
							get_query_seq_string().substr(qPosj-1,fragment_size);
					pose::Pose pose;
					fragment::make_pose_from_frags( pose, sequence, fragdatapair, false );

					// output silent file
					std::stringstream tag;  // put candidate id and fragment start positions of fragment pair in tag
					// tag format example:  12asA_1_17_28_132  DO NOT CHANGE THIS FORMAT SINCE THE SCORING APP USES THIS
					tag << thread_pairs[j][k]->get_candidate_i().first->get_pdb_id() << thread_pairs[j][k]->get_candidate_i().first->get_chain_id() << // candidate id
							"_" << qPosi << "_" << qPosj << "_" <<  // fragment query positions
							thread_pairs[j][k]->get_candidate_i().first->get_residue(1)->resi() << "_" << thread_pairs[j][k]->get_candidate_j().first->get_residue(1)->resi(); // fragment template positions
					if (!pose::is_ideal_pose(pose)) {
						tr.Warning << "skipping " << tag.str() << ": non-ideal pose from VALL" << std::endl;
						continue;
					}

					// calculate rms to native
					if (has_native) {
						// get pose CA coords
						std::vector< core::Vector > pose_coords;
						std::vector< core::Vector > native_pose_coords;
						for (Size i=1; i<=pose.total_residue(); i++)
							pose_coords.push_back( pose.residue(i).xyz("CA") );
						for (Size i=0; i<fragment_size; i++) {
							// get native CA coords for frag i
							Size respos = qPosi+i;
							native_pose_coords.push_back( nativePose->residue(respos).xyz("CA") );
						}
						for (Size i=0; i<fragment_size; i++) {
							// get native CA coords for frag j
							Size respos = qPosj+i;
							native_pose_coords.push_back( nativePose->residue(respos).xyz("CA") );
						}
						int const natoms = pose_coords.size();
						FArray2D< core::Real > p1a( 3, natoms );  // orig pose
						FArray2D< core::Real > p2a( 3, natoms );  // native pose
						for ( int i = 0; i < natoms; ++i ) {
							for ( int l = 0; l < 3; ++l ) { // l = X, Y and Z
								p1a(l+1,i+1) = pose_coords[i][l];
								p2a(l+1,i+1) = native_pose_coords[i][l];
							}
						}
						// calculate rms of native to original pose
						core::Real rms_orig_native = numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
						core::pose::setPoseExtraScore( pose, "frms", rms_orig_native );
					}

					// save fragment i and fragment j data
					core::pose::add_score_line_string( pose, "qpos_i", string_of(qPosi) ); // frag query i start position
					core::pose::add_score_line_string( pose, "qpos_j", string_of(qPosj) ); // frag query j start position
					core::pose::add_score_line_string( pose, "tpos_i", string_of(thread_pairs[j][k]->get_candidate_i().first->get_residue(1)->resi()) ); // frag template i start position
					core::pose::add_score_line_string( pose, "tpos_j", string_of(thread_pairs[j][k]->get_candidate_j().first->get_residue(1)->resi()) ); // frag template j start position
					core::pose::add_score_line_string( pose, "rank_i", string_of(thread_pairs[j][k]->get_candidate_i_rank()) ); // frag i rank
					core::pose::add_score_line_string( pose, "rank_j", string_of(thread_pairs[j][k]->get_candidate_j_rank()) ); // frag j rank
					core::pose::setPoseExtraScore( pose, "fscore_i", ms->total_score(thread_pairs[j][k]->get_candidate_i().second)); // frag i score
					core::pose::setPoseExtraScore( pose, "fscore_j", ms->total_score(thread_pairs[j][k]->get_candidate_j().second)); // frag j score
					core::pose::setPoseExtraScore( pose, "fscore", ms->total_score(thread_pairs[j][k]->get_candidate_i().second) + ms->total_score(thread_pairs[j][k]->get_candidate_j().second)); // frag i+j score

					// save contact counts
					for ( iter = contact_type_cnt.begin(); iter != contact_type_cnt.end(); iter++ )
						core::pose::add_score_line_string( pose, contact_name(iter->first), string_of(iter->second) );
					core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( pose );
					ss->fill_struct( pose, tag.str() );
					sfd.write_silent_struct( *ss, silent_out_file_name );

				}
			}


			for (Size i = 1; i <= contacts.size(); ++i) {

				// output all contacts
				if (output_all)
					// i j type dist frag_i_pos frag_j_pos frag_i_rank frag_j_rank
					contacts_output_all << contacts[i]->i() << " " << contacts[i]->j() << " " << contacts[i]->type_name() << " " <<
							format::F(5, 2, contacts[i]->dist()) << " " << qPosi << " " << qPosj << " " <<
							thread_pairs[j][k]->get_candidate_i_rank() << " " << thread_pairs[j][k]->get_candidate_j_rank() << std::endl;

				// iterate contact counts
				// sorry for the copy and paste code below
				if (contacts[i]->type() == CEN) {
					std::pair<Real,ContactType> p(0,CEN);
					std::pair<Size,Size> querypair(contacts[i]->i(), contacts[i]->j());
					contact_counts[p]->iterate(querypair);

					// iterate neighboring contact counts
					if (neighbors > 0) {
						int m_min_tmp = contacts[i]->i()-neighbors;
						int n_min_tmp = contacts[i]->j()-neighbors;
						Size m_min = (m_min_tmp >= 1) ? m_min_tmp : 1;
						Size n_min = (n_min_tmp >= 1) ? n_min_tmp : 1;
						Size m_max = contacts[i]->i()+neighbors;
						Size n_max = contacts[i]->j()+neighbors;
						// m == query position i
						for (Size m = m_min; m <= m_max; ++m) {
							if (m > size_of_query()) continue;
							// chunk_i = chunk position i
							int chunk_i = cPosi_offset + m;
							if (chunk_i < 1 || chunk_i > (int)chunki->size()) continue;
							// n == query position j
							for (Size n = n_min; n <= n_max; ++n) {
								if (n > size_of_query()) continue;
								if (m == contacts[i]->i() && n == contacts[i]->j()) continue;
								// chunk_j = chunk position j
								int chunk_j = cPosj_offset + n;
								if (chunk_j < 1 || chunk_j > (int)chunkj->size()) continue;

								// skip local contacts relative to query
								if (std::abs(int(m-n)) < (int)contacts_min_seq_sep_) continue;
								// skip local contacts relative to fragments
								if (std::abs(int( chunki->at(chunk_i)->resi() - chunkj->at(chunk_j)->resi() )) < (int)contacts_min_seq_sep_) continue;

								// contact distance
								Real dist_squared = chunki->at(chunk_i)->distance_squared(chunkj->at(chunk_j), contacts[i]->type());
								if (dist_squared <= sidechain_contact_dist_cutoff_->get_cutoff_squared( chunki->at(chunk_i)->aa(), chunkj->at(chunk_j)->aa() )) {
									std::pair<Size,Size> neighbor_pair(m, n);
									contact_counts[p]->iterate_neighbor(querypair, neighbor_pair);
								}
							}
						}
					}
				} else {

					for (Size cdi=1; cdi<=contacts_dist_cutoffs_squared_.size();++cdi) {
						if (contacts[i]->dist_squared() <= contacts_dist_cutoffs_squared_[cdi]) {
							std::pair<Real,ContactType> p(contacts_dist_cutoffs_squared_[cdi],contacts[i]->type());
							std::pair<Size,Size> querypair(contacts[i]->i(), contacts[i]->j());
							contact_counts[p]->iterate(querypair);

							// iterate neighboring contact counts
							if (neighbors > 0) {
								int m_min_tmp = contacts[i]->i()-neighbors;
								int n_min_tmp = contacts[i]->j()-neighbors;
								Size m_min = (m_min_tmp >= 1) ? m_min_tmp : 1;
								Size n_min = (n_min_tmp >= 1) ? n_min_tmp : 1;
								Size m_max = contacts[i]->i()+neighbors;
								Size n_max = contacts[i]->j()+neighbors;
								// m == query position i
								for (Size m = m_min; m <= m_max; ++m) {
									if (m > size_of_query()) continue;
									// chunk_i = chunk position i
									int chunk_i = cPosi_offset + m;
									if (chunk_i < 1 || chunk_i > (int)chunki->size()) continue;
									// n == query position j
									for (Size n = n_min; n <= n_max; ++n) {
										if (n > size_of_query()) continue;
										if (m == contacts[i]->i() && n == contacts[i]->j()) continue;
										// chunk_j = chunk position j
										int chunk_j = cPosj_offset + n;
										if (chunk_j < 1 || chunk_j > (int)chunkj->size()) continue;

										// skip local contacts relative to query
										if (std::abs(int(m-n)) < (int)contacts_min_seq_sep_) continue;
										// skip local contacts relative to fragments
										if (std::abs(int( chunki->at(chunk_i)->resi() - chunkj->at(chunk_j)->resi() )) < (int)contacts_min_seq_sep_) continue;

										// contact distance
										Real dist_squared = chunki->at(chunk_i)->distance_squared(chunkj->at(chunk_j), contacts[i]->type());
										if (dist_squared <= contacts_dist_cutoffs_squared_[cdi]) {
											 std::pair<Size,Size> neighbor_pair(m, n);
											 contact_counts[p]->iterate_neighbor(querypair, neighbor_pair);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if (output_all) contacts_output_all.close();

	// now output pair counts
	// sorry for the copy and paste code below
	for (it=contact_types_.begin(); it!=contact_types_.end(); it++) {
		if (*it == CEN) {
			std::string scale_factor = string_of(sidechain_contact_dist_cutoff_->scale_factor());
			replace( scale_factor.begin(), scale_factor.end(), '.', '_' );
			const std::string out_file_name_contacts = prefix_ + "." + contact_name(*it) + "." + string_of(contacts_min_seq_sep_) + "." + scale_factor + "."	+
					string_of(n_frags_) + "." + string_of(fragment_size) + "mers.nonlocal_pairs.contacts";
			utility::io::ozstream output_contacts(out_file_name_contacts);
			output_contacts << "# i j count";
			if (neighbors > 0) output_contacts << " neighbors_" << neighbors << "_i_j_count";
			output_contacts << std::endl;
			std::pair<Real,ContactType> p(0,*it);
			std::map<std::pair<Size,Size>, Size> query_counts = contact_counts[p]->counts();
			std::map<std::pair<Size,Size>, Size>::iterator iter;
			for ( iter = query_counts.begin(); iter != query_counts.end(); iter++ ) {
				std::pair<Size,Size> query_pair = iter->first;
				output_contacts << query_pair.first << " " << query_pair.second << " " << iter->second;
				if (neighbors > 0 && contact_counts[p]->neighbor_counts_exist(query_pair)) {
					std::map<std::pair<Size,Size>, Size> neighbor_counts = contact_counts[p]->neighbor_counts(query_pair);
					std::map<std::pair<Size,Size>, Size>::iterator neigh_iter;
					for ( neigh_iter = neighbor_counts.begin(); neigh_iter != neighbor_counts.end(); neigh_iter++ ) {
						std::pair<Size,Size> neighbor_pair = neigh_iter->first;
						output_contacts << " " << neighbor_pair.first << " " << neighbor_pair.second << " " << neigh_iter->second;
					}
				}
				output_contacts << std::endl;
			}
			output_contacts.close();
		} else {
			for (Size i=1; i<=contacts_dist_cutoffs_squared_.size();++i) {
				const std::string out_file_name_contacts = prefix_ + "." + contact_name(*it) + "." + string_of(contacts_min_seq_sep_) + "." + string_of(sqrt(contacts_dist_cutoffs_squared_[i])) +
						"." + string_of(n_frags_) + "."	+ string_of(fragment_size) + "mers.nonlocal_pairs.contacts";
				utility::io::ozstream output_contacts(out_file_name_contacts);
				output_contacts << "# i j count";
				if (neighbors > 0) output_contacts << " neighbors_" << neighbors << "_i_j_count";
				output_contacts << std::endl;
				std::pair<Real,ContactType> p(contacts_dist_cutoffs_squared_[i],*it);
				std::map<std::pair<Size,Size>, Size> query_counts = contact_counts[p]->counts();
				std::map<std::pair<Size,Size>, Size>::iterator iter;
				for ( iter = query_counts.begin(); iter != query_counts.end(); iter++ ) {
					std::pair<Size,Size> query_pair = iter->first;
					output_contacts << query_pair.first << " " << query_pair.second << " " << iter->second;
					if (neighbors > 0 && contact_counts[p]->neighbor_counts_exist(query_pair)) {
						std::map<std::pair<Size,Size>, Size> neighbor_counts = contact_counts[p]->neighbor_counts(query_pair);
						std::map<std::pair<Size,Size>, Size>::iterator neigh_iter;
						for ( neigh_iter = neighbor_counts.begin(); neigh_iter != neighbor_counts.end(); neigh_iter++ ) {
							std::pair<Size,Size> neighbor_pair = neigh_iter->first;
							output_contacts << " " << neighbor_pair.first << " " << neighbor_pair.second << " " << neigh_iter->second;
						}
					}
					output_contacts << std::endl;
				}
				output_contacts.close();
			}
		}
	}

	time_t time_end = time(NULL);

	tr.Info << "... done.  Processed " << query_positions_.size() << " positions.  Time elapsed: "
		<< (time_end - time_start) << " seconds." << std::endl;
	tr.flush();

	option[frags::write_ca_coordinates].value(orig_opt);
}


// should be thread safe
void FragmentPicker::pick_chunk_candidates(utility::vector1<VallChunkOP> const & chunks, Size const & index) {
	for (Size i=1; i<=chunks.size(); ++i) {
		VallChunkOP chunk = chunks[i];
		scores_[index]->do_caching(chunk);
		scores::FragmentScoreMapOP empty_map = scores_[index]->create_empty_map();
		for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
			Size fragment_size = frag_sizes_[iFragSize];
			if (chunk->size() < fragment_size) continue; // This fragment is too short
			CandidatesCollectorOP sink = candidates_sinks_[index][fragment_size];
			for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) { // loop over positions in a query
				Size iPos = query_positions_[iqpos];
				if ( iPos > size_of_query() - fragment_size + 1 ) continue;
				// split chunk into fragment candidates and score them
				for (Size j = 1; j <= chunk->size() - fragment_size + 1; j++) {
					FragmentCandidateOP f( new FragmentCandidate(iPos, j, chunk, fragment_size) );
					if (scores_[index]->score_fragment_from_cache(f, empty_map)) {
						std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP> p(f,empty_map);
						if(sink->add(p)) empty_map  = scores_[index]->create_empty_map();
					}
				}
			} // all query positions done
		} // all fragment sizes done
		scores_[index]->clean_up();
	} // all chunks
}

void FragmentPicker::pick_candidates() {

	PROF_START( basic::FRAGMENTPICKING );

	tr.Info << "Picking candidates..." << std::endl;
	tr.flush();

	time_t time_start = time(NULL);

#if (defined MULTI_THREADED && defined CXX11) || defined USE_BOOST_THREAD

	if (max_threads_ > 1) {
		utility::vector1<utility::vector1<VallChunkOP> > chunks_to_run( max_threads_ );
		Size valid_chunks_cnt = 0;
		for (Size i = 1; i <= chunks_->size(); ++i) { // loop over provided chunks
			VallChunkOP chunk = chunks_->at(i);
			if (!is_valid_chunk( chunk )) continue;
			valid_chunks_cnt++;
		}
		const Size chunks_per_thread = valid_chunks_cnt/max_threads_;
		Size thread = 1;
		for (Size i = 1; i <= chunks_->size(); ++i) { // loop over provided chunks
			VallChunkOP chunk = chunks_->at(i);
			if (!is_valid_chunk( chunk )) continue;
			chunks_to_run[thread].push_back( chunk );
			if (chunks_to_run[thread].size() >= chunks_per_thread && thread < max_threads_) ++thread;
		}
#if defined MULTI_THREADED && defined CXX11
		utility::vector1<std::thread> threads;
#elif defined USE_BOOST_THREAD
		boost::thread_group threads;
#endif
		tr.super_mute(true); // lets suppress tracer output when running multi threads
		for (Size j = 1; j <= max_threads_; ++j) {
			if (chunks_to_run[j].size() > 0) {
				std::cout << "thread: " << j << " - " << chunks_to_run[j].size() << " chunks" << std::endl;
#if defined MULTI_THREADED && defined CXX11
				threads.push_back(std::thread(&FragmentPicker::pick_chunk_candidates,this,chunks_to_run[j],j));
#elif defined USE_BOOST_THREAD
				threads.create_thread(boost::bind(&FragmentPicker::pick_chunk_candidates, this, boost::ref(chunks_to_run[j]), j));
#endif
			}
		}
#if defined MULTI_THREADED && defined CXX11
		for (auto& th : threads) th.join();
#elif defined USE_BOOST_THREAD
		threads.join_all();
#endif
		tr.super_mute(false);

		time_t time_end = time(NULL);
		tr.Info << "... done.  Processed " << chunks_->size() << " chunks.  Time elapsed: "
			<< (time_end - time_start) << " seconds." << std::endl;
		tr.flush();

		PROF_STOP( basic::FRAGMENTPICKING );

		return;
	}

#endif  // (defined MULTI_THREADED && defined CXX11) || defined USE_BOOST_THREAD

	scores::FragmentScoreMapOP empty_map = scores_[1]->create_empty_map();

	for (Size i = 1; i <= chunks_->size(); i++) { // loop over provided chunks
		VallChunkOP chunk = chunks_->at(i); // For each chunk from a provider...
		if (!is_valid_chunk( chunk )) continue;
		tr.Trace << "Processing sequence from vall: " << chunk->get_sequence() << std::endl;

		// cache the new chunk
		scores_[1]->do_caching(chunk);

		for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
			Size fragment_size = frag_sizes_[iFragSize];
			if (chunk->size() < fragment_size) continue; // This fragment is too short

			Size maxqpos = size_of_query() - fragment_size + 1;
			CandidatesCollectorOP sink = candidates_sinks_[1][fragment_size];
			tr.Trace << "Picking fragments of size "<<fragment_size<<
				" at "<<query_positions_.size()<<" query positions"<<std::endl;
			for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) { // loop over positions in a query
				Size iPos = query_positions_[iqpos];
				if ( iPos > maxqpos ) continue;

				// split chunk into fragment candidates and score them
				for (Size j = 1; j <= chunk->size() - fragment_size + 1; ++j) {
					FragmentCandidateOP f( new FragmentCandidate(iPos, j, chunk, fragment_size) );
					if (scores_[1]->score_fragment_from_cache(f, empty_map)) {
						std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP> p(f,empty_map);
						if(sink->add(p)) empty_map  = scores_[1]->create_empty_map();
					}
				}
			} // all query positions done
		} // all fragment sizes done
		scores_[1]->clean_up();
		tr.Trace << chunk->get_pdb_id() << " done" << std::endl;
		if ( (i*100) % (chunks_->size()/100*100) == 0 ) tr.Info << (i*100) / chunks_->size()
									<< "% done at "<< chunk->get_pdb_id() << std::endl;
	} // all chunks done

	time_t time_end = time(NULL);
	tr.Info << "... done.  Processed " << chunks_->size() << " chunks.  Time elapsed: "
		<< (time_end - time_start) << " seconds." << std::endl;
	tr.flush();

	PROF_STOP( basic::FRAGMENTPICKING );
}

double FragmentPicker::total_score(scores::FragmentScoreMapOP f, Size index) {

	utility::vector1<Real> components = f->get_score_components();
	utility::vector1<Real> weights = scores_[index]->get_weights();
	Real total = 0.0;
	for (Size i = 1; i <= components.size(); i++)
		total += components.at(i) * weights.at(i);

	return total;
}


void FragmentPicker::read_ss_files(utility::vector1<std::string> sec_str_input) {
	tr.Debug << sec_str_input.size() / 2 << " secondary structure assignment(s):\n";
	for (Size i = 1; i <= sec_str_input.size(); i += 2) {
		tr.Debug << i / 2 << " " << sec_str_input[i]
			<< " file will be loaded under \"" << sec_str_input[i + 1] << "\" name\n";
		read_ss_file(sec_str_input[i], sec_str_input[i + 1]);
	}
	tr.Debug << std::endl;
}

void FragmentPicker::read_ss_file(std::string const & file_name,
		std::string prediction_name) {

	utility::io::izstream data( file_name.c_str() );
	 if ( !data ) {
			data.close();
			utility_exit_with_message( "Can't read secondary structure file: "+file_name );
	}

	std::string line, l1, l2, l3, l4, l5;
	getline( data, line );
	data.close();

	std::istringstream line_stream( line );
	line_stream >> l1 >> l2 >> l3 >> l4 >> l5;

	if ( (l1 == "#") && (l2 == "PSIPRED") && (l3 == "VFORMAT")
			 && (l4 == "(PSIPRED") ) {
		read_psipred_ss2( file_name, prediction_name);
	} else if ( ( (l1 == "REMARK") && (l2 == "Neural") && (l3 == "network")
			 && (l4 == "secondary") && (l5 == "structure") ) || ( (l1 == "REMARK") && (l2 == "TALOS-N") ) ) {
			read_talos_ss( file_name, prediction_name);
		} else {
			utility_exit_with_message( "Can't identify secondary structure file type (needs vertical psipred_ss2 or talos+ pred.ss): "+file_name );
		}
}

void FragmentPicker::read_psipred_ss2(std::string const & file_name,
		std::string prediction_name) {

	core::fragment::SecondaryStructureOP ss_profile( new core::fragment::SecondaryStructure() );
	ss_profile->read_psipred_ss2(file_name);

	std::string query_ss_as_string;
	for (Size i = 1; i <= ss_profile->total_residue(); i++)
		query_ss_as_string += ss_profile->secstruct(i);

	query_ss_as_string_[prediction_name] = query_ss_as_string;
	query_ss_profile_[prediction_name] = ss_profile;
}

void FragmentPicker::read_talos_ss(std::string const & file_name,
		std::string prediction_name) {

	core::fragment::SecondaryStructureOP ss_profile( new core::fragment::SecondaryStructure() );
	ss_profile->read_talos_ss(file_name);
	//TALOS files can be shortened if there is no data at end of sequence. Fill up with 1/3 1/3 1/3 propensities until end is reached
	ss_profile->extend(query_profile_->length());
	for ( Size pos = ss_profile->total_residue()+1; pos <= query_profile_->length(); pos++ ) {
		ss_profile->set_fractions( pos, 1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0 );
	}

	std::string query_ss_as_string;
	for (Size i = 1; i <= ss_profile->total_residue(); i++)
		query_ss_as_string += ss_profile->secstruct(i);

	query_ss_as_string_[prediction_name] = query_ss_as_string;
	query_ss_profile_[prediction_name] = ss_profile;
}

void FragmentPicker::read_depth(std::string const & file_name) {

	utility::io::izstream data( file_name.c_str() );
	if ( !data ) {
		data.close();
		utility_exit_with_message( "Can't read DEPTH file: "+file_name );
	}
	std::string line, jnk, aathree;
	core::Real depth;
	query_residue_depth_.clear();
	getline( data, line ); // skip header
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
//# chain:residue all-atom        all-atom(stdev) MC-atom MC-atom(stdev)  SC-atom SC-atom(stdev)  SC-polar-atom   SC-polar-atom(stdev)    SC-nonpolar-atom        SC-nonpolar-atom(stdev)
		line_stream >> jnk >> aathree >> depth;
		if (aathree != "UNK")
			query_residue_depth_.push_back( depth );
		if ( line_stream.fail() )
			utility_exit_with_message( "Error reading in FragmentPicker::read_depth()!" );
	}
	data.close();
	if (query_residue_depth_.size() != size_of_query())
		utility_exit_with_message( "Error reading in FragmentPicker::read_depth(): does not match size of query!" );
}

void FragmentPicker::read_spine_x(std::string const & file_name) {

	utility::io::izstream data( file_name.c_str() );
	if ( !data ) {
		data.close();
		utility_exit_with_message( "Can't read spine-x file: "+file_name );
	}
	std::string line, jnk;
	core::Real phi, psi, asa, pkc_phi, pkc_psi;
	query_sa_prediction_.clear();
	query_phi_prediction_.clear();
	query_psi_prediction_.clear();
	getline( data, line ); // skip header
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
//#                index   AA    SS     phi1  psi1    P_E    P_C    P_H    phi0   psi0   ASA    S_pk   S_SS  pk_phi pk_psi  pkc_phi    pkc_ps
		line_stream >> jnk >> jnk >> jnk >> phi >> psi >> jnk >> jnk >> jnk >> jnk >> jnk >> asa >> jnk >> jnk >> jnk >> jnk >> pkc_phi >> pkc_psi;
		query_sa_prediction_.push_back( asa );
		query_phi_prediction_.push_back( phi );
		query_psi_prediction_.push_back( psi );
		query_phi_prediction_conf_.push_back( pkc_phi );
		query_psi_prediction_conf_.push_back( pkc_psi );
		if ( line_stream.fail() )
			utility_exit_with_message( "Error reading in FragmentPicker::read_spine_x()!" );
	}
	data.close();
	if (query_sa_prediction_.size() != size_of_query())
		utility_exit_with_message( "Error reading in FragmentPicker::read_spine_x(): does not match size of query!" );
}

void FragmentPicker::add_query_ss(std::string query_secondary,
		std::string prediction_name) {

	core::fragment::SecondaryStructureOP ss_profile( new core::fragment::SecondaryStructure() );
	ss_profile->extend(query_secondary.length());

	for (Size i = 1; i <= query_secondary.length(); ++i) {
		char ss = query_secondary[i - 1];
		if (ss == 'E')
			ss_profile->set_fractions(i, 0.0, 1.0, 0.0);
		else if (ss == 'L')
			ss_profile->set_fractions(i, 0.0, 0.0, 1.0);
		else
			ss_profile->set_fractions(i, 1.0, 0.0, 0.0);
	}
	query_ss_as_string_[prediction_name] = query_secondary;
	query_ss_profile_[prediction_name] = ss_profile;
}

void FragmentPicker::save_fragments() {
	using namespace ObjexxFCL;
	tr.Info << "Saving Fragments..." << std::endl;
	const bool skip_merge = (candidates_sinks_.size() == 1) ? true : false;
	tr.Debug << "skip_merge: " << ( skip_merge ? "true" : "false" ) << std::endl;
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
		Size fragment_size = frag_sizes_[iFragSize];
		Size maxqpos = size_of_query() - fragment_size + 1;

		utility::vector1<Candidates> final_fragments(maxqpos); // final fragments

		for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) {
			Size qPos = query_positions_[iqpos];
			if ( qPos > maxqpos) continue;
			tr.Debug << "saving " << fragment_size << "mers for position..." << qPos << std::endl;
			if (!skip_merge) {  // merge candidates
				for (Size i=1;i<=candidates_sinks_.size();++i)
					candidates_sink_[fragment_size]->insert(qPos, candidates_sinks_[i][fragment_size]);
			}
			Candidates in, out;
			if (skip_merge) {
				in = candidates_sinks_[1][fragment_size]->get_candidates(qPos);
			} else {
				in = candidates_sink_[fragment_size]->get_candidates(qPos);
			}
			if ( in.size() == 0 ) continue;
			selector_->select_fragments(in, out);
			final_fragments[qPos] = out;
		}
		tr.Debug << "call output_fragments now: " << std::endl;
		output_fragments( fragment_size, final_fragments );
		if (skip_merge) {
			//			tr.Debug << "write report..." << std::endl;
			//		candidates_sinks_[1][fragment_size]->print_report(tr.Info, get_score_manager());
		} else {
			//			candidates_sink_[fragment_size]->print_report(tr.Info, get_score_manager());
		}
	}
}

void FragmentPicker::save_candidates() {
	using namespace ObjexxFCL;
	const bool skip_merge = (candidates_sinks_.size() == 1) ? true : false;
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
		Size fragment_size = frag_sizes_[iFragSize];
		Size maxqpos = size_of_query() - fragment_size + 1;
		std::string out_file_name = prefix_ + "." + string_of(fragment_size)
				+ "mers";
		utility::io::ozstream output(out_file_name);

		for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) {
			Size qPos = query_positions_[iqpos];
			if ( qPos > maxqpos) continue;

			if (!skip_merge) {  // merge candidates
				for (Size i=1;i<=candidates_sinks_.size();++i)
					candidates_sink_[fragment_size]->insert(qPos, candidates_sinks_[i][fragment_size]);
			}
			Candidates out;
			if (skip_merge) {
				out = candidates_sinks_[1][fragment_size]->get_candidates(qPos);
			} else {
				out = candidates_sink_[fragment_size]->get_candidates(qPos);
			}
			if (out.size() == 0) continue;
			output << "position: " << I(12, qPos) << " neighbors:   " << I(10,
				out.size()) << std::endl << std::endl;
			for (Size fi = 1; fi <= out.size(); ++fi) {
				out[fi].first->print_fragment(output);
				output << std::endl;
			}
		}
		if (skip_merge) {
			candidates_sinks_[1][fragment_size]->print_report(tr.Debug, get_score_manager());
		} else {
			candidates_sink_[fragment_size]->print_report(tr.Debug, get_score_manager());
		}
		output.close();
	}
}

void FragmentPicker::pick_candidates(Size i_pos,Size frag_len) {
	if (candidates_sinks_.size() > 1)
		utility_exit_with_message( "pick_candidates(Size i_pos,Size frag_len) does not support multiple CandidateCollectors" );
	scores::FragmentScoreMapOP empty_map = scores_[1]->create_empty_map();
	for (Size i = 1; i <= chunks_->size(); i++) { // loop over provided chunks
		VallChunkOP chunk = chunks_->at(i); // For each chunk from a provider...
		if (!is_valid_chunk( frag_len, chunk )) continue;

		tr.Trace << "Processing sequence from vall: " << chunk->get_sequence() << std::endl;

		CandidatesCollectorOP sink = candidates_sinks_[1][frag_len];

		// split chunk into fragment candidates and score them
		for (Size j = 1; j <= chunk->size() - frag_len + 1; j++) {
			FragmentCandidateOP f( new FragmentCandidate(i_pos, j,chunk, frag_len) );
			if (scores_[1]->score_fragment(f, empty_map)) {
				scores::FragmentScoreMapOP new_map = empty_map->clone();
				std::pair<FragmentCandidateOP, scores::FragmentScoreMapOP> p(f, new_map);
				sink->add(p);
			}
		} // All chunk locations done
		tr.Trace << chunk->get_pdb_id() << " done" << std::endl;
		tr.Trace << sink->count_candidates()<<" candidates stored at pos. "
				<<i_pos<<", "<<sink->count_candidates()<<" in total"<< std::endl;
		tr.flush();
	} // all chunks done
}

// called in main
void FragmentPicker::parse_command_line() {

#if (defined MULTI_THREADED && defined CXX11) || defined USE_BOOST_THREAD
	//## multi-threaded?
	if (option[ frags::j ].user()) max_threads_ = option[ frags::j ]();
#endif
	// score with multiple threads
	while (max_threads_ > scores_.size())
		scores_.push_back(scores::FragmentScoreManagerOP( new scores::FragmentScoreManager() ));
	while (max_threads_ > candidates_sinks_.size()) {
		CandidatesSink storage;
		candidates_sinks_.push_back(storage);
	}

	//## -------- setup query profile
	if (option[in::file::checkpoint].user()) {
		core::sequence::SequenceProfileOP q_prof( new core::sequence::SequenceProfile );
		tr.Info << "reading a query profile from: "
			<< option[in::file::checkpoint]() << std::endl;
		q_prof->read_from_checkpoint(option[in::file::checkpoint]());
		set_query_seq(q_prof);
		tr.Info << "picking fragments for query profile: "
			<< get_query_seq_string() << std::endl;
	}
	//## -------- setup query profile : legacy blast's binary checkpoint file
	if (option[in::file::binary_chk].user()) {
		core::sequence::SequenceProfileOP q_prof( new core::sequence::SequenceProfile );
		tr.Info << "reading a query profile from: "
			<< option[in::file::binary_chk]() << std::endl;
		q_prof->read_from_binary_chk(option[in::file::binary_chk]());
		set_query_seq(q_prof);
		tr.Info << "picking fragments for query profile: "
			<< get_query_seq_string() << std::endl;
	}
	if (option[in::file::pssm].user()) {
		core::sequence::SequenceProfileOP q_prof( new core::sequence::SequenceProfile );
		tr.Info << "reading a query profile from: "
			<< option[in::file::pssm]()[1] << std::endl;
		q_prof->read_from_file(option[in::file::pssm]()[1] );
		q_prof->convert_profile_to_probs(1.0); // was previously implicit in read_from_file()
		set_query_seq(q_prof);
		tr.Info << "picking fragments for query profile: "
			<< get_query_seq_string() << std::endl;
	}

	//Fasta file trumps sequence profile as far as query_seq_string_ is concerned
	if (option[in::file::fasta].user()) {
		std::string q_seq = core::sequence::read_fasta_file(option[in::file::fasta]()[1])[1]->sequence();
		tr.Info << "reading a query sequence from: "
			<< option[in::file::fasta]()[1] << std::endl;

		set_query_seq(q_seq);
		tr.Info << "picking fragments for query sequence: "
			<< get_query_seq_string() << std::endl;
	}

	// --------- setup query secondary structure
	if (option[frags::ss_pred].user()) {
		utility::vector1<std::string> sec_str_input(option[frags::ss_pred]());
		read_ss_files(sec_str_input);
	}

	// --------- setup query phi,psi,sa predictions from spine-x file
	if (option[frags::spine_x].user()) {
		read_spine_x(option[frags::spine_x]());
	}

	// --------- setup query residue depth values from DEPTH file
	if (option[frags::depth].user()) {
		read_depth(option[frags::depth]());
	}

	//---------- setup chunk filters
	if (option[frags::allowed_pdb].user()) {
		AllowPdbIdFilterOP allow( new AllowPdbIdFilter() );
		allow->load_pdb_id_from_file(option[frags::allowed_pdb]());
		add_chunk_filter(allow);
		tr.Info << "Allowed PDB chains:\n";
		allow->show_pdb_ids(tr.Info);
	}

	if (option[frags::denied_pdb].user()) {
		DenyPdbIdFilterOP deny( new DenyPdbIdFilter() );
		deny->load_pdb_id_from_file(option[frags::denied_pdb]());
		add_chunk_filter(deny);
		tr.Info << "Excluded PDBs:\n";
		deny->show_pdb_ids(tr.Info);
	}

	// ##--------- setup VALL
	PROF_START( basic::FRAGMENTPICKING_READ_VALL );
	if (option[in::file::vall].user()) {
		read_vall(option[in::file::vall]());
	}
	PROF_STOP( basic::FRAGMENTPICKING_READ_VALL );

	// -------- fragment sizes
	if (option[frags::frag_sizes].user()) {
		utility::vector1<Size> frag_sizes_tmp = option[frags::frag_sizes]();
		for (Size i = 1; i <= frag_sizes_tmp.size(); ++i) {
				if(frag_sizes_tmp[i] > max_frag_size_)
					max_frag_size_ = frag_sizes_tmp[i];
				frag_sizes_.push_back(frag_sizes_tmp[i]);
		}
	} else {
		max_frag_size_ = 9;
		frag_sizes_.push_back(3);
		frag_sizes_.push_back(9);
	}
	tr.Info << "Will pick fragments of size:";
	for (Size i = 1; i <= frag_sizes_.size(); ++i)
		tr.Info << frag_sizes_[i] << " ";
	tr.Info << std::endl;

	//---------- setup scoring scheme
	tr.Info << "Creating fragment scoring scheme" << std::endl;
	if (option[frags::scoring::config].user()) {
		// todo:  the create scores method should be improved so the score files aren't read more than once -dk
		for (Size i = 1; i <= scores_.size(); ++i)
			scores_[i]->create_scores(option[frags::scoring::config](), get_self_ptr());
	}

	// -------- how many fragments and candidates
	n_frags_ = option[frags::n_frags]();
	n_candidates_ = option[frags::n_candidates]();

	if (n_frags_ > n_candidates_) n_candidates_ = n_frags_;

	tr.Info << "Picking " << n_frags_ << " fragments based on "
		<< n_candidates_ << " candidates" << std::endl;

	//-------- this comparator is used both for collecting and selecting fragments
	// note: The comparator is based on the first score manager so the score managers have to have the same scoring scheme! -dk
	CompareTotalScore comparator(get_score_manager());

	//---------- setup scoring scheme for the selection step
	tr.Info << "Creating fragment scoring scheme for the selection step" << std::endl;
	FragmentScoreManagerOP selection_scoring;
	if (option[frags::picking::selecting_scorefxn].user()) {
		selection_scoring = FragmentScoreManagerOP( new FragmentScoreManager() );
		selection_scoring->create_scores(option[frags::picking::selecting_scorefxn](), get_self_ptr());
		selector_ = FragmentSelectingRuleOP( new CustomScoreSelector(n_frags_, selection_scoring) );
	} else {
		// note: The selector is based on the first score manager so the score managers have to have the same scoring scheme! -dk
		selector_ = FragmentSelectingRuleOP( new BestTotalScoreSelector(n_frags_, get_score_manager()) );
	}

	//-------- collector & selector set up
	if (option[frags::quota_protocol].user() || option[frags::picking::quota_config_file].user()) {
	// This setup is a bit more complicated when user needs quota.
	// The quota version of this code was moved into a separate method
		parse_quota_command_line();
	// This setup is a bit more complicated, when user needs quota. The quota version of this code was moved into a separate method
	} else {
		if (option[frags::keep_all_protocol].user()) {
			for (Size i = 1; i <= frag_sizes_.size(); ++i) {
				CandidatesCollectorOP collector( new GrabAllCollector(size_of_query()) );
				set_candidates_collector(frag_sizes_[i], collector);
				tr.Info << "Collector for fragment size: " << frag_sizes_[i] << " set to: GrabAllCollector" << std::endl;
			}
		} else {
			for (Size i = 1; i <= frag_sizes_.size(); ++i) {
				for (Size j = 0; j <= max_threads_; ++j) {  // 0 for merged collector
					CandidatesCollectorOP collector( new BoundedCollector<CompareTotalScore> (size_of_query(), n_candidates_,
						comparator,get_score_manager()->count_components()) );
					set_candidates_collector(frag_sizes_[i], collector, j);
				}
				tr.Info << "Collector for fragment size: " << frag_sizes_[i] << " set to: BoundedCollector" << std::endl;
			}
		}
		//-------- Selecting fragments from candidates
/*			if (option[frags::picking::selecting_rule].user()) {
		std::string type = option[frags::picking::selecting_rule]();
		if (type.compare("BestTotalScoreSelector")==0) {
			selector_ = new BestTotalScoreSelector(n_frags_, selection_scoring);
			tr.Info << "Fragment selector: BestTotalScoreSelector"
					<< std::endl;
		} else {
				utility_exit_with_message("[ERROR]: unknown fragment selecting rule: " + type + "!");
		}
			} else {
		selector_ = new BestTotalScoreSelector(n_frags_, selection_scoring);
		tr.Info << "Fragment selector: BestTotalScoreSelector" << std::endl;
					}*/
	}
	// # ---------- output file prefix:
	if (option[out::file::frag_prefix].user()) {
		prefix_ = option[out::file::frag_prefix]();
	}

	if (option[frags::picking::query_pos].user()) {
		set_picked_positions( option[frags::picking::query_pos]() );
	}

	// non-local contacts options
	nonlocal_min_contacts_per_res_ = option[ frags::nonlocal::min_contacts_per_res ]();

	// frag contacts options
	contact_types_.clear();
	utility::vector1<std::string> contact_types = option[ frags::contacts::type ]();
	for (Size i = 1; i <= contact_types.size(); ++i) {
		contact_types_.insert(contact_type(contact_types[i]));
		if (contact_type(contact_types[i]) == CEN)
			sidechain_contact_dist_cutoff_ = SidechainContactDistCutoffOP( new SidechainContactDistCutoff( option[ frags::contacts::centroid_distance_scale_factor ]() ) );
	}
	// sequence separation
	contacts_min_seq_sep_ = option[ frags::contacts::min_seq_sep ](); // j>=i+contacts_min_seq_sep_
	// distance cutoffs
	utility::vector1<Real> dist_cutoffs = option[ frags::contacts::dist_cutoffs ]();
	Real max_dist = 0.0;
	for (Size i = 1; i <= dist_cutoffs.size(); ++i) {
		if (dist_cutoffs[i] > max_dist) max_dist = dist_cutoffs[i];
		contacts_dist_cutoffs_squared_.push_back( dist_cutoffs[i]*dist_cutoffs[i] );
	}
	contacts_dist_cutoff_squared_ = max_dist*max_dist;

	show_scoring_methods(tr);
	tr << std::endl;
}

/// @brief sets the query sequence
/// @detailed Well, it is a sequence profile, but the sequence can be extracted from it
void FragmentPicker::set_query_seq(core::sequence::SequenceProfileOP query_sequence) {
	query_profile_ = query_sequence;
	query_seq_as_string_ = query_profile_->sequence();
	set_picked_positions(1,query_sequence->length());
}

/// @brief sets the query sequence
void FragmentPicker::set_query_seq(std::string & query_sequence) {
	if (query_profile_ == 0) {
		query_profile_ = core::sequence::SequenceProfileOP( new core::sequence::SequenceProfile() );
		tr.Warning << "CAUTION: No sequence profile supplied. Profile-dependant options/scoring will not work." << std::endl;
	}
	query_profile_->sequence(query_sequence);
	query_seq_as_string_ = query_sequence;
	set_picked_positions(1,query_profile_->length());
}

void FragmentPicker::set_up_ss_abego_quota() {

	std::string quota_config_file("UNKNOWN-QUOTA-CONFIG_FILE");
	if (option[frags::picking::quota_config_file].user())
		quota_config_file = option[frags::picking::quota_config_file]();
	quota::ABEGO_SS_Config q_config(quota_config_file);

	utility::vector1<Size> components;
	utility::vector1<Real> weights;
	utility::vector1<Real> scoring_weights = scores_[1]->get_weights();
	for (Size i = 1; i <= scores_[1]->count_components(); ++i) {
		ABEGO_SS_ScoreOP s0 =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::ABEGO_SS_Score > ( scores_[1]->get_component(i) );
		if (s0 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s0->get_id()] );
		}
		ProfileScoreL1OP s1 =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::ProfileScoreL1 > ( scores_[1]->get_component(i) );
		if (s1 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s1->get_id()] );
		}

		RamaScoreOP s2 =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::RamaScore > ( scores_[1]->get_component(i) );
		if (s2 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s2->get_id()] );
		}

		CSScoreOP s3 =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::CSScore > ( scores_[1]->get_component(i) );
		if (s3 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s3->get_id()] );
		}

		SecondarySimilarityOP s4 =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::SecondarySimilarity > ( scores_[1]->get_component(i) );
		if (s4 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s4->get_id()] );
		}
	}

	tr.Debug<<"Scoring scheme for ABEGO_SS quota pool sorting is:";
	for(Size l=1;l<=weights.size();l++) {
			tr.Debug<<"\n\t"<<components[l]<<"\t"<<weights[l];
	}
	tr.Debug<<std::endl;
	Size buffer_factor = 5;
	for(Size f=1;f<=frag_sizes_.size();f++) {
		for (Size j = 0; j <= max_threads_; ++j) { // 0 for the merged collector
			quota::QuotaCollectorOP collector( new quota::QuotaCollector( size_of_query(), frag_sizes_[f] ) );
			set_candidates_collector(frag_sizes_[f],collector, j);
		}
		Size middle = frag_sizes_[f] / 2 + 1;
		assert( size_of_query() == q_config.size() ); // Test if the abego-ss table has the same size as the query sequence
		for(Size j=1;j<=size_of_query()-frag_sizes_[f]+1;j++) {
			tr.Debug<<"Creating "<<q_config.n_columns()<<" quota pools at pos "<<j<<std::endl;
			for(Size i=1;i<=q_config.n_columns();i++) {
				Real prob = q_config.probability(j+middle-1,i);
				for (Size k = 0; k <= max_threads_; ++k) {  // 0 for the merged collector
					CandidatesCollectorOP storage = get_candidates_collector(frag_sizes_[f], k);
					quota::QuotaCollectorOP collector = utility::pointer::dynamic_pointer_cast< quota::QuotaCollector > ( storage );
					quota::QuotaPoolOP p( new quota::ABEGO_SS_Pool(n_candidates_,q_config.get_pool_name(i),
							q_config.get_pool_bins((i)),components,weights,prob,scores_[1]->count_components(),buffer_factor) );
					collector->add_pool(j,p);
				}
			}
		}
	}
}

void FragmentPicker::set_up_quota_nnmake_style() {
	std::string quota_config_file("UNKNOWN-QUOTA-CONFIG_FILE");
	if (option[frags::picking::quota_config_file].user())
		quota_config_file = option[frags::picking::quota_config_file]();
	quota::QuotaConfig q_config(quota_config_file);

	utility::vector1<Size> components;
	utility::vector1<Real> weights;
	components.push_back( 0 );		// the free entry in the vector is for secondary structure score (only one for each pool)
	weights.push_back( 0.0 );		// score weight for SecondarySimilarity; will be changed later
	components.push_back( 0 );		// this free entry in the vector is for RamaScore
	weights.push_back( 0.0 );		// score weight for RamaScore; will be changed later
	utility::vector1<Real> scoring_weights = scores_[1]->get_weights();
	for (Size i = 1; i <= scores_[1]->count_components(); ++i) {
		ProfileScoreL1OP s1 =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::ProfileScoreL1 > ( scores_[1]->get_component(i) );
		if (s1 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s1->get_id()] );
		}

/****** RamaScore is a special case, dispatched below
		RamaScore *s2 =
			dynamic_cast<RamaScore*> (scores_->get_component(i).get());
		if (s2 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s2->get_id()] );
		}
*********/
		CSScoreOP s3 =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::CSScore > ( scores_[1]->get_component(i) );
		if (s3 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s3->get_id()] );
		}
		ABEGO_SS_ScoreOP s4 =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::ABEGO_SS_Score > ( scores_[1]->get_component(i) );
		if (s4 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s4->get_id()] );
		}
		TorsionBinSimilarityOP s5 =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::TorsionBinSimilarity > ( scores_[1]->get_component(i) );
		if (s5 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s5->get_id()] );
		}
	}

	utility::vector1<core::fragment::SecondaryStructureOP> predictions;
	utility::vector1<Real> ss_weights;
	std::map<std::string, core::fragment::SecondaryStructureOP>::iterator it;
	Real weight = 1.0 / ((Real) query_ss_profile_.size());
	for ( it=query_ss_profile_.begin() ; it != query_ss_profile_.end(); it++ ) {
			predictions.push_back((*it).second);
			ss_weights.push_back(weight);
	}
//	core::fragment::SecondaryStructureOP avg_ss = new core::fragment::SecondaryStructure(predictions,ss_weights);

	for(Size f=1;f<=frag_sizes_.size();f++) {
		for (Size j = 0; j <= max_threads_; ++j) { // 0 for the merged collector
			quota::QuotaCollectorOP collector( new quota::QuotaCollector( size_of_query(), frag_sizes_[f] ) );
			set_candidates_collector(frag_sizes_[f], collector, j);
		}
// --------- This part puts RamaScore into quota scoring; each Rama is based on a certain SS prediction and this part of the code
// --------- dispatches each Rama into a proper pool
		for (Size i = 1; i <= scores_[1]->count_components(); ++i) {
			RamaScoreOP sr = utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::RamaScore > ( scores_[1]->get_component(i) );
			if (sr != 0) {
				std::string & name = sr->get_prediction_name();
				if( ! q_config.is_valid_quota_pool_name( name ) ) continue;
				components[2] = i;
				weights[2] = scoring_weights[sr->get_id()];
				tr.Warning<<"RamaScore with ID "<<sr->get_id()<<" named "<<name<<
					" has been attached to its quota pool with weight "<<weights[2]<<std::endl;
			}
		}

// ---------- end of RamaScore dispatch

// Create secondary structure pools (if any)
		for (Size i = 1; i <= scores_[1]->count_components(); ++i) {
			SecondarySimilarityOP ss = utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::SecondarySimilarity > ( scores_[1]->get_component(i) );

		//PartialSecondarySimilarity is a variant of SecondarySimilarity, this means they're not compatible
		//So what is it compatible with???
//		if (s == 0) {
//			PartialSecondarySimilarity *s =
//				dynamic_cast<PartialSecondarySimilarity*> (scores_->get_component(i).get());
//		}

			if (ss != 0) {
				std::string & name = ss->get_prediction_name();
				if( ! q_config.is_valid_quota_pool_name( name ) ) continue;
				components[1] = i;
				weights[1] = scoring_weights[ss->get_id()];
				Size size = (Size)(q_config.get_fraction( name ) * n_candidates_);
				if( size == 0 ) {
					tr.Warning<<"Config file couldn't provide quota fraction for the pool named "
						<<name<<". Skipping the pool"<<std::endl;
					continue;
				}

				for (Size j = 0; j <= max_threads_; ++j) { // 0 for the merged collector
					CandidatesCollectorOP storage = get_candidates_collector(frag_sizes_[f], j);
					quota::QuotaCollectorOP collector = utility::pointer::dynamic_pointer_cast< quota::QuotaCollector > ( storage );

					core::fragment::SecondaryStructureOP ss_prediction( get_query_ss( name ) );
					if( ! ss_prediction ) {
						utility_exit_with_message("Unable to get secondary structure prediction for " + name);
					}
					collector->attach_secondary_structure_pools(q_config.get_fraction( name ) ,
					ss_prediction,name,n_candidates_,components,weights,scores_[1]->count_components());
//					avg_ss,name,n_candidates_,components,weights,scores_->count_components());
				}
			}
		}
	}
}

void FragmentPicker::parse_quota_command_line() {

	set_up_quota_nnmake_style();

	if (option[frags::picking::query_pos].user()) {
			set_picked_positions( option[frags::picking::query_pos]() );
	}
}

void FragmentPicker::read_vall( utility::vector1< std::string > const & fns ) {
	chunks_ = VallProviderOP( new VallProvider() );
	chunks_->vallChunksFromLibraries(fns);
}

void FragmentPicker::read_vall( std::string const & fn ) {
	chunks_ = VallProviderOP( new VallProvider() );
	chunks_->vallChunksFromLibrary(fn);
}

void FragmentPicker::set_picked_positions(Size from,Size to) {
	query_positions_.clear();
	for(Size i=from;i<=to;i++)
		query_positions_.push_back( i );
}

void FragmentPicker::set_picked_positions(utility::vector1<Size> q_positions) {
	query_positions_.clear();
	for(Size i=1;i<=q_positions.size();i++)
		query_positions_.push_back( q_positions[i] );
}


Size QuotaDebug::max_pools() {
	return tags_.size();
}

void QuotaDebug::write_summary() {
	tr<< "Quota report: difference between the total expected and total picked foreach pool"<<std::endl;
	tr<< "This table is for first "<<nFrags_<<" fragments"<<std::endl;
	tr<< "Negative value says that was picked more that expected."<<std::endl;
	tr<< this->str()<<std::endl;
	this->str("");
}

void QuotaDebug::log(Size frag_len,Size q_pos,utility::vector1<Real> data) {
	*this  << std::setw(4)<<q_pos<< std::setw(4)<<frag_len;
	for(Size i=1;i<=data.size();i++)
		if(data[i]<1000000)
			*this  << std::setw(10)<<std::setprecision(3)<<data[i];
		else
			*this  << std::setw(10)<<"  --- ";
	*this  << std::endl;
}

void QuotaDebug::setup_summary(
	quota::QuotaCollector const & collector
)
{
	Size last_tag = 0;
	for(Size i=1;i<=collector.query_length();++i) {
		for(Size j=1;j<=collector.count_pools(i);++j) {
			if( tag_map_.find(collector.get_pool(i,j)->get_pool_name())==tag_map_.end() ) {
				tags_.push_back(collector.get_pool(i,j)->get_pool_name());
				last_tag++;
				tag_map_[collector.get_pool(i,j)->get_pool_name()] = last_tag;
			}
		}
	}

	*this <<"\n#len pos ";
	for(Size i=1;i<=tags_.size();i++) {
		*this << std::setw(10)<<tags_[i];
	}
	*this<<std::endl;
}

utility::vector1<ConstantLengthFragSetOP> FragmentPicker::getFragSet(int residueInPose_){

	// Storage for the resulting FragSets
	utility::vector1<ConstantLengthFragSetOP> result;

	// Loop over various sizes of fragments
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) {

		ConstantLengthFragSetOP myFragSet( new ConstantLengthFragSet() );

		Size fragment_size = frag_sizes_[iFragSize];
		CandidatesCollectorOP storage = get_candidates_collector(fragment_size);
		for (Size qPos = 1; qPos <= size_of_query(); ++qPos) {
			if(storage->get_candidates(qPos).size() == 0) continue;

			Candidates out;
			selector_->select_fragments(storage->get_candidates(qPos), out);

			//FrameOP frame = new Frame(out[1].first->get_residue(1)->resi()); // start pos = residue id from the fragment residue (sequence from original pdb file?) or is this some internal index?

			/*
				start pos = residue id from the fragment residue (sequence from original pdb file?) or is this some internal index?

			 In  ConstantLengthFragSet, this is set from insertion_pos of fragment file.

			 Therefore,  I think this is the connection between Pose and Fragment File. In design mode,
			 we do not have this, therefore we just add an incremented index for each qPos
			*/
			FrameOP frame( new Frame(residueInPose_++) );

			for (Size fi = 1; fi <= out.size(); ++fi) {

						FragDataOP current_fragment( NULL );

						for (Size i = 1; i <= out[1].first->get_length(); ++i) {
							VallResidueOP r   =  out[fi].first->get_residue(i);
							string pdbid      = out[fi].first->get_pdb_id();
							//char chainid      = out[fi].first->get_chain_id();
							Size index        = r->resi();
							char aa           = toupper(r->aa());
							char ss           = r->ss();
							Real phi          = r->phi();
							Real psi          = r->psi();
							Real omega        = r->omega();

							if (i == 1){
								current_fragment = FragDataOP( new AnnotatedFragData( pdbid, index ) );
							}
							utility::pointer::shared_ptr< BBTorsionSRFD > res_torsions( new BBTorsionSRFD(3,ss,aa) ); // 3 protein torsions
							res_torsions->set_torsion   ( 1, phi   ); // ugly numbers 1-3, but pose.set_phi also uses explicit numbers
							res_torsions->set_torsion   ( 2, psi   );
							res_torsions->set_torsion   ( 3, omega );
							res_torsions->set_secstruct ( ss );

							// Add residue to fragment
							current_fragment->add_residue( res_torsions );

						} // End VallResidue loop

						if (current_fragment) { // != NULL) {
							current_fragment->set_valid(); //it actually containts data

							// Add fragment to frame
							if (!frame->add_fragment(current_fragment)){
								cerr << "ERROR Bad fragment : "<<endl;
								current_fragment->show(cout);
								exit(1111);
							}
						}
			} // End FragmentCandidate loop

			// Add frame to myFragSet.
			myFragSet->add(frame);

		} // End size of query

		// For each size fragment add to vector of ConstantLengthFragSet
		result.push_back(myFragSet);

	} // End size of frags


	return (result);

}

bool FragmentPicker::is_valid_chunk( VallChunkOP chunk ) {
	bool flag = true;
	for (Size iFilter = 1; iFilter <= filters_.size(); iFilter++) {
		if ((flag = filters_[iFilter]->test_chunk(chunk)) == false) {
			tr.Debug << "Chunk: " << chunk->get_pdb_id()
				<< " didn't pass a filter" << std::endl;
			break;
		}
	}
	return flag;
}

bool FragmentPicker::is_valid_chunk( Size const frag_len, VallChunkOP chunk ) {
	if (chunk->size() < frag_len) return false; // This fragment is too short
	return is_valid_chunk( chunk );
}

// Output fragments
void FragmentPicker::output_fragments( Size const fragment_size, utility::vector1<Candidates> const & final_fragments ) {
	using namespace ObjexxFCL;

	// find and output nonlocal pairs
	if (option[frags::nonlocal_pairs].user())
		nonlocal_pairs( fragment_size, final_fragments );

	// find and output fragment contacts
	if (option[frags::fragment_contacts].user())
		fragment_contacts( fragment_size, final_fragments );


	std::string out_file_name = prefix_ + "." + string_of(n_frags_) + "." + string_of(fragment_size) + "mers";
	std::string silent_out_file_name = out_file_name + ".out";
	utility::io::ozstream output_file(out_file_name);
	utility::io::ozstream output_info_file;
	if (option[frags::describe_fragments].user()) {
		std::string describe_name = option[frags::describe_fragments]()+"." + string_of(n_frags_) +"."+string_of(fragment_size)+"mers";
		output_info_file.open(describe_name.c_str());
	}

	FragmentScoreManagerOP ms = get_score_manager();
	Size maxqpos = size_of_query() - fragment_size + 1;
	for ( Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos ) {
		Size qPos = query_positions_[iqpos];
		if ( qPos > maxqpos) continue;
		//		std::cerr << "h1 " << qPos << std::endl;
		output_file << "position: " << I(12, qPos) << " neighbors:   " << I(10, final_fragments[qPos].size()) << std::endl << std::endl;
		for (Size fi = 1; fi <= final_fragments[qPos].size(); ++fi) {
			if (option[frags::write_sequence_only]()) {
				final_fragments[qPos][fi].first->print_fragment_seq(output_file);
			} else {
				if ( !final_fragments[qPos][fi].first || !final_fragments[qPos][fi].second ) {
					tr.Warning << "final_frag candidate " << fi << " at position " << qPos << " is corrupted. skipping... " << std::endl;
					continue;
				}
				final_fragments[qPos][fi].first->print_fragment(output_file, final_fragments[qPos][fi].second, ms);
			}
			output_file << std::endl;
		}
		if ( ms->if_late_scoring_for_zeros() )  {
			//			std::cerr << "h3" << std::endl;
			for ( Size fi = 1; fi <= final_fragments[qPos].size(); ++fi ) {
				if ( !final_fragments[qPos][fi].first || !final_fragments[qPos][fi].second ) {
					tr.Warning << "final_frag candidate " << fi << " at position " << qPos << " is corrupted. skipping... " << std::endl;
					continue;
				}
				ms->score_zero_scores(final_fragments[qPos][fi].first,final_fragments[qPos][fi].second);
			}
		}
		if ( option[frags::describe_fragments].user() ) {
			//			std::cerr << "h4" << std::endl;
			ms->describe_fragments(final_fragments[qPos], output_info_file);
		}
	}
	output_file.close();
	output_info_file.close();

	// silent file output
	if (option[frags::output_silent]() || option[frags::score_output_silent]()) {
		core::io::silent::SilentFileData sfd;
		for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) {
			Size qPos = query_positions_[iqpos];
			if ( qPos > maxqpos) continue;
			std::string const & sequence = get_query_seq_string().substr(qPos-1,fragment_size);
			for (Size fi = 1; fi <= final_fragments[qPos].size(); ++fi) {
				std::string tag = "frag_" + ObjexxFCL::lead_zero_string_of(qPos,6) + "_" + ObjexxFCL::lead_zero_string_of(fi,6);
				final_fragments[qPos][fi].first->output_silent( sfd, sequence, silent_out_file_name, tag, final_fragments[qPos][fi].second, ms );
			}
		}
	}

}

QuotaDebug FragmentPicker::log_25_(25);
QuotaDebug FragmentPicker::log_200_(200);

} // frag_picker
} // protocols
