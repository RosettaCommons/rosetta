// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/legacy_sewing/conformation/Assembly.cc
///
/// @author Tim Jacobs
/// @author Benfeard Williams

//Unit headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>
#include <protocols/legacy_sewing/sampling/SewGraph.fwd.hh>

//Package headers
#include <protocols/legacy_sewing/util/io.hh>
#include <protocols/legacy_sewing/util/util.hh>

//Protocol headers
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/AddResiduesRotamerSetOperation.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <core/conformation/Residue.functions.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/types.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>
#include <core/pose/motif/reference_frames.hh>

#include <protocols/loops/Loop.hh>
#include <core/select/util/SelectResiduesByLayer.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/legacy_sewing.OptionKeys.gen.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/xyzTransform.io.hh>

#include <utility/LexicographicalIterator.hh>
#include <utility/vector1.hh>


//temporary
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>


namespace protocols {
namespace legacy_sewing  {

static basic::Tracer TR("protocols.legacy_sewing.Assembly");

Assembly::Assembly()= default;

///@brief Only used to add the first model to the Assembly, subsequent models
///should be added using the follow_edge function
void
Assembly::add_model(
	SewGraphCOP graph,
	Model const & model,
	bool available /*= true*/
){
	core::Size pre_segments_size = segments_.size();

	for ( core::Size i=1; i<=model.segments_.size(); ++i ) {
		segments_.push_back(model.segments_[i]);
		utility::vector1<SewSegment> seg_vec;
		seg_vec.push_back(model.segments_[i]);
		all_segments_.push_back(seg_vec);
	}

	//copy connections
	std::map<core::Size, core::Size> forward_connections = model.segments_.forward_connections();
	std::map<core::Size, core::Size>::const_iterator f_it = forward_connections.begin();
	std::map<core::Size, core::Size>::const_iterator f_end = forward_connections.end();
	for ( ; f_it != f_end; ++f_it ) {
		segments_.add_connection(f_it->first + pre_segments_size, f_it->second + pre_segments_size);
	}

	if ( available ) {
		std::set<core::Size> node_indices = graph->get_node_indices_from_model_id(model.model_id_);
		available_nodes_.insert(node_indices.begin(), node_indices.end());
	}
}


ModelNode const *
Assembly::starting_node() const {
	return edge_path_.front().get<0>();
}

ModelNode const *
Assembly::ending_node() const {
	return edge_path_.back().get<1>();
}

utility::vector1< boost::tuple<ModelNode const *, ModelNode const *, HashEdge const *> >
Assembly::edges() const{
	return edge_path_;
}


void
Assembly::update_coords_from_pose(
	core::pose::Pose const & pose
) {

	ModelIterator<SewSegment> it = assembly_begin();
	ModelIterator<SewSegment> it_end = assembly_end();
	for ( core::Size i=1; i<=pose.size(); ++i ) {

		it.residue()->chi_angles_ = pose.residue(i).chi();
		for ( Size ii=1; ii <= 4; ++ii ) {
			it.atom()->coords_ = pose.residue(i).xyz(ii);
			++it;
		}

	}
	runtime_assert(it == it_end);
}

ScoreResult
Assembly::get_score_result(
	ModelNode const * model_1_node,
	HashEdge const * const cur_edge,
	SewGraphCOP graph_
) const {
	Model const reference_model = regenerate_model(model_1_node->model().model_id_);
	SewResidue const reference_residue = reference_model.get_residue( cur_edge->model_resnum( reference_model.model_id_ ) );

	//Get the model to be aligned from the graph
	Model mobile_model( graph_->get_model_node(cur_edge->get_other_ind(model_1_node->get_node_index()))->model() );
	SewResidue mobile_residue = mobile_model.get_residue( cur_edge->model_resnum( mobile_model.model_id_ ) );

	Hasher hasher;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::legacy_sewing;
	core::Size box_length = basic::options::option[ basic::options::OptionKeys::legacy_sewing::box_length ];

	ScoreResult edge_score = hasher.score_one(reference_model, reference_residue, mobile_model, mobile_residue, box_length);

	return edge_score;
}

std::map<SewSegment,SewSegment>
Assembly::get_matching_model_segments(
	Model const & mobile_model,
	ScoreResult const & edge_score
) const {
	//Somewhat confusing logic below. All that is happening is that we are adding the
	//new model's segments to the all_segments_ vector at the correct positions (the
	//same positions as the segments they are matching).

	std::map< SewSegment, SewSegment > model_segment_pairs;
	//std::set<SewSegment> matched_model_segments;
	int ref_model_id = edge_score.first.first.model_id;
	int mobile_model_id = edge_score.first.second.model_id;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Appending new model " << mobile_model_id << " using reference model " << ref_model_id << std::endl;
	}

	std::map< SegmentPair, AtomMap > segment_map = edge_score.second.segment_matches;
	std::map< SegmentPair, AtomMap >::const_iterator it = segment_map.begin();
	std::map< SegmentPair, AtomMap >::const_iterator it_end = segment_map.end();
	for ( ; it != it_end; ++it ) {
		for ( core::Size i=1; i<=all_segments_.size(); ++i ) {
			for ( core::Size j=1; j<=all_segments_[i].size(); ++j ) {
				if ( ref_model_id == all_segments_[i][j].model_id_ && it->first.first == all_segments_[i][j].segment_id_ ) {
					for ( core::Size k=1; k<=mobile_model.segments_.size(); ++k ) {
						if ( it->first.second == mobile_model.segments_[k].segment_id_ ) {
							model_segment_pairs.insert(std::make_pair(all_segments_[i][j],mobile_model.segments_[k]));
							if ( TR.Debug.visible() ) {
								TR.Debug << "Matching new segment " << mobile_model.segments_[k].segment_id_ << " to reference segment " << all_segments_[i][j].segment_id_ << std::endl;
							}
						}
					}
				}
			}
		}
	}
	runtime_assert(model_segment_pairs.size() > 0);
	return model_segment_pairs;
}



SewSegment
Assembly::create_chimera_segment(
	SewSegment const & reference_segment,
	SewSegment const & mobile_segment,
	AtomMap const & atom_map,
	bool reference_is_nter
) const {

	SewSegment chimera;
	chimera.model_id_ = CHIMERA_SEGMENT;
	chimera.chimera_ = true;
	if ( reference_segment.dssp_ != mobile_segment.dssp_ ) {
		//utility_exit_with_message("Mobile and reference segments have different dssp codes!");
		//A match between two segments with different DSSP codes tends to happen quite often when matching beta-strands, so
		//instead of failing, simply use whichever segment isn't being identified as a loop
		if ( reference_segment.dssp_ == 'L' ) {
			chimera.dssp_ = mobile_segment.dssp_;
		} else {
			chimera.dssp_ = reference_segment.dssp_;
		}
	} else {
		chimera.dssp_ = reference_segment.dssp_;
	}

	std::map<core::id::AtomID, core::id::AtomID> largest_continuous =
		largest_continuous_atom_map(atom_map);
	std::map<core::id::AtomID, core::id::AtomID> current_stretch;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Creating chimera segment " << reference_segment.model_id_ << ", " << reference_segment.segment_id_
			<< " -> " << mobile_segment.model_id_ << ", " << mobile_segment.segment_id_ << std::endl;

		TR.Debug << "Full Atom Map: " << std::endl;
		auto it = atom_map.begin();
		auto it_end = atom_map.end();
		for ( ; it != it_end; ++it ) {
			TR.Debug << "\t" << it->first.rsd() << ":" << it->first.atomno() << " -> "
				<< it->second.rsd() << ":" << it->second.atomno() << std::endl;
		}

		TR.Debug << "Using Atom Map for chimera: " << std::endl;
		it = largest_continuous.begin();
		it_end = largest_continuous.end();
		for ( ; it != it_end; ++it ) {
			TR.Debug << "\t" << it->first.rsd() << ":" << it->first.atomno() << " -> "
				<< it->second.rsd() << ":" << it->second.atomno() << std::endl;
		}
	}

	//Now pick a position in the largest continuous section to use
	//as the 'switch' position from one segment to the next
	std::pair<core::id::AtomID, core::id::AtomID> switch_position;

	//Check to see if either segment comes from is a pose model, if so, keep as much as possible
	//in the chimera. Otherwise, split the chimera in the middle of the largest segment
	if ( reference_segment.model_id_ < 0 ) {
		switch_position = *largest_continuous.rbegin();
		chimera.model_id_ = reference_segment.model_id_;
	} else if ( mobile_segment.model_id_ < 0 ) {
		switch_position = *largest_continuous.begin();
		chimera.model_id_ = mobile_segment.model_id_;
	} else {
		std::map<core::id::AtomID, core::id::AtomID>::const_iterator middle_it = largest_continuous.begin();
		std::advance(middle_it, largest_continuous.size()/2);
		switch_position = *middle_it;
	}

	//Figure out which segment is the nterm segment and which is the cterm segment
	//in the new chimera.
	SewSegment nterm_segment, cterm_segment;
	core::Size nterm_switch_res, cterm_switch_res;
	core::Size nterm_start_res, cterm_start_res;
	if ( reference_is_nter ) {
		nterm_segment = reference_segment;
		nterm_switch_res = switch_position.first.rsd();

		cterm_segment = mobile_segment;
		cterm_switch_res = switch_position.second.rsd();

		nterm_start_res = largest_continuous.begin()->first.rsd();
		cterm_start_res = largest_continuous.begin()->second.rsd();
	} else {
		nterm_segment = mobile_segment;
		nterm_switch_res = switch_position.second.rsd();

		cterm_segment = reference_segment;
		cterm_switch_res = switch_position.first.rsd();

		nterm_start_res = largest_continuous.begin()->second.rsd();
		cterm_start_res = largest_continuous.begin()->first.rsd();
	}

	//Finally, create the chimera by first adding residues from the
	//nterminal segment
	core::Size overlap_counter = 0;
	for ( core::Size i = 1; i <= nterm_segment.residues_.size(); ++i ) {
		chimera.residues_.push_back(nterm_segment.residues_[i]);
		if ( nterm_segment.residues_[i].resnum_ >= nterm_start_res ) {
			for ( core::Size j = 1; j <= cterm_segment.residues_.size(); ++j ) {
				if ( cterm_segment.residues_[j].resnum_ == cterm_start_res+overlap_counter ) {
					chimera.residues_.back().matched_residues_.push_back(cterm_segment.residues_[j]);
				}
			}
			++overlap_counter;
		}
		if ( nterm_segment.residues_[i].resnum_ == nterm_switch_res ) {
			break;
		}
	}

	//Now add in the residues from the cterminal segment
	bool found=false;
	for ( core::Size i = 1; i <= cterm_segment.residues_.size(); ++i ) {
		if ( found ) {
			chimera.residues_.push_back(cterm_segment.residues_[i]);
			for ( core::Size j = 1; j <= nterm_segment.residues_.size(); ++j ) {
				if ( nterm_segment.residues_[j].resnum_ == nterm_start_res+overlap_counter ) {
					chimera.residues_.back().matched_residues_.push_back(nterm_segment.residues_[j]);
				}
			}
			++overlap_counter;
		}
		if ( cterm_segment.residues_[i].resnum_ == cterm_switch_res ) {
			found = true;
		}
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Generated chimera has " << chimera.residues_.size() << " residues" << std::endl;
	}
	return chimera;
}


utility::vector1<SewSegment>
Assembly::get_chimera_segments(
	std::map<SewSegment, SewSegment> const & matching_segments,
	std::map<SegmentPair, AtomMap> const & segment_matches,
	Model const & mobile_model
) {
	utility::vector1<SewSegment> chimera_segs;

	bool reference_is_nter = false;
	auto it = segment_matches.begin();
	auto it_end = segment_matches.end();
	for ( ; it != it_end; ++it ) {
		auto it2 = matching_segments.begin();
		auto it2_end = matching_segments.end();

		//it2->first = reference segment in assembly
		//it2->second = mobile segment
		for ( ; it2 != it2_end; ++it2 ) {
			if ( it2->first.segment_id_ == it->first.first &&
					it2->second.segment_id_ == it->first.second ) {
				AtomMap atom_map = it->second;

				//By default, use the reference model to determine whether the reference is nterminal or cterminal to the
				//mobile model. However, if the reference is only a single segment then use the mobile model to get this info
				if ( segments_.size() > 1 ) {
					for ( core::Size i = 1; i <= segments_.size(); ++i ) {
						if ( it2->first == segments_[i] && !segments_.has_previous(i) ) {
							reference_is_nter = false;
							chimera_segs.push_back(create_chimera_segment(it2->first, it2->second,
								atom_map, reference_is_nter));
							continue;
						}
						if ( it2->first == segments_[i] && !segments_.has_next(i) ) {
							reference_is_nter = true;
							chimera_segs.push_back(create_chimera_segment(it2->first, it2->second,
								atom_map, reference_is_nter));
							continue;
						}
					}
				} else {
					for ( core::Size i = 1; i <= mobile_model.segments_.size(); ++i ) {
						if ( it2->second == mobile_model.segments_[i] && !mobile_model.segments_.has_previous(i) ) {
							reference_is_nter = true;
							chimera_segs.push_back(create_chimera_segment(it2->first, it2->second,
								atom_map, reference_is_nter));
							continue;
						}
						if ( it2->second == mobile_model.segments_[i] && !mobile_model.segments_.has_next(i) ) {
							reference_is_nter = false;
							chimera_segs.push_back(create_chimera_segment(it2->first, it2->second,
								atom_map, reference_is_nter));
							continue;
						}
					}
				}
			}
		}
	}

	runtime_assert(chimera_segs.size() > 0);
	return chimera_segs;
}


utility::vector1<SewSegment>
Assembly::segments() const {
	return segments_;
}

utility::vector1< utility::vector1<SewSegment> >
Assembly::all_segments() const {
	return all_segments_;
}



utility::vector1<SewSegment>
Assembly::get_model_segments(
	int model_id
) const {
	utility::vector1<SewSegment> model_segments;
	for ( core::Size i=1; i<=all_segments_.size(); ++i ) {
		for ( core::Size j=1; j<=all_segments_[i].size(); ++j ) {
			if ( all_segments_[i][j].model_id_ == model_id ) {
				model_segments.push_back(all_segments_[i][j]);
			}
		}
	}
	return model_segments;
}



std::set<core::Size>
Assembly::model_ids() const {
	std::set<core::Size> model_ids;

	for ( core::Size i=1; i<=all_segments_.size(); ++i ) {
		for ( core::Size j=1; j<=all_segments_[i].size(); ++j ) {
			model_ids.insert(all_segments_[i][j].model_id_);
		}
	}
	return model_ids;
}



Model
Assembly::regenerate_model(
	int model_id
) const {
	utility::vector1<SewSegment> assembly_segments =
		get_model_segments(model_id);

	Model regenerated_model;
	for ( utility::vector1<SewSegment>::const_iterator it = assembly_segments.begin(); it != assembly_segments.end(); ++it ) {
		regenerated_model.model_id_ = it->model_id_;
		regenerated_model.segments_.push_back(*it);
	}
	return regenerated_model;
}



void
Assembly::add_loop_segment(
	core::pose::Pose const & pose,
	protocols::loops::Loop loop,
	core::Size segment_index
) {
	SewSegment new_segment;
	new_segment.model_id_ = 0;
	new_segment.segment_id_ = 0;
	for ( core::Size i=loop.start(); i<=loop.stop(); ++i ) {
		SewResidue new_res;
		new_res.resnum_ = i;
		new_res.residue_type_ = pose.residue(i).type().name();
		//don't set any chi angles, yet

		for ( core::Size j=1; j<=4; ++j ) {
			SewAtom new_atom;
			new_atom.atomno_=j;
			new_atom.coords_=pose.residue(i).atom(j).xyz();
			new_res.basis_atoms_.push_back(new_atom);
		}
		new_segment.residues_.push_back(new_res);
	}
	segments_.insert_connected_segment(new_segment, segment_index);
	utility::vector1<SewSegment> loop_seg_vec;
	loop_seg_vec.push_back(new_segment);
	all_segments_.insert(all_segments_.begin()+segment_index, loop_seg_vec);
}




AtomMap
Assembly::atom_map_from_score_result(
	ScoreResult const & alignment_scores
) const {
	//insert basis atoms first
	AtomMap atom_map;
	atom_map.insert(std::make_pair(core::id::AtomID(1, alignment_scores.first.first.resnum), core::id::AtomID(1, alignment_scores.first.second.resnum)));
	atom_map.insert(std::make_pair(core::id::AtomID(2, alignment_scores.first.first.resnum), core::id::AtomID(2, alignment_scores.first.second.resnum)));
	atom_map.insert(std::make_pair(core::id::AtomID(3, alignment_scores.first.first.resnum), core::id::AtomID(3, alignment_scores.first.second.resnum)));

	//insert other alignments
	std::map< SegmentPair, AtomMap > segment_map = alignment_scores.second.segment_matches;
	std::map< SegmentPair, AtomMap >::const_iterator it = segment_map.begin();
	std::map< SegmentPair, AtomMap >::const_iterator it_end = segment_map.end();
	for ( ; it != it_end; ++it ) {
		atom_map.insert(it->second.begin(), it->second.end());
	}
	return atom_map;
}



void
Assembly::follow_edge(
	SewGraphCOP graph,
	HashEdge const * const edge,
	core::Size source_index
) {

	core::Size other_index = edge->get_other_ind(source_index);
	if ( TR.Debug.visible() ) {
		TR << "Appending node " << *graph->get_model_node(other_index) << " to assembly" << std::endl;
	}

	//Regenerate the source model from the Assembly (it can have changed coordinates from the ModelNode)
	Model const reference_model = regenerate_model(graph->get_model_node(source_index)->model().model_id_);
	SewResidue const reference_residue = reference_model.get_residue( edge->model_resnum( reference_model.model_id_ ) );

	//Get the model to be aligned from the graph
	Model mobile_model( graph->get_model_node(other_index)->model() );
	SewResidue mobile_residue = mobile_model.get_residue( edge->model_resnum( mobile_model.model_id_ ) );

	Hasher hasher;

	core::Size box_length = basic::options::option[ basic::options::OptionKeys::legacy_sewing::box_length ];

	ScoreResult edge_score = hasher.score_one(reference_model, reference_residue, mobile_model, mobile_residue, box_length);

	// std::pair< BasisPair, HashResult > ScoreResult

	//This 'check' is due to slight differences in geometric hashing before and after a
	//model transformation. For example, Model X, model Y, and model Z are scored from the initial model
	//file and model X has the requisite number of overlapping atoms with both model Y and model Z.
	//Now, during Assembly generation, we start with model Y, then superimpose model X using the
	//score result. This supimposition obviously results in a change to model X's coordinates, we'll
	//call this X'. The next step of assembly brings in model Z. In order to figure how how to superimpose
	//Z onto the Assembly, we need to rescore Z and X'. The slight differences in atomic coordinates can
	//trigger edge-effects with the geometric hashing, and thus differences in the output. In especially rare
	//cases, this causes two segments that previously had no overlapping atoms, to have overlapping atoms.
	//It's reasonable to assume that any overlaps of 2 atoms or less are the result of this happening, and therefore
	//we delete them here.
	auto it = edge_score.second.segment_matches.begin();
	auto it_end = edge_score.second.segment_matches.end();
	for ( ; it != it_end; ) {
		if ( it->second.size() <= 2 ) {
			edge_score.second.segment_matches.erase(it++);
		} else {
			++it;
		}
	}

	AtomMap atom_alignments = atom_map_from_score_result(edge_score);

	//map_residues(reference_model.model_id_, mobile_model, atom_alignments);

	align_model(mobile_model, atom_alignments, reference_model.model_id_);
	append_model(mobile_model, edge_score);

	edge_path_.push_back(boost::make_tuple(graph->get_model_node(source_index), graph->get_model_node(other_index), edge));
	path_.push_back(std::make_pair(reference_model.model_id_, mobile_model.model_id_));

	//Add other nodes from the mobile model so that they can be built on next and remove the reference node
	//so we don't built off it again
	std::set<core::Size> node_indices = graph->get_node_indices_from_model_id(mobile_model.model_id_);
	available_nodes_.insert(node_indices.begin(), node_indices.end());

	available_nodes_.erase(available_nodes_.find(source_index));
	available_nodes_.erase(available_nodes_.find(other_index));

}



///@details Use the given atom alignments to create mappings between SewResidues in the Assembly.
///These mappings will be used to generate 'Native' Rotamers for positions in the final Pose
void
Assembly::map_residues(
	int reference_model_id,
	Model mobile_model,
	AtomMap const & atom_alignments
) {
	for ( auto const & atom_alignment : atom_alignments ) {

		SewResidue * ref_res = nullptr;
		for ( core::Size seg_i=1; seg_i<=all_segments_.size(); ++seg_i ) {
			for ( core::Size seg_j=1; seg_j<=all_segments_[seg_i].size(); ++seg_j ) {
				SewSegment & cur_seg = all_segments_[seg_i][seg_j];
				for ( core::Size res_i=1; res_i<=cur_seg.residues_.size(); ++res_i ) {
					SewResidue & cur_res = cur_seg.residues_[res_i];
					for ( core::Size atom_i=1; atom_i<=cur_res.basis_atoms_.size(); ++atom_i ) {
						SewAtom const & cur_atom = cur_res.basis_atoms_[atom_i];
						if ( reference_model_id == cur_seg.model_id_ &&
								atom_alignment.first.rsd() == cur_res.resnum_ &&
								atom_alignment.first.atomno() == cur_atom.atomno_ ) {
							ref_res = &cur_res;
							break;
						}
					}
				}
			}
		}
		runtime_assert(ref_res != nullptr);

		SewResidue * mobile_res = nullptr;
		for ( core::Size seg_i=1; seg_i<=mobile_model.segments_.size(); ++seg_i ) {
			SewSegment & cur_seg = mobile_model.segments_[seg_i];
			for ( core::Size res_i=1; res_i<=cur_seg.residues_.size(); ++res_i ) {
				SewResidue & cur_res = cur_seg.residues_[res_i];
				for ( core::Size atom_i=1; atom_i<=cur_res.basis_atoms_.size(); ++atom_i ) {
					SewAtom const & cur_atom = cur_res.basis_atoms_[atom_i];
					if ( atom_alignment.second.rsd() == cur_res.resnum_ &&
							atom_alignment.second.atomno() == cur_atom.atomno_ ) {
						mobile_res = &cur_res;
						break;
					}
				}
			}
		}
		runtime_assert(mobile_res != nullptr);

		//Add the entire set of residues matched to the reference to the newly matched residue
		mobile_res->matched_residues_.push_back(*ref_res);
		mobile_res->matched_residues_.insert(mobile_res->matched_residues_.end(), ref_res->matched_residues_.begin(), ref_res->matched_residues_.end());

		//Add the newly matched residue to the matched list for the reference residue
		ref_res->matched_residues_.push_back(*mobile_res);
		if ( TR.Debug.visible() ) {
			TR.Debug << "Mobile residue has " << mobile_res->matched_residues_.size() << " matched residues" << std::endl;
			TR.Debug << "Reference residue has " << ref_res->matched_residues_.size() << " matched residues" << std::endl;
		}
	}
}


///@details create a map of sequence positions to ResidueOP objects for each of the 'natives' seen during the Assembly process. Ensure that if a residue
///is the first or last in the Assembly then we create lower and upper terminus types for that residues at tha position.
NativeRotamersMap
Assembly::generate_native_rotamers_map()
const {
	core::chemical::ResidueTypeSetCOP res_type_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	NativeRotamersMap nat_ro_map;
	core::Size counter = 0;
	for ( core::Size i=1; i<=segments_.size(); ++i ) {
		for ( core::Size j=1; j<=segments_[i].residues_.size(); ++j ) {
			++counter;
			//core::Size const seqpos = pose_num(segments_[i].model_id_, segments_[i].residues_[j].resnum_);
			core::Size const seqpos = counter;

			core::chemical::ResidueType const & base_type = res_type_set->name_map(segments_[i].residues_[j].residue_type_);
			core::conformation::ResidueOP new_residue;
			if ( seqpos == 1 ) {
				core::chemical::ResidueType const & res_type = res_type_set->get_residue_type_with_variant_added(base_type, core::chemical::LOWER_TERMINUS_VARIANT);
				new_residue = core::conformation::ResidueFactory::create_residue(res_type);
			} else if ( seqpos == total_residue() ) {
				core::chemical::ResidueType const & res_type = res_type_set->get_residue_type_with_variant_added(base_type, core::chemical::UPPER_TERMINUS_VARIANT);
				new_residue = core::conformation::ResidueFactory::create_residue(res_type);
			} else {
				core::chemical::ResidueType const & res_type1 = res_type_set->get_residue_type_with_variant_removed(base_type, core::chemical::LOWER_TERMINUS_VARIANT);
				core::chemical::ResidueType const & res_type2 = res_type_set->get_residue_type_with_variant_removed(res_type1, core::chemical::UPPER_TERMINUS_VARIANT);
				new_residue = core::conformation::ResidueFactory::create_residue(res_type2);
			}
			new_residue->seqpos(seqpos);
			//Don't save rotamers for which we have no chi angles (these are currently only residues generated by loophash fragment insertions)
			if ( segments_[i].residues_[j].chi_angles_.size() == new_residue->nchi() ) {
				new_residue->set_all_chi(segments_[i].residues_[j].chi_angles_);
				nat_ro_map[seqpos].push_back(std::make_pair(true, new_residue));
			} else {
				nat_ro_map[seqpos].push_back(std::make_pair(false, new_residue));
			}

			//Save the matched residues at each position
			for ( core::Size mapped_res_i=1; mapped_res_i <= segments_[i].residues_[j].matched_residues_.size(); ++mapped_res_i ) {
				core::chemical::ResidueType const & mapped_base_type = res_type_set->name_map(segments_[i].residues_[j].matched_residues_[mapped_res_i].residue_type_);

				core::conformation::ResidueOP new_matched_residue;
				if ( seqpos == 1 ) {
					core::chemical::ResidueType const & res_type = res_type_set->get_residue_type_with_variant_added(mapped_base_type, core::chemical::LOWER_TERMINUS_VARIANT);
					new_matched_residue = core::conformation::ResidueFactory::create_residue(res_type);
				} else if ( seqpos == total_residue() ) {
					core::chemical::ResidueType const & res_type = res_type_set->get_residue_type_with_variant_added(mapped_base_type, core::chemical::UPPER_TERMINUS_VARIANT);
					new_matched_residue = core::conformation::ResidueFactory::create_residue(res_type);
				} else {
					core::chemical::ResidueType const & res_type1 = res_type_set->get_residue_type_with_variant_removed(mapped_base_type, core::chemical::LOWER_TERMINUS_VARIANT);
					core::chemical::ResidueType const & res_type2 = res_type_set->get_residue_type_with_variant_removed(res_type1, core::chemical::UPPER_TERMINUS_VARIANT);
					new_matched_residue = core::conformation::ResidueFactory::create_residue(res_type2);
				}
				new_matched_residue->seqpos(seqpos);

				//Don't save rotamers for which we have no chi angles (these are currently only residues generated by loophash fragment insertions)
				if ( segments_[i].residues_[j].matched_residues_[mapped_res_i].chi_angles_.size() == new_matched_residue->nchi() ) {
					new_matched_residue->set_all_chi(segments_[i].residues_[j].matched_residues_[mapped_res_i].chi_angles_);
					nat_ro_map[seqpos].push_back(std::make_pair(true, new_matched_residue));
				} else {
					nat_ro_map[seqpos].push_back(std::make_pair(false, new_matched_residue));
				}
			} //foreach mapped residue
		} //foreach residue
	} //foreach segment
	return nat_ro_map;
}


///@details prepare the pose for packing by adding all the rotmaers defined by the chi-angles in SewResidue. Also,
///add residue type constraints to the Pose for each of the 'native' residue types for each position in the pose.
void
Assembly::prepare_for_packing(
	core::pose::Pose & pose,
	core::pack::task::TaskFactoryOP task_factory,
	core::Real /*base_native_bonus*/,
	core::Size /*neighbor_cutoff*/
) const {
	core::chemical::ResidueTypeSetCOP res_type_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	pose.update_residue_neighbors();

	///// layered design ////
	core::select::util::SelectResiduesByLayerOP layer_select( new core::select::util::SelectResiduesByLayer );
	layer_select->compute(pose, "");//empty string for secondary structure
	utility::vector1<core::Size> surface_residues = layer_select->selected_surface_residues();

	core::pack::task::operation::DisallowIfNonnativeOP disallow( new core::pack::task::operation::DisallowIfNonnative );
	disallow->restrict_to_residue(surface_residues);
	disallow->disallow_aas("WFY");
	task_factory->push_back(disallow);
	///// End layered design ////

	NativeRotamersMap nat_ro_map = generate_native_rotamers_map();
	// for(NativeRotamersMap::const_iterator map_it = nat_ro_map.begin();
	//   map_it != nat_ro_map.end(); ++map_it) {
	//
	//  core::Size seqpos = map_it->first;
	//  utility::vector1<core::conformation::ResidueOP> residues = map_it->second;
	//  for(core::Size i=1; i<=residues.size(); ++i) {
	//   if(pose.energies().tenA_neighbor_graph().get_node(seqpos)->num_neighbors_counting_self() >= neighbor_cutoff) {
	//    core::scoring::constraints::ResidueTypeConstraintOP matched_nat_res_constraint =
	//     new core::scoring::constraints::ResidueTypeConstraint(seqpos,
	//      residues[i]->type().name3(), residues[i]->type().name3(), base_native_bonus);
	//    pose.add_constraint(matched_nat_res_constraint);
	//    //Need to update residue neighbors every time we add a constraint. As adding a constraint clears the energies object
	//    pose.update_residue_neighbors();
	//   }
	//  }
	//
	//  core::pack::rotamer_set::AddResiduesRotamerSetOperation const nat_ro_set(residues);
	//  task_factory->push_back(new core::pack::task::operation::AppendResidueRotamerSet(seqpos, nat_ro_set.clone()) );
	// }
}



void
Assembly::align_model(
	Model & mobile_model,
	AtomMap const & atom_alignments,
	int reference_model_id
) const {
	utility::vector1< numeric::xyzVector< core::Real > > reference_coords;
	utility::vector1< numeric::xyzVector< core::Real > > mobile_coords;

	for ( auto const & atom_alignment : atom_alignments ) {

		for ( core::Size seg_i=1; seg_i<=all_segments_.size(); ++seg_i ) {
			for ( core::Size seg_j=1; seg_j<=all_segments_[seg_i].size(); ++seg_j ) {
				SewSegment const & cur_seg = all_segments_[seg_i][seg_j];
				//If this is a chimera, we shouldn't be aligning anything to it. If we don't check for this then
				//we may run into issues with PDB models (who retain their model id, it isn't set to CHIMERA_SEGMENT)
				//with multiple segments that can be built off of (the second time we build off of the pdb model we
				//have the chance to add non-pdb model residues from the chimera to the reference coords vector)
				if ( cur_seg.chimera_ ) {
					continue;
				}
				for ( core::Size res_i=1; res_i<=cur_seg.residues_.size(); ++res_i ) {
					SewResidue const & cur_res = cur_seg.residues_[res_i];
					for ( core::Size atom_i=1; atom_i<=cur_res.basis_atoms_.size(); ++atom_i ) {
						SewAtom const & cur_atom = cur_res.basis_atoms_[atom_i];
						if ( reference_model_id == cur_seg.model_id_ &&
								atom_alignment.first.rsd() == cur_res.resnum_ &&
								atom_alignment.first.atomno() == cur_atom.atomno_ ) {
							reference_coords.push_back(cur_atom.coords_);
						}
					}
				}
			}
		}

		ModelIterator<SewSegment> mobile_it = mobile_model.model_begin();
		ModelIterator<SewSegment> mobile_end = mobile_model.model_end();
		for ( ; mobile_it != mobile_end; ++mobile_it ) {
			if ( atom_alignment.second.rsd() == mobile_it.residue()->resnum_ &&
					atom_alignment.second.atomno() == mobile_it.atom()->atomno_ ) {
				mobile_coords.push_back(mobile_it.atom()->coords_);
			}
		}
	}
	if ( reference_coords.size() != mobile_coords.size() ) {
		std::stringstream err;
		err << "reference coords " << reference_coords.size();
		err << " != mobile_coords " << mobile_coords.size() << std::endl;
		err << "Alignment size " << atom_alignments.size() << std::endl;
		utility_exit_with_message(err.str());
	}
	//runtime_assert(reference_coords.size() == mobile_coords.size());

	numeric::xyzVector<core::Real> ref_com = numeric::center_of_mass(reference_coords);
	numeric::xyzVector<float> mobile_com = numeric::center_of_mass(mobile_coords);

	utility::vector1<core::Real> weights(reference_coords.size(), 1.0);
	numeric::xyzMatrix<core::Real> uu;
	core::Real sigma3;
	numeric::model_quality::findUU(reference_coords, mobile_coords, weights, (int)reference_coords.size(), uu, sigma3);

	numeric::xyzTransform<core::Real> transformer(uu,ref_com);
	ModelIterator<SewSegment> transform_it = mobile_model.model_begin();
	for ( ; transform_it != mobile_model.model_end(); ++transform_it ) {
		transform_it.atom()->coords_ = transform_it.atom()->coords_ - mobile_com;//move to origin
		transform_it.atom()->coords_ = transformer(transform_it.atom()->coords_);//rotate and move to reference
	}
}



///@details convert this assembly to a pose. This involves
///creating ideal residues and placing them on the backbone atom
///coordinates held by the assembly. This function assembles the pose
///residues in the order they are contained in the segments_ vector
///to re-order segments, use the reorder function
core::pose::Pose
Assembly::to_pose(
	std::string residue_type_set,
	bool create_cuts /*=true*/
) const {
	core::chemical::ResidueTypeSetCOP res_type_set =
		core::chemical::ChemicalManager::get_instance()->residue_type_set( residue_type_set );

	utility::vector1< std::pair< std::string, std::string > > atom_pairs;
	atom_pairs.push_back(std::make_pair("N", "N"));
	atom_pairs.push_back(std::make_pair("CA", "CA"));
	atom_pairs.push_back(std::make_pair("C", "C"));

	core::Size nres = 0;
	std::string sequence = "";
	for ( auto const & segment : segments_ ) {
		for ( core::Size res_i=1; res_i<=segment.residues_.size(); ++res_i ) {
			sequence += res_type_set->name_map(segment.residues_[res_i].residue_type_).name1();
			++nres;
		}
	}

	core::pose::Pose model_pose;
	core::pose::make_pose_from_sequence(model_pose, sequence, *res_type_set);

	core::Size counter = 1;
	for ( core::Size seg_i=1; seg_i<=segments_.size(); ++seg_i ) {
		for ( core::Size res_i=1; res_i<=segments_[seg_i].residues_.size(); ++res_i ) {
			SewResidue sew_res = segments_[seg_i].residues_[res_i];

			core::conformation::ResidueOP template_residue = model_pose.residue(counter).clone();
			template_residue->atom(1).xyz(sew_res.basis_atoms_[1].coords_);
			template_residue->atom(2).xyz(sew_res.basis_atoms_[2].coords_);
			template_residue->atom(3).xyz(sew_res.basis_atoms_[3].coords_);

			core::conformation::ResidueOP new_residue;
			core::chemical::ResidueType const & base_type = model_pose.residue(counter).type();
			if ( counter != 1 && base_type.is_lower_terminus() ) {
				core::chemical::ResidueType const & res_type = res_type_set->get_residue_type_with_variant_removed(base_type, core::chemical::LOWER_TERMINUS_VARIANT);
				new_residue = core::conformation::ResidueFactory::create_residue(res_type);
			} else if ( counter != model_pose.size() && base_type.is_upper_terminus() ) {
				core::chemical::ResidueType const & res_type = res_type_set->get_residue_type_with_variant_removed(base_type, core::chemical::UPPER_TERMINUS_VARIANT);
				new_residue = core::conformation::ResidueFactory::create_residue(res_type);
			} else {
				new_residue = core::conformation::ResidueFactory::create_residue(base_type);
			}
			//new_residue->orient_onto_residue(*template_residue, atom_pairs);
			new_residue->orient_onto_residue(*template_residue);
			if ( sew_res.chi_angles_.size() == new_residue->nchi() ) {
				new_residue->set_all_chi(sew_res.chi_angles_);
			}
			//   else if(segments_[seg_i].model_id_ != 0 && residue_type_set == core::chemical::FA_STANDARD) {
			//    TR.Warning << "No Chi angles found for model with non-zero model_id: " << segments_[seg_i].model_id_ << " " << sew_res.resnum_ << " " << sew_res.residue_type_ << std::endl;
			//   }

			//Fix backbone oxygens and hydrogens
			new_residue->atom(4).xyz(sew_res.basis_atoms_[4].coords_);

			model_pose.replace_residue(counter, *new_residue, false);
			++counter;
		}
	}

	//Fixup the hydrogens
	for ( core::Size i=1; i<=model_pose.size(); ++i ) {
		core::conformation::ResidueOP ires = model_pose.residue( i ).clone();
		core::conformation::idealize_hydrogens( *ires, model_pose.conformation() );
		model_pose.replace_residue( i, *ires, false );
	}

	//Set the fold tree - doesn't work with continuous assemblies due to
	//the chimera segments
	if ( create_cuts ) {
		core::kinematics::FoldTree ft;
		int jump_counter=0;
		for ( core::Size i=1; i<=segments_.size(); ++i ) {
			core::Size pose_seg_start = pose_num(segments_[i].model_id_, segments_[i].residues_[1].resnum_);
			core::Size pose_seg_end = pose_num(segments_[i].model_id_, segments_[i].residues_.back().resnum_);
			if ( i != 1 ) {
				++jump_counter;
				ft.add_edge(1, (int)pose_seg_start, jump_counter);
			}
			ft.add_edge((int)pose_seg_start, (int)pose_seg_end, core::kinematics::Edge::PEPTIDE);
		}
		model_pose.fold_tree(ft);
	}

	runtime_assert(model_pose.size() == nres);

	//finally, add the partner
	if ( partner_pose_ ) {
		core::pose::Pose partner_pose_copy = *partner_pose_->clone();
		if ( partner_pose_->conformation().residue_typeset_mode(false) != res_type_set->mode() ) {
			core::util::switch_to_residue_type_set( partner_pose_copy, residue_type_set, true );
		}
		core::pose::append_pose_to_pose(model_pose, partner_pose_copy, true);
	}

	return model_pose;
}



///@details for testing purposes, create a pose that includes all untrimmed
///segments. These will (or at least, should) heavily overlap with one another.
///Make each model a new chain
core::pose::Pose
Assembly::to_multichain_pose(
	std::string residue_type_set
) const {
	core::chemical::ResidueTypeSetCOP res_type_set =
		core::chemical::ChemicalManager::get_instance()->residue_type_set( residue_type_set );

	utility::vector1< std::pair< std::string, std::string > > atom_pairs;
	atom_pairs.push_back(std::make_pair("N", "N"));
	atom_pairs.push_back(std::make_pair("CA", "CA"));
	atom_pairs.push_back(std::make_pair("C", "C"));

	std::string sequence = "";
	std::set<core::Size> all_model_ids = model_ids();
	auto it = all_model_ids.begin();
	auto it_end = all_model_ids.end();
	for ( ; it != it_end; ++it ) {
		utility::vector1<SewSegment> model_segs = get_model_segments(*it);
		for ( core::Size seg_i=1; seg_i<=model_segs.size(); ++seg_i ) {
			for ( core::Size res_i=1; res_i<=model_segs[seg_i].residues_.size(); ++res_i ) {
				sequence += res_type_set->name_map(model_segs[seg_i].residues_[res_i].residue_type_).name1();
			}
		}
	}

	core::pose::Pose model_pose;
	core::pose::make_pose_from_sequence(model_pose, sequence, *res_type_set);

	core::Size counter = 1;
	utility::vector1<core::Size> chain_endings;
	it = all_model_ids.begin();
	for ( ; it != it_end; ++it ) {
		utility::vector1<SewSegment> model_segs = get_model_segments(*it);
		for ( core::Size seg_i=1; seg_i<=model_segs.size(); ++seg_i ) {
			for ( core::Size res_i=1; res_i<=model_segs[seg_i].residues_.size(); ++res_i ) {
				SewResidue sew_res = model_segs[seg_i].residues_[res_i];

				core::conformation::ResidueOP template_residue = model_pose.residue(counter).clone();
				template_residue->atom(1).xyz(sew_res.basis_atoms_[1].coords_);
				template_residue->atom(2).xyz(sew_res.basis_atoms_[2].coords_);
				template_residue->atom(3).xyz(sew_res.basis_atoms_[3].coords_);

				core::conformation::ResidueOP new_residue = core::conformation::ResidueFactory::create_residue(model_pose.residue(counter).type());
				//new_residue->orient_onto_residue(*template_residue, atom_pairs);
				new_residue->orient_onto_residue(*template_residue);
				if ( sew_res.chi_angles_.size() == new_residue->nchi() ) {
					new_residue->set_all_chi(sew_res.chi_angles_);
				}

				//Fix backbone oxygens and hydrogens
				new_residue->atom(4).xyz(sew_res.basis_atoms_[4].coords_);

				model_pose.replace_residue(counter, *new_residue, false);

				++counter;
			}
		}
		chain_endings.push_back(counter-1);
	}
	chain_endings.pop_back();
	model_pose.conformation().chain_endings(chain_endings);

	for ( core::Size i=1; i<=model_pose.size(); ++i ) {
		core::conformation::ResidueOP ires = model_pose.residue( i ).clone();
		core::conformation::idealize_hydrogens( *ires, model_pose.conformation() );
		model_pose.replace_residue( i, *ires, false );
	}

	return model_pose;
}



bool
Assembly::reorder_randomly(
	core::Real max_loop_distance
) {
	utility::vector1< utility::vector1<core::Size> > possible_orders =
		find_possible_orders(max_loop_distance);
	if ( possible_orders.size() == 0 ) {
		return false;
	}
	core::Size rand_index = numeric::random::random_range(1, (int)possible_orders.size());
	segments_.reorder(possible_orders[rand_index]);

	return true;
}



void
Assembly::reorder(
	utility::vector1<core::Size> new_order
) {
	segments_.reorder(new_order);
	if ( basic::options::option[ basic::options::OptionKeys::legacy_sewing::dump_pdbs ] ) {
		to_pose(core::chemical::CENTROID).dump_pdb("after_reorder.pdb");
	}
}


void
Assembly::append_segments(
	utility::vector1<SewSegment> const & segments
) {
	segments_.insert(segments_.end(), segments.begin(), segments.end());
	for ( core::Size i=1; i<=segments.size(); ++i ) {
		utility::vector1<SewSegment> new_vec;
		new_vec.push_back(segments[i]);
		all_segments_.push_back(new_vec);
	}
}

void
Assembly::delete_segments(
	core::Size seg_index
) {
	segments_.erase(segments_.begin() + (seg_index-1));
	all_segments_.erase(all_segments_.begin() + (seg_index-1));
}


///@details find all possible segment orders where the CA-CA
///distance of all unclosed loops is less than the max distance
///respect any current connections.
utility::vector1< utility::vector1<core::Size> >
Assembly::find_possible_orders(
	core::Real max_loop_distance
) const {

	core::Real max_distance_sq = max_loop_distance * max_loop_distance;

	//if a segment is the next_segment of any other segment then don't consider it in the optimization
	utility::vector1< core::Size > free_segments;
	for ( core::Size i=1; i<=segments_.size(); ++i ) {
		if ( !segments_.has_previous(i) ) {
			free_segments.push_back(i);
		}
	}

	core::Size n_dims = free_segments.size();
	utility::vector1<core::Size> dim_sizes =
		utility::vector1<core::Size>(free_segments.size(), free_segments.size());

	utility::vector1< utility::vector1<core::Size> > possible_orders;
	for ( utility::LexicographicalIterator lex( dim_sizes ); ! lex.at_end(); ++lex ) {

		bool valid = true;
		for ( core::Size i=1; i<=n_dims; ++i ) {
			for ( core::Size j=i+1; j<=n_dims; ++j ) {
				if ( lex[i] == lex[j] ) {
					valid=false;
					break;
				}
			}
			if ( !valid ) {
				break;
			}
		}

		if ( valid ) {
			for ( core::Size i=1; i<=n_dims-1; ++i ) {

				SewSegment n_term_segment = segments_[segments_.c_segment(free_segments[lex[i]])];
				SewResidue n_term_res = n_term_segment.residues_.back();

				SewSegment c_term_segment = segments_[free_segments[lex[i+1]]];
				SewResidue c_term_res = c_term_segment.residues_[1];
				core::Real dist_sq = n_term_res.basis_atoms_[2].coords_.distance_squared(c_term_res.basis_atoms_[2].coords_);

				if ( dist_sq > max_distance_sq ) {
					valid=false;
					break;
				}
				if ( TR.Debug.visible() ) {
					TR.Debug << "Distance between " << pose_num(n_term_segment.model_id_, n_term_res.resnum_) << " and " << pose_num(c_term_segment.model_id_, c_term_res.resnum_) << ": " << dist_sq << std::endl;
				}
			}
		}

		if ( valid ) {
			utility::vector1<core::Size> new_order;
			for ( core::Size i=1; i<=n_dims; ++i ) {
				new_order.push_back(free_segments[lex[i]]);
			}
			possible_orders.push_back(new_order);
		}
	}

	TR << "Found " << possible_orders.size() << " valid orders for assembly." << std::endl;
	return possible_orders;
}


core::Size
Assembly::get_next_reference_node(
	SewGraphOP /*graph*/
) const{
	if ( available_nodes_.size() == 0 ) {
		utility_exit_with_message("No available nodes to choose from!");
	}
	core::Size rand_index = numeric::random::random_range(0, (int)available_nodes_.size()-1);
	auto it = available_nodes_.begin();
	std::advance(it, rand_index);
	return *it;
}


core::Size
Assembly::pose_num(
	int model_id,
	core::Size resnum
) const {
	core::Size counter=1;
	for ( core::Size i=1; i<=segments_.size(); ++i ) {
		for ( core::Size j=1; j<=segments_[i].residues_.size(); ++j ) {
			if ( segments_[i].model_id_ == model_id && segments_[i].residues_[j].resnum_ == resnum ) {
				return counter;
			}
			++counter;
		}
	}
	utility_exit_with_message("No residue " + utility::to_string(resnum) + " with model ID " + utility::to_string(model_id));
	return 0;
}



core::Size
Assembly::total_residue()
const {
	core::Size sum=0;
	for ( core::Size i=1; i<=segments_.size(); ++i ) {
		sum+=segments_[i].residues_.size();
	}
	return sum;
}

///@brief Records which positions in the given pose match the sequence identify of at least
///one of the amino acids from the underlying models at that position.
std::set<core::Size>
Assembly::native_positions(
	core::pose::Pose const & pose
) const {
	core::Size residue_counter=1;

	core::chemical::ResidueTypeSetCOP res_type_set;
	if ( pose.is_fullatom() ) {
		res_type_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	} else if ( pose.is_centroid() ) {
		res_type_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
	} else {
		utility_exit_with_message("Not a centroid or fullatom pose. Panic!");
	}

	std::set<core::Size> native_positions;
	for ( core::Size i=1; i<=segments_.size(); ++i ) {
		for ( core::Size j=1; j<=segments_[i].residues_.size(); ++j ) {
			if ( pose.residue(residue_counter).type().name3() == res_type_set->name_map(segments_[i].residues_[j].residue_type_).name3() ) {
				native_positions.insert(residue_counter);
				++residue_counter;
				continue;
			}
			for ( core::Size m=1; m <= segments_[i].residues_[j].matched_residues_.size(); ++m ) {
				if ( TR.Debug.visible() ) {
					TR.Debug << "Checking for matched 'native' at position " << residue_counter << ": "
						<< res_type_set->name_map(segments_[i].residues_[j].matched_residues_[m].residue_type_).name3() << std::endl;
				}
				if ( pose.residue(residue_counter).type().name3() == res_type_set->name_map(segments_[i].residues_[j].matched_residues_[m].residue_type_).name3() ) {
					native_positions.insert(residue_counter);
					break;
				}
			}
			++residue_counter;
		}
	}
	return native_positions;
}

core::Real
Assembly::percent_native(
	core::pose::Pose const & pose
) const {
	return native_positions(pose).size() / (core::Real)pose.size();
}


std::string
Assembly::natives_select(
	core::pose::Pose const & pose,
	std::string object_name /*=""*/
) const {
	std::stringstream pymol_select;

	std::set<core::Size> native_positions_set = native_positions(pose);
	auto it = native_positions_set.begin();
	auto it_end = native_positions_set.end();
	// /object-name/segi-identifier/chain-identifier/resi-identifier/name-identifier
	pymol_select << "select " << object_name << "_natives, /" << object_name << "///";
	for ( ; it != it_end; ++it ) {
		pymol_select << *it << "+";
	}
	return pymol_select.str();
}



utility::vector1<core::Size>
Assembly::pose_loop_anchors() const {
	utility::vector1<core::Size> loop_anchors;
	for ( core::Size i=1; i<segments_.size(); ++i ) {
		if ( !segments_.has_next(i) ) {
			loop_anchors.push_back(pose_num(segments_[i].model_id_, segments_[i].residues_.back().resnum_));
		}
	}
	return loop_anchors;
}


utility::vector1<core::Size>
Assembly::disconnected_segments() const {
	utility::vector1<core::Size> disconnected_segments;
	for ( core::Size i=1; i<segments_.size(); ++i ) {
		if ( !segments_.has_next(i) ) {
			disconnected_segments.push_back(i);
		}
	}
	return disconnected_segments;
}


utility::vector1< std::pair<core::Size, core::Size> >
Assembly::path() const {
	return path_;
}

std::string
Assembly::string_path() const {
	std::stringstream name;
	for ( core::Size i=1; i<=path_.size(); ++i ) {
		name << "(" << path_[i].first << "," << path_[i].second << ")-";
	}
	return name.str();
}

std::string
Assembly::string_blosum(
	utility::vector1< std::pair<core::Real, core::Size> > blosum_history
) const {
	std::stringstream name;
	for ( core::Size i=1; i<=blosum_history.size(); ++i ) {
		name << "(" << blosum_history[i].first << "," << blosum_history[i].second << ")-";
	}
	return name.str();
}

void
Assembly::set_partner(
	core::pose::PoseOP partner_pose
){
	partner_pose_ = partner_pose;
}

core::pose::PoseOP
Assembly::get_partner() const {
	return partner_pose_;
}

std::ostream &
operator<<(std::ostream& out, Assembly const & assembly ) {
	out << "Assembly: " << assembly.segments_.size() << std::endl;
	for ( core::Size i=1; i<= assembly.segments_.size(); ++i ) {
		out << "\tSegment " << i << std::endl;
		out << "\t\tmodel: " << assembly.segments_[i].model_id_ << std::endl;
		out << "\t\tsegment id: " << assembly.segments_[i].segment_id_ << std::endl;
		if ( assembly.all_segments_[i].size() > 1 ) {
			for ( core::Size j=1; j<=assembly.all_segments_[i].size(); ++j ) {
				out << "\t\t\thidden segment " << assembly.all_segments_[i][j].model_id_ <<
					"-" << assembly.all_segments_[i][j].segment_id_ << std::endl;
			}
		}
	}
	return out;
}

//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////   Iterator Stuff   ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

ModelConstIterator<SewSegment>
Assembly::assembly_begin() const {
	auto seg_it = segments_.begin();
	utility::vector1<SewResidue>::const_iterator res_it;
	utility::vector1<SewAtom>::const_iterator atom_it;

	if ( seg_it == segments_.end() ) {
		return ModelConstIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
	}
	res_it = seg_it->residues_.begin();

	if ( res_it == seg_it->residues_.end() ) {
		return ModelConstIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
	}
	atom_it = res_it->basis_atoms_.begin();

	return ModelConstIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
}

ModelIterator<SewSegment>
Assembly::assembly_begin() {
	auto seg_it = segments_.begin();
	utility::vector1<SewResidue>::iterator res_it;
	utility::vector1<SewAtom>::iterator atom_it;

	if ( seg_it == segments_.end() ) {
		return ModelIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
	}
	res_it = seg_it->residues_.begin();

	if ( res_it == seg_it->residues_.end() ) {
		return ModelIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
	}

	atom_it = res_it->basis_atoms_.begin();
	return ModelIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
}

ModelConstIterator<SewSegment>
Assembly::assembly_end() const {
	auto seg_it = segments_.end();
	utility::vector1<SewResidue>::const_iterator res_it;
	utility::vector1<SewAtom>::const_iterator atom_it;

	if ( seg_it == segments_.begin() ) {
		return ModelConstIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
	}
	--seg_it;
	res_it = seg_it->residues_.end();

	if ( res_it == seg_it->residues_.begin() ) {
		return ModelConstIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
	}
	--res_it;
	atom_it = res_it->basis_atoms_.end();

	return ModelConstIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
}

ModelIterator<SewSegment>
Assembly::assembly_end() {
	auto seg_it = segments_.end();
	utility::vector1<SewResidue>::iterator res_it;
	utility::vector1<SewAtom>::iterator atom_it;

	if ( seg_it == segments_.begin() ) {
		return ModelIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
	}
	--seg_it;
	res_it = seg_it->residues_.end();

	if ( res_it == seg_it->residues_.begin() ) {
		return ModelIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
	}
	--res_it;
	atom_it = res_it->basis_atoms_.end();

	return ModelIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
}



} //legacy_sewing namespace
} //protocols namespace
