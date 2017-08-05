// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/conformation/Model.cc
/// @author   Tim Jacobs

//Unit headers
#include <protocols/sewing/conformation/Model.hh>
#include <protocols/sewing/conformation/Assembly.hh>
#include <protocols/sewing/conformation/DisembodiedAssembly.hh>

//Package headers
#include <basic/Tracer.hh>

//Protocol headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/database/sql_utils.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <protocols/features/ProteinSilentReport.hh>

#include <core/conformation/Residue.functions.hh>

#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/DbDataType.hh>


//Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

//Numeric headers
#include <numeric/random/random.hh>

namespace protocols {
namespace sewing  {

static THREAD_LOCAL basic::Tracer TR("protocols.sewing.Model");

///////////////////////////////////////////////////////////////////////
////////////////////   SewResidue Class    ////////////////////////////
///////////////////////////////////////////////////////////////////////

SewResidue
Model::get_residue(
	core::Size resnum
) const {
	for ( core::Size i=1; i<=segments_.size(); ++i ) {
		for ( core::Size j=1; j<=segments_[i].residues_.size(); ++j ) {
			if ( segments_[i].residues_[j].resnum_ == resnum ) {
				return segments_[i].residues_[j];
			}
		}
	}
	std::stringstream err;
	err << "Residue " << resnum << " not found in model " << model_id_;
	utility_exit_with_message(err.str());
	return segments_[1].residues_[1];
}

///////////////////////////////////////////////////////////////////////
////////////////////   SegmentGraph Class    //////////////////////////
///////////////////////////////////////////////////////////////////////

void
SegmentGraph::add_connection(
	core::Size i,
	core::Size j
){
	forward_connections_.insert(std::make_pair(i, j));
	reverse_connections_.insert(std::make_pair(j, i));
}

core::Size
SegmentGraph::c_segment(
	core::Size ind
) const {
	while ( true ) {
		std::map<core::Size, core::Size>::const_iterator it = forward_connections_.find(ind);
		if ( it == forward_connections_.end() ) {
			return ind;
		}
		ind = it->second;
	}
}

core::Size
SegmentGraph::n_segment(
	core::Size ind
) const {
	while ( true ) {
		std::map<core::Size, core::Size>::const_iterator it = reverse_connections_.find(ind);
		if ( it == reverse_connections_.end() ) {
			return ind;
		}
		ind = it->second;
	}
}

bool
SegmentGraph::has_next(
	core::Size ind
) const {
	return forward_connections_.find(ind) != forward_connections_.end();
}

core::Size
SegmentGraph::next(
	core::Size ind
) const {
	return forward_connections_.find(ind)->second;
}

bool
SegmentGraph::has_previous(
	core::Size ind
) const {
	return reverse_connections_.find(ind) != reverse_connections_.end();
}

core::Size
SegmentGraph::previous(
	core::Size ind
) const {
	return reverse_connections_.find(ind)->second;
}

void
SegmentGraph::insert_connected_segment(
	SewSegment segment,
	core::Size ind
) {

	insert(begin()+ind, segment);

	std::map<core::Size, core::Size> new_forward;
	std::map<core::Size, core::Size> new_reverse;
	std::map<core::Size, core::Size>::iterator f_it = forward_connections_.begin();
	std::map<core::Size, core::Size>::iterator f_end = forward_connections_.end();
	std::map<core::Size, core::Size>::iterator r_it = reverse_connections_.begin();
	for ( ; f_it != f_end; ++f_it, ++r_it ) {
		core::Size new_f_start = f_it->first;
		core::Size new_f_end = f_it->second;
		core::Size new_r_start = r_it->first;
		core::Size new_r_end = r_it->second;
		if ( f_it->first > ind ) ++new_f_start;
		if ( f_it->second > ind ) ++new_f_end;
		if ( r_it->first > ind ) ++new_r_start;
		if ( r_it->second > ind ) ++new_r_end;
		new_forward.insert(std::make_pair(new_f_start, new_f_end));
		new_reverse.insert(std::make_pair(new_r_start, new_r_end));
	}
	forward_connections_ = new_forward;
	reverse_connections_ = new_reverse;
	forward_connections_.insert(std::make_pair(ind, ind+1));
	reverse_connections_.insert(std::make_pair(ind+1, ind));

}

std::set<SewSegment>
SegmentGraph::erase_segments(
	std::set<SewSegment> segments_to_erase
){

	//Convert segments to indices
	std::set<core::Size> indices_to_erase;
	for ( core::Size i=1; i<=size(); ++i ) {
		if ( segments_to_erase.find(utility::vector1<SewSegment>::operator[](i)) != segments_to_erase.end() ) {
			indices_to_erase.insert(i);
		}
	}

	//Remove any segments chemically connected to a segment being deleted
	std::set<core::Size> linker_segments_to_delete;
	for ( core::Size i=1; i<=size(); ++i ) {
		if ( (has_next(i) && indices_to_erase.find(next(i)) != indices_to_erase.end()) ||
				(has_previous(i) && indices_to_erase.find(previous(i)) != indices_to_erase.end()) ) {
			linker_segments_to_delete.insert(i);
		}
	}
	indices_to_erase.insert(linker_segments_to_delete.begin(), linker_segments_to_delete.end());

	//Remove all connections to segments we are deleting
	std::set<core::Size>::const_iterator it = indices_to_erase.begin();
	std::set<core::Size>::const_iterator it_end = indices_to_erase.end();
	for ( ; it != it_end; ++it ) {
		if ( forward_connections_.find(*it) != forward_connections_.end() ) {
			reverse_connections_.erase(forward_connections_.find(*it)->second);
			forward_connections_.erase(*it);
		}
		if ( reverse_connections_.find(*it) != reverse_connections_.end() ) {
			forward_connections_.erase(reverse_connections_.find(*it)->second);
			reverse_connections_.erase(*it);
		}
	}

	//delete from vector
	std::set<SewSegment> deleted_segments;
	it = indices_to_erase.begin();
	core::Size rm_counter=0;
	for ( ; it != it_end; ++it ) {
		core::Size adjusted_index = *it - rm_counter - 1;
		deleted_segments.insert(*(begin()+adjusted_index));
		utility::vector1<SewSegment>::erase(begin() + adjusted_index);
		++rm_counter;
	}

	//update connections
	std::map<core::Size, core::Size> new_forward;
	std::map<core::Size, core::Size> new_reverse;
	std::map<core::Size, core::Size>::iterator f_it = forward_connections_.begin();
	std::map<core::Size, core::Size>::iterator f_end = forward_connections_.end();
	std::map<core::Size, core::Size>::iterator r_it = reverse_connections_.begin();
	for ( ; f_it != f_end; ++f_it, ++r_it ) {
		it = indices_to_erase.begin();
		core::Size new_f_start = f_it->first;
		core::Size new_f_end = f_it->second;
		core::Size new_r_start = r_it->first;
		core::Size new_r_end = r_it->second;
		for ( ; it != it_end; ++it ) {
			if ( f_it->first > *it ) --new_f_start;
			if ( f_it->second > *it ) --new_f_end;
			if ( r_it->first > *it ) --new_r_start;
			if ( r_it->second > *it ) --new_r_end;
		}
		new_forward.insert(std::make_pair(new_f_start, new_f_end));
		new_reverse.insert(std::make_pair(new_r_start, new_r_end));
	}
	forward_connections_ = new_forward;
	reverse_connections_ = new_reverse;

	return deleted_segments;
}

std::map<core::Size,core::Size>
SegmentGraph::forward_connections() const { return forward_connections_; }

std::map<core::Size,core::Size>
SegmentGraph::reverse_connections() const { return reverse_connections_; }

void
SegmentGraph::reorder(
	utility::vector1<core::Size> new_order
) {
	std::map<core::Size, core::Size> old_to_new;
	core::Size counter=0;
	for ( core::Size i=1; i<=new_order.size(); ++i ) {
		++counter;
		old_to_new[new_order[i]] = counter;

		core::Size temp = new_order[i];
		while ( has_next(temp) ) {
			++counter;
			temp = next(temp);
			old_to_new[temp] = counter;
		}
	}

	//swap
	utility::vector1<SewSegment> old = *this;
	std::map<core::Size, core::Size>::const_iterator it = old_to_new.begin();
	std::map<core::Size, core::Size>::const_iterator it_end = old_to_new.end();
	for ( ; it != it_end; ++it ) {
		utility::vector1<SewSegment>::operator[](it->second) = old[it->first];
	}

	//remap edges
	std::map<core::Size, core::Size> new_forward;
	std::map<core::Size, core::Size> new_reverse;

	std::map<core::Size, core::Size>::iterator f_it = forward_connections_.begin();
	std::map<core::Size, core::Size>::iterator f_end = forward_connections_.end();
	std::map<core::Size, core::Size>::iterator r_it = reverse_connections_.begin();
	for ( ; f_it != f_end; ++f_it, ++r_it ) {
		new_forward[old_to_new[f_it->first]] = old_to_new[f_it->second];
		new_reverse[old_to_new[r_it->first]] = old_to_new[r_it->second];
	}
	forward_connections_ = new_forward;
	reverse_connections_ = new_reverse;
}

void
SegmentGraph::clear_connections(
){
	forward_connections_.clear();
	reverse_connections_.clear();
}

///////////////////////////////////////////////////////////////////////
////////////////////////   Model Class    /////////////////////////////
///////////////////////////////////////////////////////////////////////

ModelIterator<SewSegment>
Model::model_begin() {
	utility::vector1<SewSegment>::iterator seg_it = segments_.begin();

	if ( seg_it == segments_.end() ) {
		utility_exit_with_message("Model " + utility::to_string(seg_it->model_id_) + " has no segments to iterator over!");
	}
	utility::vector1<SewResidue>::iterator res_it = seg_it->residues_.begin();

	if ( res_it == seg_it->residues_.end() ) {
		utility_exit_with_message("Model " + utility::to_string(seg_it->model_id_) + ", Segment " +
			utility::to_string(seg_it->segment_id_) + " has no residues to iterator over!");
	}
	utility::vector1<SewAtom>::iterator atom_it = res_it->basis_atoms_.begin();
	return ModelIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
}

ModelConstIterator<SewSegment>
Model::model_begin() const {
	utility::vector1<SewSegment>::const_iterator seg_it = segments_.begin();

	if ( seg_it == segments_.end() ) {
		utility_exit_with_message("Model " + utility::to_string(seg_it->model_id_) + " has no segments to iterator over!");
	}
	utility::vector1<SewResidue>::const_iterator res_it = seg_it->residues_.begin();

	if ( res_it == seg_it->residues_.end() ) {
		utility_exit_with_message("Model " + utility::to_string(seg_it->model_id_) + ", Segment " +
			utility::to_string(seg_it->segment_id_) + " has no residues to iterator over!");
	}
	utility::vector1<SewAtom>::const_iterator atom_it = res_it->basis_atoms_.begin();
	return ModelConstIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
}

ModelIterator<SewSegment>
Model::model_end() {
	utility::vector1<SewSegment>::iterator seg_it;
	utility::vector1<SewResidue>::iterator res_it;
	if ( segments_.size() == 0 ) {
		utility_exit_with_message("Model " + utility::to_string(seg_it->model_id_) + " has no segments to iterator over!");
	}
	seg_it = segments_.end();
	--seg_it;

	if ( segments_[segments_.size()].residues_.size() == 0 ) {
		utility_exit_with_message("Model " + utility::to_string(seg_it->model_id_) + ", Segment " +
			utility::to_string(seg_it->segment_id_) + " has no residues to iterator over!");
	}
	res_it = seg_it->residues_.end();
	--res_it;

	utility::vector1<SewAtom>::iterator atom_it =
		res_it->basis_atoms_.end();
	return ModelIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
}

ModelConstIterator<SewSegment>
Model::model_end() const {
	utility::vector1<SewSegment>::const_iterator seg_it;
	utility::vector1<SewResidue>::const_iterator res_it;
	if ( segments_.size() == 0 ) {
		utility_exit_with_message("Model " + utility::to_string(seg_it->model_id_) + " has no segments to iterator over!");
	}
	seg_it = segments_.end();
	--seg_it;

	if ( segments_[segments_.size()].residues_.size() == 0 ) {
		utility_exit_with_message("Model " + utility::to_string(seg_it->model_id_) + ", Segment " +
			utility::to_string(seg_it->segment_id_) + " has no residues to iterator over!");
	}
	res_it = seg_it->residues_.end();
	--res_it;

	utility::vector1<SewAtom>::const_iterator atom_it =
		res_it->basis_atoms_.end();
	return ModelConstIterator<SewSegment>(seg_it, segments_.end(), res_it, atom_it);
}



core::pose::Pose
Model::to_pose_from_db() const {
	core::pose::Pose test_pose;
	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );
	protocols::features::ProteinSilentReportOP protein_silent_report( new protocols::features::ProteinSilentReport() );
	protein_silent_report->load_pose(db_session, structure_id_, test_pose);

	return test_pose;
}


std::set<core::Size>
Model::segment_ids() const {
	std::set<core::Size> segment_ids;
	for ( core::Size i=1; i<=segments_.size(); ++i ) {
		segment_ids.insert(segments_[i].segment_id_);
	}
	return segment_ids;
}

///@brief Use the given segment-ranges to create a Model. Give this
///model an id that doesn't currently exist in the given map of models
Model
create_model_from_pose(
	core::pose::Pose const & pose,
	int model_id
){
	utility::vector1< std::pair<core::Size,core::Size> > segments;
	segments.push_back( std::make_pair(1, pose.size()) );
	return create_model_from_pose(pose, segments, model_id);
}

///@brief Use the given segment-ranges to create a Model. Give this
///model an id that doesn't currently exist in the given map of models
Model
create_model_from_pose(
	core::pose::Pose const & pose,
	utility::vector1< std::pair<core::Size,core::Size> > const & segments,
	int model_id
){
	Model model;
	model.model_id_ = model_id;

	for ( core::Size i=1; i<=segments.size(); ++i ) {
		SewSegment segment = SewSegment();
		segment.model_id_ = model_id;
		segment.hash_ = true;
		segment.segment_id_ = i;

		for ( core::Size j=segments[i].first; j<=segments[i].second; ++j ) {
			SewResidue residue;
			residue.resnum_ = j;
			residue.residue_type_ = pose.residue(j).type().name();
			residue.chi_angles_ = pose.residue(j).chi();

			for ( core::Size k=1; k<=4; k++ ) {
				SewAtom atom;
				atom.atomno_ = k;
				atom.coords_ = pose.residue(j).atom(k).xyz();
				residue.basis_atoms_.push_back(atom);
			}
			segment.residues_.push_back(residue);
			if ( j == (segments[i].first + segments[i].second)/2 ) {
				segment.dssp_ = pose.secstruct(j);
			}
		}
		model.segments_.push_back(segment);
	}

	for ( core::Size i=1; i<segments.size(); ++i ) {
		model.segments_.add_connection(i, i+1);
	}
	return model;
}



void
Model::trim_db_pose(
	core::pose::Pose & pose
) const {
	std::set<core::Size> resnums;
	for ( core::Size seg_i=1; seg_i<=segments_.size(); ++seg_i ) {
		for ( core::Size res_i=1; res_i<=segments_[seg_i].residues_.size(); ++res_i ) {
			resnums.insert(segments_[seg_i].residues_[res_i].resnum_);
		}
	}

	core::Size num_removed_residues=0;
	core::Size total_res = pose.size();
	for ( core::Size i=1; i<=total_res; ++i ) {
		if ( resnums.find(i) == resnums.end() ) {
			pose.conformation().delete_residue_slow(i-num_removed_residues);
			++num_removed_residues;
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////       UTILITY FUNCTIONS             //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

std::map< int, Model >
result_to_models(
	cppdb::result res
){

	std::string input_tag, res_type, dssp;
	int model_id;
	core::Size struct_id, segment_id, seqpos, atomno;
	core::Real distance, hoist, packing, meridian, chi1, chi2, chi3, chi4, x, y, z;

	std::map< int, Model > models;

	bool first = true;
	Model cur_model;
	SewSegment cur_segment = SewSegment();
	SewResidue cur_residue;
	while ( res.next() ) {
		res >> input_tag >> struct_id >> model_id >> distance >> hoist >> packing >> meridian >>
			segment_id >> dssp >>
			res_type >> chi1 >> chi2 >> chi3 >> chi4 >> seqpos >>
			atomno >> x >> y >> z;

		//new model
		if ( model_id != cur_model.model_id_ ) {
			if ( !first ) {
				cur_segment.residues_.push_back(cur_residue);
				cur_model.segments_.push_back(cur_segment);
				models[cur_model.model_id_] = cur_model;
			} else {
				first=false;
			}

			cur_model = Model();
			cur_model.model_id_ = model_id;
			cur_model.structure_id_ = struct_id;
			cur_model.pdb_code_ = input_tag;
			cur_model.distance_ = distance;
			cur_model.hoist_angle_degrees_ = hoist;
			cur_model.packing_angle_degrees_ = packing;
			cur_model.meridian_angle_degrees_ = meridian;

			cur_segment = SewSegment();
			cur_segment.model_id_ = model_id;
			cur_segment.segment_id_ = segment_id;
			cur_segment.dssp_ = dssp[0];

			cur_residue = SewResidue();
			cur_residue.resnum_ = seqpos;
			cur_residue.residue_type_ = res_type;
			if ( chi1 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi1);
			}
			if ( chi2 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi2);
			}
			if ( chi3 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi3);
			}
			if ( chi4 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi4);
			}
		} else if ( segment_id != cur_segment.segment_id_ ) {
			//new segment
			cur_segment.residues_.push_back(cur_residue);
			cur_model.segments_.push_back(cur_segment);

			cur_segment = SewSegment();
			cur_segment.model_id_ = model_id;
			cur_segment.segment_id_ = segment_id;
			cur_segment.dssp_ = dssp[0];

			cur_residue = SewResidue();
			cur_residue.resnum_ = seqpos;
			cur_residue.residue_type_ = res_type;
			if ( chi1 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi1);
			}
			if ( chi2 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi2);
			}
			if ( chi3 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi3);
			}
			if ( chi4 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi4);
			}
		} else if ( seqpos != cur_residue.resnum_ ) {
			//new residue
			cur_segment.residues_.push_back(cur_residue);

			cur_residue = SewResidue();
			cur_residue.resnum_ = seqpos;
			cur_residue.residue_type_ = res_type;
			if ( chi1 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi1);
			}
			if ( chi2 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi2);
			}
			if ( chi3 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi3);
			}
			if ( chi4 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi4);
			}
		}

		SewAtom atom;
		numeric::xyzVector<core::Real> coords(x,y,z);
		atom.atomno_ = atomno;
		atom.coords_ = coords;
		cur_residue.basis_atoms_.push_back(atom);
	}
	cur_segment.residues_.push_back(cur_residue);
	cur_model.segments_.push_back(cur_segment);
	models[cur_model.model_id_] = cur_model;
	return models;
} //result_to_models // for original model with 3 secondary_structure_segments


std::map< int, Model >
result_to_five_ss_models(
	cppdb::result res
){
	// As of 2015/10/26, this result_to_five_ss_models fn is not super-perfect, since some ~0.5% of model may have L at terminal segment.
	// TR << "result_to_five_ss_models " << std::endl; // for now, just 5 ss
	std::string input_tag, res_type, dssp;
	int model_id;
	core::Size struct_id, segment_id, seqpos, atomno;
	//core::Real distance, hoist, packing, meridian, chi1, chi2, chi3, chi4, x, y, z;
	core::Real chi1, chi2, chi3, chi4, x, y, z;

	utility::vector1<int> residue_numbers;

	std::map< int, Model > models;

	bool first = true;
	core::Size offset_of_model_id = -1; // due_to_different_struct_id or starting model_id
	Model cur_model;
	SewSegment cur_segment = SewSegment();
	SewResidue cur_residue;
	int new_model_id = 0; // initial value
	int segment_id_count_in_each_model_id = 1; // intentional meaningful initial value

	if ( TR.Debug.visible() ) {
		TR.Debug << "segment_id_count_in_each_model_id before res.next: " << segment_id_count_in_each_model_id << std::endl;
	}

	while ( res.next() ) {
		//res >> input_tag >> struct_id >> model_id >> distance >> hoist >> packing >> meridian >> segment_id >> dssp >> res_type >> chi1 >> chi2 >> chi3 >> chi4 >> seqpos >> atomno >> x >> y >> z;
		res >> input_tag >> struct_id >> model_id >> segment_id >> dssp >> res_type >> chi1 >> chi2 >> chi3 >> chi4 >> seqpos >> atomno >> x >> y >> z;

		if ( TR.Debug.visible() ) {
			TR.Debug << "struct_id: " << struct_id << ", model_id: " << model_id << ", segment_id: " << segment_id << ", dssp: " << dssp << ", seqpos: " << seqpos << ", atomno: " << atomno << std::endl;
			TR.Debug << "cur_segment.segment_id_: " << cur_segment.segment_id_ << std::endl;
			TR.Debug << "segment_id_count_in_each_model_id before if statement: " << segment_id_count_in_each_model_id << std::endl;
		}

		if ( cur_segment.segment_id_ != 0 && segment_id != cur_segment.segment_id_ ) {
			segment_id_count_in_each_model_id++;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "segment_id_count_in_each_model_id after if statement: " << segment_id_count_in_each_model_id << std::endl;
		}

		if ( dssp != "L" && segment_id_count_in_each_model_id==3 ) {
			int result = std::count( residue_numbers.begin(), residue_numbers.end(), seqpos );
			if ( TR.Debug.visible() ) {
				TR.Debug << "Number of seqpos in residue_numbers_vector: " << result << std::endl;
				TR.Debug << std::endl;
			}

			if ( result > 3 ) {
				continue;
			}
			residue_numbers.push_back(seqpos);
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << std::endl;
		}

		if ( segment_id_count_in_each_model_id==4 ) {
			residue_numbers.clear();
		}

		if ( cur_model.structure_id_ != struct_id ) {
			offset_of_model_id++;
		}

		core::Size remnant_after_model_id_divided_by_2 = model_id%2 ;

		if ( remnant_after_model_id_divided_by_2 == 1 ) {
			new_model_id = (model_id+1)/2 + offset_of_model_id;
		} else {
			new_model_id = (model_id)/2 + offset_of_model_id;
		}

		//new model
		if ( new_model_id != cur_model.model_id_ ) {
			if ( !first ) {
				segment_id_count_in_each_model_id = 1;
				cur_segment.residues_.push_back(cur_residue);
				cur_model.segments_.push_back(cur_segment);
				models[cur_model.model_id_] = cur_model;
			} else {
				first=false;
			}

			cur_model = Model();

			cur_model.model_id_ = new_model_id;

			cur_model.structure_id_ = struct_id;
			cur_model.pdb_code_ = input_tag;
			//   cur_model.distance_ = distance;
			//   cur_model.hoist_angle_degrees_ = hoist;
			//   cur_model.packing_angle_degrees_ = packing;
			//   cur_model.meridian_angle_degrees_ = meridian;

			cur_segment = SewSegment();

			cur_segment.model_id_ = new_model_id;

			cur_segment.segment_id_ = segment_id;
			cur_segment.dssp_ = dssp[0];

			cur_residue = SewResidue();
			cur_residue.resnum_ = seqpos;
			cur_residue.residue_type_ = res_type;
			if ( chi1 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi1);
			}
			if ( chi2 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi2);
			}
			if ( chi3 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi3);
			}
			if ( chi4 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi4);
			}
		} else if ( segment_id != cur_segment.segment_id_ ) {
			//new segment

			cur_segment.residues_.push_back(cur_residue);
			cur_model.segments_.push_back(cur_segment);

			cur_segment = SewSegment();

			cur_segment.model_id_ = new_model_id;

			cur_segment.segment_id_ = segment_id;
			cur_segment.dssp_ = dssp[0];

			cur_residue = SewResidue();
			cur_residue.resnum_ = seqpos;
			cur_residue.residue_type_ = res_type;
			if ( chi1 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi1);
			}
			if ( chi2 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi2);
			}
			if ( chi3 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi3);
			}
			if ( chi4 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi4);
			}
		} else if ( seqpos != cur_residue.resnum_ ) {
			//new residue
			cur_segment.residues_.push_back(cur_residue);

			cur_residue = SewResidue();
			cur_residue.resnum_ = seqpos;
			cur_residue.residue_type_ = res_type;
			if ( chi1 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi1);
			}
			if ( chi2 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi2);
			}
			if ( chi3 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi3);
			}
			if ( chi4 != 0.0 ) {
				cur_residue.chi_angles_.push_back(chi4);
			}
		}

		SewAtom atom;
		numeric::xyzVector<core::Real> coords(x,y,z);
		atom.atomno_ = atomno;
		atom.coords_ = coords;
		cur_residue.basis_atoms_.push_back(atom);
	} // while(res.next())

	cur_segment.residues_.push_back(cur_residue);
	cur_model.segments_.push_back(cur_segment);
	models[cur_model.model_id_] = cur_model;
	return models;
} //result_to_five_ss_models // for new model definition with 5 secondary_structure_segments

//std::map< int, Model >
//get_discontinuous_models_from_db(){
//
// utility::sql_database::sessionOP db_session( basic::database::get_db_session() );
//
// std::string select_string =
//  "SELECT s.input_tag, s.struct_id, sm.model_id, ss.segment_id, r.res_type, prc.chi1, prc.chi2, prc.chi3, prc.chi4, coords.seqpos, coords.atomno, coords.x, coords.y, coords.z\n"
//  "FROM sewing_models sm\n"
//  "JOIN structures s ON\n"
//  " sm.struct_id = s.struct_id\n"
//  "JOIN model_segments ms ON\n"
//  " ms.model_id = sm.model_id\n"
//  "JOIN secondary_structure_segments ss ON\n"
//  " s.struct_id = ss.struct_id AND\n"
//  " ms.segment_id = ss.segment_id\n"
//  "JOIN residues r ON\n"
//  " s.struct_id = r.struct_id AND\n"
//  " r.resnum BETWEEN ss.residue_begin AND ss.residue_end\n"
//  "JOIN protein_residue_conformation prc ON\n"
//  " s.struct_id = prc.struct_id AND\n"
//  " r.resnum = prc.seqpos\n"
//  "JOIN residue_atom_coords coords ON\n"
//  " s.struct_id = coords.struct_id AND\n"
//  " r.resnum = coords.seqpos\n"
//  "WHERE\n"
//  " coords.atomno IN (1,2,3,4)\n"
//  "ORDER BY s.struct_id, sm.model_id, ss.segment_id, coords.seqpos, coords.atomno;";
// cppdb::statement select_stmt=basic::database::safely_prepare_statement(select_string, db_session);
// cppdb::result res=basic::database::safely_read_from_database(select_stmt);
//
// TR << "Done selecting from database " << std::endl;
// return result_to_models(res);
//}

std::map< int, Model >
get_discontinuous_models_from_db(){

	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );

	std::string select_string =
		"SELECT s.input_tag, s.struct_id, sm.model_id,'1','1','1','1',ss.segment_id, ss.dssp, r.res_type, prc.chi1, prc.chi2, prc.chi3, prc.chi4, coords.seqpos, coords.atomno, coords.x, coords.y, coords.z\n"
		"FROM sewing_models sm\n"
		"JOIN structures s ON\n"
		"\tsm.struct_id = s.struct_id\n"
		"JOIN model_segments ms ON\n"
		"\tms.model_id = sm.model_id\n"
		"JOIN secondary_structure_segments ss ON\n"
		"\ts.struct_id = ss.struct_id AND\n"
		"\tms.segment_id = ss.segment_id\n"
		"JOIN residues r ON\n"
		"\ts.struct_id = r.struct_id AND\n"
		"\tr.resnum BETWEEN ss.residue_begin AND ss.residue_end\n"
		"JOIN protein_residue_conformation prc ON\n"
		"\ts.struct_id = prc.struct_id AND\n"
		"\tr.resnum = prc.seqpos\n"
		"JOIN residue_atom_coords coords ON\n"
		"\ts.struct_id = coords.struct_id AND\n"
		"\tr.resnum = coords.seqpos\n"
		"WHERE\n"
		" coords.atomno IN (1,2,3,4) AND\n"
		" (SELECT COUNT(*) FROM model_segments ms2\n"
		" JOIN secondary_structure_segments ss2 ON\n"
		"\t\ts.struct_id = ss2.struct_id AND\n"
		"\t\tms2.segment_id = ss2.segment_id AND\n"
		"\t\tss2.dssp = 'E'\n"
		"\tWHERE\n"
		"\t\tms.model_id = ms2.model_id\n"
		"\t) = 0\n"
		"ORDER BY s.struct_id, sm.model_id, ss.segment_id, coords.seqpos, coords.atomno;";
	cppdb::statement select_stmt=basic::database::safely_prepare_statement(select_string, db_session);
	cppdb::result res=basic::database::safely_read_from_database(select_stmt);

	TR << "Done selecting from database " << std::endl;
	return result_to_models(res);
}//get_discontinuous_models_from_db

std::map< int, Model >
get_strand_sew_models_from_db(){

	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );

	std::string select_string =
		"SELECT s.input_tag, s.struct_id, sm.model_id, ss.segment_id, r.res_type, prc.chi1, prc.chi2, prc.chi3, prc.chi4, coords.seqpos, coords.atomno, coords.x, coords.y, coords.z\n"
		"FROM sewing_models sm\n"
		"JOIN structures s ON\n"
		"\tsm.struct_id = s.struct_id\n"
		"JOIN model_segments ms ON\n"
		"\tms.model_id = sm.model_id\n"
		"JOIN secondary_structure_segments ss ON\n"
		"\ts.struct_id = ss.struct_id AND\n"
		"\tms.segment_id = ss.segment_id\n"
		"JOIN residues r ON\n"
		"\ts.struct_id = r.struct_id AND\n"
		"\tr.resnum BETWEEN ss.residue_begin AND ss.residue_end\n"
		"JOIN residue_atom_coords coords ON\n"
		"\ts.struct_id = coords.struct_id AND\n"
		"\tr.resnum = coords.seqpos\n"
		"JOIN protein_residue_conformation prc ON\n"
		"\ts.struct_id = prc.struct_id AND\n"
		"\tr.resnum = prc.seqpos\n"
		"WHERE\n"
		" coords.atomno IN (1,2,3,4) AND\n"
		" (SELECT COUNT(*) FROM model_segments ms2\n"
		" JOIN secondary_structure_segments ss2 ON\n"
		"\t\ts.struct_id = ss2.struct_id AND\n"
		"\t\tms2.segment_id = ss2.segment_id AND\n"
		"\t\tss2.dssp = 'H'\n"
		"\tWHERE\n"
		"\t\tms.model_id = ms2.model_id\n"
		"\t) = 0\n"
		"ORDER BY s.struct_id, sm.model_id, ss.segment_id, coords.seqpos, coords.atomno;";
	cppdb::statement select_stmt=basic::database::safely_prepare_statement(select_string, db_session);
	cppdb::result res=basic::database::safely_read_from_database(select_stmt);

	TR << "Done selecting from database " << std::endl;
	return result_to_models(res);
} //get_strand_sew_models_from_db

std::map< int, Model >
get_continuous_models_from_db(std::string hash_between){

	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );

	std::string select_string =
		"SELECT s.input_tag, s.struct_id, sm.smotif_id as model_id, sm.distance, sm.hoist, sm.packing, sm.meridian,\n"
		" ss.segment_id, ss.dssp,\n"
		" r.res_type, prc.chi1, prc.chi2, prc.chi3, prc.chi4, coords.seqpos,\n"
		" coords.atomno, coords.x, coords.y, coords.z\n"
		"FROM smotifs sm\n"
		"JOIN structures s ON\n"
		"\tsm.struct_id = s.struct_id\n"
		"JOIN secondary_structure_segments ss ON\n"
		"\ts.struct_id = ss.struct_id AND\n"
		" ss.segment_id IN (sm.secondary_struct_segment_id_1, sm.secondary_struct_segment_id_2, loop_segment_id)\n"
		"JOIN residues r ON\n"
		"\ts.struct_id = r.struct_id AND\n"
		"\tr.resnum BETWEEN ss.residue_begin AND ss.residue_end\n"
		"JOIN residue_atom_coords coords ON\n"
		"\ts.struct_id = coords.struct_id AND\n"
		"\tr.resnum = coords.seqpos\n"
		"JOIN protein_residue_conformation prc ON\n"
		"\ts.struct_id = prc.struct_id AND\n"
		"\tr.resnum = prc.seqpos\n"
		"WHERE\n"
		" coords.atomno IN (1,2,3,4) \n"
		"ORDER BY s.struct_id, sm.smotif_id, ss.segment_id, coords.seqpos, coords.atomno;";

	cppdb::statement select_stmt=basic::database::safely_prepare_statement(select_string, db_session);
	cppdb::result res=basic::database::safely_read_from_database(select_stmt);

	std::map< int, Model > models = result_to_models(res);

	//Now make sure we aren't hashing the linker segments
	std::string linker_id_select =
		"SELECT sm.smotif_id as model_id, ss.segment_id\n"
		"FROM smotifs sm\n"
		"JOIN structures s ON\n"
		"\tsm.struct_id = s.struct_id\n"
		"JOIN secondary_structure_segments ss ON\n"
		"\ts.struct_id = ss.struct_id AND\n"
		" ss.segment_id = loop_segment_id\n;";
	cppdb::statement linker_id_select_stmt=basic::database::safely_prepare_statement(linker_id_select, db_session);
	cppdb::result linker_id_res=basic::database::safely_read_from_database(linker_id_select_stmt);

	std::set< std::pair<core::Size, core::Size> > linker_segs;
	while ( linker_id_res.next() ) {
		int model_id;
		core::Size ss_id;
		linker_id_res >> model_id >> ss_id;
		linker_segs.insert(std::make_pair(model_id, ss_id));
	}
	TR << "Found " << linker_segs.size() << " linker segments" << std::endl;
	TR << "hash_between: " << hash_between << std::endl;

	std::map< int, Model >::iterator it = models.begin();
	std::map< int, Model >::iterator it_end = models.end();
	for ( ; it != it_end; ++it ) {
		Model & cur_model = it->second;
		for ( core::Size i=1; i <= cur_model.segments_.size(); ++i ) {
			if ( hash_between == "hash_between_any_HEs" ) {
				if ( linker_segs.find(std::make_pair(cur_model.model_id_, cur_model.segments_[i].segment_id_)) != linker_segs.end() ) {
					cur_model.segments_[i].hash_ = false;
				}
			} else { // (hash_between == "hash_tag_only_terminal_Es")
				if ( ((i != 1) && (i != cur_model.segments_.size()))
						|| ((cur_model.segments_[i].dssp_) != 'E') ) {
					cur_model.segments_[i].hash_ = false;
				}
			}
		}
	}

	return models;
}//get_continuous_models_from_db




std::map< int, Model >
get_5_ss_models_from_db(std::string hash_between){

	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );

	std::string select_string =
		"SELECT s.input_tag, s.struct_id, sm.smotif_id as model_id,\n"
		" ss.segment_id, ss.dssp,\n"
		" r.res_type, prc.chi1, prc.chi2, prc.chi3, prc.chi4, coords.seqpos,\n"
		" coords.atomno, coords.x, coords.y, coords.z\n"
		"FROM smotifs sm\n"
		"JOIN structures s ON\n"
		"\tsm.struct_id = s.struct_id\n"
		"JOIN secondary_structure_segments ss ON\n"
		"\ts.struct_id = ss.struct_id AND\n"
		" ss.segment_id IN (sm.secondary_struct_segment_id_1, sm.secondary_struct_segment_id_2, loop_segment_id)\n"
		"JOIN residues r ON\n"
		"\ts.struct_id = r.struct_id AND\n"
		"\tr.resnum BETWEEN ss.residue_begin AND ss.residue_end\n"
		"JOIN residue_atom_coords coords ON\n"
		"\ts.struct_id = coords.struct_id AND\n"
		"\tr.resnum = coords.seqpos\n"
		"JOIN protein_residue_conformation prc ON\n"
		"\ts.struct_id = prc.struct_id AND\n"
		"\tr.resnum = prc.seqpos\n"
		"WHERE\n"
		" coords.atomno IN (1,2,3,4) \n"
		"ORDER BY s.struct_id, sm.smotif_id, ss.segment_id, coords.seqpos, coords.atomno;";

	cppdb::statement select_stmt=basic::database::safely_prepare_statement(select_string, db_session);
	cppdb::result res=basic::database::safely_read_from_database(select_stmt);

	// for new model definition with 5 secondary_structure_segments
	std::map< int, Model > models_w_5_sss = result_to_five_ss_models(res);

	////////////// for models_w_5_sss
	//Now make sure we aren't hashing the linker segments
	std::map< int, Model >::iterator it_5 = models_w_5_sss.begin();
	std::map< int, Model >::iterator it_end_5 = models_w_5_sss.end();
	for ( ; it_5 != it_end_5; ++it_5 ) {
		Model & cur_model = it_5->second;
		for ( core::Size i=1; i <= cur_model.segments_.size(); ++i ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "hash_between: " << hash_between << std::endl;
			}
			if ( hash_between == "hash_between_any_HEs" ) {
				if ( ((i != 1) && (i != cur_model.segments_.size()))
						|| ((cur_model.segments_[i].dssp_) == 'L') ) {
					cur_model.segments_[i].hash_ = false;
				}
			} else { // (hash_between == "hash_tag_only_terminal_Es")
				if ( ((i != 1) && (i != cur_model.segments_.size()))
						|| ((cur_model.segments_[i].dssp_) != 'E') ) {
					cur_model.segments_[i].hash_ = false;
				}
			}
		}
	}

	return models_w_5_sss;
} //get_5_ss_models_from_db


///Will only work for smotif models!!!
void
remove_models_from_dssp(
	std::map< int, Model > & models,
	char dssp1,
	char dssp2
) {
	std::set<int> invalid_model_ids;
	std::map< int, Model >::iterator it = models.begin();
	std::map< int, Model >::iterator it_end = models.end();
	for ( ; it != it_end; ++it ) {
		Model const & cur_model = it->second;
		if ( cur_model.segments_.size() != 3 ) {
			utility_exit_with_message("Current model " + utility::to_string(cur_model.model_id_) + " has " + utility::to_string(cur_model.segments_.size()) + " segments! Needs 3");
		}
		if ( cur_model.segments_[1].dssp_ == dssp1 && cur_model.segments_[3].dssp_ == dssp2 ) {
			invalid_model_ids.insert(cur_model.model_id_);
		}
	}

	TR << "Removing " << invalid_model_ids.size() << " models matching DSSP pattern " << dssp1 << "-" << dssp2 << std::endl;
	std::set<int>::const_iterator remove_it = invalid_model_ids.begin();
	std::set<int>::const_iterator remove_it_end = invalid_model_ids.end();
	for ( ; remove_it != remove_it_end; ++remove_it ) {
		models.erase(models.find(*remove_it));
	}
}//remove_models_from_dssp

void
add_num_neighbors(
	std::map< int, Model > & /*models*/
){
	// core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	//
	// std::map< int, Model >::iterator it = models.begin();
	// std::map< int, Model >::iterator it_end = models.end();
	// for(; it != it_end; ++it) {
	//
	//  AssemblyOP model_assembly = new DisembodiedAssembly(it->second);
	//
	//  core::pose::Pose pose = model_assembly->to_pose(core::chemical::FA_STANDARD);
	//  scorefxn->score(pose);
	//  //pose.update_residue_neighbors();
	//
	//  Model & cur_model = it->second;
	//  for(core::Size i=1; i<=cur_model.segments_.size(); ++i) {
	//   for(core::Size j=1; j<=cur_model.segments_[i].residues_.size(); ++j) {
	//    core::Size pose_num = model_assembly->pose_num(cur_model.model_id_, cur_model.segments_[i].residues_[j].resnum_);
	//    core::Size num_neighbors = pose.energies().tenA_neighbor_graph().get_node(pose_num)->num_neighbors_counting_self();
	//    cur_model.segments_[i].residues_[j].num_neighbors_ = num_neighbors;
	//    //TR << "Neighbs: " << num_neighbors << it->second.segments_[i].residues_[j].num_neighbors_;
	//   }
	//  }



	//  for(core::Size i=model_assembly.segments_.size(); i<=pose.size(); ++i) {
	//
	//  }

	//  ModelIterator<SewSegment> model_it1 = it->second.model_begin();
	//  ModelIterator<SewSegment> model_end = it->second.model_end();
	//  for(; model_it1 != model_end; ++model_it1) {
	//   if(model_it1.atom()->atomno_ == 2) {
	//
	//    std::set<core::Size> neighbors;
	//    ModelIterator<SewSegment> model_it2 = model_it1;
	//    ++model_it2;
	//    for(; model_it2 != model_end; ++model_it2) {
	//     if(model_it2.atom()->atomno_ != 2) {
	//      core::Real dist2 = model_it1.atom()->coords_.distance_squared(model_it2.atom()->coords_);
	//      if(dist2 < dist_cutoff2){
	//       neighbors.insert(model_it2.residue()->resnum_);
	//      }
	//     }
	//    }
	//
	//    model_it1.residue()->num_neighbors_ = neighbors.size();
	//   }
	//  }
	// }
}


///@details Go through the list of models and add extra 'linker' segments
///to Models that have secondary structure segments that are separated by a
///single "linker" segments. These linker segments will not be hashed by the
///Hasher class. Segments with linkers will have their next_segment_ populated
///so that the order is maintained in the assembly
void
add_linker_segments(
	std::map< int, Model > & models
){
	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );

	std::string select_string =
		"SELECT r.res_type, prc.chi1, prc.chi2, prc.chi3, prc.chi4, coords.seqpos, coords.atomno, coords.x, coords.y, coords.z\n"
		"FROM structures s\n"
		"JOIN secondary_structure_segments ss ON\n"
		"\ts.struct_id = ss.struct_id\n"
		"JOIN residues r ON\n"
		"\ts.struct_id = r.struct_id AND\n"
		"\tr.resnum BETWEEN ss.residue_begin AND ss.residue_end\n"
		"JOIN protein_residue_conformation prc ON\n"
		"\ts.struct_id = prc.struct_id AND\n"
		"\tr.resnum = prc.seqpos\n"
		"JOIN residue_atom_coords coords ON\n"
		"\ts.struct_id = coords.struct_id AND\n"
		"\tr.resnum = coords.seqpos\n"
		"WHERE\n"
		" s.struct_id = ? AND\n"
		"\tss.segment_id = ? AND\n"
		" coords.atomno IN (1,2,3,4)\n"
		"ORDER BY coords.seqpos, coords.atomno;";
	cppdb::statement select_stmt=basic::database::safely_prepare_statement(select_string, db_session);

	std::map< int, Model >::iterator it = models.begin();
	std::map< int, Model >::iterator it_end = models.end();
	for ( ; it != it_end; ++it ) {
		Model & cur_model = it->second;
		utility::vector1<SewSegment> linker_segments;
		for ( core::Size lower=1; lower <= cur_model.segments_.size(); ++lower ) {
			for ( core::Size upper=1; upper <= cur_model.segments_.size(); ++upper ) {
				//TR << "Checking segments " << cur_model.segments_[i+1].segment_id_ << " and " << cur_model.segments_[i].segment_id_ << std::endl;
				if ( cur_model.segments_[upper].segment_id_ - cur_model.segments_[lower].segment_id_ == 2 ) {
					TR << "Getting linker segment for model " << cur_model.model_id_ << std::endl;
					//get segment
					core::Size linker_segment_id = cur_model.segments_[lower].segment_id_+1;
					select_stmt.bind(1, cur_model.structure_id_);
					select_stmt.bind(2, linker_segment_id);
					cppdb::result res=basic::database::safely_read_from_database(select_stmt);

					std::string res_type;
					core::Size seqpos, atomno;
					core::Real chi1, chi2, chi3, chi4, x, y, z;

					bool first = true;
					SewSegment cur_segment = SewSegment();
					cur_segment.model_id_ = cur_model.model_id_;
					cur_segment.segment_id_ = linker_segment_id;
					cur_segment.hash_ = false;
					SewResidue cur_residue;
					while ( res.next() ) {
						res >> res_type >> chi1 >> chi2 >> chi3 >> chi4 >> seqpos >> atomno >> x >> y >> z;

						//new residue
						if ( seqpos != cur_residue.resnum_ ) {
							if ( !first ) {
								cur_segment.residues_.push_back(cur_residue);
							} else {
								first = false;
							}

							cur_residue = SewResidue();
							cur_residue.resnum_ = seqpos;
							cur_residue.residue_type_ = res_type;
							if ( chi1 != 0.0 ) {
								cur_residue.chi_angles_.push_back(chi1);
							}
							if ( chi2 != 0.0 ) {
								cur_residue.chi_angles_.push_back(chi2);
							}
							if ( chi3 != 0.0 ) {
								cur_residue.chi_angles_.push_back(chi3);
							}
							if ( chi4 != 0.0 ) {
								cur_residue.chi_angles_.push_back(chi4);
							}
						}
						SewAtom atom;
						numeric::xyzVector<core::Real> coords(x,y,z);
						atom.atomno_ = atomno;
						atom.coords_ = coords;
						cur_residue.basis_atoms_.push_back(atom);
					}
					//Add the last residue
					cur_segment.residues_.push_back(cur_residue);

					//Add connections
					cur_model.segments_.add_connection(lower, cur_model.segments_.size());
					cur_model.segments_.add_connection(cur_model.segments_.size(), upper);
					linker_segments.push_back(cur_segment);
				}
			}
		}
		cur_model.segments_.insert(cur_model.segments_.end(), linker_segments.begin(), linker_segments.end());
	}
	TR << "Done adding extra segment" << std::endl;
}//add_linker_segments

//New function for beta-alpha-beta models -- unfinished
std::map< int, Model >
get_alpha_beta_models_from_db(){

	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );

	std::string select_string =
		"SELECT s.input_tag, s.struct_id, sm.model_id, ss.segment_id, r.res_type, prc.chi1, prc.chi2, prc.chi3, prc.chi4, coords.seqpos, coords.atomno, coords.x, coords.y, coords.z\n"
		"FROM sewing_models sm\n"
		"JOIN structures s ON\n"
		"\tsm.struct_id = s.struct_id\n"
		"JOIN secondary_structure_segments ss ON\n"
		"\ts.struct_id = ss.struct_id AND\n"
		"\tsm.ss_segment_id = ss.segment_id\n"
		"JOIN residues r ON\n"
		"\ts.struct_id = r.struct_id AND\n"
		"\tr.resnum BETWEEN ss.residue_begin AND ss.residue_end\n"
		"JOIN residue_atom_coords coords ON\n"
		"\ts.struct_id = coords.struct_id AND\n"
		"\tr.resnum = coords.seqpos\n"
		"JOIN protein_residue_conformation prc ON\n"
		"\ts.struct_id = prc.struct_id AND\n"
		"\tr.resnum = prc.seqpos\n"
		"WHERE\n"
		" coords.atomno IN (1,2,3,4)\n"
		"ORDER BY s.struct_id, sm.model_id, ss.segment_id, coords.seqpos, coords.atomno;";
	cppdb::statement select_stmt=basic::database::safely_prepare_statement(select_string, db_session);
	cppdb::result res=basic::database::safely_read_from_database(select_stmt);

	std::map< int, Model > models = result_to_models(res);

	//Now make sure we aren't hashing the linker segments
	std::string linker_id_select =
		"SELECT sm.model_id, sm.ss_segment_id\n"
		"FROM sewing_models sm\n"
		"JOIN structures s ON\n"
		"\tsm.struct_id = s.struct_id\n"
		"JOIN secondary_structure_segments ss ON\n"
		"\ts.struct_id = ss.struct_id AND\n"
		"\ts.segment_id = sm.ss_segment_id\n"
		"WHERE\n"
		"\tss.dssp IN ('L')\n";
	cppdb::statement linker_id_select_stmt=basic::database::safely_prepare_statement(linker_id_select, db_session);
	cppdb::result linker_id_res=basic::database::safely_read_from_database(linker_id_select_stmt);

	std::set< std::pair<core::Size, core::Size> > linker_segs;
	while ( linker_id_res.next() ) {
		int model_id;
		core::Size ss_id;
		linker_id_res >> model_id >> ss_id;
		linker_segs.insert(std::make_pair(model_id, ss_id));
	}
	TR << "Found " << linker_segs.size() << " linker segments" << std::endl;

	std::map< int, Model >::iterator it = models.begin();
	std::map< int, Model >::iterator it_end = models.end();
	for ( ; it != it_end; ++it ) {
		Model & cur_model = it->second;
		for ( core::Size i=1; i <= cur_model.segments_.size(); ++i ) {
			if ( linker_segs.find(std::make_pair(cur_model.model_id_, cur_model.segments_[i].segment_id_)) != linker_segs.end() ) {
				cur_model.segments_[i].hash_ = false;
			}
		}
	}

	return models;
}//get_alpha_beta_models_from_db

struct segment{
	core::Size struct_id;
	core::Size segment_id;
	std::string dssp;
};

//New function for beta-alpha-beta models -- unfinished
void
create_alpha_beta_models_table(){

	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );

	using namespace basic::database::schema_generator;
	using namespace utility;

	//******smotifs******//
	Column model_id("model_id", DbDataTypeOP( new DbInteger() ), true);
	Column struct_id_column("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column ss_segment_id("ss_segment_id", DbDataTypeOP( new DbInteger() ), false);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id_column);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	ForeignKey foreign_key1(foreign_key_columns1, "structures", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id_column);
	foreign_key_columns2.push_back(ss_segment_id);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("segment_id");
	ForeignKey foreign_key2(foreign_key_columns2, "secondary_structure_segments", reference_columns2, true);

	Schema sewing_models("sewing_models", model_id);
	sewing_models.add_column(struct_id_column);
	sewing_models.add_column(ss_segment_id);
	sewing_models.add_foreign_key(foreign_key1);
	sewing_models.add_foreign_key(foreign_key2);

	sewing_models.write(db_session);

	std::string sewing_model_insert_string =
		"INSERT INTO sewing_models (struct_id, ss_segment_id) VALUES(?,?)";

	cppdb::statement sewing_model_insert_string_stmt =
		basic::database::safely_prepare_statement(sewing_model_insert_string, db_session);

	std::string select_string =
		"SELECT s.input_tag, s.struct_id, ss.segment_id, ss.dssp\n"
		"FROM sewing_models sm\n"
		"JOIN structures s ON\n"
		"\tsm.struct_id = s.struct_id\n"
		"JOIN secondary_structure_segments ss ON\n"
		"\ts.struct_id = ss.struct_id AND\n"
		"\tsm.ss_segment_id = ss.segment_id\n"
		"ORDER BY ss.struct_id, sm.model_id, ss.segment_id;";
	cppdb::statement select_stmt=basic::database::safely_prepare_statement(select_string, db_session);
	cppdb::result res=basic::database::safely_read_from_database(select_stmt);

	core::Size input_tag, struct_id, segment_id;
	std::string dssp;

	utility::vector1<segment> segments;
	while ( res.next() ) {
		res >> input_tag >> struct_id >> segment_id >> dssp;
		segment new_segment;
		new_segment.struct_id = struct_id;
		new_segment.segment_id = segment_id;
		new_segment.dssp = dssp;

		if ( segments.size() < 5 ) {
			segments.push_back(new_segment);
		} else {
			segments.erase(segments.begin()); //erase first element
			segments.push_back(new_segment);

			runtime_assert(segments.size() == 5);
			std::set<core::Size> current_struct_ids;
			std::string sequence;
			for ( core::Size i=1; i<=segments.size(); ++i ) {
				//check to see if all segments have same struct id and the order matches secondary structure.
				//if it does match, set valid to true
				current_struct_ids.insert(segments[i].struct_id);
				sequence += segments[i].dssp;
			}
			if ( current_struct_ids.size() == 1 && sequence == "ELHLE" ) {
				for ( core::Size i=1; i<=segments.size(); ++i ) {
					sewing_model_insert_string_stmt.bind(1, segments[i].struct_id);
					sewing_model_insert_string_stmt.bind(1, segments[i].segment_id);
					basic::database::safely_write_to_database(sewing_model_insert_string_stmt);
				}
			}
		}
	}

	TR << "Done creating new database table" << std::endl;
}//create_alpha_beta_models_table


void
write_model_file(
	std::string const & comments,
	std::map< int, Model > models,
	std::string filename
){
	utility::io::ozstream file;
	file.open(filename);

	file << comments << std::endl;
	for ( std::map< int, Model >::const_iterator it = models.begin(); it != models.end(); ++it ) {
		Model cur_model = it->second;
		file << "MODEL " << cur_model.model_id_ << " " << cur_model.structure_id_ << " "
			<< cur_model.distance_ << " " << cur_model.hoist_angle_degrees_ << " " << cur_model.packing_angle_degrees_ << " "
			<< cur_model.meridian_angle_degrees_ << " " << cur_model.pdb_code_ << std::endl;

		for ( core::Size segment_i=1; segment_i <= cur_model.segments_.size(); ++segment_i ) {
			SewSegment cur_segment = cur_model.segments_[segment_i];
			file << "SEGMENT " << cur_segment.segment_id_ << " " << cur_segment.dssp_ << " " << cur_segment.hash_ << std::endl;

			for ( core::Size residue_i=1; residue_i <= cur_segment.residues_.size(); ++residue_i ) {
				SewResidue cur_residue = cur_segment.residues_[residue_i];
				file << "RESIDUE " << cur_residue.resnum_ << " " << cur_residue.num_neighbors_ << " " << cur_residue.residue_type_ << " ";
				for ( core::Size chi_i=1; chi_i<=cur_residue.chi_angles_.size(); ++chi_i ) {
					file << cur_residue.chi_angles_[chi_i];
					if ( chi_i < cur_residue.chi_angles_.size() ) {
						file << " ";
					}
				}
				file << std::endl;

				for ( core::Size atom_i=1; atom_i <= cur_residue.basis_atoms_.size(); ++atom_i ) {
					SewAtom cur_atom = cur_residue.basis_atoms_[atom_i];
					file << "ATOM " << cur_atom.atomno_ << " " << cur_atom.coords_.x() << " " << cur_atom.coords_.y() << " " << cur_atom.coords_.z() << std::endl;
				}
			}
		}
	}
	file.close();
}

std::map<int, Model>
read_model_file(
	std::string filename
){
	core::Size starttime = time(NULL);
	utility::io::izstream file(filename);
	if ( !file.good() ) {
		utility_exit_with_message("Could not find Models file with name: " + filename);
	}

	TR << "Reading models from file: " << filename << std::endl;
	std::map< int, Model > models;
	//boost::unordered_map< int, Model > models;
	std::string line;

	Model cur_model;
	cur_model.model_id_ = 0;

	SewSegment cur_segment;
	SewResidue cur_residue;
	while ( getline( file, line ) ) {

		//skip comment lines and empty lines
		utility::trim(line);
		if ( line.length() == 0 || line[0] == '#' ) {
			continue;
		}

		utility::vector1<std::string> tokens = utility::string_split(line);
		runtime_assert(tokens.size() > 0);
		if ( tokens[1]=="MODEL" ) {
			runtime_assert(tokens.size() == 8);
			if ( cur_model.model_id_ != 0 ) {
				cur_segment.residues_.push_back(cur_residue);
				cur_model.segments_.push_back(cur_segment);
				models[cur_model.model_id_] = cur_model;
			}

			//////THIS IS THE PROBLEM LINE///////
			cur_model = Model();
			cur_segment = SewSegment();
			cur_residue = SewResidue();

			cur_model.model_id_ = utility::string2int(tokens[2]);
			cur_model.structure_id_ = utility::string2int(tokens[3]);
			cur_model.distance_ = utility::string2float(tokens[4]);
			cur_model.hoist_angle_degrees_ = utility::string2float(tokens[5]);
			cur_model.packing_angle_degrees_ = utility::string2float(tokens[6]);
			cur_model.meridian_angle_degrees_ = utility::string2float(tokens[7]);
			cur_model.pdb_code_ = tokens[8];
		} else if ( tokens[1]=="SEGMENT" ) {
			runtime_assert(tokens.size() == 4);
			if ( cur_segment.segment_id_ != 0 && cur_segment.segment_id_ != (core::Size)utility::string2int(tokens[2]) ) {
				cur_segment.residues_.push_back(cur_residue);
				cur_model.segments_.push_back(cur_segment);
			}

			cur_segment = SewSegment();
			cur_residue = SewResidue();

			cur_segment.model_id_ = cur_model.model_id_;
			cur_segment.segment_id_ = utility::string2int(tokens[2]);
			cur_segment.dssp_ = tokens[3][0];
			cur_segment.hash_ = utility::string2int(tokens[4]);
		} else if ( tokens[1]=="RESIDUE" ) {
			runtime_assert(tokens.size() >= 4);
			if ( cur_residue.resnum_ != 0 && cur_residue.resnum_ != (core::Size)utility::string2int(tokens[2]) ) {
				cur_segment.residues_.push_back(cur_residue);
			}

			cur_residue = SewResidue();
			cur_residue.resnum_ = utility::string2int(tokens[2]);
			cur_residue.num_neighbors_ = utility::string2int(tokens[3]);
			cur_residue.residue_type_ = tokens[4];
			for ( core::Size i=5; i<= tokens.size(); ++i ) {
				cur_residue.chi_angles_.push_back(utility::string2float(tokens[i]));
			}
		} else if ( tokens[1]=="ATOM" ) {
			SewAtom atom;
			atom.atomno_ = utility::string2int(tokens[2]);
			core::Real x = utility::string2float(tokens[3]);
			core::Real y = utility::string2float(tokens[4]);
			core::Real z = utility::string2float(tokens[5]);
			atom.coords_ = numeric::xyzVector<core::Real>(x,y,z);
			cur_residue.basis_atoms_.push_back(atom);
		} else {
			utility_exit_with_message("Error: malformed model file!");
		}
	}
	cur_segment.residues_.push_back(cur_residue);
	cur_model.segments_.push_back(cur_segment);
	models[cur_model.model_id_] = cur_model;

	//Keep a vector of model_ids for models containing segments with 1 or fewer residues
	std::set<int> invalid_model_ids;

	//Infer linkers based on segment ids. Consecutive segments = linked
	std::map< int, Model >::iterator it = models.begin();
	std::map< int, Model >::iterator it_end = models.end();
	for ( ; it != it_end; ++it ) {
		Model & curr_model = it->second;
		for ( core::Size i=1; i<=curr_model.segments_.size(); ++i ) {
			for ( core::Size j=1; j<=curr_model.segments_.size(); ++j ) {
				if ( curr_model.segments_[i].segment_id_ == (curr_model.segments_[j].segment_id_ - 1) ) {
					curr_model.segments_.add_connection(i, j);
				}
			}

			if ( curr_model.segments_[i].residues_.size() <= 1 ) {
				// if any segment in a model is constitued with less than 2 residues, then don't care this model, if cared, it will crash generating too high atom numbers
				invalid_model_ids.insert(curr_model.model_id_);
				continue;
			}
		}
	}

	std::set<int>::const_iterator remove_it = invalid_model_ids.begin();
	std::set<int>::const_iterator remove_it_end = invalid_model_ids.end();
	for ( ; remove_it != remove_it_end; ++remove_it ) {
		models.erase(models.find(*remove_it));
	}

	core::Size endttime = time(NULL);
	TR << "Read " << models.size() << " models in " << endttime - starttime << " seconds" << std::endl;

	return models;
} //read_model_file


core::Size
Model::pose_number(
	core::Size resnum
) const {
	core::Size counter=1;
	for ( core::Size i=1; i<=segments_.size(); ++i ) {
		for ( core::Size j=1; j<=segments_[i].residues_.size(); ++j ) {
			if ( segments_[i].residues_[j].resnum_ == resnum ) {
				return counter;
			}
			++counter;
		}
	}
	utility_exit_with_message("No residue " + utility::to_string(resnum) + " in model with ID " + utility::to_string(model_id_));
	return 0;
}//pose_number

} //sewing namespace
} //protocols namespace
