// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/MergePDBMover.cc
/// @brief This class will allign & combine parts of the pdb.
/// @author TJ Brunette (tjbrunette@gmail.com)
///

// Unit headers
#include <protocols/simple_moves/MergePDBMoverCreator.hh>
#include <protocols/simple_moves/MergePDBMover.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.MergePDBMover" );

#include <utility/tag/Tag.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/xyzVector.hh>
#include <protocols/toolbox/superimpose.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {

using namespace core;

std::string MergePDBMoverCreator::keyname() const
{
	return MergePDBMoverCreator::mover_name();
}

protocols::moves::MoverOP
MergePDBMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MergePDBMover );
}

std::string
MergePDBMoverCreator::mover_name()
{
	return "MergePDB";
}

MergePDBMover::MergePDBMover()
: moves::Mover("MergePDB"),
	parent_pose_( /* NULL */ )
{
}

std::string
MergePDBMover::get_name() const {
	return MergePDBMoverCreator::mover_name();
}

moves::MoverOP
MergePDBMover::clone() const
{
	return moves::MoverOP( new MergePDBMover( *this ) );
}

moves::MoverOP
MergePDBMover::fresh_instance() const
{
	return moves::MoverOP( new MergePDBMover );
}

void MergePDBMover::determine_overlap(Pose const pose, Size & start_overlap_parent, Size & end_overlap_parent, Size & start_overlap_pose, Size & end_overlap_pose){
	using namespace core::scoring;
	if ( overlap_location_pose_== "c_term" ) {
		Size initial_start_res_parent=1;
		Size end_res_parent=overlap_max_length_;
		Size initial_start_res_pose=pose.total_residue()-overlap_max_length_+1;
		Size end_res_pose=pose.total_residue();
		Real best_rmsd = 9999;
		for ( Size ii=0; ii<overlap_range_; ++ii ) {
			utility::vector1<Size> parent_pose_positions;
			utility::vector1<Size> pose_positions;
			Size start_res_parent = initial_start_res_parent+ii;
			Size start_res_pose = initial_start_res_pose+ii;
			for ( Size jj=start_res_parent; jj<=end_res_parent; ++jj ) {
				parent_pose_positions.push_back(jj);
			}
			for ( Size jj=start_res_pose; jj<=end_res_pose; ++jj ) {
				pose_positions.push_back(jj);
			}
			pose::Pose parent_pose_slice;
			pose::Pose ref_pose_slice;
			pdbslice(parent_pose_slice,*parent_pose_,parent_pose_positions);
			pdbslice(ref_pose_slice,pose,pose_positions);
			Real rmsd = CA_rmsd(ref_pose_slice,parent_pose_slice);
			if ( rmsd < best_rmsd ) {
				best_rmsd = rmsd;
				start_overlap_parent=start_res_parent;
				end_overlap_parent=end_res_parent;
				start_overlap_pose=start_res_pose;
				end_overlap_pose=end_res_pose;
			}
		}
	}
	if ( overlap_location_pose_== "n_term" ) {
		Size start_res_parent=parent_pose_->total_residue()-overlap_max_length_+1;
		Size initial_end_res_parent=parent_pose_->total_residue();
		Size start_res_pose=1;
		Size initial_end_res_pose=overlap_max_length_;
		Real best_rmsd = 9999;
		for ( Size ii=0; ii<overlap_range_; ++ii ) {
			utility::vector1<Size> parent_pose_positions;
			utility::vector1<Size> pose_positions;
			Size end_res_parent = initial_end_res_parent-ii;
			Size end_res_pose = initial_end_res_pose-ii;
			for ( Size jj=start_res_parent; jj<=end_res_parent; ++jj ) {
				parent_pose_positions.push_back(jj);
			}
			for ( Size jj=start_res_pose; jj<=end_res_pose; ++jj ) {
				pose_positions.push_back(jj);
			}
			pose::Pose parent_pose_slice;
			pose::Pose ref_pose_slice;
			pdbslice(parent_pose_slice,*parent_pose_,parent_pose_positions);
			pdbslice(ref_pose_slice,pose,pose_positions);
			Real rmsd = CA_rmsd(ref_pose_slice,parent_pose_slice);
			if ( rmsd < best_rmsd ) {
				best_rmsd = rmsd;
				start_overlap_parent=start_res_parent;
				end_overlap_parent=end_res_parent;
				start_overlap_pose=start_res_pose;
				end_overlap_pose=end_res_pose;

			}
		}
	}
}

void MergePDBMover::generate_overlap(Pose & pose, Size start_overlap_parent, Size end_overlap_parent, Size start_overlap_pose, Size end_overlap_pose){
	using namespace core::id;
	using namespace core::scoring;
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, BOGUS_ATOM_ID );
	for ( Size ii=0; ii<=end_overlap_pose-start_overlap_pose; ++ii ) {
		core::id::AtomID const id1(pose.residue(start_overlap_pose+ii).atom_index("CA"),start_overlap_pose+ii);
		core::id::AtomID const id2(parent_pose_->residue(start_overlap_parent+ii).atom_index("CA"), start_overlap_parent+ii );
		atom_map[id1]=id2;
	}
	superimpose_pose(pose,*parent_pose_,atom_map);
	//create subpose of pose
	utility::vector1<Size> pose_positions;
	utility::vector1<Size> parent_pose_positions;
	if ( overlap_location_pose_== "n_term" ) {
		for ( Size ii=end_overlap_pose+1; ii<=pose.total_residue(); ++ii ) {
			pose_positions.push_back(ii);
		}
		for ( Size ii=1; ii<=end_overlap_parent; ++ii ) {
			parent_pose_positions.push_back(ii);
		}
	}
	if ( overlap_location_pose_== "c_term" ) {
		for ( Size ii=1; ii<start_overlap_pose; ++ii ) {
			pose_positions.push_back(ii);
		}
		for ( Size ii=start_overlap_parent; ii<=parent_pose_->total_residue(); ++ii ) {
			parent_pose_positions.push_back(ii);
		}
	}
	pose::Pose ref_pose_slice;
	pose::Pose parent_pose_slice;

	pdbslice(ref_pose_slice,pose,pose_positions);
	pdbslice(parent_pose_slice,*parent_pose_,parent_pose_positions);
	if ( overlap_location_pose_== "n_term" ) {
		remove_upper_terminus_type_from_pose_residue(parent_pose_slice,parent_pose_slice.total_residue());
		remove_lower_terminus_type_from_pose_residue(ref_pose_slice,1);
		append_pose_to_pose(parent_pose_slice,ref_pose_slice,false);
		pose = parent_pose_slice;
	}
	if ( overlap_location_pose_== "c_term" ) {
		remove_upper_terminus_type_from_pose_residue(ref_pose_slice,ref_pose_slice.total_residue());
		remove_lower_terminus_type_from_pose_residue(parent_pose_slice,1);
		append_pose_to_pose(ref_pose_slice,parent_pose_slice,false);
		pose = ref_pose_slice;
	}
	renumber_pdbinfo_based_on_conf_chains(pose);
}


void MergePDBMover::apply( Pose & pose )
{
	//Step 1 get location
	Size start_overlap_parent =0;
	Size end_overlap_parent=0;
	Size start_overlap_pose=0;
	Size end_overlap_pose=0;
	determine_overlap(pose, start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
	//Step 2 combine
	generate_overlap(pose,start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
}


void MergePDBMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	overlap_location_pose_ = tag->getOption< std::string >( "overlap_location_pose" ,"c_term" );
	overlap_range_ =  tag->getOption< core::Size >( "overlap_range", 5 );
	overlap_max_length_ = tag->getOption< core::Size >( "overlap_max_length", 40);
	std::string fname( tag->getOption< std::string >( "parent_pdb" ) );
	core::import_pose::pose_from_pdbstring(*parent_pose_, fname );
}

} // simple_moves
} // protocols

