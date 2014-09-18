// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/devel/sewing/BridgeFragmentMover.cc
/// @brief 
/// @author Tim Jacobs

//Unit
#include <devel/sewing/BridgeFragmentMover.hh>

//Core
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/rms_util.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>

//Basic
#include <basic/Tracer.hh>

//Utility
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>

namespace devel {
namespace sewing {

static thread_local basic::Tracer TR( "BridgeFragmentMover" );

BridgeFragmentMover::BridgeFragmentMover(utility::vector1<core::fragment::FragSetOP> frag_sets):
frag_sets_(frag_sets),
num_helical_residues_(5),
size_helical_window_(3)
{}

BridgeFragmentMover::~BridgeFragmentMover(){}

void BridgeFragmentMover::apply(core::pose::Pose & pose){

	using namespace std;
	bool reverse=false; //bool to determine whether the helix was built off of the n or c terminus (reverse true means it was built from the N)

	core::Size num_jumps(pose.fold_tree().num_jump());
	TR << "There are " << num_jumps << " jumps" << endl;
	core::Size total_low_rms_counter(0);
	utility::vector1<core::pose::Pose> min_rms_loop_poses;
	for(core::Size cur_jump_edge=1; cur_jump_edge <= num_jumps; ++cur_jump_edge){

		core::Size low_rms_counter(0);
		core::Size downstream_jump_residue(pose.fold_tree().downstream_jump_residue(cur_jump_edge));
		core::Size pose_frag_start(downstream_jump_residue-size_helical_window_);
		core::Size pose_frag_end(downstream_jump_residue+size_helical_window_-1);

		core::Real jump_distance(pose.residue(downstream_jump_residue).atom("CA").xyz().distance(pose.residue(downstream_jump_residue-1).atom("CA").xyz()));
		TR << "jump distance is: " << jump_distance << endl;

		TR << "FRAG START: " << pose_frag_start << endl;
		TR << "FRAG END: " << pose_frag_end << endl;

		core::Size total_size(0);
		for(core::Size frag_sets_index=1; frag_sets_index<=frag_sets_.size(); ++frag_sets_index){
			total_size+=frag_sets_[frag_sets_index]->size();
		}
		TR << "Frag set size: " << total_size << endl;

		core::Real min_rms(10.0);
		core::fragment::AnnotatedFragData min_frag_data("dummy", 1);
		Pose min_frag_pose;
		core::Size min_frag_start_window_offset;
		core::Size min_frag_end_window_offset;
		for(core::Size frag_sets_index=1; frag_sets_index<=frag_sets_.size(); ++frag_sets_index){
			core::fragment::FragSetOP cur_frag_set=frag_sets_[frag_sets_index];

			for(core::fragment::FrameIterator it = cur_frag_set->begin(); it!=cur_frag_set->end(); ++it){

				for(core::Size i = 1; i<=it->nr_frags(); ++i){
					core::fragment::AnnotatedFragData cur_frag =
						dynamic_cast< const core::fragment::AnnotatedFragData & > (it->fragment(i));

					std::string p2a_test("");
					ObjexxFCL::FArray2D< core::Real > pose_helical_coords( 3, size_helical_window_*2 );

					//get coordinates for the helical residues before the jump
					for(core::Size j=0; j<size_helical_window_;j++){
						numeric::xyzVector<core::Real> res_xyz(pose.residue(pose_frag_start+j).atom("CA").xyz());

						pose_helical_coords(1,j+1)=res_xyz.x();
						pose_helical_coords(2,j+1)=res_xyz.y();
						pose_helical_coords(3,j+1)=res_xyz.z();

						p2a_test += utility::to_string(pose.residue(pose_frag_start+j).name1()) + " " +
								utility::to_string(res_xyz.x()) + " " +
								utility::to_string(res_xyz.y()) + " " +
								utility::to_string(res_xyz.z()) + "\n";
					}

					//get coordinates for the helical residues after the jump
					for(int j=size_helical_window_-1; j>=0; j--){
						numeric::xyzVector<core::Real> res_xyz(pose.residue(pose_frag_end-j).atom("CA").xyz());

						pose_helical_coords(1,size_helical_window_*2-j)=res_xyz.x();
						pose_helical_coords(2,size_helical_window_*2-j)=res_xyz.y();
						pose_helical_coords(3,size_helical_window_*2-j)=res_xyz.z();

						p2a_test += utility::to_string(pose.residue(pose_frag_start+j).name1()) + " " +
								utility::to_string(res_xyz.x()) + " " +
								utility::to_string(res_xyz.y()) + " " +
								utility::to_string(res_xyz.z()) + "\n";
					}

					//Check every 'size_helical_window' number of helical residues at the beginning and end of the bridge fragment.
					//If the helical residue window used does come directly before or directly after the loop residues then the extra
					//helical residues will also be added to the pose. The goal of using this window was to ensure that the helix ends
					//are in the optimal position for loop closure without shortening the helices.
					for(core::Size start_window_offset=0; start_window_offset<num_helical_residues_-size_helical_window_+1;
							++start_window_offset){

						ObjexxFCL::FArray2D< core::Real > frag_helical_coords( 3, size_helical_window_*2 );
						std::string p1a_test("");

						//get coordinates for the helical residues on the beginning of the bridge fragment
						for(core::Size j=1; j<=size_helical_window_; j++){
							core::Size cur_resnum = j+start_window_offset;
							core::fragment::BBTorsionSRFDCOP fragment_residue =
									dynamic_cast< const core::fragment::BBTorsionSRFD* > (cur_frag.get_residue(cur_resnum)());

							frag_helical_coords(1,j)=fragment_residue->x();
							frag_helical_coords(2,j)=fragment_residue->y();
							frag_helical_coords(3,j)=fragment_residue->z();

							p1a_test += utility::to_string(fragment_residue->sequence()) + " " +
									utility::to_string(fragment_residue->x()) + " " +
									utility::to_string(fragment_residue->y()) + " " +
									utility::to_string(fragment_residue->z()) + "\n";
						}

						for(core::Size end_window_offset=0; end_window_offset<num_helical_residues_-size_helical_window_+1;
								++end_window_offset){

							//get coordinates for the helical residues on the end of the bridge fragment
							for(core::Size j=0; j<size_helical_window_; j++){
								core::Size cur_resnum = cur_frag.size()-(size_helical_window_-1-j)-end_window_offset;
								core::fragment::BBTorsionSRFDCOP fragment_residue =
														dynamic_cast< const core::fragment::BBTorsionSRFD* > (cur_frag.get_residue(cur_resnum)());

								frag_helical_coords(1,size_helical_window_+j+1)=fragment_residue->x();
								frag_helical_coords(2,size_helical_window_+j+1)=fragment_residue->y();
								frag_helical_coords(3,size_helical_window_+j+1)=fragment_residue->z();

								p1a_test += utility::to_string(fragment_residue->sequence()) + " " +
										utility::to_string(fragment_residue->x()) + " " +
										utility::to_string(fragment_residue->y()) + " " +
										utility::to_string(fragment_residue->z()) + "\n";
							}

							core::Real rms(numeric::model_quality::rms_wrapper( size_helical_window_*2, frag_helical_coords, pose_helical_coords ));

							if(rms<0.8){
								low_rms_counter++;
							}

							//TR << "Fragment RMSD: " << rms << endl;
							if(rms<min_rms){
								it->fragment_as_pose(i,min_frag_pose,core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ));

								min_rms=rms;
								min_frag_data=cur_frag;
								min_frag_start_window_offset=start_window_offset;
								min_frag_end_window_offset=end_window_offset;

								//Remove unused helical elements from the frag_pose
								for(core::Size i=1; i<=min_frag_start_window_offset; ++i){
									min_frag_pose.conformation().delete_residue_slow(1);
								}
								for(core::Size i=1; i<=min_frag_end_window_offset; ++i){
									min_frag_pose.conformation().delete_residue_slow(min_frag_pose.total_residue());
								}

							}
						}
					}
				}
			}
		}

		TR << "Min RMS: " << min_rms << "  " << min_frag_data.size() << "-residue fragment from pdb " << min_frag_data.pdbid() << endl;
		TR << "Start offset" << min_frag_start_window_offset << endl;
		TR << "End offset" << min_frag_end_window_offset << endl;
		TR << low_rms_counter << " " << pose.pdb_info()->name() << " " << cur_jump_edge << " bridge fragments < 0.8"
				<< endl;
		total_low_rms_counter+=low_rms_counter;

		std::string sequence(min_frag_data.sequence());
		TR << "Sequence: " << sequence << endl;

		core::id::AtomID_Map< core::id::AtomID > atom_map;
		// maps every atomid to bogus atom. Atoms in the query_structure that are left mapped to the bogus atoms will not be used in the superimposition
		atom_map.clear();
		core::pose::initialize_atomid_map( atom_map, min_frag_pose, core::id::BOGUS_ATOM_ID );

		//superimpose the found two-helix pose onto the query structure
		for(Size j=0; j< size_helical_window_; j++){

			TR << "superimpose residues " <<  j+1 << " of frag to residue " << j+pose_frag_start << " of pose." << endl;

			core::id::AtomID const id1( min_frag_pose.residue(j+1).atom_index("CA"),j+1 );
			core::id::AtomID const id2( pose.residue(j+pose_frag_start).atom_index("CA"), j+pose_frag_start);
			atom_map.set(id1, id2);

//			core::id::AtomID const id3( pose.residue(j+pose_frag_start).atom_index("C"), j+pose_frag_start);
//			core::id::AtomID const id4( min_frag_pose.residue(j+1).atom_index("C"), j+1 );
//			atom_map[ id4 ] = id3;
//
//			core::id::AtomID const id5( pose.residue(j+pose_frag_start).atom_index("N"), j+pose_frag_start);
//			core::id::AtomID const id6( min_frag_pose.residue(j+1).atom_index("N"), j+1 );
//			atom_map[ id6 ] = id5;
//
//			core::id::AtomID const id7( pose.residue(j+pose_frag_start).atom_index("O"), j+pose_frag_start);
//			core::id::AtomID const id8( min_frag_pose.residue(j+1).atom_index("O"), j+1 );
//			atom_map[ id8 ] = id7;
		}
		for(int j=size_helical_window_-1; j>=0; j--){
			//Sloppy way of adding all bb atoms, must be a better way...

			TR << "superimpose residues " <<  min_frag_pose.total_residue()-j << " of frag to residue " << pose_frag_end-j << " of pose." << endl;

			core::id::AtomID const id1( min_frag_pose.residue(min_frag_pose.total_residue()-j).atom_index("CA"),
					min_frag_pose.total_residue()-j );
			core::id::AtomID const id2( pose.residue(pose_frag_end-j).atom_index("CA"), pose_frag_end-j);
			atom_map.set(id1, id2);

//			core::id::AtomID const id3( pose.residue(pose_frag_end-j).atom_index("C"), pose_frag_end-j);
//			core::id::AtomID const id4( min_frag_pose.residue(min_frag_pose.total_residue()-j).atom_index("C"), min_frag_pose.total_residue()-j );
//			atom_map[ id4 ] = id3;
//
//			core::id::AtomID const id5( pose.residue(pose_frag_end-j).atom_index("N"), pose_frag_end-j);
//			core::id::AtomID const id6( min_frag_pose.residue(min_frag_pose.total_residue()-j).atom_index("N"), min_frag_pose.total_residue()-j );
//			atom_map[ id6 ] = id5;
//
//			core::id::AtomID const id7( pose.residue(pose_frag_end-j).atom_index("O"), pose_frag_end-j);
//			core::id::AtomID const id8( min_frag_pose.residue(min_frag_pose.total_residue()-j).atom_index("O"), min_frag_pose.total_residue()-j );
//			atom_map[ id8 ] = id7;
		}

		//do superimpose
		core::Real test(core::scoring::superimpose_pose(min_frag_pose, pose, atom_map));

		utility::file::FileName pdb_name (pose.pdb_info()->name());
		std::string output_name = pdb_name.base() + "_loop_" + utility::to_string(cur_jump_edge) + ".pdb";
		min_frag_pose.dump_pdb(output_name);
		min_rms_loop_poses.push_back(min_frag_pose);

//		core::chemical::ResidueTypeSetCAP fa_set =
//				core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
//
//		core::Size anchor_res(pose_frag_start+num_helical_atoms-1);
//		TR << "Pose after apply: " << pose.phi(anchor_res+1) << " " << pose.psi(anchor_res+1) << " " << pose.omega(anchor_res+1) << endl;
//		//change torsions of last elements of helices
//		min_frag_data.get_residue(num_helical_atoms)->apply(pose, anchor_res);
//		for ( int i = num_helical_atoms+1; i <= min_frag_data.size()-num_helical_atoms; i++)
//		{
//			TR << "sequence[" << i << "]: " << sequence[i-1] << endl;
//			core::chemical::ResidueTypeCOPs residue_types = fa_set->aa_map(core::chemical::aa_from_oneletter_code(sequence[i-1]));
//			TR << "restypes size: " << residue_types.size() << endl;
//			core::conformation::Residue new_residue((*residue_types[1]), true/*dummy*/);
//			TR << "appending: " << residue_types[1]->name() << " to residue " << anchor_res << endl;
//
//			//remover upper terminus notation from cap of helix so we can another loop residue
//			core::conformation::remove_upper_terminus_type_from_conformation_residue(pose.conformation(), anchor_res);
//
//			//Append idealized loop residue
//			pose.append_polymer_residue_after_seqpos(new_residue, anchor_res, true);
//
//			anchor_res++;
//		}

		//Generate sub fragment that is a single helical residue, followed by all loop residues
//		core::fragment::FragDataOP loop_frag(min_frag_data.generate_sub_fragment(num_helical_atoms, min_frag_data.size()-num_helical_atoms));

		//Apply new sub fragment to the last residue before the jump and the newly added loop residues
//		loop_frag->apply(pose, pose_frag_start+num_helical_atoms, pose_frag_start+num_helical_atoms+loop_frag->size());

//		TR << "Pose after apply: " << endl;
//		for(core::Size foo=pose_frag_start+num_helical_atoms; foo<=pose_frag_start+num_helical_atoms+loop_frag->size()-1; foo++){
//			TR << "Residue " << foo << " " << pose.phi(foo) << " " << pose.psi(foo) << " " << pose.omega(foo) << endl;
//		}
	}
	TR << total_low_rms_counter << " " << pose.pdb_info()->name() << " " << " bridge fragments < 0.8 for ALL jumps"
							<< endl;

	core::pose::Pose closed_bundle = constructClosedBundle(pose, min_rms_loop_poses);

	utility::file::FileName pdb_name (pose.pdb_info()->name());
	std::string output_name = pdb_name.base() + "_closed.pdb";
	closed_bundle.dump_pdb(output_name);

	pose=closed_bundle;
}

core::pose::Pose BridgeFragmentMover::constructClosedBundle(core::pose::Pose & pose, utility::vector1<core::pose::Pose> & loop_poses){
	utility::vector1< core::pose::PoseOP > chain_poses = pose.split_by_chain();

	core::pose::Pose new_pose;
	for(core::Size i=1; i<=chain_poses.size(); ++i){
		for(core::Size j=1; j<=chain_poses[i]->total_residue(); ++j){
			if((i==1 && j<= chain_poses[i]->total_residue()-size_helical_window_) ||
				(i==chain_poses.size() && j>size_helical_window_) ||
				(j<= chain_poses[i]->total_residue()-size_helical_window_ && j>size_helical_window_)){
				core::conformation::remove_upper_terminus_type_from_conformation_residue(chain_poses[i]->conformation(), j);
				core::conformation::remove_lower_terminus_type_from_conformation_residue(chain_poses[i]->conformation(), j);
				new_pose.append_residue_by_bond(chain_poses[i]->residue(j));
			}
		}
		if(i!=chain_poses.size()){

			for(core::Size j=1; j<=loop_poses[i].total_residue(); ++j){
				core::conformation::remove_upper_terminus_type_from_conformation_residue(loop_poses[i].conformation(), j);
				core::conformation::remove_lower_terminus_type_from_conformation_residue(loop_poses[i].conformation(), j);
				new_pose.append_residue_by_bond(loop_poses[i].residue(j));
//				std::cout << "Added loop residue - new pose size" << new_pose.total_residue() << std::endl;
			}
		}
	}
	return new_pose;
}

} //sewing namespace
} //devel namespace
