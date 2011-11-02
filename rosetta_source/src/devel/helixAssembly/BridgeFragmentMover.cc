// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/devel/helixAssembly/BridgeFragmentEvaluator.ccBridgeFragmentEvaluator.cc
/// @brief 
/// @author Tim Jacobs

#include <devel/helixAssembly/BridgeFragmentMover.hh>

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

#include <core/pose/util.hh>
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

//Utility
#include <utility/string_util.hh>

static basic::Tracer TR("BridgeFragmentMover");

BridgeFragmentMover::BridgeFragmentMover(utility::vector1<core::fragment::FragSetOP> frag_sets):
frag_sets_(frag_sets)
{}

BridgeFragmentMover::~BridgeFragmentMover(){}

void BridgeFragmentMover::apply(core::pose::Pose & pose){

	using namespace std;

	//Numer of CA atoms we are using from each helix for RMSD calculations to bridge fragments
	core::Size num_helical_residues(2);

	core::Size num_jumps(pose.fold_tree().num_jump());
	TR << "There are " << num_jumps << " jumps" << endl;
	for(core::Size cur_jump_edge=1; cur_jump_edge <= num_jumps; ++cur_jump_edge){

		core::Size low_rms_counter(0);
		core::Size downstream_jump_residue(pose.fold_tree().downstream_jump_residue(cur_jump_edge));
		core::Size pose_frag_start(downstream_jump_residue-num_helical_residues);
		core::Size pose_frag_end(downstream_jump_residue+num_helical_residues-1);

		core::Size total_size(0);
		for(core::Size frag_sets_index=1; frag_sets_index<=frag_sets_.size(); ++frag_sets_index){
			total_size+=frag_sets_[frag_sets_index]->size();
		}
		TR << "Frag set size: " << total_size << endl;

		core::Real min_rms(10.0);
		core::fragment::AnnotatedFragData min_frag_data("dummy", 1);
		Pose min_frag_pose;
		for(core::Size frag_sets_index=1; frag_sets_index<=frag_sets_.size(); ++frag_sets_index){
			core::fragment::FragSetOP cur_frag_set=frag_sets_[frag_sets_index];

			for(core::fragment::FrameIterator it = cur_frag_set->begin(); it!=cur_frag_set->end(); ++it){

				for(core::Size i = 1; i<=it->nr_frags(); ++i){
					core::fragment::AnnotatedFragData cur_frag =
						dynamic_cast< const core::fragment::AnnotatedFragData & > (it->fragment(i));

					ObjexxFCL::FArray2D< core::Real > p1a( 3, num_helical_residues*2 );
					ObjexxFCL::FArray2D< core::Real > p2a( 3, num_helical_residues*2 );

					std::string p1a_test("");
					std::string p2a_test("");
					for(core::Size j=1; j<=num_helical_residues; j++){
						core::fragment::BBTorsionSRFDCOP fragment_residue =
								dynamic_cast< const core::fragment::BBTorsionSRFD* > (cur_frag.get_residue(j)());

						p1a(1,j)=fragment_residue->x();
						p1a(2,j)=fragment_residue->y();
						p1a(3,j)=fragment_residue->z();

						p1a_test += utility::to_string(fragment_residue->sequence()) + " " +
								utility::to_string(fragment_residue->x()) + " " +
								utility::to_string(fragment_residue->y()) + " " +
								utility::to_string(fragment_residue->z()) + "\n";
					}
					for(core::Size j=0; j<num_helical_residues; j++){
						core::fragment::BBTorsionSRFDCOP fragment_residue =
												dynamic_cast< const core::fragment::BBTorsionSRFD* > (cur_frag.get_residue(cur_frag.size()-(num_helical_residues-1-j))());

						p1a(1,num_helical_residues+j+1)=fragment_residue->x();
						p1a(2,num_helical_residues+j+1)=fragment_residue->y();
						p1a(3,num_helical_residues+j+1)=fragment_residue->z();

						p1a_test += utility::to_string(fragment_residue->sequence()) + " " +
								utility::to_string(fragment_residue->x()) + " " +
								utility::to_string(fragment_residue->y()) + " " +
								utility::to_string(fragment_residue->z()) + "\n";
					}
					for(core::Size j=0; j<num_helical_residues*2;j++){
						numeric::xyzVector<core::Real> res_xyz(pose.residue(pose_frag_start+j).atom("CA").xyz());

						p2a(1,j+1)=res_xyz.x();
						p2a(2,j+1)=res_xyz.y();
						p2a(3,j+1)=res_xyz.z();

						p2a_test += utility::to_string(pose.residue(pose_frag_start+j).name1()) + " " +
								utility::to_string(res_xyz.x()) + " " +
								utility::to_string(res_xyz.y()) + " " +
								utility::to_string(res_xyz.z()) + "\n";
					}
					core::Real rms(numeric::model_quality::rms_wrapper( num_helical_residues*2, p1a, p2a ));

					if(rms<1){
						low_rms_counter++;
					}

					//TR << "Fragment RMSD: " << rms << endl;
					if(rms<min_rms){
						it->fragment_as_pose(i,min_frag_pose,core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ));
						min_rms=rms;
						min_frag_data=cur_frag;
					}

				}
			}
		}

		TR << "Min RMS: " << min_rms << " from pdb " << min_frag_data.pdbid() << endl;
		TR << "Number of bridge fragments < 1: " << low_rms_counter << endl;
		std::string sequence(min_frag_data.sequence());
		TR << "Sequence: " << sequence << endl;

		core::id::AtomID_Map< core::id::AtomID > atom_map;
		// maps every atomid to bogus atom. Atoms in the query_structure that are left mapped to the bogus atoms will not be used in the superimposition
		atom_map.clear();
		core::pose::initialize_atomid_map( atom_map, min_frag_pose, core::id::BOGUS_ATOM_ID );

		//superimpose the found two-helix pose onto the query structure
		for(Size j=0; j< num_helical_residues; j++){
			//Sloppy way of adding all bb atoms, must be a better way...

//			TR << "superimpose residues " <<  j+1 << " of frag to residue " << j+pose_frag_start << " of pose." << endl;

			core::id::AtomID const id1( min_frag_pose.residue(j+1).atom_index("CA"), j+1 );
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
		for(int j=num_helical_residues-1; j>=0; j--){
			//Sloppy way of adding all bb atoms, must be a better way...

//			TR << "superimpose residues " <<  min_frag_pose.total_residue()-j << " of frag to residue " << pose_frag_end-j << " of pose." << endl;

			core::id::AtomID const id1( min_frag_pose.residue(min_frag_pose.total_residue()-j).atom_index("CA"), min_frag_pose.total_residue()-j );
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
//			core::chemical::ResidueTypeCAPs residue_types = fa_set->aa_map(core::chemical::aa_from_oneletter_code(sequence[i-1]));
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
}
