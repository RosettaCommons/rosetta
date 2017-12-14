// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RepeatPropagationMover.cc
/// @brief Can repeat both a symmetric, non symmetric & broken pose
/// @author TJ Brunette

// Unit headers
#include <protocols/simple_moves/RepeatPropagationMover.hh>
#include <protocols/simple_moves/RepeatPropagationMoverCreator.hh>


// Project Headers

#include <basic/Tracer.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <protocols/idealize/idealize.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/relax/CentroidRelax.hh>
#include <protocols/relax/cst_util.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/toolbox/superimpose.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <numeric/conversions.hh>
#include <numeric/alignment/QCP_Kernel.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/xyz.functions.hh>

#include <string>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

using namespace core;
using namespace std;
using utility::vector1;

namespace protocols {
namespace simple_moves {
static basic::Tracer TR( "protocols.simple_moves.RepeatPropagationMover" );
// XRW TEMP std::string RepeatPropagationMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return RepeatPropagationMover::mover_name();
// XRW TEMP }


// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RepeatPropagationMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RepeatPropagationMover );
// XRW TEMP }


// XRW TEMP std::string RepeatPropagationMover::mover_name(){
// XRW TEMP  return "RepeatPropagation";
// XRW TEMP }

// apl note -- this appeared in master around 9/28?std::string RepeatPropagationMoverCreator::mover_name(){
// apl note -- this appeared in master around 9/28? return "RepeatPropagationMover";
// apl note -- this appeared in master around 9/28?}


RepeatPropagationMover::RepeatPropagationMover() :
	moves::Mover( mover_name() )
{}

void RepeatPropagationMover::apply(core::pose::Pose & pose) {
	if ( pose.size()== last_res_ ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "can not handle situation where last_res is being repeated (yet...)");
	}
	if ( orig_pdb_length_!=9999 ) {
		detect_last_repeat_residue(pose);
	}
	if ( repeat_without_replacing_pose_ ) {
		copy_phi_psi_omega(pose,pose);
	} else {
		core::pose::PoseOP repeat_poseOP( new core::pose::Pose() );

		initialize_repeat_pose(pose,*repeat_poseOP);
		copy_phi_psi_omega(pose,*repeat_poseOP);
		if ( maintain_cap_ ) {
			add_caps(pose,*repeat_poseOP);
		}
		if ( maintain_ligand_ ) {
			repeat_ligand(pose,*repeat_poseOP);
			fix_ligand_residues(pose,*repeat_poseOP);
			repeat_ligand_constraints(pose,*repeat_poseOP);
		}
		pose = *repeat_poseOP;
		pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) );
		if ( deal_with_ideal_loop_closure_scar_ ) {
			trim_back_repeat_to_repair_scar(pose);
		}
	}

	//rd2
	//  std::cout <<"after propogation_______________________________________________________________________________" << std::endl;
	//  for(Size ii=1; ii<= 75; ii++){
	//   Real distN_CA = pose.residue(ii).xyz("N").distance(pose.residue(ii).xyz("CA"));
	//   Real distCA_C = pose.residue(ii).xyz("CA").distance(pose.residue(ii).xyz("C"));
	//   Real distC_N = pose.residue(ii).xyz("C").distance(pose.residue(ii+1).xyz("N"));
	//   conformation::Residue const & rsd1( pose.residue( ii ) );
	//   conformation::Residue const & rsd2( pose.residue( ii+1 ) );
	//   AtomID const n1( rsd1.atom_index("N" ), ii );
	//   AtomID const ca1( rsd1.atom_index("CA"), ii );
	//   AtomID const c1( rsd1.atom_index("C" ), ii );
	//   AtomID const n2( rsd2.atom_index("N" ), ii+1 );
	//   AtomID const ca2( rsd2.atom_index("CA"), ii+1 );
	//   AtomID const c2( rsd2.atom_index("C" ), ii+1 );
	//   Real angle1 = pose.atom_tree().bond_angle(n1,ca1,c1);
	//   Real angle2 = pose.atom_tree().bond_angle(ca1,c1,n2);
	//   Real angle3 = pose.atom_tree().bond_angle(c1,n2,ca2);
	//   std::cout << ii << ",a:" << pose.phi(ii) <<"," <<  pose.psi(ii) << "," << pose.omega(ii) << "," <<distN_CA << ","<< distCA_C << "," << distC_N <<"||" <<angle1 <<"," << angle2 <<"," << angle3 << std::endl;
	//  }
}

void RepeatPropagationMover::initialize_repeat_pose( Pose & pose, Pose & repeat_pose){
	duplicate_residues_by_type(pose,repeat_pose);
}

void RepeatPropagationMover::add_caps(Pose & pose, Pose & repeat_pose){
	//create c-term subpose first
	Size repeat_size = last_res_-first_res_+1;
	if ( cTerm_cap_size_>0 ) {
		core::pose::Pose cCap_pose;
		vector1<Size> cTerm_pose_positions;
		for ( Size res=pose.size()-cTerm_cap_size_-repeat_size; res<=pose.size(); ++res ) {
			cTerm_pose_positions.push_back(res);
		}
		core::kinematics::FoldTree mod_ft(cTerm_pose_positions.size());
		core::pose::create_subpose(pose,cTerm_pose_positions,mod_ft,cCap_pose);
		Size start_overlap_parent = 0;
		Size end_overlap_parent = 0;
		Size start_overlap_pose = 0;
		Size end_overlap_pose = 0;
		determine_overlap(repeat_pose,cCap_pose,repeat_size,5, "c_term", start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
		generate_overlap(repeat_pose,cCap_pose, "c_term",start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
	}
	if ( nTerm_cap_size_>0 ) {
		core::pose::Pose nCap_pose;
		vector1<Size> nTerm_pose_positions;
		for ( Size res=1; res<=nTerm_cap_size_+repeat_size; ++res ) {
			nTerm_pose_positions.push_back(res);
		}
		core::kinematics::FoldTree mod_ft(nTerm_pose_positions.size());
		core::pose::create_subpose(pose,nTerm_pose_positions,mod_ft,nCap_pose);
		Size start_overlap_parent = 0;
		Size end_overlap_parent = 0;
		Size start_overlap_pose = 0;
		Size end_overlap_pose = 0;
		determine_overlap(repeat_pose,nCap_pose,repeat_size,5, "n_term",start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
		generate_overlap(repeat_pose,nCap_pose, "n_term",start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
	}
}

vector<Real> RepeatPropagationMover::get_center_of_mass(const Real* coordinates, int number_of_atoms){
	vector<Real> center;
	center.push_back(0);
	center.push_back(0);
	center.push_back(0);

	for ( int n = 0; n < number_of_atoms * 3; n += 3 ) {
		center[0] += coordinates[n + 0];
		center[1] += coordinates[n + 1];
		center[2] += coordinates[n + 2];
	}

	center[0] /= number_of_atoms;
	center[1] /= number_of_atoms;
	center[2] /= number_of_atoms;

	return(center);
}

void RepeatPropagationMover::repeat_ligand(Pose & pose, Pose & repeat_pose){
	using namespace chemical;
	using Matrix = numeric::xyzMatrix<Real>;
	Size ligand_residue = 0;
	for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( ! pose.residue( i ).is_protein() ) {
			ligand_residue = i;
			TR << "Residue " << i << " (type=" << pose.residue(i).name3() << ") is probably chosen ligand" << std::endl;
		}
	}
	core::conformation::ResidueOP new_ligand = pose.residue(ligand_residue).clone();
	repeat_pose.append_residue_by_jump(*new_ligand,1,"","",true);

	if ( ligand_residue == 0 ) {
		utility_exit_with_message("ligand not found in pdb");
	}
	Size repeat_length = last_res_ - first_res_+1;
	std::vector< numeric::xyzVector<numeric::Real> > ligand1_coordinates;
	std::vector< numeric::xyzVector<numeric::Real> > ligand1_com_coordinates;
	std::vector< numeric::xyzVector<numeric::Real> > repeat1_coordinates;
	std::vector< numeric::xyzVector<numeric::Real> > repeat1_com_coordinates;
	numeric::alignment::QCP_Kernel<core::Real> qcp;
	for ( Size ii = first_res_;  ii <=last_res_; ++ii ) {
		repeat1_coordinates.push_back(repeat_pose.residue(ii).xyz("CA"));
		repeat1_com_coordinates.push_back(repeat_pose.residue(ii).xyz("CA"));
	}
	Size n_atoms = new_ligand->natoms();
	for ( Size atom_id = 1; atom_id <= n_atoms; ++atom_id ) {
		ligand1_coordinates.push_back(new_ligand->xyz(atom_id));
		ligand1_com_coordinates.push_back(new_ligand->xyz(atom_id));
	}
	numeric::alignment::QCP_Kernel<core::Real>::remove_center_of_mass( &ligand1_com_coordinates.front().x() , ligand1_com_coordinates.size());
	numeric::alignment::QCP_Kernel<core::Real>::remove_center_of_mass( &repeat1_com_coordinates.front().x() , repeat1_com_coordinates.size());
	for ( Size repeat_id = 2; repeat_id <= numb_repeats_-1; ++repeat_id ) { //switched to 2
		//generate new ligand------------
		core::conformation::ResidueOP new_ligand = pose.residue(ligand_residue).clone();
		//get xyz for ligand-----------
		vector1<Real> rot_vector;
		for ( Size ii = 0;  ii < 9; ++ii ) {
			rot_vector.push_back(0);
		}
		std::vector< numeric::xyzVector<numeric::Real> > repeat2_coordinates;
		std::vector< numeric::xyzVector<numeric::Real> > repeat2_com_coordinates;
		for ( Size ii = first_res_+repeat_length*(repeat_id-1);  ii <=first_res_+repeat_length*repeat_id-1; ++ii ) {
			repeat2_coordinates.push_back(repeat_pose.residue(ii).xyz("CA"));
			repeat2_com_coordinates.push_back(repeat_pose.residue(ii).xyz("CA"));
		}
		numeric::alignment::QCP_Kernel<core::Real>::remove_center_of_mass( &repeat2_com_coordinates.front().x() , repeat2_com_coordinates.size());
		qcp.calc_centered_coordinate_rmsd( &repeat2_com_coordinates.front().x(), &repeat1_com_coordinates.front().x(), repeat1_com_coordinates.size(), &rot_vector[1]);
		//TR << "tmp_rmsd" << tmp_rmsd << std::endl;
		Matrix rot_matrix = numeric::xyzMatrix<Real>::rows(rot_vector[1],rot_vector[2],rot_vector[3],rot_vector[4],rot_vector[5],rot_vector[6],rot_vector[7],rot_vector[8],rot_vector[9]);
		Size n_atoms = new_ligand->natoms();
		std::vector<Real> repeat1_com_xyz = get_center_of_mass(&repeat1_coordinates.front().x(),repeat1_coordinates.size());
		std::vector<Real> repeat2_com_xyz = get_center_of_mass(&repeat2_coordinates.front().x(),repeat2_coordinates.size());
		std::vector<Real> ligand1_com_xyz = get_center_of_mass(&ligand1_coordinates.front().x(),ligand1_coordinates.size());
		numeric::xyzVector<numeric::Real> ligand2_com_xyz;
		numeric::xyzVector<numeric::Real> ligand2_tmp_vector;
		ligand2_tmp_vector.x(ligand1_com_xyz[0]-repeat1_com_xyz[0]);
		ligand2_tmp_vector.y(ligand1_com_xyz[1]-repeat1_com_xyz[1]);
		ligand2_tmp_vector.z(ligand1_com_xyz[2]-repeat1_com_xyz[2]);
		//get the distance then apply the rotation.
		numeric::xyzVector<numeric::Real> ligand_tmp_vector_rot = rot_matrix*ligand2_tmp_vector;
		ligand2_com_xyz.x(repeat2_com_xyz[0] +ligand_tmp_vector_rot[0]);
		ligand2_com_xyz.y(repeat2_com_xyz[1] +ligand_tmp_vector_rot[1]);
		ligand2_com_xyz.z(repeat2_com_xyz[2] +ligand_tmp_vector_rot[2]);
		//note numeric xyz runs from 0 and ligand atom count starts at 1
		for ( Size atom_id = 0; atom_id < n_atoms; ++atom_id ) {
			numeric::xyzVector<numeric::Real> ligand_rot = rot_matrix*ligand1_com_coordinates[atom_id];
			numeric::xyzVector<numeric::Real> tmp_xyz;
			tmp_xyz.x(ligand_rot.x()+ligand2_com_xyz[0]);
			tmp_xyz.y(ligand_rot.y()+ligand2_com_xyz[1]);
			tmp_xyz.z(ligand_rot.z()+ligand2_com_xyz[2]);
			new_ligand->atom(atom_id+1).xyz(tmp_xyz);
		}
		//append ligand
		repeat_pose.append_residue_by_jump(*new_ligand,1,"","",true);
	}
}


vector1<Size> RepeatPropagationMover::initial_constrained_residues(const Pose & pose){
	using namespace core::scoring::constraints;
	vector1<Size> constrained_residues;
	ConstraintCOPs constraints = pose.constraint_set()->get_all_constraints();
	for ( ConstraintCOP const c : constraints ) {
		if ( c->type() == "MultiConstraint" ) {
			MultiConstraintCOP multi_cst( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::MultiConstraint const > ( c ) );
			ConstraintCOPs constraints_l2=multi_cst->member_constraints();
			for ( ConstraintCOP const c_l2 : constraints_l2 ) {
				if ( c_l2->type() == "AtomPair" ) {
					AtomPairConstraintCOP cst( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::AtomPairConstraint const > ( c_l2 ) );
					Size res1  =  cst->atom1().rsd();
					Size res2  =  cst->atom2().rsd();
					if ( pose.residue( res1 ).is_protein() ) {
						constrained_residues.push_back(res1);
					}
					if ( pose.residue( res2 ).is_protein() ) {
						constrained_residues.push_back(res2);
					}
				}
			}
		} else {
			std::cout << "*************currently only implemented for MultiConstraint" << std::endl;
		}
	}
	return(constrained_residues);
}


void RepeatPropagationMover::fix_ligand_residues(Pose & pose, Pose & repeat_pose){
	// std::cout << "showing definition" << std::endl;
	// pose.constraint_set()->show_definition( std::cout, pose );
	Size repeat_length = last_res_ - first_res_+1;
	vector1<Size> constrained_residues = initial_constrained_residues(pose);
	for ( Size ii=1; ii<=constrained_residues.size(); ++ii ) {
		if ( pose.residue(constrained_residues[ii]).is_protein() ) {
			Size residue_in_repeat = constrained_residues[ii]%repeat_length;
			for ( Size jj=1; jj<=numb_repeats_; ++jj ) {
				Size residue_to_fix = residue_in_repeat+repeat_length*(jj-1);
				repeat_pose.replace_residue(residue_to_fix, pose.residue(constrained_residues[ii]), true );
			}
		}
	}
}

core::id::SequenceMapping RepeatPropagationMover::setup_repeat_seqmap(Size repeat_number,Size ligand_in_pose, Size ligand_in_repeat){
	Size repeat_length = last_res_ - first_res_+1;
	core::id::SequenceMapping seqmap;
	for ( Size ii=1; ii<=repeat_length*2; ++ii ) {
		seqmap.insert_aligned_residue_safe(ii,ii+((int)repeat_number-1)*repeat_length);
	}
	seqmap.insert_aligned_residue_safe(ligand_in_pose,ligand_in_repeat);
	return(seqmap);
}

void RepeatPropagationMover::repeat_ligand_constraints(Pose & pose, Pose & repeat_pose){
	using namespace core::scoring::constraints;
	using namespace toolbox::match_enzdes_util;
	// std::cout << "showing definition" << std::endl;
	// pose.constraint_set()->show_definition( std::cout, pose );
	repeat_pose.remove_constraints();
	Size repeat_length = last_res_ - first_res_+1;
	ConstraintCOPs constraints = pose.constraint_set()->get_all_constraints();
	for ( Size ii=1; ii<=numb_repeats_-1; ++ii ) {
		for ( ConstraintCOP const c : constraints ) {
			if ( c->type() == "MultiConstraint" ) {
				MultiConstraintCOP multi_cst_old( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::MultiConstraint const > ( c ) );
				ConstraintCOPs constraints_l2=multi_cst_old->member_constraints();
				for ( ConstraintCOP const c_l2 : constraints_l2 ) {
					if ( c_l2->type() == "AtomPair" ) { //it seemed like the wrong residue numbering was stored in the top level multi-constraint
						AtomPairConstraintCOP old_cst( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::AtomPairConstraint const > ( c_l2 ) );
						core::id::SequenceMapping seq_map;
						Size res1  =  old_cst->atom1().rsd();
						Size res2  =  old_cst->atom2().rsd();
						if ( pose.residue( res1 ).is_protein() ) {
							res1 = res1+repeat_length*(ii-1);
						} else {
							res1 = res1+ii-1;//ligand location
						}
						if ( pose.residue( res2 ).is_protein() ) {
							res2 = res2+repeat_length*(ii-1);
						} else {
							res2 = res2+ii-1;//ligand location
						}
						seq_map.insert_aligned_residue_safe(old_cst->atom1().rsd(),res1);
						seq_map.insert_aligned_residue_safe(old_cst->atom2().rsd(),res2);
						MultiConstraintOP multi_cst_new(utility::pointer::dynamic_pointer_cast< core::scoring::constraints::MultiConstraint> (multi_cst_old->remap_resid(seq_map) ));
						repeat_pose.add_constraint(multi_cst_new);
					}
				}
			}
		}
	}
	//repeat the ligand EnzdesCacheableObserverOP but ignoring cloud because I'm unsure how that would propogate
	//setup observer
	EnzdesCacheableObserverOP enz_obs( get_enzdes_observer( repeat_pose ) );
	for ( core::Size ii = 1; ii <= repeat_pose.total_residue(); ii++ ) {
		if ( repeat_pose.residue( ii ).is_ligand() ) {
			Size ligand_residue_pos = ii;
			core::conformation::ResidueCOP ligand_residue = repeat_pose.residue( ii ).get_self_ptr();
			utility::vector1< core::conformation::ResidueCOP > tmp_residue_cloud;
			tmp_residue_cloud.push_back(ligand_residue);
			enz_obs->set_rigid_body_confs_for_lig(ligand_residue_pos,tmp_residue_cloud);
		}
	}
	//initialize repeat_cst_cache and io--------------
	EnzdesCstCacheOP repeat_cst_cache = get_enzdes_observer( repeat_pose )->cst_cache();
	EnzdesCstCacheOP pose_cst_cache = get_enzdes_observer( pose )->cst_cache();
	Size n_csts_repeat_pose = pose_cst_cache->ncsts()*(numb_repeats_-1);
	Size n_csts_pose = pose_cst_cache->ncsts();
	EnzConstraintIOCOP pose_cstio( pose_cst_cache->enzcst_io()); //should be able to use the pose_io but modify the cst_cache
	if ( !repeat_cst_cache ) {
		toolbox::match_enzdes_util::get_enzdes_observer( repeat_pose )->set_cst_cache( toolbox::match_enzdes_util::EnzdesCstCacheOP( new EnzdesCstCache( pose_cstio, n_csts_repeat_pose+1 ) ) ); //+1 is for the positions to hold fix
	}
	repeat_cst_cache = get_enzdes_observer( repeat_pose )->cst_cache();
	//fill repeat_cst_cache
	//setup cst_cache
	Size pose_ligand_position = 999;
	for ( Size ii=1; ii<=pose.total_residue() && pose_ligand_position==999; ++ii ) {
		if ( repeat_pose.residue( ii ).is_ligand() ) {
			pose_ligand_position = ii;
		}
	}
	if ( !pose_cst_cache->contains_position(pose_ligand_position) ) {
		utility_exit_with_message("ligand not constrainted the code below will need to be adjusted to accomodate");
	}
	for ( Size ii=1; ii<=numb_repeats_-1; ++ii ) {
		Size repeat_ligand_position = repeat_length*numb_repeats_+ii; //protein length + number of repeats.
		core::id::SequenceMapping seqmap = setup_repeat_seqmap(ii,pose_ligand_position,repeat_ligand_position);
		std::cout << "seqmap about to be applied" << std::endl;
		seqmap.show(std::cout);
		for ( Size jj=1; jj<=n_csts_pose; ++jj ) {
			Size cst_location_in_cst_cache = ((int)ii-1)*n_csts_pose+jj;
			EnzdesCstParamCacheOP tmp_param_cache =pose_cst_cache->param_cache(jj);
			EnzdesCstParamCacheOP copy_of_param_cache = EnzdesCstParamCacheOP( new EnzdesCstParamCache(*tmp_param_cache));
			copy_of_param_cache->remap_resid(seqmap);
			repeat_cst_cache->set_param_cache(cst_location_in_cst_cache,copy_of_param_cache);
		}
	}
	//----------ligand in first and last repeat
	std::map<Size,std::string> last_repeat_constrained_res; //I use a map so if there's multiple constraints on a residue I only add 1 pose_constraint
	std::map<Size,std::string> first_repeat_constrained_res;
	vector1<Size> pose_constrained_positions = pose_cst_cache->ordered_constrained_positions(pose);
	for ( Size jj=1; jj<=pose_constrained_positions.size(); ++jj ) {
		Size residue_numb = pose_constrained_positions[jj];
		if ( !pose.residue( jj ).is_ligand() ) {
			if ( residue_numb<=repeat_length ) {
				Size location_in_repeat = residue_numb + (numb_repeats_-1)*repeat_length; //place in first repeat + 3 repeats
				last_repeat_constrained_res.insert(std::pair<Size,std::string>(location_in_repeat,pose.residue(residue_numb).name3() ) );
			}
			if ( residue_numb>repeat_length && residue_numb <=2*repeat_length ) {
				Size location_in_repeat = residue_numb-repeat_length; //place in first repeat + 3 repeats
				first_repeat_constrained_res.insert(std::pair<Size,std::string>(location_in_repeat,pose.residue(residue_numb).name3() ) );
			}
		}
	}
	//debug print out first and last residues to constrain
	EnzdesCstParamCacheOP extra_res_param_cache = EnzdesCstParamCacheOP( new EnzdesCstParamCache());
	for ( auto & first_repeat_constrained_re : first_repeat_constrained_res ) {
		extra_res_param_cache->template_res_cache( 1 )->add_position_in_pose( first_repeat_constrained_re.first );
		TR <<"Constraining additional residue" <<  first_repeat_constrained_re.first << std::endl;

	}
	for ( auto & last_repeat_constrained_re : last_repeat_constrained_res ) {
		extra_res_param_cache->template_res_cache( 1 )->add_position_in_pose( last_repeat_constrained_re.first );
		TR <<"Constraining additional residue" <<  last_repeat_constrained_re.first << std::endl;
	}
	repeat_cst_cache->set_param_cache(n_csts_repeat_pose+1 , extra_res_param_cache);
	// EnzdesCstParamCacheOP tmp_param_cache = EnzdesCstParamCacheOP( new EnzdesCstParamCache());
	// std::cout << "hereA" << std::endl;
	// tmp_param_cache->set_position_for_missing_res(10);
	// std::cout << "hereB" << std::endl;
	// repeat_cst_cache->set_param_cache(n_csts_repeat_pose+1,tmp_param_cache);
	// std::cout << "hereC" << std::endl;

	//debug code
	// std::cout << "all fixed constraints" << std::endl;
	// for(Size ii=1; ii<=repeat_cst_cache->ncsts(); ++ii){
	//  std::cout << "cst index" << ii << std::endl;
	//  EnzdesCstParamCacheOP tmp_param_cache = repeat_cst_cache->param_cache(ii);
	//  utility::vector1< ConstraintCOP > & cur_active_constraints = tmp_param_cache->active_pose_constraints();
	//  for ( utility::vector1< ConstraintCOP >::iterator cst_it = cur_active_constraints.begin(); cst_it != cur_active_constraints.end(); ++cst_it ) {
	//   (*cst_it)->show(std::cout);
	//  }
	// }
}



void RepeatPropagationMover::duplicate_residues_by_type( Pose & pose, Pose & repeat_pose){
	//Repeat-----------------------
	if ( first_res_ == 1 ) {
		remove_lower_terminus_type_from_pose_residue(pose, 1);
	}
	if ( last_res_ == pose.size() ) {
		remove_upper_terminus_type_from_pose_residue(pose,pose.size());
	}
	Size tmp_numb_repeats = numb_repeats_;
	for ( Size rep=1; rep<=tmp_numb_repeats; ++rep ) {
		for ( Size res=first_res_; res<=last_res_; ++res ) {
			repeat_pose.append_residue_by_bond(pose.residue(res),true/*ideal bonds*/);

		}
	}
	if ( first_res_==1 ) {
		add_lower_terminus_type_to_pose_residue(pose,1);
	}
	if ( last_res_ == pose.size() ) {
		add_upper_terminus_type_to_pose_residue(pose,pose.size());
	}
	//teminal residues
	add_lower_terminus_type_to_pose_residue(repeat_pose,1);
	add_upper_terminus_type_to_pose_residue(repeat_pose,repeat_pose.size());
}

void RepeatPropagationMover::copy_phi_psi_omega(Pose & pose, Pose & repeat_pose){
	//Repeat region
	Real loop_phi = 0;
	Real loop_psi = 0;
	Real loop_omega = 0;
	char loop_secstruct = 'H';
	Size segment_length = last_res_ - first_res_+1;
	for ( Size ii=1; ii<=segment_length; ++ii ) {
		Size pose_res = first_res_-1+ii;
		Size repeat_pose_res = ii;
		if ( pose_res == 1 ) { //if the first and last positions are involved, loop around
			loop_phi = pose.phi(pose_res+segment_length);
			loop_psi = pose.psi(pose_res+segment_length);
			loop_omega = pose.omega(pose_res+segment_length);
			loop_secstruct = pose.secstruct(pose_res+segment_length);
		} else {
			loop_phi = pose.phi(pose_res);
			loop_psi = pose.psi(pose_res);
			loop_omega =pose.omega(pose_res);
			loop_secstruct = pose.secstruct(pose_res);
		}
		for ( int rep=0; rep<(int)numb_repeats_; ++rep ) {
			repeat_pose.set_phi(repeat_pose_res+( segment_length*rep), loop_phi );
			repeat_pose.set_psi(repeat_pose_res+( segment_length*rep), loop_psi );
			repeat_pose.set_omega( repeat_pose_res+(segment_length*rep),loop_omega );
			repeat_pose.set_secstruct( repeat_pose_res+(segment_length*rep),loop_secstruct );
		}
	}
}

void RepeatPropagationMover::determine_overlap(const Pose & pose, Pose & parent_pose,Size overlap_max_length,Size overlap_range, std::string overlap_location_pose,Size & start_overlap_parent, Size & end_overlap_parent, Size & start_overlap_pose, Size & end_overlap_pose){
	using namespace core::scoring;
	if ( overlap_location_pose== "c_term" ) {
		Size initial_start_res_parent=1;
		Size end_res_parent=overlap_max_length;
		Size initial_start_res_pose=pose.size()-overlap_max_length+1;
		Size end_res_pose=pose.size();
		Real best_rmsd = 9999;
		for ( Size ii=0; ii<overlap_range; ++ii ) {
			utility::vector1<Size> parent_posepositions;
			utility::vector1<Size> pose_positions;
			Size start_res_parent = initial_start_res_parent+ii;
			Size start_res_pose = initial_start_res_pose+ii;
			for ( Size jj=start_res_parent; jj<=end_res_parent; ++jj ) {
				parent_posepositions.push_back(jj);
			}
			for ( Size jj=start_res_pose; jj<=end_res_pose; ++jj ) {
				pose_positions.push_back(jj);
			}
			pose::Pose parent_poseslice;
			pose::Pose ref_pose_slice;
			pdbslice(parent_poseslice,parent_pose,parent_posepositions);
			pdbslice(ref_pose_slice,pose,pose_positions);
			Real rmsd = CA_rmsd(ref_pose_slice,parent_poseslice);
			if ( rmsd < best_rmsd ) {
				best_rmsd = rmsd;
				start_overlap_parent=start_res_parent;
				end_overlap_parent=end_res_parent;
				start_overlap_pose=start_res_pose;
				end_overlap_pose=end_res_pose;
			}
		}
	}
	if ( overlap_location_pose== "n_term" ) {
		Size start_res_parent=parent_pose.size()-overlap_max_length+1;
		Size initial_end_res_parent=parent_pose.size();
		Size start_res_pose=1;
		Size initial_end_res_pose=overlap_max_length;
		Real best_rmsd = 9999;
		for ( Size ii=0; ii<overlap_range; ++ii ) {
			utility::vector1<Size> parent_posepositions;
			utility::vector1<Size> pose_positions;
			Size end_res_parent = initial_end_res_parent-ii;
			Size end_res_pose = initial_end_res_pose-ii;
			for ( Size jj=start_res_parent; jj<=end_res_parent; ++jj ) {
				parent_posepositions.push_back(jj);
			}
			for ( Size jj=start_res_pose; jj<=end_res_pose; ++jj ) {
				pose_positions.push_back(jj);
			}
			pose::Pose parent_poseslice;
			pose::Pose ref_pose_slice;
			pdbslice(parent_poseslice,parent_pose,parent_posepositions);
			pdbslice(ref_pose_slice,pose,pose_positions);
			Real rmsd = CA_rmsd(ref_pose_slice,parent_poseslice);
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

void RepeatPropagationMover::generate_overlap(Pose & pose, Pose & parent_pose, std::string overlap_location_pose,Size start_overlap_parent, Size end_overlap_parent, Size start_overlap_pose, Size end_overlap_pose){
	using namespace core::id;
	using namespace core::scoring;
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, AtomID::BOGUS_ATOM_ID() );
	for ( Size ii=0; ii<=end_overlap_pose-start_overlap_pose; ++ii ) {
		core::id::AtomID const id1(pose.residue(start_overlap_pose+ii).atom_index("CA"),start_overlap_pose+ii);
		core::id::AtomID const id2(parent_pose.residue(start_overlap_parent+ii).atom_index("CA"), start_overlap_parent+ii );
		atom_map[id1]=id2;
	}
	superimpose_pose(pose,parent_pose,atom_map);
	//create subpose of pose
	utility::vector1<Size> pose_positions;
	utility::vector1<Size> parent_posepositions;
	if ( overlap_location_pose== "n_term" ) {
		for ( Size ii=end_overlap_pose+1; ii<=pose.size(); ++ii ) {
			pose_positions.push_back(ii);
		}
		for ( Size ii=1; ii<=end_overlap_parent; ++ii ) {
			parent_posepositions.push_back(ii);
		}
	}
	if ( overlap_location_pose== "c_term" ) {
		for ( Size ii=1; ii<start_overlap_pose; ++ii ) {
			pose_positions.push_back(ii);
		}
		for ( Size ii=start_overlap_parent; ii<=parent_pose.size(); ++ii ) {
			parent_posepositions.push_back(ii);
		}
	}
	pose::Pose ref_pose_slice;
	pose::Pose parent_poseslice;

	pdbslice(ref_pose_slice,pose,pose_positions);
	pdbslice(parent_poseslice,parent_pose,parent_posepositions);
	if ( overlap_location_pose== "n_term" ) {
		remove_upper_terminus_type_from_pose_residue(parent_poseslice,parent_poseslice.size());
		remove_lower_terminus_type_from_pose_residue(ref_pose_slice,1);
		append_pose_to_pose(parent_poseslice,ref_pose_slice,false);
		pose = parent_poseslice;
	}
	if ( overlap_location_pose== "c_term" ) {
		remove_upper_terminus_type_from_pose_residue(ref_pose_slice,ref_pose_slice.size());
		remove_lower_terminus_type_from_pose_residue(parent_poseslice,1);
		append_pose_to_pose(ref_pose_slice,parent_poseslice,false);
		pose = ref_pose_slice;
	}
	renumber_pdbinfo_based_on_conf_chains(pose);
}

// XRW TEMP std::string RepeatPropagationMover::get_name() const { return "RepeatPropagationMover"; }

void RepeatPropagationMover::detect_last_repeat_residue(const Pose & pose){
	int length_change = pose.total_residue()-orig_pdb_length_;
	Size orig_repeat_length = orig_pdb_length_/orig_number_repeats_;
	Size new_repeat_length = orig_repeat_length+length_change;
	//std::cout << "length_change" << length_change << " orig_repeat_length"<< orig_repeat_length << " new_repeat_length" << new_repeat_length << std::endl;
	first_res_ = 15; //Randomly picked. to skip 1/2 the first repeat. May need to be changed
	last_res_ = first_res_+new_repeat_length-1;
	//std::cout << "first_res" << first_res_ << " last_res" << last_res_ << std::endl;
}

void RepeatPropagationMover::trim_back_repeat_to_repair_scar(Pose & pose){
	using namespace core::id;
	//erase residues 1,2 from beginning and last repeat.
	Size repeat_length = last_res_ - first_res_+1;
	Size trim_start_Nterm = 1;
	Size trim_stop_Nterm = repeat_length-first_res_+1;
	Size trim_start_Cterm =pose.total_residue()-first_res_+2; //unsure why +2 but it looks most right
	Size trim_stop_Cterm = pose.total_residue();
	std::cout << "trimming" << trim_start_Nterm << " to " << trim_stop_Nterm << " and " << trim_start_Cterm << " to " << trim_stop_Cterm<< std::endl;
	pose.conformation().delete_residue_range_slow(trim_start_Cterm,trim_stop_Cterm);
	pose.conformation().delete_residue_range_slow(trim_start_Nterm,trim_stop_Nterm);
	renumber_pdbinfo_based_on_conf_chains(pose,true,false,false,false);
}


// XRW TEMP std::string RepeatPropagationMover::get_name() const { return "RepeatPropagationMover"; }

void RepeatPropagationMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & //pose
)
{
	if ( !tag->hasOption("first_template_res") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "repeat mover requires the first residue be entered with first_res tag");
	}
	first_res_ = tag->getOption<Size>("first_template_res");
	last_res_ = tag->getOption<Size>("last_template_res",9999);
	std::string detect_repeat_length_change_str (tag->getOption< std::string >( "detect_repeat_length_change","9999,9999") ); //number of repeats,original length
	utility::vector1< std::string > detect_repeat_length_change_split( utility::string_split( detect_repeat_length_change_str , ',' ) );
	orig_number_repeats_ = atoi(detect_repeat_length_change_split[1].c_str());
	orig_pdb_length_ = atoi(detect_repeat_length_change_split[2].c_str());
	if ( !tag->hasOption("last_template_res") && orig_pdb_length_==9999 ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "repeat mover requires the last residue to be entered with last_res tag or auto detected");
	}
	if ( !tag->hasOption("numb_repeats") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "repeat mover requires the number of repeats be entered with numb_repeats tag");
	}
	numb_repeats_ = tag->getOption<Size>("numb_repeats");
	deal_with_ideal_loop_closure_scar_ = false;
	if ( orig_pdb_length_==9999 ) {
		deal_with_ideal_loop_closure_scar_ = false;
	} else {
		deal_with_ideal_loop_closure_scar_ = true;
		numb_repeats_=numb_repeats_+1;
	}
	repeat_without_replacing_pose_ = tag->getOption<bool>("repeat_without_replacing_pose",false);//for speed.
	maintain_cap_ = tag->getOption<bool>("maintain_cap",false);
	maintain_ligand_ = tag->getOption<bool>("maintain_ligand",false);
	nTerm_cap_size_=0;
	cTerm_cap_size_=0;
	if ( maintain_cap_ ) {
		if ( (!tag->hasOption("nTerm_cap_size"))||(!tag->hasOption("cTerm_cap_size")) ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "repeat mover requires cTerm_cap and nTerm_cap defined if trying to maintain cap sequence or structure");
		}
		nTerm_cap_size_ = tag->getOption<Size>("nTerm_cap_size");
		cTerm_cap_size_ = tag->getOption<Size>("cTerm_cap_size");
	}
}

std::string RepeatPropagationMover::get_name() const {
	return mover_name();
}

std::string RepeatPropagationMover::mover_name() {
	return "RepeatPropagationMover";
}

void RepeatPropagationMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "first_template_res", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "last_template_res", xsct_non_negative_integer, "XRW TO DO", "9999" )
		+ XMLSchemaAttribute::attribute_w_default( "detect_repeat_length_change", xsct_size_cs_pair, "XRW TO DO", "9999,9999" )
		+ XMLSchemaAttribute::required_attribute( "numb_repeats", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "repeat_without_replacing_pose", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "maintain_cap", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "maintain_ligand", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "nTerm_cap_size", xsct_non_negative_integer, "XRW TO DO")
		+ XMLSchemaAttribute( "cTerm_cap_size", xsct_non_negative_integer, "XRW TO DO" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string RepeatPropagationMoverCreator::keyname() const {
	return RepeatPropagationMover::mover_name();
}

protocols::moves::MoverOP
RepeatPropagationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RepeatPropagationMover );
}

void RepeatPropagationMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RepeatPropagationMover::provide_xml_schema( xsd );
}



} // moves
} // protocols
