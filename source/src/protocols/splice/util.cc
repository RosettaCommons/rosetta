// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Utility functions useful in Splice mover.
/// @author Gideon Lapidoth (glapidoth@gmail.com)

// Unit Headers
#include <core/pose/extra_pose_info_util.hh>
#include <protocols/splice/util.hh>
#include <protocols/splice/Splice.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/CutChainMover.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/protein_interface_design/movers/AddChainBreak.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <boost/algorithm/string/predicate.hpp>//for comparing string case insensitive
#include <core/scoring/ScoreFunction.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/task_operations/ProteinInterfaceDesignOperation.hh>
#include <protocols/rosetta_scripts/util.hh>










static  basic::Tracer TR( "protocols.splice.util" );
namespace protocols {
namespace splice {

//@brief This function will return the number of the pose resiude given a db entry number
core::Size nearest_to_entry_stop_on_pose (core::pose::Pose const & pose,core::pose::Pose const & template_pose,core::Size residue,std::string tail_segment, std::string protein_family,core::Size chain,std::string segment){
	if ( tail_segment=="" ) {
		//TR<<"I am here"<<std::endl;
		return protocols::rosetta_scripts::find_nearest_res(pose, template_pose, residue, 1);
	}
	if ( tail_segment=="c" ) {
		if ( protein_family=="antibodies" ) {
			if ( segment=="H3" ) {
				return pose.conformation().chain_end( chain );
			} else if ( segment=="L3" ) {
				protocols::simple_moves::CutChainMover ccm;
				core::pose::PoseOP copy_pose_ = core::pose::PoseOP( new core::pose::Pose(pose) );
				return  ccm.chain_cut(*copy_pose_);
			}
		} else { //for now I assume a genral case, might change in the future
			return pose.conformation().chain_end( chain );
		}
	} else if ( tail_segment=="n" ) {
		return protocols::rosetta_scripts::find_nearest_res(pose, template_pose, residue, 0);
	}
	return 0;//case no option found, this indicates an error
}
//@brief This function changes the db used for splicing in so it only inclueds backbone segments that don't clash. e.g. If I am splicing in segmet 2 which is in contatct with segments 1 and 3
// then only use backbone segments that are compatible with what currently is in segment 1 and 3. Gideon
void modify_dbase_with_compatible_backbones(utility::vector1< ResidueBBDofs > torsion_database, utility::vector1< core::Size > & dbase,utility::vector1 <std::string> intersect_pdbs){
	utility::vector1< core::Size > intersect_dbase;
	//TR<<"Size of dbase is: "<<dbase.size()<<std::endl;
	for ( std::vector<core::Size>::iterator it = dbase.begin(); it != dbase.end(); ++it ) {
		//TR<<torsion_database[*it].source_pdb()<<std::endl;
		std::vector<std::string>::iterator it_name;
		it_name = find (intersect_pdbs.begin(), intersect_pdbs.end(), torsion_database[*it].source_pdb());
		if ( it_name != intersect_pdbs.end() ) {
			intersect_dbase.push_back(*it);
		}
	}
	if ( intersect_dbase.size()==0 ) {
		TR << TR.Red << "WARNING! intersect of torsion database of desired delta and compatible PDBs is 0, using only compatible PDBs"<<TR.Reset << std::endl;
		for ( std::vector<ResidueBBDofs>::iterator it = torsion_database.begin(); it != torsion_database.end(); ++it ) {
			std::vector<std::string>::iterator it_name;
			it_name = find (intersect_pdbs.begin(), intersect_pdbs.end(), it->source_pdb());
			if ( it_name != intersect_pdbs.end() ) {
				intersect_dbase.push_back(it-torsion_database.begin()+1);
			}
		}
	}
	dbase = intersect_dbase;
}
void load_pdb_segments_from_pose_comments(core::pose::Pose const & pose, std::map< std::string/*which segment (L1,L2...)*/, std::string/*pdb name*/ > & pdb_segments) {
	// if(use_sequence_profiles_){
	//If we are using sequence profiles then the condition is true and function can run
	using namespace std;
	map<string, string> const comments = core::pose::get_all_comments(pose);
	if ( comments.size() < 1 ) { /// SJF changed from <3 22Jul14
		utility_exit_with_message(
			"Please check comments field in the pdb file (header= ##Begin comments##), could not find any comments");
	}
	TR << "The size of comments is: " << comments.size() << std::endl;
	for ( std::map<string, string>::const_iterator i = comments.begin(); i != comments.end(); ++i ) {
		//TR<<"the size of j is: "<<j<<std::endl;
		std::string const key(i->first);
		//TR<<"the size of j after i->first is: "<<j<<std::endl;
		std::string const val(i->second);
		//TR<<"the size of j after i->second is: "<<j<<std::endl;
		if ( key.substr(0, 7) != "segment" ) { /// the expected format is segment_??, where we're interested in ??
			continue;
		}
		std::string const short_key(key.substr(8, 1000));
		pdb_segments[short_key] = val;
		TR << "recording segment/pdb pair: " << short_key << '/' << val << std::endl;
	}
	// }
}

//@brief used to make sure that I am not adding too many coordiante constraints to the pose. Gideon 23Aug15
core::Size report_coordinate_constraints(core::pose::Pose const & pose){
	using namespace core::scoring::constraints;
	ConstraintCOPs constraints;
	constraints = pose.constraint_set()->get_all_constraints();
	// TR<<" Total number of constraints in the pose is: "<<constraints.size()<<std::endl;
	core::Size cst_num( 0 );
	for ( ConstraintCOP const c: constraints ) {
		if ( c->type() == "CoordinateConstraint" ) {
			//core::Size const seqpos( c->residues()[1]  );
			//TR<<"Residue "<<seqpos<<" has coordinate constraint"<<std::endl;
			cst_num++;
		}
	}
	TR<<"The total number of coordinate constraint is: "<<cst_num<<std::endl;
	return cst_num;
}


//@brief used to see which positions are designable/packable given a TaskOperation
void report_designable_packable_residues(core::pack::task::TaskFactoryOP const tf, core::pose::Pose pose){
	using namespace protocols::rosetta_scripts;
	TR<<"tf:"<<tf<<std::endl;
	utility::vector1<core::Size> designable_residues = residue_packer_states(pose, tf, true, false);
	TR << "Residues allowed to design: " << std::endl;
	for ( utility::vector1<core::Size>::const_iterator i(designable_residues.begin()); i != designable_residues.end(); ++i ) {
		TR << *i << "+";
	}
	TR << std::endl;
	core::pack::task::PackerTaskOP ptask = tf->create_task_and_apply_taskoperations(pose);
	TR<<"Allowed identities at each position:"<<std::endl;
	for ( utility::vector1<core::Size>::const_iterator i(designable_residues.begin()); i != designable_residues.end(); ++i ) {
		TR<<*i<<": ";
		std::list< core::chemical::ResidueTypeCOP > allowed_aas =ptask->residue_task(*i).allowed_residue_types();
		for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator restype = allowed_aas.begin(); restype != allowed_aas.end(); ++restype ) {
			TR<<(*restype )->name1()<<",";
		}
		TR<<std::endl;;
	}
	utility::vector1<core::Size> packable_residues = residue_packer_states(pose, tf, false, true);
	TR << "Residues allowed to repack: " << std::endl;
	for ( utility::vector1<core::Size>::const_iterator i(designable_residues.begin()); i != designable_residues.end(); ++i ) {
		TR << *i << "+";
	}
}

//@brief calculate rms between spliced out loop and source pdb. RMSD is calculated in one of two ways:
// 1) rmsd over the entire loop (dflt)
// 2) rmsd over ss parts and loop parts separately
bool calculate_rmsd(core::pose::Pose & pose, core::pose::Pose const & source_pose,core::Size total_residue_new,
	core::Size nearest_to_from, core::Size from_res, core::Real rms_A, core::Real rms_B  ){
	core::Real rms(0);
	core::Real rms_ss(0);//secondary structure rmsd
	core::Real rms_l(0); // loop rmsd
	core::Size num_of_ss_res(0);
	core::Size num_of_loop_res(0);
	std::string Result_filter;
	std::ostringstream convert_filter;
	TR<<"Source from res: "<<nearest_to_from<<std::endl;
	TR<<"pose from res: "<<from_res<<std::endl;
	TR<<"total_residue_new: "<<total_residue_new<<std::endl;
	TR<<"SS_rms cutoff is: "<<rms_A<<"and the loop rms is: "<<rms_B<<std::endl;
	if ( rms_B<0 ) { //if rms_B is negative (dflt value) then calcaulte rmsd over entire segment
		for ( core::Size i = 0; i <= total_residue_new - 1; ++i ) {
			core::Real const dist(pose.residue(from_res + i).xyz("CA").distance(source_pose.residue(nearest_to_from + i).xyz("CA")));
			rms += dist;
		}
		core::Real const average_rms(rms / total_residue_new);
		TR << "Average distance of spliced segment to original: " << average_rms << std::endl;
		convert_filter << average_rms;
		Result_filter = convert_filter.str();
		core::pose::add_comment(pose, "Average RMSD:", Result_filter);
		if ( average_rms >= rms_A ) {
			TR << "Failing because rmsd = " << average_rms << std::endl;
			return false;
		} else {
			return true;
		}
	} else { //fi
		core::scoring::dssp::Dssp dssp(source_pose);
		for ( core::Size i = 0; i <= total_residue_new - 1; ++i ) {
			//TR<<dssp.get_dssp_secstruct(nearest_to_from + i)<<std::endl;
			if ( dssp.get_dssp_secstruct(nearest_to_from + i)!=' ' ) { //if reside is not a loop
				core::Real const dist(
					pose.residue(from_res + i).xyz("CA").distance(source_pose.residue(nearest_to_from + i).xyz("CA")));
				//TR<<"dist between pose residue "<<from_res + i<<" and source residue "<<nearest_to_from + i<<" is "<<dist<<std::endl;
				rms_ss += dist;
				num_of_ss_res ++;
			} else { //fi// loop residue
				core::Real const dist(
					pose.residue(from_res + i).xyz("CA").distance(source_pose.residue(nearest_to_from + i).xyz("CA")));
				rms_l += dist;
				num_of_loop_res++;
				//TR<<"dist between pose residue "<<from_res + i<<" and source residue "<<nearest_to_from + i<<" is "<<dist<<std::endl;
			}
		}//for
		//TR<<"num of loop res:"<<num_of_loop_res<<std::endl;
		core::Real const average_loop_rms(rms_l / num_of_loop_res);
		core::Real const  average_ss_rms(rms_ss / num_of_ss_res);
		TR << "Average distance of secondary structure spliced segment to original: " << average_ss_rms << std::endl;
		TR << "Average distance of loop spliced segment to original: " << average_loop_rms << std::endl;
		convert_filter << average_ss_rms;
		Result_filter = convert_filter.str();
		core::pose::add_comment(pose, "Secondary structure RMSD:", Result_filter);
		convert_filter << average_loop_rms;
		Result_filter = convert_filter.str();
		core::pose::add_comment(pose, "Loop structure RMSD:", Result_filter);
		if ( (average_ss_rms >= rms_A) || (average_loop_rms >= rms_B) ) {
			TR << "Failing because ss rmsd = " << average_ss_rms<<" and loop rmsd: "<< average_loop_rms << std::endl;
			return false;
		} else {
			return true;
		}
	}//else
}//function
//@ minimize spliced segment
void min_seg(core::pose::Pose & pose, ResidueBBDofs dofs,bool debug, core::Size from_res, std::string tail_segment, core::Size cut_site, core::Size cut_vl_vh_after_llc,ResidueBBDofs tail_dofs,core::scoring::ScoreFunctionOP scorefxn,std::string segment_type){
	pose.remove_constraints();//remove coor/dih constraints before rtmin

	core::kinematics::MoveMapOP mm;
	protocols::splice::SpliceOP splice;
	mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	mm->set_chi(false);
	mm->set_bb(false);
	mm->set_jump(false);
	if ( debug ) {
		pose.dump_pdb("Before_min_seg.pdb");
	}
	TR<<"from_res()@minimization: "<<from_res<<" tail_segment: "<<tail_segment<<std::endl;
	//set movemap for segment
	for ( core::Size i = 1; i < dofs.size() - 1; ++i ) { // No flexibility to the residue after the first C and before the second C
		TR << "Allow flexibility of resi " << from_res + i << std::endl;
		core::Size const pose_resi(from_res + i);
		mm->set_chi(pose_resi, true);
		mm->set_bb(pose_resi, true);
	}
	TR<<"from_res: "<<from_res<<",to_res: "<<from_res+dofs.size()-1<<",dofs size: "<<dofs.size()<<std::endl;

	splice->add_coordinate_constraints(pose, pose, from_res,from_res+dofs.size()-1, 1);
	splice->add_coordinate_constraints(pose, pose, from_res,from_res+dofs.size()-1, 1,"O");
	splice->add_coordinate_constraints(pose, pose, from_res,from_res+dofs.size()-1, 1,"N");
	TR<<pose.conformation().chain_end(1)<<std::endl;
	TR<<dofs.size()<<std::endl;

	//runtime_assert(pose.conformation().chain_end(1) >= from_res+dofs.size()-1);
	if ( !(pose.conformation().chain_end(1) >= from_res+dofs.size()-1) ) {
		utility_exit_with_message( "ERROR: segment end is not within the size of the pose\n");
	}

	//add_dihedral_constraints(pose, pose, from_res(), to_res() + residue_diff);
	//add_dihedral_constraints(pose, pose, from_res(), from_res()+dofs.size()-1);
	if ( !(boost::iequals(tail_segment, "c")) ) { //For tail segments there are no chain breaks
		protocols::protein_interface_design::movers::AddChainBreak acb;
		acb.resnum(utility::to_string(cut_site));
		acb.find_automatically(false);
		acb.change_foldtree(false);
		acb.apply(pose);
	}
	if ( boost::iequals(tail_segment, "n") ) { //This only matters for n-ter tail because C-ter tail residue are included in the previous mm
		splice->tail_fold_tree(pose, cut_vl_vh_after_llc, cut_site,segment_type);
		//add edge to fold tree
		//add_coordinate_constraints(pose, pose, find_nearest_disulfide(pose,from_res())-tail_dofs.size(), find_nearest_disulfide(pose,from_res())-1, find_nearest_disulfide(pose,from_res()));//add coordinate constraint to tail
		//add_coordinate_constraints(pose, pose, find_nearest_disulfide(pose,from_res())-tail_dofs.size(), find_nearest_disulfide(pose,from_res())-1, find_nearest_disulfide(pose,from_res()),"O");//add coordinate constraint to tail
		for ( core::Size i =protocols::rosetta_scripts::find_nearest_disulfide(pose,from_res)-tail_dofs.size(); i<protocols::rosetta_scripts::find_nearest_disulfide(pose,from_res) - 1 ; ++i ) { // No flexibility for the residue before the C and the C.
			TR << "Allow flexibility of resi in tail " << i << std::endl;
			mm->set_chi(i, true);
			mm->set_bb(i, true);
		}
		splice->add_coordinate_constraints(pose, pose, protocols::rosetta_scripts::find_nearest_disulfide(pose,from_res)-tail_dofs.size(),protocols::rosetta_scripts::find_nearest_disulfide(pose,from_res) - 1, 1);
		//tail_fold_tree(pose, cut_vl_vh_after_llc, cut_site);
		protocols::protein_interface_design::movers::AddChainBreak acb;
		acb.resnum(utility::to_string(cut_site));
		acb.find_automatically(false);
		acb.change_foldtree(false);
		acb.apply(pose);
		//mm->show();
		//tail_fold_tree(pose, cut_vl_vh_after_llc, cut_site);
		//TR<<"pose disulfides: "<<find_nearest_disulfide(pose,to_res())<<","<<find_nearest_disulfide(pose,from_res())<<",template disulfides: "<<find_nearest_disulfide(*template_pose_,dofs.stop_loop())<<","<<find_nearest_disulfide(*template_pose_,dofs.start_loop())<<std::endl;
	} else if ( tail_segment!="" ) { //Add this statement for cases I don't use tail segments at all, gideon 11Mar15
		splice->tail_fold_tree(pose, cut_vl_vh_after_llc, 0,segment_type);

	}
	mm->show();
	TR<<"Fold tree before minimization: "<<std::endl;
	pose.fold_tree();
	core::scoring::ScoreFunctionOP scorefxn_with_chainbrk = scorefxn->clone();
	scorefxn_with_chainbrk->set_weight( core::scoring::chainbreak, 1.0 );
	protocols::minimization_packing::MinMover min_mover( mm, scorefxn_with_chainbrk, "dfpmin_armijo_nonmonotone", 0.01, true /*use_nblist*/ );

	TR << "scorefxn before min mover " << std::endl;
	scorefxn_with_chainbrk->show(pose);
	if ( debug ) {
		pose.dump_pdb("before_min_seg.pdb");
	}
	min_mover.apply(pose);
	scorefxn_with_chainbrk->show(pose);
	//remove coordinate constraints post minimization
	pose.remove_constraints();
	TR << "scorefxn after min mover " << std::endl;
	//splice->add_sequence_constraints(pose);
	if ( debug ) {
		pose.dump_pdb("after_min_seg.pdb");
	}
}
utility::vector1< numeric::xyzVector< core::Real > >
coords( core::pose::Pose const & pose, utility::vector1< core::Size > const positions ){
	utility::vector1< numeric::xyzVector< core::Real > > coords;

	coords.clear();
	for ( core::Size const pos: positions ) {
		coords.push_back( pose.residue( pos ).xyz( "N" ) );
		coords.push_back( pose.residue( pos ).xyz( "CA" ) );
		coords.push_back( pose.residue( pos ).xyz( "C" ) );
		coords.push_back( pose.residue( pos ).xyz( "O" ) );
	}
	return coords;
}
utility::vector1<core::Size> find_residues_on_chain1_inside_interface(core::pose::Pose const & pose, core::Size chainNum) {
	utility::vector1<core::Size> const chain1_interface;
	if ( pose.conformation().num_chains()==1 ) {
		return chain1_interface;//if no ligand then there are no interface residues;
	}


	using namespace protocols::task_operations;
	ProteinInterfaceDesignOperationOP pido( new ProteinInterfaceDesignOperation );
	if ( chainNum == 1 ) {
		pido->repack_chain1(true);
		pido->design_chain1(true);
		pido->repack_chain2(false);
		pido->design_chain2(false);
		pido->interface_distance_cutoff(8.0);
		core::pack::task::TaskFactoryOP tf_interface(new core::pack::task::TaskFactory);
		tf_interface->push_back(pido);
		///// FIND COMPLEMENT ////////
		utility::vector1<core::Size> const chain1_interface(
			protocols::rosetta_scripts::residue_packer_states(pose, tf_interface, true, true)); /// find packable but not designable residues; according to pido specifications above these will be on chain1 outside an 8A shell around chain2
		return chain1_interface;
	} else { //if user defined chain 2 as binder
		pido->repack_chain1(false);
		pido->design_chain1(false);
		pido->repack_chain2(true);
		pido->design_chain2(true);
		pido->interface_distance_cutoff(8.0);
		core::pack::task::TaskFactoryOP tf_interface(new core::pack::task::TaskFactory);
		tf_interface->push_back(pido);
		///// FIND COMPLEMENT ////////
		utility::vector1<core::Size> const chain1_interface(
			protocols::rosetta_scripts::residue_packer_states(pose, tf_interface, true, true)); /// find packable but not designable residues; according to pido specifications above these will be on chain1 outside an 8A shell around chain2
		return chain1_interface;
	}

}
std::string parse_pdb_code(std::string pdb_file_name) {
	const size_t last_slash_idx = pdb_file_name.find_last_of("\\/");
	if ( std::string::npos != last_slash_idx ) {
		pdb_file_name.erase(0, last_slash_idx + 1);
	}
	// Remove extension if present.
	const size_t period_idx = pdb_file_name.rfind('.');
	if ( std::string::npos != period_idx ) {
		pdb_file_name.erase(pdb_file_name.find('.'));
	}
	return pdb_file_name;

}
//@brief Sometimes the chain break residue after splice in has overlaping O and N atoms. This function erases these and create new ones.
void fix_chain_break_residue(core::pose::Pose & pose,SpliceManager splicemanager ){
	//pose.dump_pdb(splicemanager.mover_name()+ "_before_test_O_change.pdb");//
	//The purpose of this code block is to fix a problem I saw in Splice in. For some reason (probably because of using ideal residues in the loop)
	// The amide O atom is positioned wrong. This is taken care of during CCD because of the chain break score but when doing Splice in where the
	// BB dihedral angles are just imposed this problem does not get fixed and casue really bad chain break scores. To fix this I am using the
	// "build_atom_ideal" function to re-postion the O atom. Gideon, Mar15
	if ( splicemanager.debug() ) { //if debug be more verbose, Gideon mar15
		TR<<"Position of O atom of cut_site residue before correction: "
			<<"X:"<<pose.residue(splicemanager.cut_site()).xyz(pose.residue(splicemanager.cut_site()).atom_index("O"))[0]
			<<" Y:"<<pose.residue(splicemanager.cut_site()).xyz(pose.residue(splicemanager.cut_site()).atom_index("O"))[1]
			<<" Z:"<<pose.residue(splicemanager.cut_site()).xyz(pose.residue(splicemanager.cut_site()).atom_index("O"))[2]<<std::endl;
		TR<<"Get ideal xyz coordiantes:"<<std::endl;
	}
	core::Vector const O_xyz = pose.residue(splicemanager.cut_site()).build_atom_ideal(pose.residue(splicemanager.cut_site()).atom_index("O"),pose.conformation());
	core::Size O_atom_index = pose.residue(splicemanager.cut_site()).atom_index("O");
	core::conformation::Residue temp_residue=pose.residue(splicemanager.cut_site());
	temp_residue.set_xyz( O_atom_index, O_xyz );
	if ( splicemanager.debug() ) {
		TR<<"Position of the BB O atom of cut_site residue after idealization: "
			<<"X:"<<temp_residue.xyz(pose.residue(splicemanager.cut_site()).atom_index("O"))[0]
			<<" Y:"<<temp_residue.xyz(pose.residue(splicemanager.cut_site()).atom_index("O"))[1]
			<<" Z:"<<temp_residue.xyz(pose.residue(splicemanager.cut_site()).atom_index("O"))[2]<<std::endl;
	}
	core::conformation::Residue const residue(temp_residue);
	core::Size const seqpos = splicemanager.cut_site();
	pose.replace_residue(seqpos,residue,false);

	core::Vector const N_xyz = pose.residue(splicemanager.cut_site()).build_atom_ideal(pose.residue(splicemanager.cut_site()).atom_index("N"),pose.conformation());
	core::Size N_atom_index = pose.residue(splicemanager.cut_site()).atom_index("N");
	core::conformation::Residue temp_residue_N=pose.residue(splicemanager.cut_site());
	temp_residue_N.set_xyz( N_atom_index, N_xyz );
	core::conformation::Residue const residue_N(temp_residue_N);
	pose.replace_residue(seqpos,residue,false);
	//pose.dump_pdb(splicemanager.mover_name()+ "_after_test_O_change.pdb");//


}
//@brief Prints to the log the packable/designable residues given a
void report_des_pack_res(core::pose::Pose const & pose,core::pack::task::TaskFactoryOP tf){
	using namespace protocols::rosetta_scripts;

	utility::vector1<core::Size> packable_residues = residue_packer_states(pose, tf, false, true);

	TR << "Residues Allowed to Repack Before RTmin: " << std::endl;
	for ( utility::vector1<core::Size>::const_iterator i(packable_residues.begin()); i != packable_residues.end(); ++i ) {
		TR << *i << "+";
	}
	TR << std::endl;

	utility::vector1<core::Size> designable_residues = residue_packer_states(pose, tf, true, false);

	TR << "Residues Allowed to Design Before RTmin: " << std::endl;
	for ( utility::vector1<core::Size>::const_iterator i(packable_residues.begin()); i != packable_residues.end(); ++i ) {
		TR << *i << "+";
	}
	TR << std::endl;
}


} //splice
} //protocols
