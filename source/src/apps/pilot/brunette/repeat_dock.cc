// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.AlignmentCluster

/// @file   apps/pilot/brunette/repeat_dock
///
/// @brief  analyzes helix oritentation in repeat proteins

/// @usage:

/// @author TJ Brunette


// Utility Headers
#include <basic/Tracer.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>

// Core Headers
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>

#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>


#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>


#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/pose/motif/reference_frames.hh>
#include <numeric/xyzTransform.hh>


#include <core/types.hh>

#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/jumping/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <devel/init.hh>

#include <utility/string_constants.hh>
#include <utility/vector1.hh>
//basic & utility
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility/io/ozstream.hh>
#include <iostream>
#include <utility/io/ozstream.hh>

#include <ObjexxFCL/format.hh>
#include <sstream>
#include <boost/foreach.hpp>

using namespace ObjexxFCL::format;
using utility::vector1;
using core::Size;
using core::Real;

static THREAD_LOCAL basic::Tracer TR( "repeat_dock" );

struct Dock{
    Size h1_pos_start;
    Size h1_pos_end;
    Size h2_pos_start;
    Size h2_pos_end;
    Size h1_nat_start;
    Size h1_nat_end;
    Size h2_nat_start;
    Size h2_nat_end;
    Size nat_final_poseCoords_start;
    Size nat_final_poseCoords_end;
    Size repeat_final_poseCoords_start;
    Size repeat_final_poseCoords_end;
    Size redesign_final_poseCoords_start;
    Size redesign_final_poseCoords_end;
    Size nat_final_natCoords_start;
    Size nat_final_natCoords_end;
    core::pose::Pose unclosed_pose_full;
    core::pose::Pose unclosed_pose;
    core::pose::Pose closed_pose;
    core::pose::Pose final_pose;
    Real rmsd;
    bool single_helix;
    Real centroid_score;
    Real vdw_score;
    Real chainA_motif_score; //motif score for the 2 helices post attachment point
    Real overall_motif_score; //motif score for the 2 helices post attachment point in context of the rest of the structure
    Dock(Size h1_pos_start_t,Size h1_pos_end_t,Size h2_pos_start_t,Size h2_pos_end_t,Size h1_nat_start_t,Size h1_nat_end_t,Size h2_nat_start_t,Size h2_nat_end_t):h1_pos_start(h1_pos_start_t),h1_pos_end(h1_pos_end_t),h2_pos_start(h2_pos_start_t),h2_pos_end(h2_pos_end_t),h1_nat_start(h1_nat_start_t),h1_nat_end(h1_nat_end_t),h2_nat_start(h2_nat_start_t),h2_nat_end(h2_nat_end_t){single_helix=false;}
    Dock(Size h1_pos_start_t,Size h1_pos_end_t,Size h1_nat_start_t,Size h1_nat_end_t):h1_pos_start(h1_pos_start_t),h1_pos_end(h1_pos_end_t),h1_nat_start(h1_nat_start_t),h1_nat_end(h1_nat_end_t){
    single_helix=true;
    h2_pos_start = 0;
    h2_pos_end = 0;
    h2_nat_start = 0;
    h2_nat_end = 0;
    }

};


void get_helices(core::pose::Pose& pose, protocols::loops::LoopsOP helicesOP){
    using numeric::xyzVector;
	protocols::jumping::assign_ss_dssp( pose );
	char lastSecStruct = pose.secstruct(1);
	Size startHelix = 0;
	Size endHelix = 0;
	if(pose.secstruct(1) == 'H')
		startHelix = 1;
	for ( core::Size ii = 2; ii <= pose.total_residue(); ++ii ) {
		if(pose.secstruct(ii) == 'H' && lastSecStruct != 'H')
			startHelix = ii;
		if(pose.secstruct(ii) != 'H' && lastSecStruct == 'H'){
			endHelix = ii-1;
			if(endHelix-startHelix >= 2)
				helicesOP->add_loop(startHelix,endHelix);
		}
		lastSecStruct = pose.secstruct(ii);
	}
}


Real calculate_helical_tail_variance(core::pose::Pose pose, core::pose::Pose native_pose,Dock testDock){
	using namespace core::scoring;
	using namespace core::id;
    //std::cout << "pose_range" << testDock.h1_pos_start << "," << testDock.h1_pos_end  << "," <<  testDock.h2_pos_start << "," << testDock.h2_pos_end  << std::endl; 
    //std::cout << "native_pose_range" << testDock.h1_nat_start << "," << testDock.h1_nat_end  << "," <<  testDock.h2_nat_start << "," << testDock.h2_nat_end  << std::endl; 
    int h1_pos_incr;
    int h2_pos_incr;
    int h1_nat_incr;
    int h2_nat_incr;
    std::map<Size, Size> residues;
    if(testDock.h1_pos_start<testDock.h1_pos_end)
        h1_pos_incr = 1;
    else
        h1_pos_incr = -1;
    if(testDock.h2_pos_start<testDock.h2_pos_end)
        h2_pos_incr = 1;
    else
        h2_pos_incr = -1;
    if(testDock.h1_nat_start<testDock.h1_nat_end)
        h1_nat_incr = 1;
    else
        h1_nat_incr = -1;
    if(testDock.h2_nat_start<testDock.h2_nat_end)
        h2_nat_incr = 1;
    else
        h2_nat_incr = -1;
	for (int ii= 0; ii<= std::abs((int)testDock.h1_pos_start-(int)testDock.h1_pos_end); ++ii){
        Size pose1_pos = (Size)((int)ii*h1_pos_incr+testDock.h1_pos_start);
        Size pose2_pos = (Size)((int)ii*h1_nat_incr+testDock.h1_nat_start);
        residues.insert(std::pair<Size,Size>(pose1_pos,pose2_pos));
    }
    if(std::abs((int)testDock.h2_pos_start-(int)testDock.h2_pos_end)>0){
        for (int ii= 0; ii<= std::abs((int)testDock.h2_pos_start-(int)testDock.h2_pos_end); ++ii){
            Size pose1_pos = (Size)(ii*h2_pos_incr+testDock.h2_pos_start);
            Size pose2_pos = (Size)(ii*h2_nat_incr+testDock.h2_nat_start);
	        residues.insert(std::pair<Size,Size>(pose1_pos,pose2_pos));
        }
    }
    Real rmsd = core::scoring::CA_rmsd(pose,native_pose,residues);
    return(rmsd);
}

core::pose::Pose minimize_to_close(core::pose::Pose pose){
    using namespace core::optimization;
    core::scoring::ScoreFunctionOP  lowres_scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" ) );
    lowres_scorefxn->set_weight( core::scoring::cart_bonded, 1.0 );
    //lowres_scorefxn->set_weight( core::scoring::coordinate_constraint, 0.1 );
    //lowres_scorefxn->set_weight( core::scoring::dihedral_constraint, 1.0 );
    core::optimization::MinimizerOptions options_minilbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
    lowres_scorefxn->show(std::cout, pose);
    options_minilbfgs.max_iter(100);
    core::optimization::CartesianMinimizer minimizer;
    core::kinematics::MoveMap mm;
	mm.set_bb  ( true );
	mm.set_chi ( true );
	mm.set_jump( true );
    minimizer.run( pose, mm, *lowres_scorefxn, options_minilbfgs );
    return(pose);
}


void generate_combined_model(core::pose::Pose pose, core::pose::Pose native_pose,Dock & testDock){
    using namespace core::id;
    using namespace core::scoring;
    //step1:superimpose pose.
    //
    core::id::AtomID_Map< core::id::AtomID > atom_map;
    core::pose::initialize_atomid_map( atom_map, pose, BOGUS_ATOM_ID );
    int h1_pos_incr;
    int h2_pos_incr;
    int h1_nat_incr;
    int h2_nat_incr;
    if(testDock.h1_pos_start<testDock.h1_pos_end)
        h1_pos_incr = 1;
    else
        h1_pos_incr = -1;
    if(testDock.h2_pos_start<testDock.h2_pos_end)
        h2_pos_incr = 1;
    else
        h2_pos_incr = -1;
    if(testDock.h1_nat_start<testDock.h1_nat_end)
        h1_nat_incr = 1;
    else
        h1_nat_incr = -1;
    if(testDock.h2_nat_start<testDock.h2_nat_end)
        h2_nat_incr = 1;
    else
        h2_nat_incr = -1;
	for (int ii= 0; ii<= std::abs((int)testDock.h1_pos_start-(int)testDock.h1_pos_end); ++ii){
        Size pose1_pos = (Size)((int)ii*h1_pos_incr+testDock.h1_pos_start);
        Size pose2_pos = (Size)((int)ii*h1_nat_incr+testDock.h1_nat_start);
        core::id::AtomID const id1( pose.residue(pose1_pos).atom_index("CA"), pose1_pos );
        core::id::AtomID const id2(native_pose.residue(pose2_pos).atom_index("CA"),pose2_pos);
        atom_map[id1]=id2;
    }
    if(std::abs((int)testDock.h2_pos_start-(int)testDock.h2_pos_end)>0){
        for (int ii= 0; ii<= std::abs((int)testDock.h2_pos_start-(int)testDock.h2_pos_end); ++ii){
            Size pose1_pos = (Size)(ii*h2_pos_incr+testDock.h2_pos_start);
            Size pose2_pos = (Size)(ii*h2_nat_incr+testDock.h2_nat_start);
            core::id::AtomID const id1( pose.residue(pose1_pos).atom_index("CA"), pose1_pos );
            core::id::AtomID const id2(native_pose.residue(pose2_pos).atom_index("CA"),pose2_pos);
            atom_map[id1]=id2;
	    }
    }
    superimpose_pose(pose,native_pose,atom_map);
    //step3: create reference pose without first/last chunck and attach
    //assume you want the longer part of the helix to maintain. Shouldn't be close
    int dist_to_end = pose.total_residue()-numeric::max(testDock.h1_pos_start,testDock.h1_pos_end,testDock.h2_pos_start,testDock.h2_pos_end);
    int dist_to_begin = numeric::min(testDock.h1_pos_start,testDock.h1_pos_end,testDock.h2_pos_start,testDock.h2_pos_end)-1;
    if(testDock.single_helix){
        dist_to_end = pose.total_residue()-numeric::max(testDock.h1_pos_start,testDock.h1_pos_end);
        dist_to_begin = numeric::min(testDock.h1_pos_start,testDock.h1_pos_end)-1;
    }
    vector1<Size> pose_positions;
    Size firstRes = 0;
    Size lastRes = 0;
    if(dist_to_end>dist_to_begin){
        firstRes = numeric::max(testDock.h1_pos_start,testDock.h1_pos_end,testDock.h2_pos_start,testDock.h2_pos_end)+1;
        lastRes = pose.total_residue();
    }
    else{
        firstRes = 1;
        if(testDock.single_helix)
            lastRes = numeric::min(testDock.h1_pos_start,testDock.h1_pos_end)-1;
        else
            lastRes = numeric::min(testDock.h1_pos_start,testDock.h1_pos_end,testDock.h2_pos_start,testDock.h2_pos_end)-1;
    }
     for (Size ii= firstRes; ii<=lastRes; ++ii){ //append all residued in clone to this, or vice versa
		pose_positions.push_back(ii);
	}
    core::pose::Pose mod_pose;
    core::kinematics::FoldTree f_mod(pose_positions.size());
    core::pose::create_subpose(pose,pose_positions,f_mod,mod_pose);
    //chop out chain A. 
    Size chain_id =  get_chain_id_from_chain("A",native_pose);
    core::pose::PoseOP chainA_ = native_pose.split_by_chain(chain_id);
    core::pose::PoseOP chainA_clone_ = chainA_->clone();
    //append mod pose to chain A if dist_to_end>dist_to_begin, alternatively append chain A to mod pose.
    core::pose::Pose chainA_front;
    core::pose::Pose chainA_back;
    Size tmp_start;
    Size tmp_end;
    if(dist_to_end>dist_to_begin){
        vector1<Size> slice_res;
        chainA_back = mod_pose;
        if(testDock.single_helix){
            tmp_start = 1;
            tmp_end = testDock.h1_nat_end;
        }
        else{
            tmp_start = 1;
            tmp_end = testDock.h2_nat_end;
        }
        for(Size ii=tmp_start; ii<=tmp_end; ++ii)
            slice_res.push_back(ii);
        pdbslice(chainA_front,*chainA_,slice_res);
    }
    else{
        vector1<Size> slice_res;
        chainA_front = mod_pose;
        if(testDock.single_helix){
            tmp_start = testDock.h1_nat_start;
            tmp_end = chainA_->total_residue();
        }
        else{
            tmp_start = testDock.h1_nat_start;
            tmp_end = chainA_->total_residue();
        }
        for(Size ii=tmp_start; ii<=tmp_end; ++ii)
            slice_res.push_back(ii);
         pdbslice(chainA_back,*chainA_,slice_res);
    }


    remove_upper_terminus_type_from_pose_residue(chainA_front,chainA_front.total_residue());
    remove_lower_terminus_type_from_pose_residue(chainA_back,1);
    for(Size ii=1; ii<=chainA_back.total_residue(); ++ii){
        chainA_front.append_residue_by_bond(chainA_back.residue(ii),false,0,chainA_front.total_residue());
    }
    core::pose::Pose chainA_full = chainA_front;
    Size finalNativeLength = tmp_end - tmp_start+1;
    testDock.nat_final_natCoords_start = tmp_start;
    testDock.nat_final_natCoords_end = tmp_end;
    if(dist_to_end<dist_to_begin){
        testDock.redesign_final_poseCoords_start = 1;
        testDock.redesign_final_poseCoords_end = 1;
        testDock.repeat_final_poseCoords_start = 1;
        testDock.repeat_final_poseCoords_end = chainA_full.total_residue()-finalNativeLength-1;
        testDock.nat_final_poseCoords_start = chainA_full.total_residue()-finalNativeLength;
        testDock.nat_final_poseCoords_end = chainA_full.total_residue();
    }
    else{
        testDock.redesign_final_poseCoords_start = 1;
        testDock.redesign_final_poseCoords_end = 1;
        testDock.repeat_final_poseCoords_start = finalNativeLength+1;
        testDock.repeat_final_poseCoords_end = chainA_full.total_residue();
        testDock.nat_final_poseCoords_start = 1;
        testDock.nat_final_poseCoords_end = finalNativeLength;
    }
    testDock.unclosed_pose = chainA_front;
    //put rest of structure back together --------------------------------------
    std::string alphabet = "BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"; //A is excluded because of the matching
    BOOST_FOREACH( char c, alphabet) {
        Size chain_id = -1;
        try{
            chain_id = get_chain_id_from_chain(c,native_pose);
            core::pose::PoseOP tmp_chain_pose = native_pose.split_by_chain(chain_id);
            append_pose_to_pose(chainA_full,*tmp_chain_pose);
        }
        catch(utility::excn::EXCN_Base const & e ){
            //not a worry. I expect many chains to not exist.
        }
    }
    core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
    scorefxn->set_weight( core::scoring::vdw, 1.0 );
    Real score = scorefxn->score(chainA_full);
    testDock.vdw_score = score;
    core::scoring::ScoreFunctionOP  lowres_scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "score3" ) );
    score = lowres_scorefxn->score(chainA_full);
    testDock.centroid_score = score;
    testDock.unclosed_pose_full = chainA_full;
}

vector1<Dock> get_all_2helix_docks(core::pose::Pose & pose, core::pose::Pose & native_pose){
    using protocols::loops::Loops;
    Size range = 4;
    using namespace core::pose;
    vector1<Dock> docks;
    Size chain_id =  get_chain_id_from_chain("A",native_pose);
    PoseOP chainA_pose = native_pose.split_by_chain(chain_id);
    protocols::loops::LoopsOP pose_helices( new protocols::loops::Loops );
    protocols::loops::LoopsOP native_helices( new protocols::loops::Loops );
    //get helices
    get_helices(pose,pose_helices);
    vector1<Size> pose_helix_selection;
    vector1<Size> native_helix_selection;
    //option A: front of repeat
    Size helix1_length = (*pose_helices)[1].length();
    Size helix2_length = (*pose_helices)[2].length();
    pose_helix_selection.push_back((*pose_helices)[1].start());
    pose_helix_selection.push_back((*pose_helices)[1].stop());
    pose_helix_selection.push_back((*pose_helices)[2].start());
    pose_helix_selection.push_back((*pose_helices)[2].stop());
    for(Size ii=0; ii<range; ++ii){
        for(Size kk =0; kk<range; ++kk){
            native_helix_selection.push_back(1+ii);
            native_helix_selection.push_back(helix1_length+ii);
            native_helix_selection.push_back(chainA_pose->total_residue()-helix2_length+1-kk);
            native_helix_selection.push_back(chainA_pose->total_residue()-kk);
            struct Dock dock_tmp(pose_helix_selection[1],pose_helix_selection[2],pose_helix_selection[3],pose_helix_selection[4],native_helix_selection[1],native_helix_selection[2],native_helix_selection[3],native_helix_selection[4]);
            docks.push_back(dock_tmp);
            native_helix_selection.clear();

        }
    }
    pose_helix_selection.clear();
    //-------------------------------------------------------------------------
    //option B back of repeat 
    pose_helix_selection.push_back((*pose_helices)[(*pose_helices).size()-1].start());
    pose_helix_selection.push_back((*pose_helices)[(*pose_helices).size()-1].stop());
    pose_helix_selection.push_back((*pose_helices)[(*pose_helices).size()].start());
    pose_helix_selection.push_back((*pose_helices)[(*pose_helices).size()].stop());
    for(Size ii=0; ii<range; ++ii){
        for(Size kk =0; kk<range; ++kk){
            native_helix_selection.push_back(1+ii);
            native_helix_selection.push_back(helix2_length+ii);
            native_helix_selection.push_back(chainA_pose->total_residue()-helix1_length+1-kk);
            native_helix_selection.push_back(chainA_pose->total_residue()-kk);
            struct Dock dock_tmp(pose_helix_selection[1],pose_helix_selection[2],pose_helix_selection[3],pose_helix_selection[4],native_helix_selection[1],native_helix_selection[2],native_helix_selection[3],native_helix_selection[4]);
            docks.push_back(dock_tmp);
            native_helix_selection.clear();
        }
    }
    return(docks);
}


    
vector1<Dock> get_all_1helix_docks(core::pose::Pose & pose, core::pose::Pose & native_pose){
    using protocols::loops::Loops;
    Size range = 3;
    using namespace core::pose;
    vector1<Dock> docks;
    Size chain_id =  get_chain_id_from_chain("A",native_pose);
    PoseOP chainA_pose = native_pose.split_by_chain(chain_id);
    protocols::loops::LoopsOP pose_helices( new protocols::loops::Loops );
    protocols::loops::LoopsOP native_helices( new protocols::loops::Loops );
    //get helices
    get_helices(pose,pose_helices);
    vector1<Size> pose_helix_selection;
    vector1<Size> native_helix_selection;
    //option A: last helix native, first helix protein
    Size helix1_length = (*pose_helices)[1].length();
    pose_helix_selection.push_back((*pose_helices)[1].start());
    pose_helix_selection.push_back((*pose_helices)[1].stop());
    for(Size ii=0; ii<range; ++ii){
        native_helix_selection.push_back(chainA_pose->total_residue()-helix1_length+1-ii);
        native_helix_selection.push_back(chainA_pose->total_residue()-ii);
        struct Dock dock_tmp(pose_helix_selection[1],pose_helix_selection[2],native_helix_selection[1],native_helix_selection[2]);
        docks.push_back(dock_tmp);
        native_helix_selection.clear();
        }
    pose_helix_selection.clear();
    //option B: first helix native, last helix protein
    helix1_length = (*pose_helices)[(*pose_helices).size()].length();
    pose_helix_selection.push_back((*pose_helices)[(*pose_helices).size()].start());
    pose_helix_selection.push_back((*pose_helices)[(*pose_helices).size()].stop());
    for(Size ii=0; ii<range; ++ii){
        native_helix_selection.push_back(1+ii);
        native_helix_selection.push_back(helix1_length+ii);
        struct Dock dock_tmp(pose_helix_selection[1],pose_helix_selection[2],native_helix_selection[1],native_helix_selection[2]);
        docks.push_back(dock_tmp);
        native_helix_selection.clear();
    }
    pose_helix_selection.clear();
    return(docks);
}

core::Real compute_motif_score( core::pose::Pose const & pose, Size start_res, Size end_res, core::scoring::motif::MotifHashManager *mman_){
    using namespace core::scoring;
    using namespace core::scoring::motif;
    double score = 0.0;
    core::scoring::dssp::Dssp dssp( pose );
    dssp.dssp_reduced();
    for(size_t ir = start_res; ir <= end_res; ++ir){
        Xform const ibb_stub = core::pose::motif::get_backbone_reference_frame(pose,ir);
        char ss1 = dssp.get_dssp_secstruct( ir );
        char aa1 = pose.residue(ir).name1();
        for(size_t jr = ir+1; jr <= pose.n_residue(); ++jr){
            Real dist = pose.residue(ir).xyz("CA").distance(pose.residue(jr).xyz("CA"));
            if(dist < 12){
                char ss2 = dssp.get_dssp_secstruct( jr );
                char aa2 = pose.residue(jr).name1();
                //std::cout << ss1 << ss2 << " " << aa1 << aa2 << "dist" << dist <<  std::endl;
                Xform const jbb_stub = core::pose::motif::get_backbone_reference_frame(pose,jr);
                Xform const Xbb = ibb_stub.inverse() * jbb_stub;
                core::scoring::motif::XformScoreCOP xs_bb_fxn1(mman_->get_xform_score_BB_BB(ss1,ss2,aa1,aa2));
                core::scoring::motif::XformScoreCOP xs_bb_fxn2(mman_->get_xform_score_BB_BB(ss2,ss1,aa2,aa1));
				if(xs_bb_fxn1 != NULL){
						score += xs_bb_fxn1->score_of_bin(Xbb);
						score += xs_bb_fxn2->score_of_bin(Xbb.inverse());
				}
                        //std::cout << "score pos" << ir << "," << jr << "," << score << std::endl;
            }
        }
    }
    return(score);
}

void get_second_2_helices(Dock & dock, core::pose::Pose pose,Size & start_res,Size & end_res ){
    using protocols::loops::Loops;
    //get helices----------------
    protocols::loops::LoopsOP helices( new protocols::loops::Loops );
    get_helices(pose,helices);
    //Figure out which 2 helices.
    int dist_to_end = pose.total_residue()-numeric::max(dock.h1_pos_start,dock.h1_pos_end,dock.h2_pos_start,dock.h2_pos_end);

    int dist_to_begin = numeric::min(dock.h1_pos_start,dock.h1_pos_end,dock.h2_pos_start,dock.h2_pos_end)-1;
    if(dist_to_end>dist_to_begin){
        start_res = (*helices)[3].start();
        end_res = (*helices)[4].stop();
    }
    else{
        start_res = (*helices)[(*helices).size()-3].start();
        end_res = (*helices)[(*helices).size()-2].stop();
    }
}

core::pose::Pose superimpose_A_to_all(Dock & dock,core::pose::Pose native_pose){
    using namespace core::scoring;
	using namespace core::id;
    //get chain A.. superimpose native region to B,C,D then combine
    core::pose::Pose chainA = dock.closed_pose;
    core::id::AtomID_Map< core::id::AtomID > atom_map;
    core::pose::initialize_atomid_map( atom_map, chainA, BOGUS_ATOM_ID );
   	for (Size ii=1; ii<=dock.nat_final_poseCoords_end-dock.nat_final_poseCoords_start; ++ii){
        Size native_pos =dock.nat_final_natCoords_start-1+ii;  
        Size chainA_pos =dock.nat_final_poseCoords_start-1+ii;
        core::id::AtomID const id1(chainA.residue(chainA_pos).atom_index("CA"), chainA_pos );
        core::id::AtomID const id2(native_pose.residue(native_pos).atom_index("CA"),native_pos);
        atom_map[id1]=id2;
    }
    std::string const & alphabet( utility::LETTERS ); //A is excluded because of the matching
    core::pose::Pose finalPose;
    BOOST_FOREACH( char c, alphabet) {
    Size chain_id = -1;
    try{
        chain_id = get_chain_id_from_chain(c,native_pose);
        core::pose::PoseOP tmp_chain_pose = native_pose.split_by_chain(chain_id);
        superimpose_pose(chainA,*tmp_chain_pose,atom_map);
        append_pose_to_pose(finalPose,chainA);
        }
        catch(utility::excn::EXCN_Base const & e ){
            //not a worry. I expect many chains to not exist.
        }
    }
    core::scoring::ScoreFunctionOP  lowres_scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "score3" ) );
    Real score = lowres_scorefxn->score(finalPose);
    dock.centroid_score = score;
    return(finalPose);
}


void calculate_motif_scores(Dock & dock,core::scoring::motif::MotifHashManager *mman_){
   //step1---which residues to score over.
    Size start_res;
    Size end_res;
    get_second_2_helices(dock,dock.unclosed_pose,start_res,end_res);
    Real fullPoseScore = compute_motif_score(dock.unclosed_pose_full,start_res,end_res,mman_);
    Real partialPoseScore = compute_motif_score(dock.unclosed_pose,start_res,end_res,mman_);
    dock.chainA_motif_score = partialPoseScore/((Real) end_res-start_res); 
    dock.overall_motif_score =fullPoseScore/((Real) end_res-start_res);
}

int main( int argc, char * argv [] ) {
    try {
    using namespace core::chemical;
    using namespace core::import_pose::pose_stream;
    using core::import_pose::pose_from_file;
    using namespace core::scoring;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    devel::init(argc, argv);
    Real RMSD_1_HELIX_THRESH = 0.5;
    Real RMSD_2_HELIX_THRESH = 2.0;
    Real CENTROID_SCORE_THRESH = 0.0;
    Real VDW_SCORE_THRESH = 30.0;
    Real MOTIF_SCORE_THRESH = -1.0;             //2.25;
    //create vector of input poses.
    MetaPoseInputStream input = streams_from_cmd_line();
    core::pose::Pose current_pose;
    core::pose::Pose native_pose;
    ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
	pose_from_file(
			native_pose,
			*rsd_set,
			option[ in::file::native ]()
			);
    core::scoring::motif::MotifHashManager *mman_ = core::scoring::motif::MotifHashManager::get_instance();
    core::util::switch_to_residue_type_set( native_pose, core::chemical::CENTROID);
    while(input.has_another_pose()){
        input.fill_pose(current_pose,*rsd_set_from_cmd_line().lock());
        core::util::switch_to_residue_type_set( current_pose, core::chemical::CENTROID);
        std::string input_tag = core::pose::tag_from_pose(current_pose);
        vector1< std::string > pos = utility::string_split( input_tag, '/');
        vector1< std::string > pos1 = utility::string_split( pos[pos.size()],'.');
        std::string tag = pos1[1];
        Size counter = 1;
        vector1<Dock> docks = get_all_2helix_docks(current_pose,native_pose);
        vector1<Dock> tmp_docks = get_all_1helix_docks(current_pose,native_pose);
        docks.insert(docks.end(),tmp_docks.begin(),tmp_docks.end());
        core::io::silent::SilentFileData sfd;  
        for(Size ii=1; ii<=docks.size(); ++ii){
            docks[ii].rmsd =calculate_helical_tail_variance(current_pose,native_pose,docks[ii]);
            if(((docks[ii].single_helix && (docks[ii].rmsd < RMSD_1_HELIX_THRESH)))||(!docks[ii].single_helix && (docks[ii].rmsd < RMSD_2_HELIX_THRESH))){
                generate_combined_model(current_pose, native_pose,docks[ii]);
                if(docks[ii].centroid_score < CENTROID_SCORE_THRESH && docks[ii].vdw_score < VDW_SCORE_THRESH ){
                    calculate_motif_scores(docks[ii],mman_);
                    if((docks[ii].chainA_motif_score == docks[ii].overall_motif_score) && (docks[ii].chainA_motif_score > MOTIF_SCORE_THRESH)){
                        std::stringstream ss;
                        ss << tag << "_X_" << counter << ".pdb";
                        std::string s2 = ss.str();
                        counter += 1;
                        //docks[ii].closed_pose = minimize_to_close(docks[ii].unclosed_pose);
                        docks[ii].closed_pose = docks[ii].unclosed_pose;
                        docks[ii].final_pose = superimpose_A_to_all(docks[ii],native_pose);
                        //A. ChainAll score3 
                        //B. Motif_score
                        //C. Single or double helix
                        //D. Rmsd
                        //E. native_start
                        //F. native_end
                        //G. repeat_start
                        //H. repeat_end
                        //I. Redesign_start
                        //J. Redesign_end;
                        pose::setPoseExtraScore( docks[ii].final_pose, "score3", docks[ii].centroid_score );
                        pose::setPoseExtraScore( docks[ii].final_pose, "motif", docks[ii].overall_motif_score );
                        pose::setPoseExtraScore( docks[ii].final_pose, "singleHelix", docks[ii].single_helix );
                        pose::setPoseExtraScore( docks[ii].final_pose, "rmsd",docks[ii].rmsd );
                        pose::setPoseExtraScore( docks[ii].final_pose, "native_start",docks[ii].nat_final_poseCoords_start );
                        pose::setPoseExtraScore( docks[ii].final_pose, "native_end", docks[ii].nat_final_poseCoords_end );
                        pose::setPoseExtraScore( docks[ii].final_pose, "repeat_start",docks[ii].repeat_final_poseCoords_start );
                        pose::setPoseExtraScore( docks[ii].final_pose, "repeat_end",docks[ii].repeat_final_poseCoords_end );
                        pose::setPoseExtraScore( docks[ii].final_pose, "redesign_start",docks[ii].redesign_final_poseCoords_start );
                        pose::setPoseExtraScore( docks[ii].final_pose, "redesign_end",docks[ii].redesign_final_poseCoords_end );
                        core::io::silent::SilentStructOP ss_out_all( new core::io::silent::BinarySilentStruct );
                        ss_out_all->fill_struct(docks[ii].final_pose,s2);
                        sfd.write_silent_struct( *ss_out_all, option[ basic::options::OptionKeys::out::file::silent ]() );
                    }
                }
            }
        }
        //generate_combined_model(current_pose, native_pose,docks[1]);
    }
    }catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
        return -1;
    }
return 0;
}

