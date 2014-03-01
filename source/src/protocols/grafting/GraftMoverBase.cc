// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/GraftMoverBase.cc
/// @brief   Base class for GraftMovers
/// @author  Jared Adolf-Bryfogle (jadolfbr@gmail.com)

//Unit headers
#include <protocols/grafting/GraftMoverBase.hh>
#include <protocols/moves/Mover.hh>

//Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/TorsionID.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>

//Protocol Headers
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/util.hh>


//Basic Headers
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/PyAssert.hh>
#include <protocols/grafting/util.hh>



static basic::Tracer TR("protocols.grafting.GraftMoverBase");

namespace protocols {
namespace grafting {

using namespace core;
using namespace core::pose;
using namespace protocols::moves;
using namespace core::scoring;
using core::Size;
using core::id::AtomID;
using core::id::AtomID_Map;
using core::pose::initialize_atomid_map;
using protocols::loops::Loop;
using core::kinematics::MoveMapOP;
using core::kinematics::MoveMap;

GraftMoverBase::GraftMoverBase(Size start, Size end, std::string mover_name):
	moves::Mover(mover_name),
	start_(start),
	end_(end),
	original_end_(end)

{
}

///@brief copy ctor
//GraftMoverBase::GraftMoverBase( GraftMoverBase const & rhs ) :
//	Mover(rhs)
//{
//	*this = rhs;
//}

// Destructor
GraftMoverBase::~GraftMoverBase() {}


void
GraftMoverBase::set_piece(Pose & piece, Size Nter_overhang, Size Cter_overhang){
	piece_ = piece;
	Nter_overhang_ = Nter_overhang;
	Cter_overhang_=Cter_overhang;
	insertion_length_ = piece.total_residue()-Cter_overhang_-Nter_overhang_;
}

void
GraftMoverBase::set_insert_region(const Size start, const Size end){
    start_=start;
    end_=end;
    original_end_=end;
}

void
GraftMoverBase::set_cen_scorefunction(ScoreFunctionOP score){
    cen_scorefxn_=score->clone();
}

void
GraftMoverBase::set_fa_scorefunction(ScoreFunctionOP score){
	//Just in case modifications are made to it.
	fa_scorefxn_=score->clone();
}

void
GraftMoverBase::set_default_cen_scorefunction(){
    cen_scorefxn_=protocols::loops::get_cen_scorefxn();
}
void
GraftMoverBase::set_default_fa_scorefunction(){
    fa_scorefxn_=getScoreFunction();
}

void
GraftMoverBase::delete_overhang_residues(){

	//Nter
    for (core::Size i=1; i<=Nter_overhang_; ++i){
        //piece_.conformation().delete_polymer_residue(piece_res_num);
        piece_.conformation().delete_residue_slow(1);
    }

    //Cter
    core::Size res_num_start = (piece_.total_residue() - Cter_overhang_ + 1);
    for (core::Size i=1; i<=Cter_overhang_; ++i){
        //piece_.conformation().delete_polymer_residue(piece_res_num);
        piece_.conformation().delete_residue_slow(res_num_start);
    }
    set_overhang(0, 0);
}

void
GraftMoverBase::set_overhang(Size Nter_overhang, Size Cter_overhang){
	Nter_overhang_=Nter_overhang;
	Cter_overhang_=Cter_overhang;
}

void
GraftMoverBase::superimpose_overhangs_heavy(Pose const & pose, bool ca_only, bool silence_rms){
	AtomID_Map < AtomID > atoms_to_superimpose;
	initialize_atomid_map( atoms_to_superimpose, piece_, core::id::BOGUS_ATOM_ID );
	//Remove termini
	remove_lower_terminus_type_from_pose_residue(piece_, 1);
	remove_upper_terminus_type_from_pose_residue(piece_, piece_.total_residue());

	//Nter residues/atoms
	for (core::Size piece_res_num=1; piece_res_num<=Nter_overhang_; ++piece_res_num){
		core::Size pose_res_num = start_ - Nter_overhang_ + piece_res_num;
		if (ca_only){
			AtomID const atom_piece(piece_.residue(piece_res_num).atom_index("CA"), piece_res_num);
			AtomID const atom_pose(pose.residue(pose_res_num).atom_index("CA"), pose_res_num);
			atoms_to_superimpose[ atom_piece ]=atom_pose;
		}
		else{
			for (core::Size atom_num=1; atom_num <=4; atom_num++){
				AtomID const atom_piece(atom_num, piece_res_num);
				AtomID const atom_pose(atom_num, pose_res_num);
				atoms_to_superimpose[ atom_piece ]=atom_pose;
			}
		}
	}

	//Cter residues/atoms
	core::Size i = 0;
	for (core::Size piece_res_num=(piece_.total_residue() - Cter_overhang_ + 1); piece_res_num<=piece_.total_residue(); ++piece_res_num){
		core::Size pose_res_num = end_ + i;
		if (ca_only){
			AtomID const atom_piece(piece_.residue(piece_res_num).atom_index("CA"), piece_res_num);
			AtomID const atom_pose(pose.residue(pose_res_num).atom_index("CA"), pose_res_num);
			atoms_to_superimpose[ atom_piece ]=atom_pose;
		}
		else{
			for (core::Size atom_num=1; atom_num <=4; ++atom_num){
				AtomID const atom_piece(atom_num, piece_res_num);
				AtomID const atom_pose( atom_num, pose_res_num);
				atoms_to_superimpose[ atom_piece ]=atom_pose;
			}
		}
		++i;
	}

	core::Real rms(superimpose_pose(piece_, pose, atoms_to_superimpose));

	if (!silence_rms){
		TR << "RMS of piece from starting position to current: "<<rms<<std::endl;
	}
	return;
}

MoveMapOP
GraftMoverBase::combine_movemaps(MoveMap const & scaffold_mm, MoveMap const & insert_mm){
	using namespace core::kinematics;
	MoveMapOP mm = new MoveMap();
    TR<<"Combining movemaps"<<std::endl;
	//First time with typedefs  They are pretty narley!
	//typedef core::id::TorsionID TorsionID;
	typedef std::pair< Size, core::id::TorsionType > MoveMapTorsionID;
	typedef std::map< MoveMapTorsionID, bool > MoveMapTorsionID_Map;

	//Assert that a piece is set.

	//Copy Nterminal MoveMapTorsionIDs
	for (MoveMapTorsionID_Map::const_iterator it=scaffold_mm.movemap_torsion_id_begin(), it_end=scaffold_mm.movemap_torsion_id_end(); it !=it_end; ++it){
	    //Scaffold to new MM
	    if (it->first.first<=start_){
		MoveMapTorsionID new_id = MoveMapTorsionID(it->first.first, it->first.second);
		mm->set(new_id, it->second);
	    }
	}

	//Set insert residues
	for (MoveMapTorsionID_Map::const_iterator it=insert_mm.movemap_torsion_id_begin(), it_end=insert_mm.movemap_torsion_id_end(); it !=it_end; ++it){
	    //Add from start_
	    Size new_resnum = start_ + it->first.first;
	    MoveMapTorsionID new_id = MoveMapTorsionID(new_resnum, it->first.second);
	    mm->set(new_id, it->second);
	}

	//Set Cterminal residues after insert. We may have a deletion then an insertion that we need to change the numbers for.
	Size deleted_residues = original_end_-start_-1;
	for (MoveMapTorsionID_Map::const_iterator it=scaffold_mm.movemap_torsion_id_begin(), it_end=scaffold_mm.movemap_torsion_id_end(); it !=it_end; ++it){
	    //Check if residue exists in movemap.  Copy that info no matter what it is to the new movemap.
	    if (it->first.first >=original_end_) {
	    //Add insertion length
		Size new_resnum = it->first.first - deleted_residues+insertion_length_;
		MoveMapTorsionID new_id = MoveMapTorsionID(new_resnum, it->first.second);
		mm->set(new_id, it->second);
	    }
	}
	mm->show(TR);
	return mm;
}

Pose
GraftMoverBase::insert_piece(Pose const & pose){

	//Delete overhang if necessary.
	delete_overhang_residues();

	//strip termini variants from insert if necessary
	core::pose::remove_variant_type_from_pose_residue(piece_, core::chemical::LOWER_TERMINUS, 1);
	core::pose::remove_variant_type_from_pose_residue(piece_, core::chemical::UPPER_TERMINUS, insertion_length_);

    //Delete residues for pose.
	Pose final_pose(pose);
	delete_region(final_pose, start_+1, end_-1);
	//Pose insert(piece_); Does the piece still exist afterward?
	final_pose = insert_pose_into_pose(final_pose, piece_, start_, start_+1);

	//Update end residue number.
	end_ = start_+insertion_length_+1;
    TR <<"Insertion complete."<<std::endl;
	return final_pose;

}



core::Real
GraftMoverBase::perturb_backbone_for_test(Pose& pose, core::kinematics::MoveMapOP mm){
	Pose pre_test(pose);
	TR <<"Randomizing residues in movemap." <<std::endl;
	protocols::simple_moves::SmallMover randomize(mm, 10000, 500); //Huge KT to randomize the residues
	randomize.angle_max( 'H', 180.0 );
	randomize.angle_max( 'E', 180.0 );
	randomize.angle_max( 'L', 180.0 );
	randomize.apply(pose);

	Real RMS = core::scoring::bb_rmsd_including_O(pre_test, pose);
	TR <<"All BB Heavy Atom RMSD "<< RMS <<std::endl;
	return RMS;
	}


//////////////////////////////////////////////////////////////////////////////////////////////////////
///MOVEMAP and REGION SETUP
///Note: Only these functions interact with class variables defined here. To be used optionally in apply.
//////////////////////////////////////////////////////////////////////////////////////////////////////

void
GraftMoverBase::setup_movemap_and_regions(Pose & pose){
    if (!use_default_movemap_){
        movemap_ = combine_movemaps(*scaffold_movemap_, *insert_movemap_);
        set_regions_from_movemap(pose);

    }
    else{
    	set_regions_from_flexibility();
	    set_default_movemap();
    }

    PyAssert((Nter_loop_start_!=Cter_loop_end_), "Start and end of remodeling region the same.  Please check movemap");
}

void
GraftMoverBase::set_default_movemap(){
    TR <<"Setting default movemap"<<std::endl;
	movemap_ = new core::kinematics::MoveMap();
	for (Size i=Nter_loop_start_; i<=Nter_loop_end_; ++i){
		movemap_->set_bb(i, true);
		movemap_->set_chi(i, true);
	}
	for (Size i=Cter_loop_start_; i<=Cter_loop_end_; ++i){
		movemap_->set_bb(i, true);
		movemap_->set_chi(i, true);
	}
}

void
GraftMoverBase::set_regions_from_flexibility(){
	Nter_loop_start_ = start_-Nter_scaffold_flexibility_+1;//(loop_start/insert_loop_start)First Flexible loop residue
    Nter_loop_end_ = start_+Nter_insert_flexibility_;

    Cter_loop_start_=start_+insertion_length_+1-Cter_insert_flexibility_;
    Cter_loop_end_ = start_+insertion_length_+Cter_scaffold_flexibility_;//(loop_end) Last Flexible loop residue

    //TR <<"Nter_loop_start: "<<Nter_loop_start_<<std::endl;
    //TR <<"Nter_loop_end: "<<Nter_loop_end_<<std::endl;
    //TR <<"Cter_loop_start:  "<<Cter_loop_start_<<std::endl;
    //TR <<"Cter_loop_end:  "<<Cter_loop_end_<<std::endl;
}


void
GraftMoverBase::set_regions_from_movemap(Pose & pose){
    //N terminal
    for (Size i=1; i<=start_; ++i){
        if (movemap_->get_bb(i)){
            Nter_loop_start_=i;
            break;
        }
    }

    //C terminal end
    for (Size i=pose.total_residue(); i>=Nter_loop_start_; --i){
        if (movemap_->get_bb(i)){
        	TR <<i<<std::endl;
            Cter_loop_end_=i;
            break;
        }
    }

    //Insertion default flexiblity.  Will not effect single_loop foldtrees.
    Nter_insert_flexibility_=0;
    for (Size i = start_+1; i<=start_+insertion_length_; ++i){
        if (movemap_->get_bb(i)){
            ++Nter_insert_flexibility_;
        }
        else{
            Nter_loop_end_=start_+Nter_insert_flexibility_;

            break;
        }
    }

    Cter_insert_flexibility_=0;
    for (Size i = start_+insertion_length_; i>=1; --i){
        if (movemap_->get_bb(i)){
            ++Cter_insert_flexibility_;
        }
        else{
            Cter_loop_start_=start_+insertion_length_-Cter_insert_flexibility_;
            break;
        }
    }
    //TR <<"Nter_insert_flex: "<<Nter_insert_flexibility_<<std::endl;
    //TR <<"Nter_loop_start: "<<Nter_loop_start_<<std::endl;
    //TR <<"Nter_loop_end: "<<Nter_loop_end_<<std::endl;
    //TR <<"Cter_loop_start:  "<<Cter_loop_start_<<std::endl;
    //TR <<"Cter_loop_end:  "<<Cter_loop_end_<<std::endl;
}





/////////////////////////////////////////////////////////////////////////////////////////////////
///FOLDTREE SETUP.   options depending on how you want your graft algorithm to work!
///---Indicates Flexible regions, | indicates cutpoint. Arrows are direction of ARMs used to close the loop in conjunction with algorithm. (CCD, KIC)
///
/////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
/// @brief ****Nter_loop_start---->Piece----> | Cter_loop_end****
/// Insert will move in cartesian space
/// @params lower_cutpoint for CCD and loops is Cter_loop_end-1
///
void
GraftMoverBase::setup_single_loop_single_arm_remodeling_foldtree(Pose & pose, Size const Nter_loop_start, Size const Cter_loop_end, bool loop_modeling){
	using core::kinematics::Edge;

	//Setup the offset so the anchor of the loop is correct for either CCD loop closure or vanilla loop modeling.
	Size anchor_offset = 1;
	if (loop_modeling){anchor_offset=2;}

	core::Size const cutpoint_lower(Cter_loop_end-1); //cutpoint at the end of the loop
	core::Size const cutpoint_upper(Cter_loop_end); //cutpoint at the end of the loop
	core::Size const loop_start_foldtree_anchor(Nter_loop_start-anchor_offset); //this is the N-terminal jump anchor for the loop
	core::Size const loop_end_foldtree_anchor(Cter_loop_end+anchor_offset); //C-terminal jump anchor

	TR << "loop_start_foldtree_anchor " << loop_start_foldtree_anchor << std::endl;
	TR << "cutpoint_lower " << cutpoint_lower << std::endl;
	TR << "cutpoint_upper " << cutpoint_upper << std::endl;
	TR << "loop_end_foldtree_anchor " << loop_end_foldtree_anchor << std::endl;

	core::kinematics::FoldTree remodeling_tree(pose.total_residue());
	remodeling_tree.clear();
	remodeling_tree.add_edge(Edge(1, loop_start_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, cutpoint_lower, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(cutpoint_upper, loop_end_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_end_foldtree_anchor, pose.total_residue(), Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, loop_end_foldtree_anchor, 1));
	remodeling_tree.reorder(1);
	TR << remodeling_tree << std::endl;
	pose.fold_tree(remodeling_tree);
}


//////////////////////////////////////////////////////////////////
/// @brief ****Nter_loop_start---->Piece | <----Nter_loop_end****
/// Insert will move in cartesian space
/// @params lower_cutpoint for CCD and loops is end_-1
///
void
GraftMoverBase::setup_single_loop_double_arm_remodeling_foldtree(Pose & pose, Size const Nter_loop_start, Size const Cter_loop_end, bool loop_modeling){
	//Note: Differs from single arm by moving the cutpoint.
	TR <<"Setting up single loop, double arm remodeling foldtree"<<std::endl;
	using core::kinematics::Edge;
	insertion_length_ = piece_.total_residue();//Double check.

	//Setup the offset so the anchor of the loop is correct for either CCD loop closure or vanilla loop modeling.
	Size anchor_offset = 1;
	if (loop_modeling){anchor_offset=2;}

	core::Size const cutpoint_lower(end_-1);//Cutpoint is at the end of the insert
	core::Size const cutpoint_upper(end_);
	core::Size const loop_start_foldtree_anchor(Nter_loop_start-anchor_offset); //this is the N-terminal jump anchor for the loop
	core::Size const loop_end_foldtree_anchor(Cter_loop_end+anchor_offset); //C-terminal jump anchor

	TR << "loop_start_foldtree_anchor " << loop_start_foldtree_anchor << std::endl;
	TR << "cutpoint_lower " << cutpoint_lower << std::endl;
	TR << "cutpoint_upper " << cutpoint_upper << std::endl;
	TR << "loop_end_foldtree_anchor " << loop_end_foldtree_anchor << std::endl;

	core::kinematics::FoldTree remodeling_tree(pose.total_residue());
	remodeling_tree.clear();
	remodeling_tree.add_edge(Edge(1, loop_start_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, cutpoint_lower, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(cutpoint_upper, loop_end_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_end_foldtree_anchor, pose.total_residue(), Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, loop_end_foldtree_anchor, 1));
	remodeling_tree.reorder(1);
	TR << remodeling_tree << std::endl;
	pose.fold_tree(remodeling_tree);
}

// Mover methods
std::string
GraftMoverBase::get_name() const
{
	return "GraftMoverBase";
}



}  // namespace grafting
}  // namespace protocols
