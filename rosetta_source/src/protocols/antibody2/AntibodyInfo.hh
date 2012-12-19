// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody2/AntibodyInfo.hh
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

#ifndef INCLUDED_protocols_antibody2_AntibodyInfo_hh
#define INCLUDED_protocols_antibody2_AntibodyInfo_hh

#include <protocols/antibody2/AntibodyInfo.fwd.hh>

// Rosetta Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
#include <protocols/docking/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/FoldTree.hh>


///////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace utility;
namespace protocols {
namespace antibody2 {
    
enum AntibodyCDRNameEnum{
    start_cdr_loop = 1,
    h1 = start_cdr_loop,
    h2,
    h3,
    H_chain_last_loop = h3,
    camelid_last_loop = h3,
    l1,
    l2,
    l3,
    L_chain_last_loop = l3,
    num_cdr_loops = l3
};
    
enum AntibodyNumberingEnum{
    Aroop = 1,
    Chothia,
    Kabat,
    Enhanced_Chothia,
    AHO,
    IMGT
};
	
enum BeginEndEnum{
    Begin = 1,
	End,
	Pack_Angle_Begin,
	Pack_Angle_End
};
    
enum H3BaseTypeEnum{
    Kinked = 1,
    Extended,
    Neutral
};
    
    
struct FrameWork{
    Size    start;
    Size    stop;
    char chain_name;
};
    
    

    
class AntibodyInfo : public pointer::ReferenceCount {

public:
	AntibodyInfo( pose::Pose const & pose,
                 AntibodyNumberingEnum const & numbering_scheme = Aroop,
				 bool const & cdr_pdb_numbered = true);


public:
	
	/// @brief: get the current numbeirng scheme being used
	std::string get_Current_AntibodyNumberingScheme() const {
		return get_string_numbering_scheme()[numbering_scheme_];
	}
	
	/// @brief return the loop of a certain loop type
    loops::LoopsOP get_CDR_in_loopsop( AntibodyCDRNameEnum const & cdr_name ) const {
		return vector1_loopsop_having_cdr_[cdr_name];
	}
	
	/// @brief return the loop of a certain loop type
	loops::Loop get_CDR_loop( AntibodyCDRNameEnum const & cdr_name ) const{
		return (*get_CDR_in_loopsop(cdr_name))[1];
	}
	
	/// @brief return a LoopsOP object, which saves all the CDR Loop object
	loops::LoopsOP get_AllCDRs_in_loopsop() const {
        return loopsop_having_allcdrs_;
    }
	
	/// @brief return the sequence of a particular CDR loop
	vector1<char> get_CDR_Sequence_with_Stem( AntibodyCDRNameEnum const & cdr_name,
											 Size left_stem = 0,
											 Size right_stem = 0) const;
	
	/// @brief return the antibody sequence of LH or just H for camelid
	vector1<char> const & get_Ab_Sequence() const {
		return ab_sequence_;
	}
    
    /// @brief intput an enum, and get a string for it
	std::string get_CDR_Name(AntibodyCDRNameEnum const & cdr_name) const {
        return get_string_cdr_name()[cdr_name];
   }

    /// @brief return this antibody is camelid or not
    bool is_Camelid()  const {
        return is_camelid_;
    }
    
	/// @brief return the framework numbering information 
   	vector1< vector1<FrameWork> > get_AntibodyFrameworkInfo() const {
		return framework_info_;
	}
	
    // FoldTrees //TODO: find a way to remove setup_simple_fold_tree
	kinematics::FoldTreeCOP setup_simple_fold_tree(Size const & jumppoint1,
												   Size const & cutpoint,
												   Size const & jumppoint2,
												   pose::Pose const & pose ) const;
	
    kinematics::FoldTreeCOP get_FoldTree_AllCDRs_LHDock( pose::Pose & pose) const;
	kinematics::FoldTreeCOP get_FoldTree_AllCDRs(pose::Pose const & pose) const;
	
	/// @brief SnugDock foldtrees
    kinematics::FoldTree get_FoldTree_LH_A( pose::Pose const & pose ) const;
    kinematics::FoldTree get_FoldTree_L_HA( pose::Pose const & pose ) const;
    kinematics::FoldTree get_FoldTree_LA_H( pose::Pose const & pose ) const;
	
	/// TODO: this should be a standard utility for loops?
	/// @brief get a movemap for loops
	kinematics::MoveMap get_MoveMap_for_Loops(pose::Pose const & pose,
												loops::Loops const & the_loops,
												 bool const & bb_only = false,
												 bool const & include_nb_sc = false,
												 Real const & nb_dist = 10.0) const;
	
	/// @brief get a movemap for loops and set the first jump be true
	kinematics::MoveMap get_MoveMap_for_LoopsandDock(pose::Pose const & pose,
											  loops::Loops const & the_loops,
											  bool const & bb_only = false,
											  bool const & include_nb_sc = false,
											  Real const & nb_dist = 10.0) const;
	
	/// @brief TaskFactory
    pack::task::TaskFactoryOP get_TaskFactory_AllCDRs(pose::Pose & pose) const;
	pack::task::TaskFactoryOP get_TaskFactory_OneCDR(pose::Pose & pose, AntibodyCDRNameEnum const & cdr_name) const;
    
    /// @brief return num of cdr loops, 3 (nanobody) or 6 (regular antibody)
    AntibodyCDRNameEnum get_TotalNumCDRs() const {
        return total_cdr_loops_;
    }

    /// @brief return whether this pose has antigen or not
    bool get_PoseHasAntigen() const {
        return InputPose_has_antigen_;
    }

    /// @brief get H3 cterminal kink/extended conformation (predicted by constructor)
    H3BaseTypeEnum get_Predicted_H3BaseType() const {
        return predicted_H3_base_type_;
    }
	
	/// @brief get residues used to calculate VL/VH packing angle
	vector1< Size > get_PackingAngleResidues() const {
		return packing_angle_residues_;
	}

    /// @brief use the H3 cterm coordinates in the pose to calculate the cterminal type
    //std::string calculate_H3_base_by_coordinates(pose::Pose const & pose) const;
	
	void show( std::ostream & out=std::cout );
    friend std::ostream & operator<<(std::ostream& out, const AntibodyInfo & ab_info ) ;
	
	
// private functions
private:
	/////////////////////////////////////////////////////////////////////////////////////////
	/// 								all the setters									  ///
	/////////////////////////////////////////////////////////////////////////////////////////
	///
    void set_default();
    
    /// @brief check the input pose is nanobody, antibody or wrong
    void identify_antibody(pose::Pose const & pose);
    
    /// @brief initialization 
    void init(pose::Pose const & pose);
    
	/// @brief setup the CDR loops objects based on the input numbering scheme
    void setup_CDRsInfo( pose::Pose const & pose );
	
    /// @brief setup the framework information based on the input numbering scheme
    void setup_FrameWorkInfo(pose::Pose const & pose );
	
	/// @brief setup the residues used to calculate VL/VH packing angle
	void setup_VL_VH_packing_angle( pose::Pose const & pose );
    
    /// @brief predict H3 cterminus base as Kinked or Extended
    void predict_H3_base_type( pose::Pose const & pose ) ;
    void detect_and_set_camelid_CDR_H3_stem_type( pose::Pose const & pose );
	void detect_and_set_regular_CDR_H3_stem_type( pose::Pose const & pose );
    void detect_and_set_regular_CDR_H3_stem_type_new_rule( pose::Pose const & pose );
	
	/// @brief identify CDRs on L or H sequence
	void identify_CDR_from_a_sequence(std::string const & querychain);
	///																					  ///
	/////////////////////////////////////////////////////////////////////////////////////////
	
	
	/// @brief return the numbering scheme: e.g.   numbering[Begin][h1]
	vector1< vector1<Size> > get_CDR_NumberingInfo(AntibodyNumberingEnum const & numbering_scheme) const;

	
	/// @brief copy 
	//void init_for_equal_operator_and_copy_constructor( AntibodyInfo & lhs, AntibodyInfo const & rhs);
	
	/// @brief all the static functions
	static vector1<std::string> const & get_string_cdr_name(void);
	static vector1<std::string> const & get_string_h3_base_type(void);
	static vector1<std::string> const & get_string_numbering_scheme(void);
	static core::scoring::ScoreFunctionCOP get_Pack_ScoreFxn(void);
	static core::scoring::ScoreFunctionCOP get_Dock_ScoreFxn(void);
	static core::scoring::ScoreFunctionCOP get_LoopCentral_ScoreFxn(void);
	static core::scoring::ScoreFunctionCOP get_LoopHighRes_ScoreFxn(void);
	
/// private members
private:
	
	/// the information of the antibody pose 
	bool is_camelid_;
    bool InputPose_has_antigen_;
	bool cdr_pdb_numbered_;
	
	/// the CDR and Framework information
	vector1<loops::LoopsOP> vector1_loopsop_having_cdr_; // each Loops contains one CDR
	loops::LoopsOP loopsop_having_allcdrs_;  // one Loops object containing the set of Loop objects for all CDRs
	vector1< vector1<FrameWork> > framework_info_ ;
	vector1<char> ab_sequence_;
	vector1< Size > packing_angle_residues_;
    
	/// Antibody properties
    AntibodyNumberingEnum numbering_scheme_;
	H3BaseTypeEnum predicted_H3_base_type_;
	AntibodyCDRNameEnum total_cdr_loops_;

};


} //namespace antibody2
} //namespace protocols


#endif //INCLUDED_protocols_loops_AntibodyInfo_HH
