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
	End
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
                 AntibodyNumberingEnum const & numbering_scheme = Aroop);


public:
	
	/// @brief return the loop of a certain loop type
    loops::LoopsOP get_CDR_loop( AntibodyCDRNameEnum const & cdr_name ) const {
		return vector1_of_all_cdr_loopsOP_[cdr_name];
	}
	
	/// @brief return the loop of a certain loop type
    loops::Loop get_one_cdr_loop_object( AntibodyCDRNameEnum const & cdr_name ) const{
		return (*get_CDR_loop(cdr_name))[1];
	}
	
	/// @brief return the sequence of a particular CDR loop
	vector1<char> get_CDR_Sequence_with_Stem( AntibodyCDRNameEnum const & cdr_name,
											 Size left_stem = 0,
											 Size right_stem = 0) const;
	
	/// @brief return the antibody sequence of LH or just H
	vector1<char> get_Ab_Sequence() const{
		return ab_sequence_;
	}
    
    /// @brief intput an enum, and get a string for it
    std::string get_CDR_Name(AntibodyCDRNameEnum const & cdr_name) const {
        return cdr_name_[cdr_name];
    }
	/// FIXME: this is redundent
	/// @brief return a LoopsOP object, which saves all the CDR Loop object 
	loops::LoopsOP get_all_cdr_loops() const {
        return all_cdr_loops_;
    }
    
    /// @brief return this antibody is camelid or not
    bool is_camelid()  const {
        return is_camelid_;
    }
    
	/// @brief return the framework numbering information 
   	vector1< vector1<FrameWork> > get_AntibodyFrameworkInfo() const{
		return framework_info_;
	}
	
    // JQX: all the FoldTree Stuff
    void all_cdr_fold_tree( pose::Pose & pose );
    //	void cdr_h3_fold_tree( pose::Pose & pose );
    void all_cdr_VL_VH_fold_tree( pose::Pose & pose_in);
    kinematics::FoldTree get_foldtree_LH_A( pose::Pose const & pose ) const;
    kinematics::FoldTree get_foldtree_L_HA( pose::Pose const & pose ) const;
    kinematics::FoldTree get_foldtree_LA_H( pose::Pose const & pose ) const;
    
    
    // JQX: all the MoveMap Stuff
    kinematics::MoveMapOP get_movemap_allCDRbb(pose::Pose const & pose);
    kinematics::MoveMapOP get_movemap_OneCDRbb(pose::Pose const & pose,
											   AntibodyCDRNameEnum const & cdr_loop);
	
	/// @brief TaskFactory
    pack::task::TaskFactoryOP get_taskfctory_allCDRs(pose::Pose  & pose);
    
    /// @brief return num of cdr loops, 3 (nanobody) or 6 (regular antibody)
    AntibodyCDRNameEnum get_total_num_cdr_loops() const {
        return tot_cdr_loops_enum_;
    }

    /// @brief return whether this pose has antigen or not
    bool get_pose_has_antigen() const {
        return InputPose_has_antigen_;
    }

    /// @brief predict H3 cterminal kink/extended conformation
    H3BaseTypeEnum get_predicted_H3_base_type() const {
        return H3_base_type_predicted_;
    }

    /// @brief use the H3 cterm coordinates in the pose to calculate the cterminal type
    std::string calculate_H3_base_by_coordinates(pose::Pose const & pose) const;
	
	void show( std::ostream & out=std::cout );
    friend std::ostream & operator<<(std::ostream& out, const AntibodyInfo & ab_info ) ;
    
	
// private functions
private:

    void set_default();
    
    /// @brief check the input pose is nanobody, antibody or wrong
    void check_AntibodyRelatedPose(pose::Pose const & pose);
    
    /// @brief initialization 
    void init(pose::Pose const & pose);
    
	/// @brief setup the CDR loops objects based on the input numbering scheme
    void setup_CDRsInfo( pose::Pose const & pose );
	
    /// @brief setup the framework information based on the input numbering scheme
    void setup_FrameWorkInfo(pose::Pose const & pose );
    
    /// @brief predict H3 cterminus base as Kinked or Extended
    void predict_H3_base_type( pose::Pose const & pose ) ;
    void detect_and_set_camelid_CDR_H3_stem_type(pose::Pose const & pose );
	void detect_and_set_regular_CDR_H3_stem_type( pose::Pose const & pose );
    void detect_and_set_regular_CDR_H3_stem_type_new_rule( pose::Pose const & pose );
	
	/// @brief return the numbering scheme: e.g.   numbering[Begin][h1]
	vector1< vector1<int> > get_CDR_NumberingInfo(AntibodyNumberingEnum const & numbering_scheme) const;
	
	
	/// @brief: get the current numbeirng scheme being used
	AntibodyNumberingEnum get_current_AntibodyNumberingScheme() const {
		return numbering_scheme_;
	}

	/// @brief identify CDRs on L or H sequence
	void identify_CDR_from_a_sequence(std::string & querychain);
	

	
/// private members
private:      
	/// some Enum types
    AntibodyNumberingEnum numbering_scheme_;
	H3BaseTypeEnum H3_base_type_predicted_;
	AntibodyCDRNameEnum tot_cdr_loops_enum_;
	
	/// the information of the antibody pose 
	bool is_camelid_;
    bool InputPose_has_antigen_;
	
	/// the CDR and Framework information
	vector1<loops::LoopsOP> vector1_of_all_cdr_loopsOP_;
	loops::LoopsOP all_cdr_loops_;
	vector1< vector1<FrameWork> > framework_info_ ;
	vector1<char> ab_sequence_;
    
	/// use vector1 to map Enum to string, for some display purpose
    vector1<std::string> cdr_name_;
	vector1<std::string> h3_base_type_;

};


} //namespace antibody2
} //namespace protocols


#endif //INCLUDED_protocols_loops_AntibodyInfo_HH
