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
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
#include <protocols/docking/types.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <boost/lexical_cast.hpp>


///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace antibody2 {

    /*
     class Antibody_Loop : public loops::Loop{
     public:
     Antibody_Loop(core::pose::Pose const & pose,
     AntibodyCDRNameEnum const & cdr_name,
     AntibodyNumberingEnum const & numbering);
     
     std::string get_Sequence() {return cdr_sequence_;}
     std::string get_ChainID()  {return chain_id_;}
     std::string get_CDRLoopName() {return cdr_name_;}
     std::string get_NumberingInfo() { return numbering_type_ + " numbering: " + chain_id_ + " " +
     boost::lexical_cast<std::string>(numbering_start_) + "-" +
     boost::lexical_cast<std::string>(numbering_stop_) ;}
     
     private:
     AntibodyNumberingEnum numbering_;
     std::string chain_id_; // this CDR is on heavy or light chain
     std::string cdr_sequence_;
     std::string cdr_name_;
     std::string numbering_type_;
     core::Size numbering_start_;
     core::Size numbering_stop_;
     
     };
     */
    
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
    
    
struct FrameWork{
public:
    core::Size start(){return start_;}
    core::Size stop() {return stop_;}
    std::string chain_name(){return chain_name_;}
    void set_start (core::Size num){start_=num;}
    void set_stop  (core::Size num){stop_ =num;}
    void set_chain_name (std::string name){chain_name_ = name;}
private:
    core::Size    start_;
    core::Size    stop_;
    std::string chain_name_;
};
    
    

    

/// antibody2 definition
class AntibodyInfo : public utility::pointer::ReferenceCount {

public:

	/// default constructor
	AntibodyInfo();

	/// constructor with arguments
	AntibodyInfo( core::pose::Pose const & pose );
	AntibodyInfo( core::pose::Pose const & pose,
                 AntibodyNumberingEnum const & numbering_scheme);


    void init(core::pose::Pose const & pose);
    
    /// @brief check the input pose is nanobody, antibody or wrong
    void check_AntibodyRelatedPose(core::pose::Pose const & pose);

    /// @brief setup the CDR loops objects based on the input numbering scheme
    void setup_CDRsInfo( core::pose::Pose const & pose );
    
    /// @brief setup the framework information based on the input numbering scheme
    void setup_FrameWorkInfo(core::pose::Pose const & pose );
    
    /// @brief predict kinked/extended information based on bioinformatics rules
    std::string Predict_K_E_CDRH3( core::pose::Pose const & pose );
    void detect_and_set_camelid_CDR_H3_stem_type();
	void detect_and_set_regular_CDR_H3_stem_type( core::pose::Pose const & pose );
    void detect_and_set_regular_CDR_H3_stem_type_new_rule( core::pose::Pose const & pose );

	/// @brief return the loop of a certain loop type
    loops::LoopOP get_CDR_loop( AntibodyCDRNameEnum const & cdr_name ) const;
    
    std::string get_CDR_Sequence( AntibodyCDRNameEnum const & cdr_name) const {
        return cdr_sequence_[cdr_name];
    };
    
    loops::LoopsOP get_all_cdr_loops(){
        return all_cdr_loops_;
    }
    
    std::string get_CDR_Name(AntibodyCDRNameEnum const & cdr_name) const{
        return cdr_name_[cdr_name];
    }

	// return kinked/extended
	bool is_kinked()   { return kinked_H3_;   }
	bool is_extended() { return extended_H3_; }
    bool is_camelid()  { return is_camelid_;  }
    utility::vector1<char> get_Fv_sequence() { return Fv_sequence_;}


    
    void load_CDR_query_info_to_check();

    void show( std::ostream & out=std::cout );
    friend std::ostream & operator<<(std::ostream& out, const AntibodyInfo & ab_info );

    
    
    // Start coordinates of active loop
	core::Size current_start;
	// End coordinates of active loop
	core::Size current_end;
    

    
    utility::vector1< utility::vector1<int> >  get_AntibodyNumberingScheme(AntibodyNumberingEnum const & numbering_scheme);
    utility::vector1<FrameWork> get_ab_framework(){return Ab_framework_;}
    utility::vector1<FrameWork> get_Hfr(){return Hfr_;}
    utility::vector1<FrameWork> get_Lfr(){return Lfr_;}

	void set_default( );
    void identify_CDR_from_a_sequence(std::string & querychain);


    docking::DockJumps LH_dock_jump(){ return LH_dock_jumps_;}

    std::string LH_dock_partner(){return LH_dock_partners_;}
    
    

    
    
    // JQX: all the FoldTree Stuff
    void all_cdr_fold_tree( core::pose::Pose & pose );
    //	void cdr_h3_fold_tree( core::pose::Pose & pose );
    void all_cdr_VL_VH_fold_tree( core::pose::Pose & pose_in);
    core::kinematics::FoldTreeOP get_foldtree_LH_A( core::pose::Pose const & pose );
    core::kinematics::FoldTreeOP get_foldtree_L_HA( core::pose::Pose const & pose );
    core::kinematics::FoldTreeOP get_foldtree_LA_H( core::pose::Pose const & pose );
    // JQX: all the MoveMap Stuff
    core::kinematics::MoveMapOP get_movemap_allCDRbb(core::pose::Pose const & pose);
    core::kinematics::MoveMapOP get_movemap_OneCDRbb(core::pose::Pose const & pose, AntibodyCDRNameEnum const & cdr_loop);
    core::pack::task::TaskFactoryOP get_taskfctory_allCDRs(core::pose::Pose  & pose);
    
    /// @brief return the num of cdr loops, can be 3 (nanobody) or 6 (regular antibody)
    AntibodyCDRNameEnum get_total_num_cdr_loops(){return tot_cdr_loops_enum_;}

    bool get_antigen_existence_input_pose(){return InputPose_has_antigen_;}
    
private:

    AntibodyNumberingEnum numbering_scheme_;
    std::string CDRH3CterminalPredicition_;

    
    /// @brief a "LoopsOP" object, save cdr "Loop" one by one
    loops::LoopsOP all_cdr_loops_;
    
    /// @brief a "vector1" of "LoopsOP", each "LoopsOP" only has one cdr "Loop"
    utility::vector1<loops::LoopsOP> vector1_of_all_cdr_loopsOP_;
    
    /// @brief a "vector1" of "LoopOP", each "LoopOP" itself is a cdr "LoopOP"
    utility::vector1<loops::LoopOP> vector1_of_all_cdr_loopOP_;
    
	std::string L1_seq_, L2_seq_, L3_seq_;
	std::string H1_seq_,H2_seq_,H3_seq_;

    utility::vector1<std::string> cdr_sequence_;
    utility::vector1<std::string> cdr_name_;
    
	bool is_camelid_;
	bool kinked_H3_;
	bool extended_H3_;
    bool InputPose_has_antigen_;

	utility::vector1< char > Fv_sequence_;


	utility::vector1<FrameWork> Lfr_, Hfr_, Ab_framework_ ;

	docking::DockJumps LH_dock_jumps_;
	std::string LH_dock_partners_;
    
    
	core::kinematics::FoldTreeOP LH_A_foldtree_;
	core::kinematics::FoldTreeOP L_HA_foldtree_;
	core::kinematics::FoldTreeOP LA_H_foldtree_;
    
    core::pose::Pose the_definition_pose_;  // save the pose that was used to define ab_info object
    
    utility::vector1< utility::vector1<int> > numbering_info_;
    
    utility::vector1<char> Chain_IDs_for_CDRs_;
    
    AntibodyCDRNameEnum tot_cdr_loops_enum_;
    
};


} //namespace antibody2
} //namespace protocols


#endif //INCLUDED_protocols_loops_AntibodyInfo_HH




