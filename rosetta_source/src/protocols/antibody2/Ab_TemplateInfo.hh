// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocls/antibody2/Ab_TemplateInfo.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody2_Ab_TemplateInfo_hh
#define INCLUDED_protocols_antibody2_Ab_TemplateInfo_hh


#include <map>
#include <core/pose/Pose.hh>


///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace antibody2 {
using namespace core;
    
/// antibody2 definition
class Ab_TemplateInfo {
        
        
        
public:
    typedef std::map <std::string, core::pose::Pose> TemplatePoseMap;
        
    /// default constructor
    Ab_TemplateInfo();
        
        
    pose::Pose get_one_template_pose(std::string template_name);
        
    
	//void show( std::ostream & out=std::cout );
	void show( std::ostream & out );
	friend std::ostream & operator<<(std::ostream& out, const Ab_TemplateInfo & ab_t_info );

        
        
private:
    TemplatePoseMap templates_poses_;
    pose::Pose L1_t_pose_, L2_t_pose_, L3_t_pose_, 
               H1_t_pose_, H2_t_pose_, H3_t_pose_;
//    pose::Pose Lfr_t_pose_, Hfr_t_pose_, LightHeavy_t_pose_;
  

    
    std::string LightHeavy_t_name_, Lfr_t_name_, Hfr_t_name_;
    std::string L1_t_name_, L2_t_name_, L3_t_name_,
                H1_t_name_, H2_t_name_, H3_t_name_;
    
    void obtain_templates_names();
    void load_templates_from_pdbs();







};
    
    
    

    



} //namespace antibody2
} //namespace protocols

     
     


#endif
