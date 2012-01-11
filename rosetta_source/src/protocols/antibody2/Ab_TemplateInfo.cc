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



#include <protocols/antibody2/Ab_TemplateInfo.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <fstream>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBInfo.hh>

static basic::Tracer TR("antibody2.Ab_TemplateInfo");

namespace protocols{
namespace antibody2{
    using namespace core;
    
    
    
Ab_TemplateInfo::Ab_TemplateInfo(){
    obtain_templates_names();
    load_templates_from_pdbs();

}
    
    
    
void Ab_TemplateInfo::obtain_templates_names(){
    std::ifstream inf;
    inf.open("query.matches");
    std::string temp,tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    LightHeavy_t_name_ = tttt;
    
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    L1_t_name_ = tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    L2_t_name_ = tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    L3_t_name_ = tttt;
    
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    H1_t_name_ = tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    H2_t_name_ = tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    H3_t_name_ = tttt;
}
    

    
    
void Ab_TemplateInfo::load_templates_from_pdbs(){    
    

        
    import_pose::pose_from_pdb( L1_t_pose_, "l1.pdb" );
    import_pose::pose_from_pdb( L2_t_pose_, "l2.pdb" );
    import_pose::pose_from_pdb( L3_t_pose_, "l3.pdb" );
    import_pose::pose_from_pdb( H1_t_pose_, "h1.pdb" );
    import_pose::pose_from_pdb( H2_t_pose_, "h2.pdb" );
    import_pose::pose_from_pdb( H3_t_pose_, "h3.pdb" );
        
    
    templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("l1", L1_t_pose_) );
    templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("l2", L2_t_pose_) );
    templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("l3", L3_t_pose_) );
    
    templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("h1", H1_t_pose_) );
    templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("h2", H2_t_pose_) );
    templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("h3", H3_t_pose_) );
        
}
    
    
    
    
pose::Pose Ab_TemplateInfo::get_one_template_pose(std::string template_name){
    
    TemplatePoseMap::iterator iter = templates_poses_.begin();
	iter = templates_poses_.find(template_name);
	if ( iter != templates_poses_.end() ) {return iter->second;}
    
}
    
    
    
    
    
    
    
    
/// @details  Show the complete setup of the Ab_TemplateInfo
void
Ab_TemplateInfo::show( std::ostream & out ) {
    //	if ( !flags_and_objects_are_in_sync_ ){
    //		sync_objects_with_flags();
    //	}
    out << *this;
}
    
std::ostream & operator<<(std::ostream& out, const Ab_TemplateInfo & ab_t_info )
{
    using namespace ObjexxFCL::fmt;
    // All output will be 80 characters - 80 is a nice number, don't you think?
    std::string line_marker = "///";
    out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
        out << line_marker << A( 47, "Rosetta Antibody Template Info" ) << space( 27 ) << line_marker << std::endl;
        out << line_marker << space( 74 ) << line_marker << std::endl;
        // Display the movable jumps that will be used in docking
        out << line_marker << " L1 template: "<<std::endl;
        out << line_marker << "   template_name:  "<<ab_t_info.L1_t_name_<<std::endl;
        out << line_marker << ab_t_info.L1_t_pose_ << std::endl;

        
        out << line_marker << " L2 template: "<<std::endl;
        out << line_marker << "   template_name:  "<<ab_t_info.L2_t_name_<<std::endl;
        out << line_marker << ab_t_info.L2_t_pose_<<std::endl;

        
        out << line_marker << " L3 template: "<<std::endl;
        out << line_marker << "   template_name:  "<<ab_t_info.L3_t_name_<<std::endl;
        out << line_marker << ab_t_info.L3_t_pose_<<std::endl;

        
        
        out << line_marker << " H1 template: "<<std::endl;
        out << line_marker << "   template_name:  "<<ab_t_info.H1_t_name_<<std::endl;
        out << line_marker << ab_t_info.H1_t_pose_<<std::endl;

        
        out << line_marker << " H2 template: "<<std::endl;
        out << line_marker << "   template_name:  "<<ab_t_info.H2_t_name_<<std::endl;
        out << line_marker << ab_t_info.H2_t_pose_<<std::endl;

        
        out << line_marker << " H3 template: "<<std::endl;
        out << line_marker << "   template_name:  "<<ab_t_info.H3_t_name_<<std::endl;
        out << line_marker << ab_t_info.H3_t_pose_<<std::endl;

        
        // Close the box I have drawn
        out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
        return out;
    }

    
    
    
    
    
} //namespace antibody2
} //namespace protocols
