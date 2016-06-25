// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocls/antibody/Ab_TemplateInfo.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @details
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody_Ab_TemplateInfo_hh
#define INCLUDED_protocols_antibody_Ab_TemplateInfo_hh


#include <protocols/antibody/Ab_TemplateInfo.fwd.hh>
#include <map>
#include <core/pose/Pose.hh>
#include <iostream>


///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace antibody {
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

/// @brief Specifically for AntibodyModeling protocol templates. Not for general use.
class Ab_TemplateInfo  : public utility::pointer::ReferenceCount {

public:
	typedef std::map <std::string, core::pose::Pose> TemplatePoseMap;

	/// constructors
	Ab_TemplateInfo();
	Ab_TemplateInfo(bool load_L1, bool load_L2, bool load_L3,
		bool load_H1, bool load_H2, bool load_H3) ;

	Ab_TemplateInfo(bool load_L1, bool load_L2, bool load_L3,
		bool load_H1, bool load_H2, bool load_H3, bool camelid);

	void load_templates_from_pdbs(bool load_L1, bool load_L2, bool load_L3,
		bool load_H1, bool load_H2, bool load_H3, bool camelid);

	pose::Pose get_one_template_pose(std::string cdr_name);
	bool is_camelid()  {
		return camelid_;
	}

	void show( std::ostream & out=std::cout );
	// void show( std::ostream & out );
	friend std::ostream & operator<<(std::ostream& out, const Ab_TemplateInfo & ab_t_info );


private:
	bool load_L1_, load_L2_, load_L3_,
		load_H1_, load_H2_, load_H3_;

	TemplatePoseMap templates_poses_;
	pose::Pose L1_t_pose_, L2_t_pose_, L3_t_pose_,
		H1_t_pose_, H2_t_pose_, H3_t_pose_;
	//    pose::Pose Lfr_t_pose_, Hfr_t_pose_, LightHeavy_t_pose_;

	void set_default( bool camelid );

	std::string LightHeavy_t_name_, Lfr_t_name_, Hfr_t_name_;
	std::string L1_t_name_, L2_t_name_, L3_t_name_,
		H1_t_name_, H2_t_name_, H3_t_name_;

	void obtain_templates_names();


	bool camelid_;


};


} //namespace antibody
} //namespace protocols


#endif


