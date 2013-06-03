// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/jianqing/antibody_H3.cc
/// @brief Code to model and antibody CDR H3 loop, simultaneously optimizing surrounding CDR loops and VH/VL orientation
/// @author Jianqing Xu (xubest@gmail.com), Daisuke Kuroda, Brian Weitzer, Jeff Gray
/// 09/09/2011 - 2013




#include <core/pose/Pose.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyUtil.hh>

//#include <protocols/jd2/JobDistributor.hh>
//#include <protocols/jd2/util.hh>

#include <devel/init.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// option key includes
#include <basic/options/option.hh>
#include <string>
#include <basic/Tracer.hh>

// non JD2 includes
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>



static basic::Tracer TR("protocols.antibody_metrics");


//namespaces
using namespace core;
using namespace core::pose;
using namespace utility;
using namespace protocols::antibody;


int
main( int argc, char* argv[] ){

	devel::init(argc, argv);

	//-s read in PDB
	core::pose::Pose pose;
	std::string pdbname(basic::options::option[ basic::options::OptionKeys::in::file::s ].value()[1]);
	core::import_pose::pose_from_pdb( pose, pdbname );

	// Build atom subsets for computing SASA
	AntibodyInfo ab_info = AntibodyInfo(pose);
	TR << ab_info;

	core::Real bbHbond = kink_bb_Hbond(pose, ab_info);
	TR << "bbHbond: " << bbHbond << std::endl;

	core::Real Hbond = kink_Hbond(pose, ab_info);
	TR << "Hbond: " << Hbond << std::endl;
	
	core::Real WHbond = kink_Trp_Hbond(pose, ab_info);
	TR << "WHbond: " << WHbond << std::endl;
	
	std::pair<core::Real,core::Real> q = kink_dihedral(pose, ab_info);
	TR << "q: " << q.first << " " << q.second << std::endl;

	return 0;
}


//
//
//int
//main( int argc, char * argv [] )
//{
//	try {
//
//	using namespace basic::options;
//	using namespace protocols::antibody;
//	using namespace protocols::jd2;
//
//	AntibodyModelerProtocol::register_options();
//	protocols::jd2::register_options();
//	// initialize core
//	devel::init(argc, argv);
//
//
//	AntibodyModelerProtocolOP ab_m_h3 = new AntibodyModelerProtocol();
//	TR<<*ab_m_h3<<std::endl;
//	JobDistributor::get_instance()->go(ab_m_h3);
//
//	 } catch ( utility::excn::EXCN_Base const & e ) { 
//		 std::cout << "caught exception " << e.msg() << std::endl;
//	}
//}
//
//
//
