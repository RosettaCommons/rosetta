// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/jianqing/antibody_assemble_CDRs.cc
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)
/// 09/09/2011



#include <protocols/antibody2/AntibodyModelerProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>

#include <devel/init.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
//#include <utility/tools/make_vector1.hh>

// option key includes
#include <basic/options/option.hh>
#include <string>
#include <basic/Tracer.hh>


static basic::Tracer TR("protocols.antibody2");


int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace protocols::antibody2;
	using namespace protocols::jd2;

	AntibodyModelerProtocol::register_options();
	protocols::jd2::register_options();
	// initialize core
	devel::init(argc, argv);


	AntibodyModelerProtocolOP ab_m_h3 = new AntibodyModelerProtocol();
	TR<<*ab_m_h3<<std::endl;
//    exit(-1);

	JobDistributor::get_instance()->go(ab_m_h3);

}



