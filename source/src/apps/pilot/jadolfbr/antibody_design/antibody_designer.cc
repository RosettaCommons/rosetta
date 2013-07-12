// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/antibody_designer.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/antibody/design/AntibodyDesignMover.hh>

#include <utility/excn/Exceptions.hh>
#include <devel/init.hh>


using namespace protocols::antibody::design;

//Description:  This application will become the Rosetta Antibody Designer.  Main code is handled through AntibodyDesignMover
int main(int argc, char* argv[])
{
	try{
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go(new AntibodyDesignMover);
	}catch(utility::excn::EXCN_Base & excn){
		std::cout << "Exception: "<<std::endl;
		excn.show(std::cerr);
	}
	
	return(0);
}
