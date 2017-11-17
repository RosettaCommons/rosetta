// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief app to clean PDBs - the input is guarunteed to be Rosetta-readable and indexed from 1 (great for writing loopsfiles, etc)
/// @author Steven Lewis

// Unit Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/moves/Mover.hh>

// Utility Headers
#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>

class ForgetMover : public protocols::moves::Mover {
public:
	ForgetMover() : Mover() {}

	virtual
	void
	apply(core::pose::Pose & pose ){
		pose.pdb_info(NULL); //NUKE IT FROM ORBIT.  IT'S THE ONLY WAY TO BE SURE
		return;
	}

	virtual
	std::string
	get_name() const { return "ForgetMover"; }

};

int
main( int argc, char* argv[] )
{

	try {

	using basic::options::option;
	devel::init(argc, argv);

	protocols::jd2::JobDistributor::get_instance()->go(new ForgetMover);

	basic::T("done") << "************************d**o**n**e**************************************" << std::endl;

	return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
