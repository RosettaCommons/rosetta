// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		src/apps/pilot/ralford/membrane_relax.cc
///
/// @brief 		Membrane Protein-Protein Docking Application
/// @details    Docking application for protein partners in membranes. Currently
///				using a naive application which just plugs in scoring components and makes
///				no other modifications
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 6/8/14

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/CreateMembranePoseMover.hh> 
#include <protocols/docking/DockingProtocol.hh>

#include <protocols/jd2/util.hh> 

#include <utility/pointer/owning_ptr.hh> 

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>

using namespace protocols;

using utility::vector1;

/// @brief Membrane Relax Mover: Create a Membrane Pose and Apply Relax Protocol
class MembraneDockingxMover : public moves::Mover {

public: 

	/// @brief Get Mover Name
	std::string get_name() const { return "MembraneDockingMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace protocols::membrane;
		using namespace protocols::moves; 
		using namespace protocols::relax; 
		using namespace core::scoring;

		// Create a membrane pose
		CreateMembranePoseMoverOP cpm = new CreateMembranePoseMover();
		cpm->apply(pose);

		pose.conformation().detect_disulfides();

		if ( !pose.conformation().is_membrane() ) {
			std::cout << "Warning, trying to drown a soluble protein in a sea of lipids" << std::endl;
		}

		// Setup custom membrane energy funcitons
		ScoreFunctionOP mpdocking_low = ScoreFunctionFactory::create_score_function("mpframework_cen_2014");
		ScoreFunctionOP mpdocking_high = ScoreFunctionFactory::create_score_function("mpframework_fa_2014");

		// Setup Relax Base Protocol
		DockingProtocolOP docking = DockingProtocol( 1, false, true, true, mpdocking_low, mpdocking_high ); 
		docking->apply(pose); 

	}
};

typedef utility::pointer::owning_ptr< MembraneDockingMover > MembraneDockingMoverOP; 

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols;
		using namespace scoring;
		using namespace basic::options;

		docking::register_options();
		jd2::register_options();

		devel::init(argc, argv);

		// Create and kick off a new relax mover
		jd2::JobDistributor::get_instance()->go( new MembraneDockingMover() );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
