// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//Headers are generally organized by either what they do or where they come from.  This organization is first core library headers, then protocols library, then utility stuff.


// Project Headers
#include <core/pose/Pose.fwd.hh>




#include <core/select/residue_selector/ResidueIndexSelector.hh>


#include <protocols/jd2/JobDistributor.hh>

// Mover headers
#include <protocols/ncbb/oop/OopCreatorMover.hh>


//Basic headers

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>

#include <basic/options/keys/OptionKeys.hh> // AUTO IWYU For IntegerVectorOptionKey, BooleanOptionKey

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::ncbb::oop;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

//kdrew: this app adds oop patches to the given pdb strucure

// application specific options
namespace oop_creator {
// pert options
IntegerVectorOptionKey const oop_plus_positions( "oop_creator::oop_plus_positions" );
IntegerVectorOptionKey const oop_minus_positions( "oop_creator::oop_minus_positions" );
IntegerVectorOptionKey const oop_d_plus_positions( "oop_creator::oop_d_plus_positions" );
IntegerVectorOptionKey const oop_d_minus_positions( "oop_creator::oop_d_minus_positions" );
IntegerVectorOptionKey const oop_low_e_puck_positions( "oop_creator::oop_low_e_puck_positions" );
IntegerOptionKey const prepend_n_residues( "oop_creator::prepend_n_residues" );
IntegerOptionKey const append_n_residues( "oop_creator::append_n_residues" );
BooleanOptionKey const final_repack( "oop_creator::final_repack" );
BooleanOptionKey const final_minimize( "oop_creator::final_minimize" );
BooleanOptionKey const final_mc ( "oop_creator::final_mc" );
BooleanOptionKey const correct_oop_post ( "oop_creator::correct_oop_post" );

}

int
main( int argc, char* argv[] )
{
	try {
		utility::vector1< Size > empty_vector(0);
		option.add( oop_creator::oop_plus_positions, "The positions of the first residues of plus oop rings" ).def( empty_vector );
		option.add( oop_creator::oop_minus_positions, "The positions of the first residues of minus oop rings" ).def( empty_vector );
		option.add( oop_creator::oop_d_plus_positions, "The positions of the first residues of chiral d plus oop rings" ).def( empty_vector );
		option.add( oop_creator::oop_d_minus_positions, "The positions of the first residues of chiral d minus oop rings" ).def( empty_vector );
		option.add( oop_creator::oop_low_e_puck_positions, "The positions of the first oop residues, pucker will change to low e conformation" ).def( empty_vector );
		option.add( oop_creator::prepend_n_residues, "Number of residues to prepend" ).def( 0 );
		option.add( oop_creator::append_n_residues, "Number of residues to append" ).def( 0 );
		option.add( oop_creator::final_repack, "Do a final repack. Default false" ).def(false);
		option.add( oop_creator::final_minimize, "Do a final minimization. Default false" ).def(false);
		option.add( oop_creator::final_mc, "Do a final monte carlo on oop. Default false" ).def(false);
		option.add( oop_creator::correct_oop_post, "Correct oop post phi/psi to low energy well. Default false" ).def(false);

		// init command line options

		devel::init(argc, argv);

		using namespace core::select::residue_selector;

		//create mover instance
		utility::vector1< core::Size > oop_plus_positions( option[oop_creator::oop_plus_positions].value() );
		utility::vector1< core::Size > oop_minus_positions( option[oop_creator::oop_minus_positions].value() );
		utility::vector1< core::Size > oop_d_plus_positions( option[oop_creator::oop_d_plus_positions].value() );
		utility::vector1< core::Size > oop_d_minus_positions( option[oop_creator::oop_d_minus_positions].value() );
		utility::vector1< core::Size > oop_low_e_puck_positions( option[oop_creator::oop_low_e_puck_positions].value() );
		ResidueSelectorOP oop_plus_selector( nullptr );
		if ( ! oop_plus_positions.empty() ) {
			oop_plus_selector = utility::pointer::make_shared< ResidueIndexSelector >( oop_plus_positions );
		}
		ResidueSelectorOP oop_minus_selector( nullptr );
		if ( ! oop_minus_positions.empty() ) {
			oop_minus_selector = utility::pointer::make_shared< ResidueIndexSelector >( oop_minus_positions );
		}
		ResidueSelectorOP oop_d_plus_selector( nullptr );
		if ( ! oop_d_plus_positions.empty() ) {
			oop_d_plus_selector = utility::pointer::make_shared< ResidueIndexSelector >( oop_d_plus_positions );
		}
		ResidueSelectorOP oop_d_minus_selector( nullptr );
		if ( ! oop_d_minus_positions.empty() ) {
			oop_d_minus_selector = utility::pointer::make_shared< ResidueIndexSelector >( oop_d_minus_positions );
		}
		ResidueSelectorOP oop_low_e_puck_selector( nullptr );
		if ( ! oop_low_e_puck_positions.empty() ) {
			oop_low_e_puck_selector = utility::pointer::make_shared< ResidueIndexSelector >( oop_low_e_puck_positions );
		}



		OopCreatorMoverOP OC_mover( new OopCreatorMover(
			oop_plus_selector,
			oop_minus_selector,
			oop_d_plus_selector,
			oop_d_minus_selector,
			oop_low_e_puck_selector,
			option[oop_creator::prepend_n_residues].value(),
			option[oop_creator::append_n_residues].value(),
			option[oop_creator::final_repack].value(),
			option[oop_creator::final_minimize].value(),
			option[oop_creator::final_mc].value(),
			option[oop_creator::correct_oop_post].value()) );

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( OC_mover );

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}//main

