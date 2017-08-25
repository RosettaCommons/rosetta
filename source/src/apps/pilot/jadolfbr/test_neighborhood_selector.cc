// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/jadolfbr/test_neighborhood_selector.cc
/// @brief Test neighborhood selector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>


// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

#include <utility/string_util.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

static THREAD_LOCAL basic::Tracer TR("test_neighborhood_selector");


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	

	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );

}

core::Size count_subset( utility::vector1< bool > subset ){
	core::Size counts = 0;
	for (bool const & b : subset){
		if (b){
			counts+=1;
		}
	}
	return counts;
}


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::select::residue_selector;
		devel::init( argc, argv );
		register_options();

		
//this won't compile until you fill in brief and default yourself

		core::pose::PoseOP pose = core::import_pose::pose_from_file( "/Users/jadolfbr/2j88.pdb");
		
		utility::vector1< bool > focus( pose->total_residue(), false);
		focus[10] = true;
		
		TR << "Testing: " << utility::to_string(focus) << std::endl;
		
		NeighborhoodResidueSelector selector = NeighborhoodResidueSelector();
		core::scoring::ScoreFunctionOP score = core::scoring::get_score_function();
		score->score(*pose);
		
		selector.set_focus(focus);
		selector.set_include_focus_in_subset( true );
		
		utility::vector1< core::Real > distances;
		distances.push_back(2.0);
		distances.push_back(4.0);
		distances.push_back(6.0);
		distances.push_back(8.0);
		distances.push_back(10.0);
		distances.push_back(12.0);
		

		
		for (core::Real d : distances){
			selector.set_distance(d);
			utility::vector1< bool > subset = selector.apply(*pose);
			TR << "Length: " << subset.size() << std::endl;
			
			core::Size counts = count_subset( subset );
			
			//std::cout << "Distance: " << d << " " << utility::to_string( subset ) << " " << subset << std::endl;
			std::cout << "Distance: " << d << " " << "TotalTrue: " << counts << " " << std::endl;
		}
		
		selector.set_include_focus_in_subset( false );
		TR << "NOT keeping residue in subset! " << std::endl;
		for (core::Real d : distances){
			selector.set_distance(d);
			utility::vector1< bool > subset = selector.apply(*pose);
			TR << "Length: " << subset.size() << std::endl;
			
			core::Size counts = count_subset( subset );
			
			//std::cout << "Distance: " << d << " " << utility::to_string( subset ) << " " << subset << std::endl;
			std::cout << "Distance: " << d << " " << "TotalTrue: " << counts << " " << std::endl;
		}
		
		
		std::cout << "complete" << std::endl;
		
		//none::TestNBRSelectorOP mover_protocol( new none::TestNBRSelector() );

		//protocols::jd2::JobDistributor::get_instance()->go( mover_protocol );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
