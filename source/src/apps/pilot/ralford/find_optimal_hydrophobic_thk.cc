// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/ralford/find_optimal_hydrophobic_thk.cc
/// @brief An application for finding the minimum energy hydrophobic thickness of a membrane protein
/// @author Rebecca Alford (rfalford12@gmail.com)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/membrane/AddMembraneMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/jd2/JobDistributor.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>

// basic headers
#include <basic/options/keys/OptionKeys.hh>

#include <utility/string_util.hh> // AUTO IWYU For string_split

static basic::Tracer TR("find_optimal_hydrophobic_thk");

class FindOptimalHydrophobicThk : public protocols::moves::Mover {

public:

	FindOptimalHydrophobicThk() {}
	~FindOptimalHydrophobicThk() {}

	std::string get_name() const override { return "FindOptimalHydrophobicThk"; }

	void apply( core::pose::Pose & pose ) override {

		using namespace protocols::membrane;

		// Initialize the membrane framework
		AddMembraneMoverOP add_memb = AddMembraneMoverOP( new AddMembraneMover() );
		add_memb->apply( pose );

		// Initialize a new scoring function
		core::scoring::ScoreFunctionOP sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "franklin2019" );

		// Create a pack rotamers mover
		using namespace protocols::minimization_packing;
		using namespace core::pack::task;
		PackerTaskOP pack_task( TaskFactory::create_packer_task( pose ) );
		PackRotamersMoverOP pack_mover( new PackRotamersMover( sfxn, pack_task ) );

		// Create a min mover
		using namespace core::kinematics;
		MoveMapOP movemap( new MoveMap() );
		MinMoverOP min_mover = utility::pointer::make_shared< MinMover >();
		min_mover->movemap(movemap);
		min_mover->score_function(sfxn);

		// Apply minimization and packing steps
		pack_mover->apply( pose );
		min_mover->apply( pose );

		// Configure the output file and write the header
		utility::vector1< std::string > temp( utility::string_split( pose.pdb_info()->name(), '/') );
		std::string tempstr = temp[ temp.size() ].substr(0, temp[ temp.size() ].size()-4 );
		std::string filename( tempstr +  "_" + sfxn->get_name() + "_optimal_thickness.dat" );
		utility::io::ozstream output( filename );
		output << "thickness score" << std::endl;

		// For each thickness, calculate the steepness
		core::Real steepness = pose.conformation().membrane_info()->implicit_lipids()->water_steepness();
		for ( core::Real t = 0; t <= 40; t += 0.25 ) {

			// Calculate and set the pseudo-thickness
			core::Real tau( 1 / exp( -steepness*t ) );
			pose.conformation().membrane_info()->implicit_lipids()->water_pseudo_thickness( tau );

			// Score the pose
			core::Real score( sfxn->score( pose ) );

			// Append the current thickness and score to an output file
			output << t << " " << score << std::endl;
		}
	}
};

using FindOptimalHydrophobicThkOP = utility::pointer::shared_ptr<FindOptimalHydrophobicThk>;
using FindOptimalHydrophobicThkCOP = utility::pointer::shared_ptr<const FindOptimalHydrophobicThk>;

///////////////////////////////////////////////////////////////////////////////


int
main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );
		FindOptimalHydrophobicThkOP mover_protocol( new FindOptimalHydrophobicThk() );
		protocols::jd2::JobDistributor::get_instance()->go( mover_protocol );

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
