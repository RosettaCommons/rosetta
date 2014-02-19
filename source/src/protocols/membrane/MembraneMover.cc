// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/MembraneMover.cc
///
/// @brief      Top-Level Unit Test for the Membrane Protein Factory (Mover Class)
/// @details    The purpose of this application is to test the membrane protein factory
///             initialization code including all external dependencies which cannot be tested
///             in JD2, Resource Manager, and the Pose cache. This can also serve as the integration
///             test for memrane protein initialization.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (1/3/14)

#ifndef INCLUDED_protocols_membrane_MembraneMover_cc
#define INCLUDED_protocols_membrane_MembraneMover_cc

// Unit Headers
#include <protocols/membrane/MembraneMover.hh>

// Project Headers
#include <core/membrane/MembraneProteinFactory.hh>
#include <core/membrane/util/Exceptions.hh>
#include <core/membrane/util/definitions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

// Protocol Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("apps.pilot.ralford.MembraneMover");

namespace protocols {
namespace membrane {

/// Constructors /////////////////////////////

/// @brief Standard Constructor
MembraneMover::MembraneMover() :
		protocols::moves::Mover(),
		membrane_pose_(NULL),
		fullatom_(true),
		chains_("")
{}

/// @brief Copy Constructor
MembraneMover::MembraneMover( MembraneMover const & src ) :
		protocols::moves::Mover(src),
		membrane_pose_(NULL),
		fullatom_(true),
		chains_("")
{
	membrane_pose_ = src.membrane_pose_;
}

/// @brief Standard Destructor
MembraneMover::~MembraneMover() {}

//// Options System //////////

/// @brief Register Relevant Options
void MembraneMover::register_options() {

	using namespace basic::options;

	option.add_relevant( OptionKeys::in::file::membrane_chains );
	option.add_relevant( OptionKeys::in::file::fullatom );
}

/// @brief Initialize Options from the Command line
void MembraneMover::init_from_cmd() {

	using namespace basic::options;

	// Check for chains specified
	if ( option[ OptionKeys::in::file::membrane_chains].user() ) {
		chains_ = option[ OptionKeys::in::file::membrane_chains ]();
	}

	// Check for fullatom option
	if ( option[ OptionKeys::in::file::fullatom ].user() ) {
		fullatom_ = option[ OptionKeys::in::file::fullatom ]();
	}
}

/// Mover Inherited Methods ////////////////////////

/// @brief Return your poses
core::pose::PoseOP MembraneMover::get_membrane_pose() { return membrane_pose_; }
core::pose::PoseOP MembraneMover::get_non_membrane_pose() { return non_mp_pose_; }

// Mover Methods
/// @brief Return the name of the Mover.
std::string MembraneMover::get_name() const { return "MembraneMover"; }

/// @brief Clone Method
protocols::moves::MoverOP MembraneMover::clone() const { return new MembraneMover(*this); }

/// @brief Fresh Instance
protocols::moves::MoverOP MembraneMover::fresh_instance() const { return new MembraneMover(); }

///// Non membrane Pose Construction methods ////////

/// @brief Initialize Chains
/// @details Initialize Chains List from Resource Options
std::map< core::Size, std::string >
MembraneMover::initialize_chains() {

	using namespace basic::options;
	using namespace basic::resource_manager;

	// Create a map of chain descriptions and poses
	std::map< core::Size, std::string > chain_descriptions;

	// Grab file and create stream
	std::string infile = option[ OptionKeys::in::file::membrane_chains ]();
	utility::io::izstream stream (infile);

	std::string line;
	std::string desc;
	stream.open(infile);

	if (stream) {

		// Grab the first line
		getline(stream, line);

		core::Size i = 1;
		while ( !stream.eof() ) {

			std::istringstream l(line);
			l >> desc;
			chain_descriptions.insert( std::pair< core::Size, std::string >( i, desc ) );
			getline(stream, line);
			i++;
		}

	} else {
		throw new EXCN_Illegal_Arguments("Cannot open file " + infile );
	}

	return chain_descriptions;
}

/// @brief Make a Non membrane Pose
/// @details Load in the non memrbane counterpart of a membrane pose
core::pose::PoseOP
MembraneMover::initialize_non_membrane_pose( std::map< core::Size, std::string > chain_descriptions ) {

	using namespace basic::resource_manager;

	// Create list of chains
	utility::vector1< core::pose::PoseOP > chains;
	core::Size nchains = chain_descriptions.size();
	chains.resize( nchains );

	// Create a new final pose to return
	core::pose::PoseOP pose = new core::pose::Pose();

	// Read in chains by resource description
	for ( core::Size i = 1; i <= nchains; i++ ) {

		// Prefix Based Resource Tags
		std::string base_desc = chain_descriptions.at(i);

		// Get pose from resource manager
		if ( ! ResourceManager::get_instance()->has_resource_with_description( base_desc ) )
		{
			throw EXCN_Resource_Definition( "Cannot load chain with description " + base_desc );
			TR << "Cannot load any of my resource description things" << std::endl;
		}
		chains[i] = basic::resource_manager::get_resource< core::pose::Pose >( base_desc );
	}

	// Append the poses together
	// Create multi chain pose and write teh chain map
	for ( core::Size i = 1; i <= nchains; i++ ) {

		// Create a given pose from multiple chains
		core::pose::PoseOP chain_pose = chains[i];
		append_pose_to_pose(*pose, *chain_pose, true);
	}

	return pose;
}

/// Mover Apply

/// @brief Apply the corresponding move to the pose
void MembraneMover::apply( core::pose::Pose & )
{
	using namespace core::membrane;
	using namespace basic::resource_manager;
	using namespace core::membrane::util;

	TR << "Loading in a membrane pose" << std::endl;
	// Create an instance of the Membrane Protein Factory
	std::clock_t mp_start;
	double mp_duration;
	mp_start = std::clock();
	MembraneProteinFactoryOP mpf = new MembraneProteinFactory( false, chains_, fullatom_ );
	membrane_pose_ = mpf->create_membrane_pose();
	mp_duration = ( std::clock() - mp_start ) / (double) CLOCKS_PER_SEC;

	TR << "Loading in a non membrane pose" << std::endl;
	// Load in a non membrane pose
	std::clock_t non_mp_start;
	double non_mp_duration;
	non_mp_start = std::clock();
	std::map< core::Size, std::string > chain_map = initialize_chains();
	non_mp_pose_ = initialize_non_membrane_pose( chain_map );
	non_mp_duration = ( std::clock() - non_mp_start ) / (double) CLOCKS_PER_SEC;

	TR << "Printing Results of Timing" << std::endl;
	TR << "Non membrane pose multi chain load time: " << non_mp_duration << std::endl;
	TR << "Membrane pose load time: " << mp_duration << std::endl;
	TR << "Done! - Exiting mover method" << std::endl;
}

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembraneMover_cc

