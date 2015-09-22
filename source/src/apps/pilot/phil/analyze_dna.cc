// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers

#include <protocols/viewer/viewers.hh>

#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <utility/excn/Exceptions.hh>
#include <core/scoring/func/Func.hh>


#include <core/pose/Pose.hh>

#include <basic/options/util.hh>


#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>


#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>


// // C++ headers

//silly using/typedef


#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/dna.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/constraints/Constraint.hh>


using basic::T;
using basic::Error;
using basic::Warning;


using namespace core;
//using namespace protocols;

using utility::vector1;
using std::string;


static THREAD_LOCAL basic::Tracer tt( "demo.phil.analyze_dna", basic::t_trace );
static THREAD_LOCAL basic::Tracer td( "demo.phil.analyze_dna", basic::t_debug );
static THREAD_LOCAL basic::Tracer ti( "demo.phil.analyze_dna", basic::t_info );


///////////////////////////////////////////////////////////////////////////////
void
dna_geometry()
{

	vector1< string > files( basic::options::start_files() );
	for ( Size n=1; n<= files.size(); ++n ) {
		pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, files[n] );
		std::cout << "DNA_GEOMETRY: " << files[n] << std::endl;
		scoring::dna::set_base_partner( pose );
		scoring::dna::show_dna_geometry( pose, std::cout );
	}

}


///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string const mode( option[ dna::specificity::mode ].value() );

	if ( mode == "geometry" ) {
		dna_geometry();
		exit(0);
	}

	exit(0); // add new mode strings
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
// 	{ // add any new options
// 		using namespace utility::options;
// 		BooleanOptionKey const myopt = BooleanOptionKey("phil:dof_constraint_weight");
// 		option.add(myopt, "Does nothing, nothing at all");
// 	}

	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
