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
// AUTO-REMOVED #include <protocols/moves/TrialMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>

#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <utility/excn/Exceptions.hh>
// AUTO-REMOVED #include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/func/Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AngleConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>


// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>


// AUTO-REMOVED #include <core/kinematics/util.hh>
// AUTO-REMOVED #include <core/kinematics/visualize.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>

// AUTO-REMOVED #include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>

// AUTO-REMOVED #include <basic/prof.hh> // profiling

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>



#include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/random/random_permutation.hh>

#include <ObjexxFCL/string.functions.hh>


// // C++ headers
// AUTO-REMOVED #include <fstream>

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

static numeric::random::RandomGenerator RG(54323); // <- Magic number, do not change it!!!

using namespace core;
//using namespace protocols;

using utility::vector1;
using std::string;


static basic::Tracer tt( "demo.phil.analyze_dna", basic::t_trace );
static basic::Tracer td( "demo.phil.analyze_dna", basic::t_debug );
static basic::Tracer ti( "demo.phil.analyze_dna", basic::t_info  );



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
