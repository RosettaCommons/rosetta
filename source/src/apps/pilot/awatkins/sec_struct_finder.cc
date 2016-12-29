// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/ncbb/SecStructFinder.hh>

#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

//Basic headers
#include <basic/basic.hh>
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>

// C++ headers
#include <string>
#include <sstream>
#include <cmath>

// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace func;
using namespace constraints;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;


// tracer - used to replace cout
static basic::Tracer TR( "SecStructFinder" );

// application specific options
namespace ss_finder {
// pert options
StringOptionKey const residue ( "ss_finder::residue" );
IntegerOptionKey const min_length ( "ss_finder::min_length" );
RealOptionKey const bin_size ( "ss_finder::bin_size");
IntegerOptionKey const max_length ( "ss_finder::max_length" );
RealOptionKey const dump_threshold ( "ss_finder::dump_threshold");
StringOptionKey const dihedral_pattern ( "ss_finder::dihedral_pattern" );
StringOptionKey const alpha_beta_pattern ( "ss_finder::alpha_beta_pattern" );
BooleanOptionKey const min_everything ( "ss_finder::min_everything" );
RealOptionKey const dihedral_min ( "ss_finder::dihedral_min");
RealOptionKey const dihedral_max ( "ss_finder::dihedral_max");
BooleanOptionKey const cart ( "ss_finder::cart");
BooleanOptionKey const constrain ( "ss_finder::constrain");
RealOptionKey const dissimilarity ( "ss_finder::dissimilarity");
}

int
main( int argc, char* argv[] )
{
	try {

		option.add( ss_finder::residue, "Residue to form a homopolymer from. Default alanine. Specify by three letter code.").def("ALA");
		option.add( ss_finder::min_length, "Minimum length to detect key SS property/periodicity. Default 5." ).def(5);
		option.add( ss_finder::bin_size, "Dihedral bin size. Default 10." ).def(10);
		option.add( ss_finder::max_length, "Maximum length to detect key SS property/periodicity. Default 5." ).def(0);
		option.add( ss_finder::dump_threshold, "Energy below which to dump. Default -100000, i.e. dump just about nothing." ).def(-100000);
		option.add( ss_finder::dihedral_pattern, "Series of letters indicating distinct dihedral sets, e.g. 'A' or 'AB' or 'AAB'. Default A." ).def("A");
		option.add( ss_finder::alpha_beta_pattern, "Series of letters indicating alpha beta pattern sets, e.g. 'A' or 'B' or 'AB'. Default A." ).def("A");
		option.add( ss_finder::min_everything, "Do we minimize everything first (true) or only the threshold-beaters? Default false" ).def(false);
		option.add( ss_finder::dihedral_min, "Dihedral min, default -180" ).def(-180);
		option.add( ss_finder::dihedral_max, "Dihedral max, default 180" ).def(180);
		option.add( ss_finder::cart, "cart min? default false" ).def(false);
		option.add( ss_finder::constrain, "in minimization, apply dihedral constraints of approximate stdev of bin_size?" ).def(false);
		option.add( ss_finder::dissimilarity, "If A and B are given dihedrals closer than this to each other, they'll be skipped. Default 10 degrees or bin size, whichever greater" ).def(10);


		devel::init(argc, argv);
		Size min_length = option[ ss_finder::min_length ]();
		Size max_length = ( option[ ss_finder::max_length ]() == 0 ) ? min_length : option[ ss_finder::max_length ]();
		Real bin_size = option[ ss_finder::bin_size ]();
		Real dump_threshold = option[ ss_finder::dump_threshold ]();
		std::string residue = option[ ss_finder::residue ]();
		std::string dihedral_pattern = option[ ss_finder::dihedral_pattern ]();
		std::string alpha_beta_pattern = option[ ss_finder::alpha_beta_pattern ]();
		bool min_everything = option[ ss_finder::min_everything ]();
		Real dihedral_min = option[ ss_finder::dihedral_min ]();
		Real dihedral_max = option[ ss_finder::dihedral_max ]();
		bool cart = option[ ss_finder::cart ]();
		bool constrain = option[ ss_finder::constrain ]();
		Real dissimilarity = option[ ss_finder::dissimilarity ]();
		if ( dissimilarity < bin_size ) dissimilarity = bin_size;

		//Pose pose;

		protocols::ncbb::SecStructFinderOP ssf( new protocols::ncbb::SecStructFinder( residue, min_length, max_length, bin_size, dissimilarity, dihedral_min, dihedral_max, dump_threshold, dihedral_pattern, alpha_beta_pattern, min_everything, cart, constrain ) );

		protocols::jd2::JobDistributor::get_instance()->go( ssf );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
