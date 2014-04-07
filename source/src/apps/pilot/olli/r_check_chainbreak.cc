// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///
/// @author Oliver Lange

// keep these headers first for compilation with Visual Studio C++
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>


// Unit Headers
#include <protocols/abinitio/AbrelaxApplication.hh>


// Package Headers
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/abinitio/JumpingFoldConstraints.hh>
#include <protocols/relax_protocols.hh>

#include <protocols/jumping/SheetBuilder.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <protocols/jumping/SecondaryStructure.hh>
#include <utility/excn/Exceptions.hh>

// Project Headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>


#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/util.hh>

#include <core/chemical/ChemicalManager.hh>


//#include <core/conformation/ResidueFactory.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/rms_util.hh>

#include <devel/init.hh>

#include <core/sequence/util.hh>

#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <protocols/evaluation/JumpEvaluator.hh>
#include <protocols/evaluation/PCA.hh>

#include <protocols/jumping/JumpSetup.hh>

//#include <protocols/simple_moves/MinMover.hh>

//numeric headers
#include <numeric/random/random.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

//Hacky headers
#include <core/scoring/constraints/BoundConstraint.hh> //needed for IO registration ... grr


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>





namespace protocols {
namespace abinitio {

static basic::Tracer tr("main");
static numeric::random::RandomGenerator RG(82397823);  // <- Magic number, do not change it!
using namespace core;
class Application {
public:
  Application();
  static void register_options();
  void setup();
  void fix_chainbreaks( pose::Pose &pose );
  void copy_native_structure( core::pose::Pose &extended_pose );
  void run();
private:
  pose::PoseOP native_pose_;
};
}
}
//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( File, fold_tree )
///@details registering of options that are relevant for Application
void protocols::abinitio::Application::register_options() {
    OPT( in::file::native );
  NEW_OPT( fold_tree, "a fold-tree","fold_tree.dat");
}


namespace protocols {
namespace abinitio {

using core::Size;
using namespace core;
using namespace protocols;
using namespace fragment;
using namespace abinitio;
using namespace jumping;
using namespace evaluation;
using namespace basic::options;
//using namespace basic::options::OptionKeys;



///@detail c'stor - nothing special
Application::Application() :
  native_pose_( NULL )
{}


///@details setup of Application data that is used for both, fold() and run()
/// this is mainly stuff for scoring and evaluation ( process_decoys(), evaluator_ )
void Application::setup() {
  using namespace basic::options::OptionKeys;

  // read native pose
  if ( option[ in::file::native ].user() ) {
    native_pose_ = new pose::Pose;
    core::import_pose::pose_from_pdb( *native_pose_, option[ in::file::native ]() );
    core::util::switch_to_residue_type_set( *native_pose_, chemical::CENTROID ); //so that in do_rerun the native pose is the same as the other poses
    pose::set_ss_from_phipsi( *native_pose_ );
  }
}


///@detail called by setup_fold() if option[ start_native ] is active
/// the routine defines a fragment of the length of the structure
/// steals the fragment from the native and applies it to the decoy
/// native needs to be idealized!
void Application::copy_native_structure( core::pose::Pose &extended_pose ) {
  // requires that the sequences match at the beginning (1..nmatch_res) -- > use sequence alignment later
  tr.Info << " *** use native structure as starting template -- NEEDS TO BE IDEALIZED !!! *** \n";
  // determine length of segment to copy from native
  Size seg_len = std::min(extended_pose.total_residue(), native_pose_->total_residue() );
  fragment::Frame long_frame(1, seg_len);

	//create apropriate length FragData object
  FragData frag; //there should be some kind of factory to do this.
  Size nbb ( 3 ); //3 backbone torsions to steal
  for ( Size pos = 1; pos<= seg_len; pos++ ) {
    frag.add_residue( new BBTorsionSRFD( nbb, native_pose_->secstruct( pos ), 'X' ) );
  };
  // get torsion angles from native pose
  frag.steal( *native_pose_, long_frame);

  // apply native torsions to extended structue
  frag.apply(extended_pose, long_frame);
}


void Application::run() {
  using namespace basic::options::OptionKeys;
  using namespace pose;
  using namespace jumping;
  using namespace scoring;
  setup(); //read native

  utility::io::izstream file( option[ fold_tree ]() );
  core::kinematics::FoldTree f;
  file >> f;
  JumpSample jump_setup = JumpSample( f );

  Pose pose = *native_pose_;
  pose.fold_tree( f );
  jump_setup.add_chainbreaks( pose );

  ScoreFunction score;
  score.set_weight( scoring::chainbreak, 1.0 );

  tr.Info << "just a pose with cuts" << std::endl;
  score.show( std::cout, pose );
  pose.dump_pdb("pose1.pdb");

  Size seg_len = native_pose_->total_residue();
  fragment::Frame long_frame(1, seg_len);

  //create apropriate length FragData object
  FragData frag; //there should be some kind of factory to do this.
  Size nbb ( 3 ); //3 backbone torsions to steal
  for ( Size pos = 1; pos<= seg_len; pos++ ) {
    frag.add_residue( new BBTorsionSRFD( nbb, native_pose_->secstruct( pos ), 'X' ) );
  };
  // get torsion angles from native pose
  frag.steal( *native_pose_, long_frame);


  Pose pose2 = *native_pose_;
  pose2.fold_tree( f );

 // make  chain
  for ( Size pos = 1; pos <= pose2.total_residue(); pos++ ) {
    pose2.set_phi( pos, -150 );
    pose2.set_psi( pos, 150);
    pose2.set_omega( pos, 180 );
  }

  jump_setup.add_chainbreaks( pose2 );

  tr.Info << "just a extended-pose" << std::endl;
  score.show( std::cout, pose2 );
  pose2.dump_pdb("pose2.pdb");

  frag.apply( pose2, long_frame );

  tr.Info << " after re-application of frags" << std::endl;
  score.show( std::cout, pose2 );
  pose2.dump_pdb("pose3.pdb");
}

 }
}

using namespace protocols::abinitio;
int main( int argc, char** argv ) {
	try{
  Application::register_options();
  devel::init( argc, argv );
  basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "~/minirosetta_database");

  Application app;
  app.run();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
  return 0;
}

