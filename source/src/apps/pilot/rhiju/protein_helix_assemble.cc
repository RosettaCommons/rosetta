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

#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/scoring/dunbrack/DunbrackRotamer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Ramachandran.hh>
#include <protocols/farna/util.hh>

#include <protocols/viewer/viewers.hh>

//Mmmm.. constraints.
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/RT.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>

//StepWiseProtein!
#include <protocols/stepwise/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/StepWiseProteinPoseSetup.hh>
#include <protocols/stepwise/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <protocols/stepwise/StepWiseDoNothingSampleGenerator.hh>
#include <protocols/stepwise/StepWiseCombineSampleGenerator.hh>
#include <protocols/stepwise/StepWiseIdentitySampleGenerator.hh>
#include <protocols/stepwise/StepWisePoseCombineSampleGenerator.hh>
#include <protocols/stepwise/protein/StepWiseBetaAntiParallelJumpSampleGenerator.hh>
#include <protocols/stepwise/protein/util.hh>
#include <protocols/stepwise/protein/StepWiseProteinFilterer.hh>
#include <protocols/stepwise/protein/StepWiseProteinLoopBridger.hh>
#include <protocols/stepwise/protein/StepWiseProteinPoseMinimizer.hh>
#include <protocols/stepwise/protein/StepWiseProteinScreener.hh>
#include <protocols/stepwise/protein/util.hh>
//#include <protocols/stepwise/protein/StepWiseProteinConnectionSampler.hh>
#include <protocols/stepwise/protein/StepWiseProteinPacker.hh>
#include <protocols/stepwise/protein/StepWiseProteinFragmentSampleGenerator.hh>
#include <protocols/stepwise/protein/StepWiseProteinJumpSampleGenerator.hh>
#include <protocols/stepwise/protein/StepWiseProteinMainChainSampleGenerator.hh>
#include <protocols/stepwise/protein/StepWiseProteinPacker.hh>
#include <protocols/stepwise/protein/MainChainTorsionSet.hh>
#include <protocols/stepwise/InputStreamWithResidueInfo.hh>

//clustering
#include <protocols/cluster/cluster.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/optimizeH.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>
#include <core/options/option_macros.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/pose_stream/ExtendedPoseInputStream.hh>
#include <core/io/pose_stream/PoseInputStream.fwd.hh>
#include <core/io/pose_stream/PDBPoseInputStream.hh>
#include <core/io/pose_stream/SilentFilePoseInputStream.hh>
#include <core/util/datacache/BasicDataCache.hh>
#include <core/util/datacache/CacheableString.hh>
#include <core/util/basic.hh>
#include <core/io/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray3D.hh>
//RNA stuff.
//#include <protocols/farna/RNA_FragmentsClasses.hh>
//#include <protocols/farna/RNA_DeNovoProtocol.hh>
//#include <protocols/farna/RNA_StructureParameters.hh>

//Job dsitributor
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <deque>
#include <vector>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/swa.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/cluster.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;
//typedef std::map< std::string, core::pose::PoseOP > PoseList;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( String, seq )


///////////////////////////////////////////////////////////////
void
build_helix_test(){

  // Read in pdb.
  using namespace core::options;
  using namespace core::options::OptionKeys;
  using namespace core::chemical;
  using namespace core::conformation;
  using namespace core::optimization;
  using namespace core::pose;
  using namespace core::id;
  using namespace core::scoring;
  using namespace core::kinematics;
  using namespace protocols::farna;

  ResidueTypeSetCAP rsd_set;
  rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

  // REad in desired sequence
  std::string const sequence = option[ seq ]();

  // Make pose with the first residue
  Pose pose;
  make_pose_from_sequence( pose, sequence.substr(0,1),  *rsd_set );
  protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );
  remove_upper_terminus_type_from_pose_residue( pose, 1 );


  // Main loop
  for (Size n = 2; n <= sequence.size(); ++n ) {

    // Append one residue
    /////////////////////////////////////
    ResidueOP rsd1( ResidueFactory::create_residue( *(rsd_set->aa_map( aa_from_oneletter_code( sequence[n-1] ) )[1] ) ) );
    pose.append_polymer_residue_after_seqpos( *rsd1, n - 1, true /*build_ideal_geometry*/ );

    pose.set_psi( n-1, -40.0 );
    pose.set_omega( n-1, 180.0 );
    pose.set_phi  ( n, -64.0  );

    // Pack rotamers
    pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
    task->initialize_from_command_line();
    Size const nres( pose.total_residue() );
    for ( Size i = 1; i <= nres; ++i ) {
	task->nonconst_residue_task( i ).restrict_to_repacking();
	task->nonconst_residue_task(i).or_ex1( true );
	task->nonconst_residue_task(i).or_ex2( true );
	task->nonconst_residue_task(i).or_ex3( true );
	task->nonconst_residue_task(i).or_ex4( true );
	task->nonconst_residue_task(i).or_ex1aro( true );
	task->nonconst_residue_task(i).or_ex1aro_sample_level( pack::task::EX_THREE_THIRD_STEP_STDDEVS );
	task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
    }

    ScoreFunctionOP scorefxn = get_score_function();
    pack::rotamer_trials( pose, *scorefxn, task );

    // minimize
    MoveMap mm;
    mm.set_bb( true );
    mm.set_chi( true );
    mm.set_jump( true );


    MinimizerOptions options( "dfpmin_armijo_nonmonotone", 0.01, true /*use_nblist*/ );
    AtomTreeMinimizer minimizer;
    minimizer.run( pose, mm, *scorefxn, options );

  }

  // dump pdb
  pose.dump_pdb( sequence+".pdb" );

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

  using namespace core::options;

  build_helix_test();


  protocols::viewer::clear_conformation_viewers();
  exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	//Uh, options?
	NEW_OPT( seq, "sequence", "" );


	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
