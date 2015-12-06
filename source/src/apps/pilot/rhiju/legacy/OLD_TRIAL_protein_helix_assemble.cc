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
#include <core/scoring/rms_util.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/RNA_SilentStruct.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>

#include <core/options/option_macros.hh>
#include <protocols/idealize/idealize.hh>
#include <protocols/viewer/viewers.hh>

#include <protocols/farna/setup/RNA_DeNovoPoseSetup.fwd.hh>
#include <protocols/farna/setup/RNA_DeNovoPoseSetup.hh>
#include <protocols/farna/libraries/RNA_ChunkLibrary.hh>
#include <protocols/farna/libraries/RNA_ChunkLibrary.fwd.hh>
#include <protocols/farna/util.hh>
#include <protocols/stepwise/modeler/rna/helix/RNA_HelixAssembler.hh>
#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <protocols/farna/movers/RNA_Minimizer.hh>

//Minimizer stuff
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/util/basic.hh>
#include <core/io/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap

#include <core/scoring/Energies.hh> //for EnergyMap


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( String,  seq )


//////////////////////////////////////////////////////
void
set_helix_torsions( pose::Pose & pose, Size const & n )
{
  pose.set_psi  ( n, -30.0  );
  pose.set_omega( n, 180.0  );
  pose.set_phi  ( n+1, -70.0  );
}

//////////////////////////////////////////////////////
void
build_on_helix( pose::Pose & pose, Size const n, char const & seq1  ){

  using namespace core::chemical;
  using namespace core::conformation;

  static const ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

  ResidueOP rsd1( ResidueFactory::create_residue( *(rsd_set->aa_map( aa_from_oneletter_code( seq1 ) )[1] ) ) );
  pose.append_polymer_residue_after_seqpos(   *rsd1, n - 1, true /*build_ideal_geometry*/ );
  set_helix_torsions( pose, n-1 );

}

//////////////////////////////////////////////////////
void
put_constraints_on_helix( pose::Pose & pose, Size const n ){

  using namespace core::scoring::constraints;
  using namespace core::id;

  Real const helix_distance( 3.0 );
  Real const distance_stddev( 0.25 ); //Hmm. Maybe try linear instead?
  FuncOP const distance_func( new HarmonicFunc( helix_distance, distance_stddev ) );

  if ( n > 3 ) {
    pose.add_constraint( new AtomPairConstraint(
						id::AtomID( id::NamedAtomID( " O  ", n-3), pose ),
						id::AtomID( id::NamedAtomID( " N  ", n  ), pose ),
						distance_func ) );
  }

}

/////////////////////////////////////////////////////////////////////////////
void
minimize_helix( pose::Pose & pose, Size const n ){

  using namespace core::scoring;
  using namespace core::optimization;

  AtomTreeMinimizer minimizer;
  float const dummy_tol( 0.0000025);
  bool const use_nblist( true );
  MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
  options.nblist_auto_update( true );

  kinematics::MoveMap mm;
  //  mm.set_bb( false );
  //  mm.set_chi( false );
  //  mm.set_jump( false );
  //  for (Size i = n-1; i <= n+2; i++ ) {
  mm.set_bb(  true );
  mm.set_chi(  true );
    //  }
  //mm.set_jump( true );

  std::cout << "minimizing" << std::endl;
  ScoreFunctionOP scorefxn =  core::scoring::get_score_function();
  minimizer.run( pose, mm, *scorefxn, options );
  (*scorefxn)( pose );

}

///////////////////////////////////////////////////////////////
void
pack_sidechains( pose::Pose & pose ){

  using namespace core::scoring;
  using namespace core::pack::task;

  PackerTaskOP pack_task_ = pack::task::TaskFactory::create_packer_task( pose );
  pack_task_->restrict_to_repacking();
  for (Size i = 1; i <= pose.total_residue(); i++) {
    if ( !pose.residue(i).is_protein() ) continue;
    pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
    pack_task_->nonconst_residue_task(i).or_ex1( true );
    pack_task_->nonconst_residue_task(i).or_ex2( true );
    pack_task_->nonconst_residue_task(i).or_ex3( true );
    pack_task_->nonconst_residue_task(i).or_ex4( true );
    pack_task_->nonconst_residue_task(i).or_include_current( true );
    if ( pose.residue(i).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) {
    	pack_task_->nonconst_residue_task(i).prevent_repacking();
    }
  }


  ScoreFunctionOP scorefxn_ =  core::scoring::get_score_function();
  pack::rotamer_trials( pose, *scorefxn_, pack_task_ );

}


///////////////////////////////////////////////////////////////
void
protein_helix_test(){

  using namespace core::scoring;
  using namespace core::chemical;
  using namespace core::options;
  using namespace core::options::OptionKeys;
  using namespace core::pose;
  using namespace core::kinematics;
  using namespace core::io::silent;
  using namespace protocols::farna;

  // What is the sequence?
  std::string full_sequence;
  full_sequence = option[ seq ]();

  ResidueTypeSetCAP rsd_set;
  rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

  // Starting two-residues.
  Pose pose;
  make_pose_from_sequence( pose, full_sequence.substr(0,2), *rsd_set );
  remove_upper_terminus_type_from_pose_residue( pose, 2 );

  protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );


  set_helix_torsions( pose, 1 );

  Size const numres = full_sequence.size();
  for ( Size n = 3; n <= numres; n++ ) {

    std::cout << "Building on residue: " << n << std::endl;
    build_on_helix( pose, n, full_sequence[n-1] );
    put_constraints_on_helix( pose, n );
    pack_sidechains( pose );
    minimize_helix( pose, n );

  }

  pose.dump_pdb( full_sequence + ".pdb" );

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

  using namespace core::options;

  protein_helix_test();

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
  NEW_OPT( seq, "Input sequence", "" );
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
