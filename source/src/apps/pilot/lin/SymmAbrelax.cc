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
#include <core/conformation/symmetry/util.hh>

#include <protocols/fibril/fibril_util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/loops/loops_main.hh>

#include <basic/options/util.hh>//option.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

// packing
#include <protocols/moves/Mover.hh>

// Auto-header: duplicate removed #include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>

#include <core/chemical/ChemicalManager.hh>


// Auto-header: duplicate removed #include <core/sequence/util.hh>
// Auto-header: duplicate removed #include <core/sequence/Sequence.hh>
// Auto-header: duplicate removed #include <core/id/SequenceMapping.hh>

#include <core/fragment/FragmentIO.hh>


//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jumping/JumpSetup.fwd.hh>
#include <utility/io/mpistream.hh>
#include <utility/excn/Exceptions.hh>

////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace protocols;
using namespace id;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace optimization;
namespace OK = OptionKeys;
using namespace basic::options::OptionKeys;
using namespace fragment;
using namespace protocols::abinitio;
using namespace kinematics;
using namespace protocols::jumping;

using utility::vector1;
using std::string;
using core::import_pose::pose_from_file;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.lin.symmabrelax" );


void
setup_for_folding( core::pose::Pose &pose, ProtocolOP& prot_ptr ) // stolen from abrelax
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  // ==========================================================================
  ///  --------- fold()-specific setup --------------------------
  // ==========================================================================

  // FRAGMENTS:
  // stores in FragSetOP fragset_large_, fragset_small_;

  //evaluator_ = new MetaPoseEvaluator;

  //probably 9mer fragments
  core::fragment::FragSetOP fragset_large_;

  //probably 3mer fragments
  core::fragment::FragSetOP fragset_small_;

  std::string frag_large_file, frag_small_file;
  if (option[ in::file::fragA ].user()) {
    frag_large_file  = option[ in::file::fragA ]();
  } else {
    frag_large_file  = option[ in::file::frag9 ]();
  }

  if (option[ in::file::fragB ].user()) {
    frag_small_file  = option[ in::file::fragB ]();
  } else {
    frag_small_file  = option[ in::file::frag3 ]();
  }

  fragset_large_ = FragmentIO(
			      option[ OptionKeys::abinitio::number_9mer_frags ](),
			      1,
			      false
			      ).read_data( frag_large_file );

  fragset_small_ = FragmentIO(
			      option[ OptionKeys::abinitio::number_3mer_frags ],
			      1, //nr_copies
			      false
			      ).read_data( frag_small_file );


  // make a MoveMap
  kinematics::MoveMapOP movemap = new kinematics::MoveMap;
  // allow bb moves
  movemap->set_bb( false );

  for ( Size i = 81; i <= 120; ++i ) {
    //score_multiply_vector.push_back(0);
    movemap->set_bb(i,true);
  }

  ClassicAbinitio::register_options();
  /// no constraints ---> ClassicAbinitio
  //tr.Info << "run ClassicAbinitio....." << std::endl;
  prot_ptr = new ClassicAbinitio( fragset_small_, fragset_large_, movemap );

  Protocol& abinitio_protocol( *prot_ptr ); // hide the fact that protocol is a pointer

  /// initialize protocol
  abinitio_protocol.init( pose );
  abinitio_protocol.return_centroid( true );
  //if ( evaluator_->size() ) abinitio_protocol.set_evaluation( evaluator_ );
}

///////////////////////////////////////////////////////////////////////////////
core::kinematics::Jump get_jump( Pose const &native_pose, core::Size const &res1, core::Size const &res2 )
{
  core::kinematics::Jump native_jump;

  // work out the stubID
  chemical::ResidueType const& rt1 ( native_pose.residue_type ( res1 ) );
  chemical::ResidueType const& rt2 ( native_pose.residue_type ( res2 ) );

  id::AtomID a1( rt1.atom_index ("N") , res1 );
  id::AtomID a2( rt1.atom_index ("CA") , res1 );
  id::AtomID a3( rt1.atom_index ("C") , res1 );
  id::StubID down_stub_  ( a1, a2, a3 );

  id::AtomID b1( rt2.atom_index ("N") , res2 );
  id::AtomID b2( rt2.atom_index ("CA") , res2 );
  id::AtomID b3( rt2.atom_index ("C") , res2 );
  id::StubID up_stub_ (b1, b2, b3 );

  Stub native_up_ = native_pose.conformation().atom_tree().stub_from_id( up_stub_ );
  Stub native_down_ = native_pose.conformation().atom_tree().stub_from_id( down_stub_ );

  native_jump.from_stubs( native_up_, native_down_ );

  return native_jump;
}

///////////////////////////////////////////////////////////////////////////////
void
SymmAbRelaxTest()
{
  using namespace conformation::symmetry;
  using namespace chemical;
  using namespace protocols::moves;
  using namespace protocols::simple_moves::symmetry;

  // make extended chain
  Pose extended_pose;
  ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( "score3" ) );
  if( option[ in::file::fasta ].user() ) {
    std::string sequence = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]->sequence();
    core::pose::make_pose_from_sequence( extended_pose, sequence,
			    *( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ))  );
    for ( Size pos = 1; pos <= extended_pose.total_residue(); pos++ ) {
      extended_pose.set_phi( pos, -150 );
      extended_pose.set_psi( pos, 150 );
      extended_pose.set_omega( pos, 180 );
    }
  } else {
    core::import_pose::pose_from_file( extended_pose, start_file() , core::import_pose::PDB_file);
  }

  if ( option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
    protocols::loops::set_secstruct_from_psipred_ss2(extended_pose);
    TR << "Pose Secondary Structure Set by Psipred " << extended_pose.secstruct() << std::endl;
  }

  scorefxn->show( std::cout, extended_pose );
  extended_pose.dump_pdb("extended_pose.pdb");

  // set up the symmetric pose
  Pose start_pose = extended_pose;
  protocols::fibril::make_symmetric_fibril( start_pose );
  TR << "Pose Secondary Structure " << start_pose.secstruct() << std::endl;

  // set up the centroid stuff
  ScoreFunctionOP scorefxn_sym = new core::scoring::symmetry::SymmetricScoreFunction(*scorefxn);
  MoverOP cendock = new protocols::symmetric_docking::SymDockingLowRes(scorefxn_sym);
  //cendock->apply(start_pose);

  scorefxn_sym->show( std::cout, start_pose );
  start_pose.dump_pdb("start.pdb");

  return;

  Pose jump_pose;
  core::import_pose::pose_from_file( jump_pose, "fibril_jump.pdb" , core::import_pose::PDB_file);
  //jump_pose.dump_pdb("fibril_jump_dump.pdb");

  Size jump_res1 = 3;
  Size jump_res2 = 9;

  Jump beta_jump = get_jump( jump_pose, jump_res1, jump_res2 );

  //set_fibril_jumps( beta_jump, jump_res1, jump_res2, start_pose );

  /*
  core::fragment::JumpingFrameOP the_jump = generate_jump_frame( jump_res1, jump_res2, true );
  the_jump->steal( jump_pose );

  the_jump->apply( start_pose );


    kinematics::Stub up, down;
      up = extended_pose.conformation().atom_tree().stub_from_id( up_stub_ );
      down = extended_pose.conformation().atom_tree().stub_from_id( down_stub_ );
   kinematics::RT rt(up, down);
  kinematics::Stub test_down;
  rt.make_jump( native_up_, test_down );

  Real rms( 0.0 );
  for ( Size i=1; i<=3; i++ ) {
          Vector tv = test_down.build_fake_xyz( i );
          Vector nv = native_down_.build_fake_xyz( i );
          Vector d = nv-tv;
          rms += d.length();
  }
  ///rms should be 0.0 between native_pose   and extended pose
  */


  if( false ) {
    // beta_jump
    Pose jump_pose;
    core::import_pose::pose_from_file( jump_pose, start_file() , core::import_pose::PDB_file);
    core::kinematics::Jump beta_jump;

    ProtocolOP prot_ptr;
    setup_for_folding( start_pose, prot_ptr );
    Protocol& abinitio_protocol( *prot_ptr );
    abinitio_protocol.set_fullatom_scorefxn( core::scoring::get_score_function() );
    abinitio_protocol.set_centroid_scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "score3" ));

    for (Size i = 1; i<= 200; i++ ) {
      Pose abinitio_pose = start_pose;
      bool success = true;
      abinitio_protocol.apply( abinitio_pose );

      if( success ) {
	std::string outname;
	std::stringstream moreout;
	moreout << "test." << i << ".pdb";
	outname = moreout.str();
	abinitio_pose.dump_pdb(outname);
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
  using namespace basic::options;

  SymmAbRelaxTest();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {

  // initialize option and random number system
  devel::init( argc, argv );

  protocols::viewer::viewer_main( my_main );

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
