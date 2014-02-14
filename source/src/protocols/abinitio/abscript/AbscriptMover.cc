// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/AbscriptMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/AbscriptMover.hh>

// Package headers
#include <protocols/abinitio/abscript/AbscriptMoverCreator.hh>
#include <protocols/abinitio/abscript/FragmentCM.hh>
#include <protocols/abinitio/abscript/FragmentJumpCM.hh>
#include <protocols/abinitio/abscript/AbscriptStageMover.hh>
#include <protocols/abinitio/abscript/StagePreparer.hh>
#include <protocols/abinitio/abscript/AbscriptLoopCloserCM.hh>

#include <protocols/environment/claims/EnvClaim.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>

#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>

#include <core/pose/Pose.hh>

#include <protocols/moves/MonteCarlo.hh>

#include <protocols/simple_moves/SmoothFragmentMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/GunnCost.hh>

#include <protocols/canonical_sampling/mc_convergence_checks/util.hh>

#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/JumpSample.hh>

#include <protocols/checkpoint/CheckPointer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>

// Basic Headers
#include <basic/prof.hh>
#include <basic/Tracer.hh>

//option includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/fold_cst.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.movers.AbscriptMover", basic::t_info);

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace protocols::environment;
using namespace protocols::environment::claims;

// If the temperature needs to be altered, the obvious way to do it is through
// a command line option, which could/will replace this. For some reason, though
// in the FragmentSampler it's not set that way...
core::Real const DEFAULT_TEMP = 2.0;

// I would really like to define this in StageID.hh, but that gives a linker error
// and I'm not sure how to fix it without adding a .cc file specifically for this.

StageID& increment_stageid( StageID& id ) {
  id = static_cast< StageID >( static_cast<int>( id ) + 1 );
  assert( !( id > END ) );
  return id;
}


class AbscriptMover::StageTracker {

public:
  StageTracker( moves::MonteCarloOP mc );

  void begin_stage( std::string const& stagename, core::pose::Pose& );

  void do_step( core::pose::Pose&, AbscriptStageMover&, core::Real const& progress );

  void end_stage( core::pose::Pose& pose );

private:
  clock_t starttime_;
  std::string stagename_;
  moves::MonteCarloOP mc_;
  core::Size iterations_;
  checkpoint::CheckPointer checkpointer_;
};


// creator
std::string
AbscriptMoverCreator::keyname() const {
  return AbscriptMoverCreator::mover_name();
}

protocols::moves::MoverOP
AbscriptMoverCreator::create_mover() const {
  return new AbscriptMover;
}

std::string
AbscriptMoverCreator::mover_name() {
  return "AbscriptMover";
}

core::scoring::ScoreFunctionOP setup_score( std::string const& scorename,
                                            basic::options::FileVectorOptionKey optkey ) {
  using namespace core::scoring;
  using namespace basic::options;

  core::scoring::ScoreFunctionOP scorefn;

  if ( option[ optkey ].user() ) {
    scorefn = ScoreFunctionFactory::create_score_function( scorename, option[ optkey ]() );
  } else {
    scorefn = ScoreFunctionFactory::create_score_function( scorename );
  }

  scorefn->set_weight( core::scoring::atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ] );

  assert( scorefn );
  return scorefn;
}

AbscriptMover::AbscriptMover():
  ClaimingMover(),
  id_map_(),
  stage_movers_(),
  mc_()
{
  using namespace basic::options::OptionKeys;
  using namespace basic::options;

  //I don't love this. This should be const static data.
  id_map_["I"] = I; id_map_["II"] = II;
  id_map_["IIIa"] = IIIa; id_map_["IIIb"] = IIIb;
  id_map_["IVa"] = IIIa; id_map_["IVb"] = IIIb;

  mc_ = new moves::MonteCarlo( *setup_score( "score3", OptionKeys::abinitio::stage4_patch ), DEFAULT_TEMP );
  canonical_sampling::mc_convergence_checks::setup_convergence_checks_from_cmdline( *mc_ );

  core::scoring::ScoreFunctionOP tmp_score = setup_score( "score0", OptionKeys::abinitio::stage1_patch );
  tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage1 ]() );
  stage_movers_[ I ] = new AbscriptStageMover( I, mc_, tmp_score, 2000 );

  tmp_score = setup_score( "score1", OptionKeys::abinitio::stage2_patch );
  tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage2 ]() );
  stage_movers_[ II ] = new AbscriptStageMover( II, mc_, tmp_score, 2000 );

  tmp_score = setup_score( "score2", OptionKeys::abinitio::stage3a_patch );
  tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage3 ]() );
  stage_movers_[ IIIa ] = new AbscriptStageMover( IIIa, mc_, tmp_score, 2000 );

  tmp_score = setup_score( "score5", OptionKeys::abinitio::stage3b_patch );
  tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage3 ]() );
  stage_movers_[ IIIb ] = new AbscriptStageMover( IIIb, mc_, tmp_score, 2000 );

  tmp_score = setup_score( "score3", OptionKeys::abinitio::stage4_patch );
  tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage4 ]() );
  stage_movers_[ IVa ] = new AbscriptStageMover( IVa, mc_, tmp_score, 4000 );

  tmp_score = setup_score( "score3", OptionKeys::abinitio::stage4_patch );
  tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage4 ]() );
  stage_movers_[ IVb ] = new AbscriptStageMover( IVb, mc_, tmp_score, 4000 );

  core::Real seqsep_stage1 = 0.15;
  core::Real seqsep_stage3 = 1.0; //have basically all but the furthest chainbreaks active at end of stage3.
  core::Real seqsep_stage4 = 1.2; //make sure that all chainbreaks are switched on in the end

  if ( option[ OptionKeys::fold_cst::seq_sep_stages ].user() ) {
    if ( option[ OptionKeys::fold_cst::seq_sep_stages ]().size() != 3 ) {
      throw utility::excn::EXCN_BadInput("Option seq_sep_stages requires exactly 3 values.");
    }
    seqsep_stage1 = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 1 ]; //default 15
    seqsep_stage3 = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 2 ]; //default 15
    seqsep_stage4 = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 3 ]; //default 15
  }

  stage_movers_[ I    ]->set_seq_sep_ramping( 0.0, seqsep_stage1 );
  stage_movers_[ II   ]->set_seq_sep_ramping( 0.0, seqsep_stage1 );
  stage_movers_[ IIIa ]->set_seq_sep_ramping( seqsep_stage3 - seqsep_stage1 , seqsep_stage1 );
  stage_movers_[ IIIb ]->set_seq_sep_ramping( seqsep_stage3 - seqsep_stage1 , seqsep_stage1 );
  stage_movers_[ IVa  ]->set_seq_sep_ramping( seqsep_stage4 - seqsep_stage3 , seqsep_stage3 );
  stage_movers_[ IVb  ]->set_seq_sep_ramping( seqsep_stage4 - seqsep_stage3 , seqsep_stage3 );

  stage_movers_[ I    ]->set_chainbreak_ramping( 0.0       , 0.00       );
  stage_movers_[ II   ]->set_chainbreak_ramping( 0.0       , 0.25 / 3.0 );
  stage_movers_[ IIIa ]->set_chainbreak_ramping( 2.5 / 3.0 , 0.00       );
  stage_movers_[ IIIb ]->set_chainbreak_ramping( 0.5 / 3.0 , 0.00       );
  stage_movers_[ IVa  ]->set_chainbreak_ramping( 1.5 / 3.0 , 2.50 / 3.0 );
  stage_movers_[ IVb  ]->set_chainbreak_ramping( 1.5 / 3.0 , 2.50 / 3.0 );
}

AbscriptMover::AbscriptMover( AbscriptMover const& src ):
  ClaimingMover( src ),
  id_map_( src.id_map_ ),
  stage_movers_( src.stage_movers_ ),
  mc_( src.mc_ )
{}

moves::MoverOP AbscriptMover::fresh_instance() const {
  return new AbscriptMover();
}

moves::MoverOP AbscriptMover::clone() const {
  return new AbscriptMover( *this );
}

void AbscriptMover::apply( core::pose::Pose& pose ){

  using namespace basic::options;

  if( ! pose.is_centroid() ){
    throw utility::excn::EXCN_BadInput( "AbscriptMover recieved a non-centroid pose--only centroid poses are accepted." );
  }
  StageTracker tracker( mc_ );

  std::map< core::Size, core::Size > ITERATIONS = calculate_iterations( pose );

  // Stage I ---------------------------------------------------------------------
  mc_->clear_poses();
  mc_->reset( pose );

  tracker.begin_stage( "I", pose );
  for( core::Size i = 1; i <= ITERATIONS[1]; ++i ){
    core::Real progress = core::Real( i - 1 )/ITERATIONS[1];
    tracker.do_step( pose,
                     *stage_movers_[I],
                     progress );
  }
  tracker.end_stage( pose );

  // Stage II --------------------------------------------------------------------
  tracker.begin_stage( "II", pose );
  tracker.do_step( pose, *stage_movers_[II], 0.0 );
  tracker.end_stage( pose );

  // Stage III -------------------------------------------------------------------
  tracker.begin_stage( "III", pose );
  for( core::Size i = 1; i <= ITERATIONS[3]; ++i ){
    // sequence is IIIb, IIIa, IIIb, IIIa, IIIb, IIIa, IIIb, IIIa, IIIa, IIIa
    core::Real progress = core::Real( i - 1 )/ITERATIONS[3];
    tracker.do_step( pose,
                     *stage_movers_[ ( ( i % 2 ) == 0 || ( i == 9 ) ) ? IIIa : IIIb ],
                     progress );
  }
  tracker.end_stage( pose );

  // Stage IV ---------------------------------------------------------------------
  tracker.begin_stage( "IV", pose );
  for( core::Size i = 1; i <= ITERATIONS[4]; ++i ){
    // sequence is IVa, IVb, IVb
    core::Real progress = core::Real( i - 1 )/ITERATIONS[4];
    tracker.do_step( pose,
                     *stage_movers_[ ( i == 1 ) ? IVa : IVb ],
                     progress );
  }
  tracker.end_stage( pose );

  mc_->recover_low( pose );
}

void AbscriptMover::yield_submovers( MoverSet& movers_out ) const {
  for( std::map< StageID, AbscriptStageMoverOP >::const_iterator it = stage_movers_.begin();
       it != stage_movers_.end(); ++it ) {
    it->second->yield_submovers( movers_out );
  }
  movers_out.insert( closer_ );
}

void
AbscriptMover::parse_my_tag(
  utility::tag::TagCOP const tag,
  basic::datacache::DataMap & ,
  protocols::filters::Filters_map const & ,
  protocols::moves::Movers_map const & movers,
  core::pose::Pose const & ) {

  using namespace basic::options;

  //cycles control
  core::Real cycles_prefactor = 1.0;

  if( tag->hasOption( "cycles" ) ){
    cycles_prefactor = tag->getOption< core::Real >( "cycles", 1.0 );
    if( option[ OptionKeys::abinitio::increase_cycles ].user() ){
      throw utility::excn::EXCN_RosettaScriptsOption("Flag increase_cycles and AbscriptMover 'cycles' option are incompatible.");
    }
  } else if ( option[ OptionKeys::abinitio::increase_cycles ].user() ){
    cycles_prefactor = option[ OptionKeys::abinitio::increase_cycles ];
  }

  for( std::map< StageID, AbscriptStageMoverOP >::const_iterator it = stage_movers_.begin();
      it != stage_movers_.end(); ++it ){
    it->second->set_cycles_adjust( cycles_prefactor );
  }


  //standard fragments load/don't load
  if( tag->getOption<bool>( "std_frags", true ) ){
    add_default_frags( tag->getOption<std::string>("small_frags", option[ OptionKeys::in::file::frag3 ] ),
                   tag->getOption<std::string>("large_frags", option[ OptionKeys::in::file::frag9 ] ) );
  }

  //Load submovers using subtags.
  typedef utility::vector0< utility::tag::TagCOP > TagVector;
  TagVector const& subtags = tag->getTags();

  for( TagVector::const_iterator subtag_it = subtags.begin();
       subtag_it != subtags.end(); ++subtag_it ) {
    TagCOP stagetag = *subtag_it;

    if( stagetag->getName() == "Stage" ){
      std::string const& id_str = stagetag->getOption<std::string>("ids");
      StageIDs stages = parse_stage_id( id_str );

      TagVector const& movertags = stagetag->getTags();

      for( TagVector::const_iterator movertag_it = movertags.begin();
           movertag_it != movertags.end(); ++movertag_it ){
        TagCOP movertag = *movertag_it;

        if( movertag->getName() == "Mover" ||
            movertag->getName() == "Preparer" ){
          std::string name = movertag->getOption<std::string>( "name", "null" );
          if ( movers.find( name ) == movers.end() ) {
            tr.Error << "Mover not found for XML tag:\n" << movertag << std::endl;
            throw utility::excn::EXCN_RosettaScriptsOption("Mover not found for XML tag.");
          }

          if( movertag->getName() == "Mover" ){
            core::Real weight = movertag->getOption< core::Real >( "weight", 1.0 );
            register_submover( movers.find( name )->second, stages, weight );
          } else if ( movertag->getName() == "Preparer" ) {
            register_preparer( movers.find( name )->second, stages );
          } else { assert( false ); } //no other options
        }
      }
    } else {
      tr.Error << "AbscriptMover recieved illegal tag (Stage is acceptable):\n" << stagetag << std::endl;
      throw utility::excn::EXCN_RosettaScriptsOption("Illegal AbscriptMover subtag.");
    }
  }

}

void AbscriptMover::add_default_frags( std::string const& small_fragfile, std::string const& large_fragfile ){
  using namespace core::fragment;
  using namespace basic::options;
  using namespace protocols::simple_moves;

  FragmentIO small_frag_io( option[ OptionKeys::abinitio::number_3mer_frags ](), 1,
                            option[ OptionKeys::frags::annotate ]() );

  FragmentIO big_frag_io( option[ OptionKeys::abinitio::number_9mer_frags ](), 1,
                          option[ OptionKeys::frags::annotate ]() );

  FragSetOP frags_small = small_frag_io.read_data( small_fragfile );
  FragSetOP frags_large = big_frag_io.read_data( large_fragfile );

  //embed standard movers in movers that will claim and unlock for them.
  FragmentCMOP claim_large  = new FragmentCM( new ClassicFragmentMover( frags_large ), "BASE" );
  FragmentCMOP claim_small  = new FragmentCM( new ClassicFragmentMover( frags_small ), "BASE" );
  FragmentCMOP claim_smooth = new FragmentCM( new SmoothFragmentMover( frags_small, new GunnCost() ),
                                              "BASE" );

  //apply to appropriate stages
  for( StageID id = I; id <= IIIb; increment_stageid(id) ){
    stage_movers_[id]->add_submover( claim_large, 1.0 );
  }

  stage_movers_[ IVa ]->add_submover( claim_small, 1.0 );
  stage_movers_[ IVb ]->add_submover( claim_smooth, 1.0 );

  closer_ = new AbscriptLoopCloserCM( frags_small,
                                      stage_movers_[IVa]->scorefxn()->clone() );
}

std::string AbscriptMover::get_name() const {
  return "AbscriptMover";
}

StageIDs AbscriptMover::parse_stage_id( std::string const& id_str ) const {
  typedef utility::vector0< std::string > Strings;
  Strings comma_delin = utility::string_split( id_str, ',' );

  StageIDs ids;

  for( Strings::const_iterator range = comma_delin.begin();
       range != comma_delin.end(); ++range ){
    Strings hyphen_delin = utility::string_split( *range, '-' );
    if( hyphen_delin.size() == 1 ){
      ids.push_back( id_map_.at( hyphen_delin[0] ) );
    } else if( hyphen_delin.size() == 2 ){
      for( StageID id = id_map_.at(hyphen_delin[0]); id <= id_map_.at( hyphen_delin[1] ); increment_stageid(id) ){
        ids.push_back( id );
      }
    } else {
      throw utility::excn::EXCN_RosettaScriptsOption("Invalid stage ID: " + *range );
    }
  }

  return ids;
}

void AbscriptMover::register_submover( protocols::moves::MoverOP mover_in,
                                       StageIDs const& ids,
                                       core::Real weight ){
  ClaimingMoverOP mover = dynamic_cast< ClaimingMover* >( mover_in.get() );
  if( !mover ){
    tr.Error << "The mover " << mover_in->get_name()
             << " cannot be attached to the Abscript mover because it is not a ClaimingMover "
             << std::endl;
    throw utility::excn::EXCN_RosettaScriptsOption("Abscript mover takes only ClaimingMovers.");
  }

  for( StageIDs::const_iterator id = ids.begin(); id != ids.end(); ++id ){
    stage_movers_.at( *id )->add_submover( mover, weight );
  }
}

std::map< core::Size, core::Size > AbscriptMover::calculate_iterations( core::pose::Pose const& pose ){

  std::map< core::Size, core::Size > ITERATIONS;

  //Stage1 has FT-dependent seqsep ramping ---------------------------------------
  core::kinematics::ShortestPathInFoldTree shortestpath( pose.fold_tree() );
  assert( shortestpath.max_dist() <= pose.total_residue() );

  core::Size end_seqsep = stage_movers_[ I ]->seq_sep_intercept() * shortestpath.max_dist();
  core::Size begin_seqsep = std::min( end_seqsep, Size(3) );

  core::Real end_seqsep_factor   = core::Real( end_seqsep   ) / shortestpath.max_dist();
  core::Real begin_seqsep_factor = core::Real( begin_seqsep ) / shortestpath.max_dist();

  stage_movers_[ I ]->set_seq_sep_ramping( 0.0, begin_seqsep_factor );

  // Calculate number of interations for each phase ------------------------------
  ITERATIONS[ 1 ] = 1;
  if( pose.constraint_set()->has_residue_pair_constraints() ){
    stage_movers_[ I ]->set_seq_sep_ramping( end_seqsep_factor - begin_seqsep_factor, begin_seqsep_factor );
    ITERATIONS[ 1 ] += (end_seqsep - begin_seqsep) / 2;
  }

  ITERATIONS[ 2 ] = 1;
  ITERATIONS[ 3 ] = 10;
  ITERATIONS[ 4 ] = 3;

  return ITERATIONS;
}

EnvClaims AbscriptMover::yield_claims( core::pose::Pose& ) {
  EnvClaims claims; return claims;
}

void AbscriptMover::register_preparer( protocols::moves::MoverOP mover, StageIDs const& ids ){
  StagePreparerOP preparer = dynamic_cast< StagePreparer* >( mover.get() );
  if( !preparer ){
    tr.Error << "The mover '" << mover->get_name() << "' is not a preparer." << std::endl;
    throw utility::excn::EXCN_RosettaScriptsOption("Non-preparer mover added as preparer to Abscript.");
  }

  for( StageIDs::const_iterator id = ids.begin(); id != ids.end(); ++id ){
    stage_movers_.at( *id )->add_preparer( preparer );
  }

}

AbscriptMover::StageTracker::StageTracker( moves::MonteCarloOP mc ):
mc_( mc ),
checkpointer_( "Abscript" )
{}

void AbscriptMover::StageTracker::begin_stage( std::string const& stagename,
                                              core::pose::Pose& pose ){
  stagename_ = stagename;
  iterations_ = 0;

  if( stagename_ == "I" &&
     basic::options::option[ basic::options::OptionKeys::run::intermediate_structures ]() ){
    jd2::output_intermediate_pose( pose, "stage0" );
  }

  // part X ----------------------------------------
  tr.Info <<  "\n===================================================================\n";
  tr.Info <<  "   Stage " << stagename << "\n";
  tr.Info <<  "--------------------------------------------------------------------" << std::endl;
  // PROF_START( stagename() );
  starttime_ = clock();
}

void AbscriptMover::StageTracker::do_step( core::pose::Pose& pose,
                                           AbscriptStageMover& mover,
                                           core::Real const& progress ){
  iterations_ += 1;

  std::string checkpoint_id = "stage_"+stagename_+"_it_"+utility::to_string( iterations_ );

  //  Checkpointing goes here, but doesn't work with ProtectedConformation yet, and so isn't implemented.
  //  if ( !checkpointer_.recover_checkpoint( pose, jd2::current_output_name(), checkpoint_id,
  //                                          false /* fullatom */, true /* fold_tree */ ) ) {
  if( mover.setup_stage( pose, progress ) ){
    mover.apply( pose );
  } else {
    tr.Info << "Stage " << stagename_ << " step at " << progress << "% skipped." << std::endl;
  }
  //    checkpointer_.checkpoint( pose, jd2::current_output_name(), checkpoint_id, true /*fold tree */ );
  //  }
  //  checkpointer_.debug( jd2::current_output_name(), checkpoint_id, pose.energies().total_energy() );
}

void AbscriptMover::StageTracker::end_stage( core::pose::Pose& pose ){

  if ( tr.Info.visible() ){
    tr.Info << "\n";
    mc_->score_function().show( tr.Info, pose );
  }
  tr.Info << std::endl;

  mc_->show_counters();
  mc_->reset_counters();

  clock_t endtime = clock();
  // PROF_STOP( stagename() );
  if ( basic::options::option[ basic::options::OptionKeys::run::profile ] ) basic::prof_show();
  if ( basic::options::option[ basic::options::OptionKeys::abinitio::debug ]() ) {
    tr.Info << "Timeperstep: " << ( double(endtime) - starttime_ )/( CLOCKS_PER_SEC ) << std::endl;
  }

  if ( basic::options::option[ basic::options::OptionKeys::run::intermediate_structures ]() ) {
    jd2::output_intermediate_pose( pose, "stage"+stagename_ );
  }
}

} // abscript
} // abinitio
} // protocols
