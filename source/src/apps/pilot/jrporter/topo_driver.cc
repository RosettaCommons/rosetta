// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file topo_driver.cc
/// @brief  Sandbox executable for the Environment and associated classes.
/// @author Justin Porter

#include <basic/Tracer.hh>

#include <core/id/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/environment/DofPassport.hh>

#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvExcn.hh>

#include <protocols/environment/claims/TorsionClaim.hh>

#include <protocols/environment/movers/FragmentCM.hh>

#include <protocols/simple_moves/FragmentMover.hh>

#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

static basic::Tracer tr("main");

int main( int argc, char** argv ){

  devel::init( argc, argv );

  core::pose::Pose pose;
  core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");

  using namespace protocols::environment;

  using namespace core::fragment;
  using namespace basic::options;

  core::fragment::FragmentIO frag_io( option[ OptionKeys::abinitio::number_3mer_frags ](),
                                      1, option[ OptionKeys::frags::annotate ]() );
  core::fragment::FragSetOP frags_small = frag_io.read_data( "/Users/jrporter/boinc_aacs_hr1958_03_05.200_v1_3.gz" );
  core::fragment::FragSetOP frags_large = frag_io.read_data( "/Users/jrporter/boinc_aacs_hr1958_09_05.200_v1_3.gz" );

  protocols::simple_moves::FragmentMoverOP mover_small = new protocols::simple_moves::ClassicFragmentMover( frags_small );
  protocols::simple_moves::FragmentMoverOP mover_large = new protocols::simple_moves::ClassicFragmentMover( frags_large );

  FragmentCMOP claim_small = new FragmentCM( mover_small, "BASE" );
  FragmentCMOP claim_large = new FragmentCM( mover_large, "BASE" );

  Environment env( "env" );
  env.register_mover( claim_large );
  env.register_mover( claim_small );


  {
    pose.dump_pdb( "pre.pdb" );
    core::pose::Pose protected_pose = env.start( pose );
    protected_pose.dump_pdb( "init.pdb" );
    claim_small->apply( protected_pose );
    claim_small->apply( protected_pose );
    claim_small->apply( protected_pose );
    claim_small->apply( protected_pose );
    protected_pose.dump_pdb("apply.pdb" );
    env.end( protected_pose );
  }

}