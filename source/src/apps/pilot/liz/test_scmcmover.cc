// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <protocols/moves/MonteCarlo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <devel/init.hh>
#include <basic/prof.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace basic;


OPT_1GRP_KEY(Boolean, test, default_mc)
OPT_1GRP_KEY(Boolean,test, detailed_balance )
OPT_1GRP_KEY(Integer,test,ntrials)
OPT_1GRP_KEY(Integer,test,output_interval)
OPT_1GRP_KEY(Real,test,unif)
OPT_1GRP_KEY(Real,test,pert)
OPT_1GRP_KEY(Real,test,within)

int
main(int argc, char* argv []){
    try {
	using namespace basic::options::OptionKeys;
	using namespace basic::options;

	NEW_OPT(test::default_mc, "", false );
	NEW_OPT(test::detailed_balance,"",false);
	NEW_OPT(test::ntrials,"",1000);
	NEW_OPT(test::output_interval,"",1);
	NEW_OPT(test::unif,"",0.1);
	NEW_OPT(test::pert,"",0);
	NEW_OPT(test::within,"",0);

	devel::init(argc,argv);
	bool sidechainmover = option[test::default_mc];
	bool detailed_balance = option[ test::detailed_balance ];
	core::Size ntrial = option[ test::ntrials ];
	core::Size output_interval = option[ test::output_interval ];

	//probabilities
	core::Real unif = option[ test::unif ];
	core::Real pert = option[ test::pert ];
	core::Real within = option[ test::within ];

	pose::Pose pose;
	core::import_pose::pose_from_file( pose , "test.pdb" , core::import_pose::PDB_file);
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	//SilentStructOP ssin = SilentStructFactory::get_silent_struct_in();
	if( !sidechainmover ){ //38.0 accorind to util.prof
	  std::cout << "Selection: FAST-SC mover with detailed balance " << detailed_balance << std::endl;
	  protocols::simple_moves::sidechain_moves::SidechainMCMover scmc;
	  pack::task::PackerTaskOP pt = core::pack::task::TaskFactory::create_packer_task( pose );
	  pt->initialize_from_command_line();
	  if ( option[ basic::options::OptionKeys::packing::resfile ].user() )
			parse_refile(pose, *pt);
	  pt->restrict_to_repacking();
	  scmc.set_task( pt );
	  scmc.init_task( pose );
	  scmc.set_ntrials( output_interval );
	  //std::cerr << "ntrials are " << scmc.ntrials() << std::endl;
	  scmc.set_prob_uniform( unif );
	  scmc.set_prob_withinrot(within );
	  scmc.set_prob_random_pert_current( pert );
	  scmc.set_preserve_detailed_balance( detailed_balance );
	  scmc.set_temperature( 1 );
	  scmc.setup(sfxn);
	  scmc.set_sampling_temperature(1);
	  //scmc.setup(pose);
	  //scmc.set_scorefunction( sfxn );
	  //PROF_START( SIDECHAINMCMOVER );
	  //std::cerr << "starting sc-trial" << std::endl;
	  for ( Size ii = 1; ii <= (ntrial / output_interval); ++ii ) {
	    scmc.apply(pose);
	    std::cout << "FAST-SC: " << (*sfxn)( pose ) << " " << scmc.last_proposal_density_ratio() << std::endl;

	  }
	  //std::cerr << "ending sc-trial " << std::endl;
	  PROF_STOP( SIDECHAINMCMOVER );
	  basic::prof_show();
	  pose.dump_pdb( "end_test.pdb");
	  (*sfxn)( pose );
	  sfxn->show( pose );

	}else{
	  std::cout << "Selection: Regular SC mover with detailed balance " << detailed_balance;

	  protocols::simple_moves::sidechain_moves::SidechainMover sc;
	  pack::task::PackerTaskOP pt = core::pack::task::TaskFactory::create_packer_task( pose );
	  protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( *sfxn, 1.0);
	  mc->reset( pose );
	  core::pose::Pose last_accepted = pose;
	  pt->initialize_from_command_line();
	  if ( option[ basic::options::OptionKeys::packing::resfile ].user() )
			parse_resfile(pose, *pt);
	  pt->restrict_to_repacking();
	  sc.set_task( pt );
	  sc.set_sampling_temperature( 1.0);
	  sc.init_task( pose );
	  sc.set_prob_uniform( unif );
	  sc.set_prob_withinrot( within );
	  sc.set_prob_random_pert_current( pert );
	  sc.set_preserve_detailed_balance( detailed_balance );
	  //std::cerr << "sidechain sampling temp: " << sc.sampling_temperature() << " mc-temp: 1.0" << std::endl;
	  //std::cerr << "starting sc-trial" << std::endl;
	  //PROF_START( SIDECHAINMCMOVER );
	  for( Size i = 1; i <= ntrial; i++){
	    //PROF_START( SIDECHAINMOVER );
	    core::Real prev_score = (*sfxn)(pose);
	    sc.apply( pose );
	    //std::cout << "score-after-move: " << (*sfxn)(pose) << " delta-score " << ((*sfxn)(pose) - prev_score) << " " << sc.type() <<std::endl;
	    //PROF_STOP( SIDECHAINMOVER );
	    mc->boltzmann( pose, "scmove", sc.last_proposal_density_ratio());
	    if( i % output_interval == 0 ) {
	      last_accepted = mc->last_accepted_pose();
	      std::cout << "SC: " << (*sfxn)( last_accepted ) << " " << sc.last_proposal_density_ratio() << std::endl;
	      //std::cout << "SC: " << (*sfxn)( pose ) << " " << sc.last_proposal_density_ratio() << std::endl;
	    }
	  }
	  //std::cerr << "ending sc-trial " << std::endl;
	  //PROF_STOP( SIDECHAINMCMOVER );
	  mc->show_counters();

	  basic::prof_show(); // 68.0 according to util.prof
	  pose.dump_pdb( "end_sc_test.pdb");
	  (*sfxn)( pose );
	  sfxn->show( pose );
	}
    } catch (utility::excn::Exception const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
}

