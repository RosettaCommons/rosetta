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
/// @author
/// @author

//graphix
#include <protocols/viewer/viewers.hh>

//
#include <protocols/canonical_sampling/CanonicalSamplingApplication.hh>
//#include <devel/Gaussian/BBGaussianMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/BBConRotMover.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
// AUTO-REMOVED #include <protocols/jd2/MPIFileBufJobDistributor.hh>
#include <protocols/canonical_sampling/CanonicalSamplingMover.hh>
#include <protocols/canonical_sampling/CanonicalSamplingMover.fwd.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>

// AUTO-REMOVED #include <protocols/canonical_sampling/mc_convergence_checks/MPIHPool_ConvergenceCheck.hh>

// AUTO-REMOVED #include <core/pack/task/ResfileReader.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/bbg.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/canonical_sampling.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>




#ifdef USEMPI
#include <mpi.h>
#endif



using namespace basic::options;



namespace protocols {
namespace canonical_sampling {

int
canonical_sampling_main(){

	basic::Tracer tr("src.apps.pilot.liz.test_canonical_mover");

  using namespace protocols::moves;
	using namespace basic;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;
	
	if( tr.visible() ){
		tr << "resetting clock" << std::endl;
	}
	//basic::prof_reset();

	
	//bool MPI_synchronize_pools = options::option[options::OptionKeys::canonical_sampling::probabilities::MPI_sync_pools];
	bool MPI_bcast = options::option[ options::OptionKeys::canonical_sampling::probabilities::MPI_bcast ];
	bool use_fast_sc_moves = options::option[ options::OptionKeys::canonical_sampling::probabilities::fast_sc_moves ];
	bool no_jd2_output = options::option[options::OptionKeys::canonical_sampling::probabilities::no_jd2_output];
	bool use_hierarchy = options::option[options::OptionKeys::canonical_sampling::probabilities::use_hierarchical_clustering];

	//get rid of unused variable warning
	if( MPI_bcast ) {}
	if( use_hierarchy ) {}

  core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
  //protocols::canonical_sampling::CanonicalSamplingMoverOP csm(new protocols::canonical_sampling::CanonicalSamplingMover(sfxn,pool_ptr,1000));
	protocols::canonical_sampling::CanonicalSamplingMoverOP csm(new CanonicalSamplingMover);
		csm->set_scorefunction(sfxn);

	if( !use_fast_sc_moves ){
		//setup sidechain mover
		protocols::simple_moves::sidechain_moves::SidechainMoverOP scm(new protocols::simple_moves::sidechain_moves::SidechainMover());
		TaskFactoryOP main_task_factory = new TaskFactory;
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory->push_back( rtrop );
		main_task_factory->push_back( new operation::ReadResfile );
		scm->set_task_factory(main_task_factory);

		if( options::option[ OptionKeys::canonical_sampling::probabilities::sc_prob_uniform ].user() ) {
			scm->set_prob_uniform(options::option[ options::OptionKeys::canonical_sampling::probabilities::sc_prob_uniform ]);
		}
		if( options::option[ options::OptionKeys::canonical_sampling::probabilities::sc_prob_withinrot ].user() ){
			scm->set_prob_withinrot(options::option[ options::OptionKeys::canonical_sampling::probabilities::sc_prob_withinrot ]);
		}
		if( options::option[ options::OptionKeys::canonical_sampling::probabilities::sc_prob_perturbcurrent ].user() ){
			scm->set_prob_random_pert_current(options::option[ options::OptionKeys::canonical_sampling::probabilities::sc_prob_perturbcurrent ]);
		}
		scm->set_preserve_detailed_balance( csm->detailed_balance() );

		csm->add_mover(scm,options::option[options::OptionKeys::canonical_sampling::probabilities::sc]);
	}else{
		//setup sidechain mover
		protocols::simple_moves::sidechain_moves::SidechainMCMoverOP scm(new protocols::simple_moves::sidechain_moves::SidechainMCMover());
		TaskFactoryOP main_task_factory = new TaskFactory;
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory->push_back( rtrop );
		scm->set_temperature( csm->get_temp() );
		scm->set_scorefunction(  *sfxn );
		scm->setup( sfxn );
		scm->set_task_factory(main_task_factory);
		scm->set_prob_uniform(options::option[ options::OptionKeys::canonical_sampling::probabilities::sc_prob_uniform ]);
		scm->set_prob_withinrot(options::option[ options::OptionKeys::canonical_sampling::probabilities::sc_prob_withinrot ]);
		scm->set_prob_random_pert_current(options::option[ options::OptionKeys::canonical_sampling::probabilities::sc_prob_perturbcurrent ]);
		scm->set_preserve_detailed_balance( csm->detailed_balance() );
		scm->set_ntrials( options::option[ options::OptionKeys::canonical_sampling::probabilities::fast_sc_moves_ntrials ] );

		csm->add_mover(scm,options::option[options::OptionKeys::canonical_sampling::probabilities::sc]);

	}
	//add move-map to bbg-mover
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	if( options::option[in::file::movemap].user() ){
		movemap->init_from_file(options::option[in::file::movemap]);
	}else{
		movemap->set_bb( true );
		movemap->set_chi( true );
	}
	simple_moves::BBG8T3AMoverOP bbg8t3mover = new simple_moves::BBG8T3AMover();
	bbg8t3mover->movemap(movemap);
	tr << "probability of executing bbg move: " << (options::option[options::OptionKeys::canonical_sampling::probabilities::localbb]*
																									(1-options::option[options::OptionKeys::canonical_sampling::probabilities::backrub]-options::option[options::OptionKeys::canonical_sampling::probabilities::conrot]())) << std::endl;
	runtime_assert(1-options::option[options::OptionKeys::canonical_sampling::probabilities::backrub]-options::option[options::OptionKeys::canonical_sampling::probabilities::conrot] >= 0);
  csm->add_mover( bbg8t3mover, options::option[options::OptionKeys::canonical_sampling::probabilities::localbb]*(1-options::option[options::OptionKeys::canonical_sampling::probabilities::backrub]-options::option[options::OptionKeys::canonical_sampling::probabilities::conrot]()) );

	//setup conrot mover
	if( options::option[options::OptionKeys::canonical_sampling::probabilities::conrot]() > 0 ) {
		simple_moves::BBConRotMoverOP conrotmover = new simple_moves::BBConRotMover();
		conrotmover->movemap( movemap );
		csm->add_mover( conrotmover, options::option[options::OptionKeys::canonical_sampling::probabilities::localbb]*options::option[options::OptionKeys::canonical_sampling::probabilities::conrot]());
	}


	//setup backrub mover
	if( options::option[options::OptionKeys::canonical_sampling::probabilities::backrub]() > 0 ) {
		protocols::backrub::BackrubMoverOP backrubmover = new protocols::backrub::BackrubMover;
		core::scoring::methods::EnergyMethodOptions energymethodoptions(sfxn->energy_method_options());
    if (energymethodoptions.bond_angle_residue_type_param_set()) {
      backrubmover->branchopt().bond_angle_residue_type_param_set(energymethodoptions.bond_angle_residue_type_param_set());
    }
    backrubmover->set_preserve_detailed_balance( csm->detailed_balance() );
		backrubmover->init_with_options();
		//backrubmover->add_mainchain_segments_from_options();
		//if( !option[in::file::fullatom] ) {
    //  backrubmover->optimize_branch_angles(pose);
    //}
		tr << "probability of executing backrub move: " << (options::option[options::OptionKeys::canonical_sampling::probabilities::localbb]*(options::option[options::OptionKeys::canonical_sampling::probabilities::backrub])) << std::endl;
		csm->add_mover((MoverOP)(backrubmover),(core::Real)(options::option[options::OptionKeys::canonical_sampling::probabilities::localbb]*options::option[options::OptionKeys::canonical_sampling::probabilities::backrub]));
	}


	tr.Info << "using regular pool-rmsd" << std::endl;
	//protocols::jd2::JobDistributor* jd2 = protocols::jd2::JobDistributor::get_instance();
	mc_convergence_checks::Pool_RMSD_OP pool_ptr =
		new mc_convergence_checks::Pool_RMSD(options::option[mc::known_structures]);
	csm->set_poolrmsd(pool_ptr);
	if( no_jd2_output ) {
		protocols::jd2::JobDistributor::get_instance()->go( csm, new protocols::jd2::NoOutputJobOutputter );
	} else {
		protocols::jd2::JobDistributor::get_instance()->go( csm );
	}
  return 0;
}

}
}
