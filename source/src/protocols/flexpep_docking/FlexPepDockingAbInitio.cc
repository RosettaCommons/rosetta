// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (C) 199x-2008 Hebrew University, Jerusalem
//
/// @file   FlexPepDockingAbInitio.hh
///
/// @brief low-resolution part of docking protocol
/// @date August 5, 2008
/// @author Barak Raveh


#include <protocols/flexpep_docking/FlexPepDockingFlags.hh>
#include <protocols/flexpep_docking/FlexPepDockingAbInitio.hh>

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/flexPepDocking.OptionKeys.gen.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <utility/exit.hh>
#include <string>
#include  <math.h>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using namespace protocols::flexpep_docking;

static THREAD_LOCAL basic::Tracer TR( "protocols.flexPepDockingAbInitio" );


//////////////////////////////////////////////
/// @brief
/// constructor for low resolution flexpible peptide docking
//
// @param[in] scorefxn_in
//            The scoring function used for optimization
// @param[in] rb_jump
//            The FoldTree rigid body jump over
//            which rigid-body pertrubations are made
FlexPepDockingAbInitio::FlexPepDockingAbInitio
( FlexPepDockingFlags flags_in,
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::kinematics::MoveMapOP movemap_in,
	Size const rb_jump_in )
: flags_(flags_in),
	movemap_(movemap_in),
	rb_jump_(rb_jump_in),
	fragset3mer_(/* NULL */),
	fragset9mer_(/* NULL */),
	fragset5mer_(/* NULL */)
{
	using namespace basic::options;

	// TODO: create a set_defaults() function
	scorefxn_ = scorefxn_in->clone();

	// Loop modeling options
	// NOTE: most LoopRelax options are initiated automatically from cmd-line
	// TODO: LoopRelaxMover is a wrapper, perhaps user the LoopModel class explicitly
	loop_relax_mover_ = protocols::comparative_modeling::LoopRelaxMoverOP( new protocols::comparative_modeling::LoopRelaxMover() );
	// loop_relax_mover_->centroid_scorefxn(scorefxn_); // TODO: we need a chain brteak score here, so let's leave it for modeller default?
	loop_relax_mover_->refine("no"); // centroid modeling only
	loop_relax_mover_->relax("no"); // centroid modeling only
	if ( option[ OptionKeys::loops::frag_files ].user() ) {
		// these protocols optionally take a fragment set .. only load if
		// specified
		utility::vector1< core::fragment::FragSetOP > frag_libs;
		protocols::loops::read_loop_fragments( frag_libs );
		loop_relax_mover_->frag_libs( frag_libs );
	}
	// TODO: choose frag files through FlexPepDockingFlags
	// 3mer fragments
	if ( option[ basic::options::OptionKeys::in::file::frag3].user() ) {
		fragset3mer_ = core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( 3 ) );
		fragset3mer_->read_fragment_file
			( option[ basic::options::OptionKeys::in::file::frag3].value() );
	}
	// 9mer fragments
	if ( option[ basic::options::OptionKeys::in::file::frag9].user() ) {
		fragset9mer_ = core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( 9 ) );
		fragset9mer_->read_fragment_file
			( option[ basic::options::OptionKeys::in::file::frag9].value() );
	}
	//5mer fragments
	if ( option[ basic::options::OptionKeys::flexPepDocking::frag5].user() ) {
		fragset5mer_ = core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( 5 ) );
		fragset5mer_->read_fragment_file
			( option[ basic::options::OptionKeys::flexPepDocking::frag5].value() );
	}

}


// empty destructor - for good inclusion of OP clasesses
FlexPepDockingAbInitio::~FlexPepDockingAbInitio()
{}


///////////////////////////////////////////////
/// @brief initial setup for apply
void
FlexPepDockingAbInitio::setup_for_apply( core::pose::Pose& pose )
{
	double temperature = 0.8;
	mc_ = moves::MonteCarloOP( new moves::MonteCarlo( pose, *scorefxn_, temperature ) );
	// setup minimizer
	std::string min_type = "lbfgs_armijo_atol"; // armijo_nonmonotone? different tolerance?
	double min_func_tol = 0.1;
	minimizer_ = protocols::simple_moves::MinMoverOP( new protocols::simple_moves::MinMover(
		movemap_, scorefxn_, min_type, min_func_tol, true /*nb_list accel.*/ ) );
}


///////////////////////////////////////////////
// switch pose to centroid mode, if not already there
void
FlexPepDockingAbInitio::to_centroid
( core::pose::Pose & pose ) const
{
	if ( !pose.is_fullatom() ) {
		return;
	}
	TR.Debug << "Switching to centroid" << std::endl;
	protocols::simple_moves::SwitchResidueTypeSetMover
		to_centroid_mover( core::chemical::CENTROID );
	to_centroid_mover.apply(pose);
}


//////////////////////////////////////////////
// switch pose to full-atom mode, if not already
// using side-chains of referencePose
void
FlexPepDockingAbInitio::to_allatom
( core::pose::Pose & pose, core::pose::Pose& referencePose ) const
{
	runtime_assert(referencePose.is_fullatom());
	//  protocols::simple_moves::SwitchResidueTypeSetMover
	//to_all_atom_mover( core::chemical::FA_STANDARD );
	protocols::simple_moves::ReturnSidechainMover
		recover_sidechains( referencePose );
	recover_sidechains.apply( pose );
}


////////////////////////////////////////////
void
FlexPepDockingAbInitio::torsions_monte_carlo
( core::pose::Pose & pose,
	const int cycles,
	double& acceptance_rate )
{
	using namespace protocols::moves;
	using namespace protocols::simple_moves;

	// setup sub-moves
	simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( movemap_, mc_->temperature() /*temp*/, 5 /*nmoves ???*/ ) );
	small_mover->angle_max('L',180 /*angle - TODO: parametrize and slowly ramp down */);
	small_mover->angle_max('H',180 /*angle - TODO: parametrize and slowly ramp down */);
	small_mover->angle_max('E',180 /*angle - TODO: parametrize and slowly ramp down */);
	protocols::simple_moves::ShearMoverOP shear_mover( new protocols::simple_moves::ShearMover( movemap_, mc_->temperature() /*temp*/, 5 /*nmoves ???*/ ) );
	shear_mover->angle_max('L',180 /*angle - TODO: parametrize, slowly ramp down */);
	shear_mover->angle_max('H',180 /*angle - TODO: parametrize, slowly ramp down */);
	shear_mover->angle_max('E',180 /*angle - TODO: parametrize, slowly ramp down */);
	ClassicFragmentMoverOP
		frag3_mover = NULL,
		frag9_mover = NULL,
		frag5_mover = NULL;
	if ( fragset3mer_ ) { //if we have fragments
		frag3_mover = ClassicFragmentMoverOP( new ClassicFragmentMover(fragset3mer_, movemap_) );
		frag3_mover->enable_end_bias_check(false);
	}
	if ( fragset9mer_ ) { //if we have fragments
		frag9_mover = ClassicFragmentMoverOP( new ClassicFragmentMover(fragset9mer_, movemap_) );
		frag9_mover->enable_end_bias_check(false);
	}
	if ( fragset5mer_ ) { //if we have fragments
		frag5_mover = ClassicFragmentMoverOP( new ClassicFragmentMover(fragset5mer_, movemap_) );
		frag5_mover->enable_end_bias_check(false);
	}

	// setup cycle of sub-moves
	RandomMoverOP random_mover( new protocols::moves::RandomMover() );
	random_mover->add_mover(small_mover,1.0 /*weight*/);
	if ( frag3_mover ) {
		random_mover->add_mover(frag3_mover,flags_.frag3_weight /*weight*/);
	}
	if ( frag9_mover ) {
		random_mover->add_mover(frag9_mover,flags_.frag9_weight /*weight*/);
	}
	if ( frag5_mover ) {
		random_mover->add_mover(frag5_mover, flags_.frag5_weight /*weight*/);
	}
	random_mover->add_mover(shear_mover, 1.0);

	// wrap with monte-carlo trial mover
	TrialMoverOP mc_trial( new TrialMover( random_mover, mc_ ) );
	mc_trial->keep_stats_type( accept_reject ); // track stats (for acceptance rate)

	// Do initial forced perturbation // TODO: do we want to keep this part?
	//  if(frag9_mover){
	//    frag9_mover->apply(pose);
	//}
	//  small_mover->apply(pose);
	//  shear_mover->apply(pose);
	// run Monte-Carlo
	//  TR << "start MC " <<std::endl;
	for ( int i=1; i<=cycles; ++i ) {
		mc_trial->apply( pose );
	}
	//  TR << "finsihed MC " <<std::endl;
	// extract best pose and return statistics
	pose = mc_->lowest_score_pose();
	mc_->reset( pose );
	acceptance_rate = mc_trial->acceptance_rate();
}


///////////////////////////////////////////
void
FlexPepDockingAbInitio::loopclosure_monte_carlo
( core::pose::Pose & pose )
{
	using namespace protocols::moves;

	if ( flags_.peptide_nres() < 5 ) { // TODO: 3 is minimum for KIC loop closure, and flanks should be excluded
		return;
	}
	// set up and model a random loop
	Size first_res = flags_.peptide_first_res() + 1;
	Size last_res = flags_.peptide_last_res() - 1;
	protocols::loops::LoopsOP loops( new protocols::loops::Loops() );
	loops->add_loop(first_res, last_res); // TODO: cut defaults to zero, is this a random cut?
	for ( Size i = first_res; i <= last_res ; i++ ) {
		runtime_assert( movemap_->get_bb(i) ); // verify loop is movable, this should have been always true for the peptide
	}
	loop_relax_mover_->loops( loops );
	loop_relax_mover_->apply( pose );

}


///////////////////////////////////////////
// returns the acceptance rate from the RB monte-carlo
void
FlexPepDockingAbInitio::rigidbody_monte_carlo
( core::pose::Pose & pose,
	const int cycles,
	const float trans_magnitude,
	const float rot_magnitude,
	double& acceptance_rate
)
{
	using namespace protocols::moves;

	// set up monte-carlo trial moves
	rigid::RigidBodyPerturbNoCenterMoverOP rb_mover( new rigid::RigidBodyPerturbNoCenterMover(
		rb_jump_, rot_magnitude, trans_magnitude ) );
	TrialMoverOP mc_trial( new TrialMover( rb_mover, mc_ ) );
	mc_trial->keep_stats_type( accept_reject ); // track stats (for acceptance rate)

	// run monte-carlo
	for ( int i = 1; i <= cycles ; i++ ) {
		mc_trial->apply(pose);
	}

	// extract best pose and return statistics
	pose = mc_->lowest_score_pose();
	mc_->reset( pose );
	acceptance_rate = mc_trial->acceptance_rate();
}


///////////////////////////////////////////////
void
FlexPepDockingAbInitio::apply( core::pose::Pose & pose )
{
	using namespace core;
	// variable declarations: //
	pose::Pose startPose = pose;
	int outer_cycles=10; // TODO: runtime param?
	int inner_cycles_rb = 50; // TODO: runtime param?
	int inner_cycles_torsions = 50; // TODO: runtime param?
	double trans_mag = 1; // Angstrom // TODO: runtime param?
	double rot_mag = 10; // Degrees // TODO: runtime param?
	double rb_acceptance, torsions_acceptance; // MC acceptance rates
	// MC temperature params
	double init_MC_temp = 2;
	double final_MC_temp = 0.6;
	double gamma = pow( (final_MC_temp / init_MC_temp) , 1.0 / (outer_cycles - 1.0) );
	std::cout << "gamma = " << gamma << std::endl;
	double MC_temp = init_MC_temp;
	// real stuff: //
	std::set<int> pSer_positions;
	if ( flags_.pSer2Asp_centroid ) {
		convertPSERtoASP(pose, pSer_positions); // for low-res
	} else if ( flags_.pSer2Glu_centroid ) {
		convertPSERtoGLU(pose, pSer_positions); // for low-res
	}
	to_centroid(pose);
	setup_for_apply(pose);
	TR.flush();
	TR.flush();
	for ( int i=1; i<=outer_cycles; ++i ) {
		mc_->set_temperature( MC_temp);
		std::cout << "temperature = " << mc_->temperature() << std::endl;
		if ( flags_.rbMCM && ! flags_.pep_fold_only ) {
			rigidbody_monte_carlo
				( pose , inner_cycles_rb, trans_mag, rot_mag, rb_acceptance);
			minimizer_->apply(pose);
		}
		if ( flags_.torsionsMCM ) {
			torsions_monte_carlo
				( pose , inner_cycles_torsions, torsions_acceptance);
			minimizer_->apply(pose);
		}
		if ( flags_.peptide_loop_model ) {
			loopclosure_monte_carlo( pose );
			minimizer_->apply(pose);
		}
		MC_temp *= gamma;
	}
	TR.flush();
	if ( flags_.pSer2Asp_centroid || flags_.pSer2Glu_centroid ) {
		restorePSER(pose, pSer_positions); // restore pSers before hi-res
	}
	to_allatom(pose, startPose /*reference pose*/);
}

std::string
FlexPepDockingAbInitio::get_name() const {
	return "FlexPepDockingAbInitio";
}

// convert all pSer residues to Asp for low-res part, and save list of pSers to pSer_positions
void
FlexPepDockingAbInitio::convertPSERtoASP(core::pose::Pose& pose, std::set<int>& pSer_positions)
{
	using namespace core::chemical;
	using namespace core::conformation;
	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	for ( int resid = 1; resid <= (int)pose.total_residue(); resid++ ) {
		if ( pose.residue_type(resid).has_variant_type(core::chemical::PHOSPHORYLATION) &&
				pose.residue_type(resid).name3() == "SER" ) {
			pSer_positions.insert(resid);
			ResidueOP asp( ResidueFactory::create_residue( rsd_set->name_map("ASP") ));
			pose.replace_residue( resid, *asp, true );
			TR << "Replaced pSer residue " << resid << " to ASP for ab-initio centroid mode" << std::endl;
		}
	}
}


// convert all pSer residues to Glu for low-res part, and save list of pSers to pSer_positions
void
FlexPepDockingAbInitio::convertPSERtoGLU(core::pose::Pose& pose, std::set<int>& pSer_positions)
{
	using namespace core::chemical;
	using namespace core::conformation;
	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	for ( int resid = 1; resid <= (int)pose.total_residue(); resid++ ) {
		if ( pose.residue_type(resid).has_variant_type(core::chemical::PHOSPHORYLATION) &&
				pose.residue_type(resid).name3() == "SER" ) {
			pSer_positions.insert(resid);
			ResidueOP glu( ResidueFactory::create_residue( rsd_set->name_map("GLU") ));
			pose.replace_residue( resid, *glu, true );
			TR << "Replaced pSer residue " << resid << " to GLU for ab-initio centroid mode" << std::endl;
		}
	}
}

// convert all Asp residues to pSer, if they were originally pSer positions (based on pSer_positions)
void
FlexPepDockingAbInitio::restorePSER(core::pose::Pose& pose, std::set<int> const& pSer_positions)
{
	using namespace core::chemical;
	using namespace core::conformation;
	ResidueTypeSetCOP centroid_set( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );
	std::set<int>::const_iterator iter, end;
	for ( iter = pSer_positions.begin(), end = pSer_positions.end(); iter != end; ++iter ) {
		int resid = *iter;
		// Residue const & asp( pose.residue( resid ) ); // Unused variable causes warning.
		ResidueOP pSer( ResidueFactory::create_residue( centroid_set->name_map("SER") ) );
		pose.replace_residue( resid, *pSer, true );
		core::pose::add_variant_type_to_pose_residue( pose , core::chemical::PHOSPHORYLATION, resid );
	}
}
