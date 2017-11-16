// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  /src/apps/public/scenarios/beta_strand_homodimer_design/homodimer_design.cc
/// @brief  Symmetric homodimer design


// Unit headers
//none

// devel headers
#include <devel/init.hh>

//core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/ScoreType.hh>

//basic
#include <basic/MetricValue.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/Tracer.hh>

//protocols
#include <protocols/moves/Mover.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperation.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.hh>
#include <protocols/simple_moves/ddG.hh>
//#include <protocols/moves/PymolMover.hh>
//symmetry
#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

//JD2
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>


// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ headers
#include <sstream>
#include <iostream>
#include <string>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "apps.public.beta_strand_homodimer_design.homodimer_design" );

using namespace core;
using namespace utility;
using namespace protocols;
using namespace protocols::moves;
using namespace core::pack::task;
using namespace core::pack::task::operation;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::jd2;
using namespace core::pose::metrics;

//  Index for option[ OptionKeys::symmetry::perturb_rigid_body_dofs ] values
enum Index{
	trans = 1,
	rot = 2
};

// application specific options
basic::options::BooleanOptionKey const make_ala_interface( "make_ala_interface" );
basic::options::IntegerOptionKey const pack_min_runs( "pack_min_runs" );
basic::options::BooleanOptionKey const find_bb_hbond_E( "find_bb_hbond_E" );
basic::options::BooleanOptionKey const skip_hd_docking( "skip_hd_docking" );
basic::options::StringOptionKey const disallow_res( "disallow_res" );
//basic::options::BooleanOptionKey const pymolreport( "pymolreport" );

// mover deffinition
class HDdesignMover : public Mover {
public:

	HDdesignMover();

	virtual void apply( core::pose::Pose& pose );

	virtual ~HDdesignMover(){};

	void cloak_and_setup(core::pose::Pose & pose);

	void task_constraint_setup(core::pose::Pose & pose );

	void sym_repack_minimize(core::pose::Pose & pose );

	void register_calculators();

	void ala_interface(core::pose::Pose & pose);

	core::Real calc_bb_E(core::pose::Pose & pose,
		core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn);

	virtual MoverOP clone() const {
		return MoverOP( new HDdesignMover( *this ) );
	}

	virtual MoverOP fresh_instance() const {
		return MoverOP( new HDdesignMover );
	}
	virtual
	std::string
	get_name() const {
		return "HDdesignMover";
	}

private:
	//init on new input control
	//bool init_new_input_;
	core::scoring::ScoreFunctionOP scorefxn_a;
	core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn_;
	TaskFactoryOP  tf_design_;
	//kinematics::MoveMapOP movemap_;
	pack::task::PackerTaskOP task_design_;
	Size monomer_nres_;
	bool ala_interface_, find_bb_binding_E_, skip_hd_docking_;// pymolreport_;
	int n_pack_min_runs_;
	std::string disallow_res_;


	//alanine interface movers
	protocols::protein_interface_design::movers::BuildAlaPoseOP build_ala_mover_;
	protocols::protein_interface_design::movers::SaveAndRetrieveSidechainsOP get_sidechains_mover_;
	//calculators
	/// @brief InterfaceNeighborDefinition calculator name string
	std::string InterfaceNeighborDefinition_;
	//other movers
	//devel::anchored_design::InterfaceAnalyzerMoverOP interface_mover_;

};

HDdesignMover::HDdesignMover() {
	// variable definitions
	scorefxn_a = core::scoring::get_score_function();
	scorefxn_ = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn_a );
	tf_design_ = TaskFactoryOP( new TaskFactory() );
	//movemap_ = new core::kinematics::MoveMap();
	//options
	ala_interface_ = option[ make_ala_interface];
	n_pack_min_runs_ = option[ pack_min_runs ];
	find_bb_binding_E_ = option[ find_bb_hbond_E ];
	disallow_res_ =  option[ disallow_res ];
	skip_hd_docking_ = option[ skip_hd_docking ];
	//pymolreport_ = option[ pymolreport ];
	//EM options for bb-bb hbond output
	//scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn_->energy_method_options() );
	//energymethodoptions.hbond_options()->decompose_bb_hb_into_pair_energies(true);
	//scorefxn_->set_energy_method_options( energymethodoptions );

}


void HDdesignMover::cloak_and_setup( pose::Pose & pose ){

	//cloak symmetry::perturb_rigid_body_dofs
	if ( basic::options::option[ OptionKeys::symmetry::perturb_rigid_body_dofs ].user() ) {
		TR<< "Cloaking option[ OptionKeys::symmetry::perturb_rigid_body_dofs ] before setup." <<std::endl;
		utility::vector1< Real > pert_mags = basic::options::option[ OptionKeys::symmetry::perturb_rigid_body_dofs ]();
		TR << "Input symmetry::perturb_rigid_body_dofs rot=" << pert_mags[rot] << "  trans=" << pert_mags[trans] << std::endl;
		vector1<Real> zero_vector;
		zero_vector.push_back(0.0);
		zero_vector.push_back(0.0);
		basic::options::option[ OptionKeys::symmetry::perturb_rigid_body_dofs ].value(zero_vector);
		//setup apply
		protocols::simple_moves::symmetry::SetupForSymmetryMoverOP setup_mover( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
		setup_mover->apply(pose);
		//uncloak
		basic::options::option[ OptionKeys::symmetry::perturb_rigid_body_dofs ].value(pert_mags);
	} else {
		//setup apply
		protocols::simple_moves::symmetry::SetupForSymmetryMoverOP setup_mover( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
		setup_mover->apply(pose);
	}


}//end cloak_and_setup

void HDdesignMover::register_calculators(){
	Size const chain1(1), chain2(2);
	//Interface calculator  should only have 2 chains here.
	InterfaceNeighborDefinition_ = "InterfaceNeighborCalc";
	if ( CalculatorFactory::Instance().check_calculator_exists( InterfaceNeighborDefinition_ ) ) {
		TR << "In InterfaceNeighborCalc, calculator " << InterfaceNeighborDefinition_
			<< " already exists." << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator(
			InterfaceNeighborDefinition_,
			PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator(chain1, chain2) ));
		//TR<<"Registering calculator " << InterfaceNeighborDefinition_ << std::endl;
	}
}//end register_calculators

void HDdesignMover::task_constraint_setup( pose::Pose & pose ){

	//allowed_aas_[ chemical::aa_cys ] = false;
	//task factory setup
	tf_design_->clear();
	tf_design_->push_back(  TaskOperationCOP( new InitializeFromCommandline() ) );
	tf_design_->push_back( TaskOperationCOP( new operation::IncludeCurrent ) );
	//if using a resfile ignore all other task restrictions
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		TR << "Using resfile, ignoring all other task info" << std::endl;
		tf_design_->push_back( TaskOperationCOP( new operation::ReadResfile() ) );
	} else {
		if ( option[ disallow_res ].user() ) {
			TR << "Not allowing residues: " << disallow_res_ << " unless native" << std::endl;
			DisallowIfNonnativeOP disallow_op( new DisallowIfNonnative() );
			disallow_op->disallow_aas(disallow_res_);
			tf_design_->push_back(disallow_op);
		}

		//restrict to interface
		tf_design_->push_back( TaskOperationCOP( new protocols::toolbox::task_operations::RestrictToInterfaceOperation( InterfaceNeighborDefinition_ ) ) );
	}


	//apply any constraints
	if ( option[OptionKeys::enzdes::favor_native_res].user() ) {
		utility::vector1< core::scoring::constraints::ConstraintCOP > favor_native_constraints;
		using namespace core::scoring::constraints;
		core::Real bonus = option[OptionKeys::enzdes::favor_native_res].value();
		TR << "favor_native_res: adding a bonus of " << bonus << " for native residues to pose." << std::endl;
		//   //safety check first
		//   if( favor_native_constraints.size() != 0 ){
		//    TR.Warning << "when setting up favor native constraints, there might already be some previously generated favor_native constraints in the pose, trying to remove these first." << std::endl;
		//   favor_native_constraints.clear();
		//   }
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			ConstraintOP resconstraint( new ResidueTypeConstraint( pose, i, bonus ) );
			favor_native_constraints.push_back( resconstraint );
		}
		//adds to pose and scorefxn
		pose.add_constraints( favor_native_constraints );
		scorefxn_->set_weight( core::scoring::res_type_constraint, 1.0 ); //weight of 1.0 means that all is controled by favor_native_res flag
	}


	//now make the movemap symmetric so it will work and not break the pose
	//shouldn't need to do this...
	//core::pose::symmetry::make_symmetric_movemap( pose, *movemap_ );

}//end task_constraint_setup

void HDdesignMover::sym_repack_minimize( pose::Pose & pose ){
	//need to setup the movemap her to correspond to what is at the interface currently
	//calc interface for
	basic::MetricValue< std::set< Size > > interface_mv;
	pose.metric( InterfaceNeighborDefinition_,  "interface_residues", interface_mv);

#ifndef NDEBUG
	std::set< Size > interface_set = ( interface_mv.value() );
	TR<< "Interface residues are: \n";
	for ( std::set< core::Size >::const_iterator it(interface_set.begin()), end(interface_set.end());
			it != end; ++it ) {
		TR << *it << ", ";
	}
	TR << std::endl;
	TR<< "Fold tree for pose: \n" << pose.fold_tree() << std::endl;
#endif

	//std::set< Size > interface_set = ( interface_mv.value() );
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( true );
	//  make_symmetric_movemap( pose, *mm );

	//   SymAtomTreeMinimizer minimizer;
	//   MinimizerOptionsOP min_options = new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.00001, true, false, false );
	//   minimizer.run( pose, *mm, *scorefxn, *min_options );


	//   //setup move map
	//   //set move map to allow bb and sc minimization for interface residues
	//  movemap_->clear();
	//  movemap_->set_jump(true);
	//  //movemap_->set_bb(false); //this works while setting individual doesn't
	//  for( Size ii=1; ii<= pose.size(); ++ii){
	//   if( interface_set.count(ii) ){
	//    movemap_->set_bb(ii, true);
	//    movemap_->set_chi(ii, true);
	//   }
	//   else{
	//    movemap_->set_bb(ii, false);
	//    movemap_->set_chi(ii, false);
	//   }
	//  } //end the movemap creation


	protocols::simple_moves::symmetry::SymMinMoverOP sym_minmover( new protocols::simple_moves::symmetry::SymMinMover(mm, scorefxn_, option[ OptionKeys::run::min_type ].value(), 0.001, true /*use_nblist*/ ) );

	task_design_ = tf_design_->create_task_and_apply_taskoperations( pose );
	protocols::simple_moves::symmetry::SymPackRotamersMoverOP sym_pack_design( new protocols::simple_moves::symmetry::SymPackRotamersMover(scorefxn_, task_design_) );

	TR<< "Monomer total residues: "<< monomer_nres_ << " Repacked/Designed residues: "
		<< task_design_->num_to_be_packed() / 2 << std::endl;

#ifndef NDEBUG
	TR<< "DESIGN Packer Task after setup: " << *(task_design_) <<std::endl;
#endif

	TR<< "Number of repack/minimize runs to do: " << n_pack_min_runs_ << std::endl;
	TR << "Minimizing with: " << option[ OptionKeys::run::min_type ].value() << std::endl;
	for ( int ii = 1; ii<=  n_pack_min_runs_; ++ii ) {
		sym_pack_design->apply(pose);
		sym_minmover->apply(pose);
		TR << "Run " << ii << " of " << n_pack_min_runs_ << "   SCORE:"
			<<  (* scorefxn_ )(pose) << std::endl;
	}
}

//mutate the interface to all alanine if needed
void HDdesignMover::ala_interface(core::pose::Pose & pose ){
	build_ala_mover_ = protocols::protein_interface_design::movers::BuildAlaPoseOP( new protocols::protein_interface_design::movers::BuildAlaPose( true, true, 8.0) ) ;
	get_sidechains_mover_ = protocols::protein_interface_design::movers::SaveAndRetrieveSidechainsOP( new protocols::protein_interface_design::movers::SaveAndRetrieveSidechains(pose) );
	build_ala_mover_->apply( pose );
}  //end ala _interface

///////////////////////////////////////////////////////////////
//finds the bb-bb hbonding energy between chains and returns it
core::Real HDdesignMover::calc_bb_E(core::pose::Pose & pose,
	core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn){

	using namespace core::scoring::hbonds;

	//make copies to avoid screwing up real score and pose
	core::pose::Pose pose_copy ( pose );
	core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn_copy (scorefxn);
	//set up dssp info
	core::scoring::dssp::Dssp dssp( pose_copy );
	dssp.insert_ss_into_pose( pose_copy );

	//EM options for bb-bb hbond output
	scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn_copy->energy_method_options() );
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
	scorefxn_copy->set_energy_method_options( energymethodoptions );

	//now score with everything set.
	(*scorefxn_copy)(pose_copy);

	//figure out energy statistics
	core::scoring::hbonds::HBondSet hbond_set;
	Real bb_score(0.0);
	Size n_hbonds (0);
	//Real hb_wt( (*scorefxn_copy).get_weight(core::scoring::hbond_lr_bb) );
	//find backbone hbonds for pdb
	pose_copy.update_residue_neighbors();
	fill_hbond_set( pose_copy,
		false /*calc_deriv*/,
		hbond_set,
		true /*bb only*/ );
	//call to try to resize bb_don/accept arrays
	//need this for everything to work right
	hbond_set.setup_for_residue_pair_energies(pose_copy);

	//itterate through all hbonds and figure out which ones are bb-bb betas
	for ( Size ii=1; ii <= hbond_set.nhbonds(); ++ii ) {
		HBond hbond ( hbond_set.hbond(ii) );
		//now filter based on what we want
		if ( hbond.don_hatm_is_backbone() && hbond.acc_atm_is_backbone() ) {
			if ( pose_copy.chain( hbond.don_res() ) != pose_copy.chain( hbond.acc_res() ) ) {
				if ( pose_copy.secstruct( hbond.don_res() ) == 'E' && pose_copy.secstruct( hbond.acc_res() ) =='E' ) {
					bb_score += ( hbond.weight() * hbond.energy() );
					++n_hbonds;
				}
			}
		}
	}
	TR << "Design has: "<< n_hbonds << " bb-bb hbonds with total evergy: " << bb_score << std::endl;

	//reset some things to prevent bad scoring outside of this
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(false);
	scorefxn->set_energy_method_options( energymethodoptions );

	return bb_score;
}//end calc_bb_E

/////////////////////////////////////////////////
// Actual mover apply
/////////////////////////////////////////////////
void HDdesignMover::apply (pose::Pose & pose ) {

	// //for pymol viewing
	// if( pymolreport_ ){
	//  protocols::moves::PyMOLObserverOP pymol_ob =  protocols::moves::AddPyMOLObserver(pose, false, 0.1);
	// }

	TR << "Homodimer Design start."<<std::endl;
	//get job info
	protocols::jd2::JobOP const job_me( JobDistributor::get_instance()->current_job() );
	//  std::string job_name (JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );

	monomer_nres_ = pose.size();

	//setup the symmetric pose and cloak it from option[ symmetry::perturb_rigid_body_dofs ]
	//if need be.
	cloak_and_setup( pose );

	//Factory & contraint setup, do based on the structure should be
	register_calculators();
	task_constraint_setup(pose);

	//debugging checks
	//JobDistributor::get_instance()->job_outputter()->other_pose( job_me, pose, "setup_");

	//mutate interface residues to ALA if needed
	if ( ala_interface_ ) {
		TR << "Building Alanine at interface" <<std::endl;
		ala_interface(pose);
	}


	////////////////////////////////////////////////
	//APPLY MOVERS
	////////////////////////////////////////////////

	//debugging
	//kinematics::FoldTree ft (pose.fold_tree());
	//TR << "Fold Tree for pose: \n" << ft << std::endl;

	//Docking first unless skipped
	if ( !skip_hd_docking_ ) {
		//using namespace core::scoring;
		symmetric_docking::SymDockProtocolOP dock_mover( new symmetric_docking::SymDockProtocol );
		dock_mover->apply(pose);
		scoring::ScoreFunctionOP dock_score_apple = scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
		core::scoring::symmetry::SymmetricScoreFunctionOP sym_dock_score = core::scoring::symmetry::symmetrize_scorefunction( *dock_score_apple );
		TR << "Docking SCORE final: " << (*sym_dock_score)(pose) << std::endl;
		TR << "Default SCORE after docking: " << (*scorefxn_)(pose) << std::endl;
		//debugging checks
		//JobDistributor::get_instance()->job_outputter()->other_pose( job_me, pose, "postdock_");

		//Symmetric dock messes something up when it does centroid mode, repack whole protein if need be
		//if( !option[ OptionKeys::docking::docking_local_refine ]() ){
		TaskFactoryOP  tf_nataa( new TaskFactory() );
		tf_nataa->push_back(  TaskOperationCOP( new InitializeFromCommandline() ) );
		tf_nataa->push_back( TaskOperationCOP( new operation::IncludeCurrent ) );
		//want to just repack the wt pose, NO design
		RestrictResidueToRepackingOP repack_op( new RestrictResidueToRepacking() );
		for ( Size ii = 1; ii<= pose.size(); ++ii ) {
			repack_op->include_residue( ii );
		}
		//fill task factory with these restrictions
		tf_nataa->push_back( repack_op );
		PackerTaskOP task_nataa = tf_nataa->create_task_and_apply_taskoperations( pose );
		protocols::simple_moves::symmetry::SymPackRotamersMoverOP sym_pack_nataa( new protocols::simple_moves::symmetry::SymPackRotamersMover(scorefxn_, task_nataa) );
		sym_pack_nataa->apply( pose );
		TR << "Default SCORE after all NATAA repack: " << (*scorefxn_)(pose) << std::endl;
		//JobDistributor::get_instance()->job_outputter()->other_pose( job_me, pose, "nataarepack_");
		// }
	}

	//retrieve if needed
	if ( ala_interface_ ) {
		//TR<< "Packing interface before recovering native residue" <<std::endl;
		//sym_repack_minimize(pose);
		TR<< "Replacing native sidechains." << std::endl;
		get_sidechains_mover_->apply(pose);
	}


	//now repack/minimize
	sym_repack_minimize(pose);
	//debugging checks
	//JobDistributor::get_instance()->job_outputter()->other_pose( job_me, pose, "postpackmin_");

	//final minimization step
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( true );
	protocols::simple_moves::symmetry::SymMinMoverOP sym_minmover_final( new protocols::simple_moves::symmetry::SymMinMover( mm, scorefxn_, option[ OptionKeys::run::min_type ].value(), 0.01, true /*use_nblist*/ ) );
	sym_minmover_final->apply(pose);
	TR << "Final minimization  SCORE:" <<  (* scorefxn_ )(pose) << std::endl;

	//output RMSD
	if ( basic::options::option[ in::file::native ].user() ) {
		Real rms(0.0);
		pose::Pose native_pose; //native pose should be symmetry mate dimer!
		core::import_pose::pose_from_file( native_pose, basic::options::option[ in::file::native ](), core::import_pose::PDB_file);
		// allow superposition because RB min is allowed
		rms = scoring::CA_rmsd( native_pose, pose /*, 1, monomer_nres_ */ ) ;
		job_me->add_string_real_pair("rms_sym", rms );
	}

	//find ddg
	protocols::simple_moves::ddGOP ddG_mover( new protocols::simple_moves::ddG( scorefxn_, 1 /*jump*/ /* , true */ ) ); //ddG autodetects symmetry now
	ddG_mover->calculate(pose);
	core::Real ddgvalue = ddG_mover->sum_ddG();
	//some dirty filtering
	if ( ddgvalue >= 0 ) {
		set_last_move_status(protocols::moves::FAIL_RETRY);
	} else {
		job_me->add_string_real_pair("dGbind", ddgvalue);
	}


	//find bb-bb hbond E
	if ( find_bb_binding_E_ ) {
		Real  bb_score (calc_bb_E( pose,scorefxn_));
		job_me->add_string_real_pair("beta_int_E", bb_score );
		//some dirty filtering
		// if (bb_score >= 0)
		//set_last_move_status(protocols::moves::FAIL_RETRY);
		//else

	}
	//score to be safe
	(* scorefxn_ )(pose);

}//end mover apply

////////////////////////////////////////////////
// Main
////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		option.add( make_ala_interface, "Make interface residues ALA." ).def(false);
		option.add( pack_min_runs, "Number of runs of repack/minimize" ).def(1);
		option.add( find_bb_hbond_E, "Find the energy of bb-bb interactions").def(false);
		option.add( skip_hd_docking, "skips docking step and just does design").def(false);
		//option.add( pymolreport, "Report to pymol observer").def(false);
		option.add( disallow_res, "String of residues not allowed unless native").def("");

		// init
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new HDdesignMover ) );

		std::cout << "Done! -------------------------------"<< std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
} //end main
