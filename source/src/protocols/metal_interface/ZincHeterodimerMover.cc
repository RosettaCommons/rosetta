// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/metal_interface/ZincHeterodimerMover.cc
/// @brief ZincHeterodimerMover methods implemented - see apps/pilot/rjha/README for details
/// @author Steven Lewis


// Unit Headers
#include <protocols/metal_interface/ZincHeterodimerMover.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh> //CENTROID, FA_STANDARD
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>

#include <core/kinematics/FoldTree.hh>
//#include <core/kinematics/Stub.hh>
//#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/rigid/RotateJumpAxisMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/MoverContainer.hh> //Random, Sequence Mover
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/OutputMovers.hh> //PDBDump for movie
#include <protocols/moves/DualMonteCarlo.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperation.hh>

// Numeric Headers


// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <basic/options/option.hh>

// C++ Headers
#include <string>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>




using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.metal_interface.ZincHeterodimerMover" );

namespace protocols {
namespace metal_interface {


///@brief local helper function
void dump_chis( core::conformation::Residue const & res ){
	for(core::Size i(1); i<=res.nchi(); ++i){TR << " chi" << i << " " << res.chi()[i];}
	TR << std::endl;
}

void ZincHeterodimerMover::apply( core::pose::Pose & pose ){

	//extract which residue is which from input edges
	bool const mtm_start_metal(pose.residue_type(metal_to_mobile_.start()).is_protein());
	core::Size const metal_res(mtm_start_metal ? metal_to_mobile_.stop() : metal_to_mobile_.start());
	core::Size const ligand_res(mtm_start_metal ? metal_to_mobile_.start() : metal_to_mobile_.stop());
	core::Size const fixed_res( core::Size(fixed_to_metal_.start()) == metal_res ? //if
															fixed_to_metal_.stop() : //then
															fixed_to_metal_.start()); //else

	/////////////////////////////make fold tree for fullatom pose////////////////////////////////////////////////
	//this code assumes the pose is set up fixed partner1 - metal - mobile partner2, with only the three chains
	//this assumption is duplicated in the private generate_factory() function so be sure to generalize it in both places
	assert(pose.conformation().num_chains() == 3);
	assert(pose.conformation().chain_end(2) == metal_res && pose.conformation().chain_begin(2) == metal_res );
	assert(fixed_res < ligand_res);
	TR << "fixed_res: " << fixed_res << "  ligand_res: " << ligand_res << std::endl;
	using core::kinematics::Edge;
	core::kinematics::FoldTree main_tree(pose.total_residue());
	main_tree.add_edge(Edge(1, fixed_res, Edge::PEPTIDE)); //peptide edge for fixed, lower partner
	main_tree.add_edge(Edge(fixed_res, metal_res-1, Edge::PEPTIDE)); //peptide edge for fixed, lower partner
	main_tree.add_edge(Edge(metal_res+1, ligand_res, Edge::PEPTIDE)); //peptide edge for upper partner
	main_tree.add_edge(Edge(ligand_res, pose.total_residue(), Edge::PEPTIDE)); //peptide edge for upper partner
	main_tree.add_edge(fixed_to_metal_);
	main_tree.add_edge(metal_to_mobile_);
	main_tree.delete_unordered_edge(1, pose.total_residue(), Edge::PEPTIDE);
	main_tree.reorder(1);
	TR << main_tree << std::endl;
	pose.fold_tree(main_tree);

	//////////////////////////make fold_tree for centroid pose//////////////////////////////////////////
	core::kinematics::FoldTree centroid_tree(pose.total_residue());
	centroid_tree.add_edge(Edge(1, metal_res-1, Edge::PEPTIDE)); //peptide edge for fixed, lower partner
	centroid_tree.add_edge(Edge(metal_res+1, pose.total_residue(), Edge::PEPTIDE));//edge for upper
	centroid_tree.add_edge(Edge(//interchain jump
															metal_res+1, //from Cterm of chain 1
															metal_res-1, //to Nterm of chain 3
															1));           //jump number 1
	centroid_tree.add_edge(Edge(//interchain jump
															metal_res,
															pose.total_residue(), //to Cterm of chain 3
															2));           //jump number 2
	centroid_tree.delete_unordered_edge(1, pose.total_residue(), Edge::PEPTIDE);
	centroid_tree.reorder(1);
	TR << centroid_tree << std::endl;

	//////////////////////make RotateJumpAxisMover///////////////////////////
	//passing no angle settings means that it will make random rotations
	protocols::rigid::RotateJumpAxisMoverOP RJAmover(new protocols::rigid::RotateJumpAxisMover(metal_to_mobile_.label()));

	//////////////////////make SidechainMover and its task///////////////////////
	//hand-make a PackerTask for SidechainMover.  A task that is "design nowhere, pack at one position" won't work if
	//the metal ligand is histidine, because of the two tautomer isoforms (hydrogen on ND1 vs NE2); the PackerTask allows
	//swapping between those residue types even when only repacking, not designing.
	core::pack::task::PackerTaskOP task(core::pack::task::TaskFactory::create_packer_task((pose)));
	utility::vector1_bool packable(pose.total_residue(), false); //false = nobody is packable
	packable[ligand_res] = true;
	task->restrict_to_residues(packable); //now only one position is mobile
	task->restrict_to_repacking();
	task->nonconst_residue_task(ligand_res).or_fix_his_tautomer(true); //does nothing if not HIS(_D)
	task->nonconst_residue_task(ligand_res).or_ex1_sample_level(core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
	task->nonconst_residue_task(ligand_res).or_ex2_sample_level(core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
	task->nonconst_residue_task(ligand_res).or_ex3_sample_level(core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
	task->nonconst_residue_task(ligand_res).or_ex4_sample_level(core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
	TR << *task << std::endl;

	protocols::simple_moves::sidechain_moves::SidechainMoverOP SCmover( new protocols::simple_moves::sidechain_moves::SidechainMover() );
 	SCmover->set_task(task);
 	SCmover->set_prob_uniform(0); //we want only Dunbrack rotamers, 0 percent chance of uniform sampling
	SCmover->set_change_chi_without_replacing_residue(true); //moves everything downstream in foldtree rather than a pack_rotamers behavior

	/////////////////////////package movers into RandomMover (do either with equal chance)/////////////////
	protocols::moves::RandomMoverOP perturb_mover( new protocols::moves::RandomMover() );
	perturb_mover->add_mover(SCmover);
	//perturb_mover->add_mover(RJAmover);

	//////////////////////make centroid-izing mover and centroid pose//////////////////////////////
	core::pose::Pose centroid;
	copy_to_centroid(pose, centroid, centroid_tree, metal_res);
	TR << "centroid score before perturb: " << (*centroid_scorefunction_)(centroid) << std::endl;
	centroid_scorefunction_->show( TR, centroid );
	TR << std::flush; //show doesn't flush the buffer

	//////////////////////////////PDBDumpMover for debugging/////////////////////////////////////
	protocols::moves::PDBDumpMoverOP output_fa(new protocols::moves::PDBDumpMover("perturb_cycle_fa"));
	protocols::moves::PDBDumpMoverOP output_centroid(new protocols::moves::PDBDumpMover("perturb_cycle_centroid"));

	/////////////////////////make DualMonteCarlo///////////////////////////////////////////////////
	using basic::options::OptionKeys::AnchoredDesign::perturb_temp;
	using basic::options::option;
	protocols::moves::DualMonteCarlo DMC(
																			 pose, //pose for DMC
																			 centroid, //pose for MC
																			 *fullatom_scorefunction_, //scorefunction for DMC
																			 *centroid_scorefunction_, //scorefunction for MC
																			 option[perturb_temp].value()); //temperature for MC


	///////////////////////////////conformational search///////////////////////////////
	using basic::options::OptionKeys::AnchoredDesign::perturb_cycles;
	using basic::options::OptionKeys::AnchoredDesign::perturb_show;
	core::Size const perturb_cyc(option[ perturb_cycles ].value());
	bool const output_pdbs(option[ perturb_show ].value());
	TR << "Cycle\tDMC Current\tDMC Low\tMC Current\tMC Low" << std::endl;
	for( core::Size i(1); i <= perturb_cyc ; ++i){

		perturb_mover->apply(pose);
		TR << "pre-DMC  "; //two spaces to line up with post-DMC
		dump_chis(pose.residue(ligand_res));

		copy_to_centroid(pose, centroid, centroid_tree, metal_res);

		if( output_pdbs ){
			output_fa->apply(pose); //debugging!
			output_centroid->apply(centroid);
		}

		DMC.boltzmann(pose, centroid);
		if (!DMC.MC().mc_accepted()) TR << "DMC rejected" << std::endl;
		else TR << "DMC accepted" << std::endl;
		TR << "post-DMC ";
		dump_chis(pose.residue(ligand_res));

	}
	DMC.recover_low(pose, centroid);

	TR << "final centroid score after perturb: " << (*centroid_scorefunction_)(centroid) << std::endl;
	//centroid_scorefunction_->show( TR, centroid );
	//TR << std::flush; //show doesn't flush the buffer

	TR << "fullatom score after perturb: " << (*fullatom_scorefunction_)(pose) << std::endl;
	//fullatom_scorefunction_->show( TR, pose );
	//TR << std::flush; //show doesn't flush the buffer


	/////////////////////////////end perturb, begin fullatom design/refinement////////////////////

	////////////////////////////PackRotamersMover////////////////////////////////////////////
	protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
	pack_mover->task_factory( factory_ );
	pack_mover->score_function( fullatom_scorefunction_ );

	//////////////////////////////////generate minimizer mover/////////////////////////
	//movemap is empty, TAmin_mover will fill it
	core::kinematics::MoveMapOP map(new core::kinematics::MoveMap() );

	using protocols::simple_moves::MinMoverOP;
	using protocols::simple_moves::MinMover;
	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover(
																			map,
																			fullatom_scorefunction_,
																			option[ basic::options::OptionKeys::run::min_type ].value(),
																			0.01,
																			true /*use_nblist*/ );

	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	protocols::simple_moves::TaskAwareMinMoverOP TAmin_mover = new protocols::simple_moves::TaskAwareMinMover(min_mover, factory_);

	//let's examine the task
	//TR << *(factory_->create_task_and_apply_taskoperations( pose )) << std::endl;

	using basic::options::OptionKeys::AnchoredDesign::refine_cycles;
	//this says: use one cycle unless there is a command-line option saying otherwise
	core::Size const refine_applies( option[ refine_cycles ].user() ? option[ refine_cycles ].value() : 1 );
	for( core::Size cycle(1); cycle <= refine_applies; ++cycle){
		pack_mover->apply(pose);
		TAmin_mover->apply(pose);
	}

	TR << "fullatom score after refine: " << (*fullatom_scorefunction_)(pose) << std::endl;
	fullatom_scorefunction_->show( TR, pose );
	TR << std::flush; //show doesn't flush the buffer

	return;
}//apply

std::string
ZincHeterodimerMover::get_name() const {
	return "ZincHeterodimerMover";
}

///@details apply() needs to generate a centroid copy of the fullatom pose from time to time.  Some of the things in the fullatom pose are not centroid safe (metal atom, possible hydrogen-lacking metal ligand residues) and must be removed before centroid-ization.  This function takes care of that; unfortunately it needs a lot of help passed in.
void ZincHeterodimerMover::copy_to_centroid(
	core::pose::Pose const & pose,
	core::pose::Pose & centroid,
	core::kinematics::FoldTree const & centroid_tree,
	core::Size const)// metal_res)
{
	centroid = pose;
	centroid.fold_tree(centroid_tree); //"boring" fold tree is centroid safe (no jumps on nonexistent atoms)
	//centroid.conformation().delete_residue_slow(metal_res);
	core::util::switch_to_residue_type_set(centroid, core::chemical::CENTROID);
}


void ZincHeterodimerMover::generate_scorefunctions(){
	using namespace core::scoring;
	fullatom_scorefunction_ = core::scoring::get_score_function();
	TR << "Using default fullatom scorefunction (TALARIS_2013)\n"
		 << *fullatom_scorefunction_ << std::flush;

	centroid_scorefunction_ = new core::scoring::ScoreFunction;
	centroid_scorefunction_->set_weight( env,         2.0 );
	centroid_scorefunction_->set_weight( cbeta,       1.0 );
	centroid_scorefunction_->set_weight( vdw,         1.0 );
	centroid_scorefunction_->set_weight( pair,        1.0 );
	centroid_scorefunction_->set_weight( cenpack,     1.0 );
	//centroid_scorefunction_->set_weight( rama,        1.0 ); //won't matter with fixed backbones
	centroid_scorefunction_->set_weight( hbond_lr_bb, 1.0 );
	centroid_scorefunction_->set_weight( hbond_sr_bb, 1.0 );
	TR << "Using default centroid scorefunction\n" << *centroid_scorefunction_ << std::flush;

	return;
}//generate_scorefunctions

void ZincHeterodimerMover::generate_factory(){
	using namespace core::pack::task;
	using namespace basic::options;
	TaskFactoryOP task_factory = new TaskFactory();
	task_factory->push_back(operation::TaskOperationOP( new operation::InitializeFromCommandline() ));
	if ( option[ OptionKeys::packing::resfile ].user() ) {
		task_factory->push_back(operation::TaskOperationOP( new operation::ReadResfile ));
	}
	operation::PreventRepackingOP prop = new operation::PreventRepacking();
	for(utility::vector1< core::Size >::const_iterator it(metal_site_.begin()), end(metal_site_.end());	it != end; ++it)
		{prop->include_residue(*it);}
	task_factory->push_back(prop);
	//this assumes that the two protein partners are chains 1 and 3 - this is dangerous!!!!
	task_factory->push_back(operation::TaskOperationOP( new protocols::toolbox::task_operations::RestrictToInterfaceOperation(1, 3) ));

	TR << "using default TaskFactory (init from command line, read resfile, prevent repacking at metal site, detect interface" << std::endl;

	factory_ = task_factory; //store as COP so must wait for this
}//generate_factory

///@details constructor
ZincHeterodimerMover::ZincHeterodimerMover(
																										 utility::vector1< core::Size > const & metal_site,
																										 core::kinematics::Edge const & fixed_to_metal,
																										 core::kinematics::Edge const & metal_to_mobile )
	: Mover(), centroid_scorefunction_(NULL), fullatom_scorefunction_(NULL), factory_(NULL),
		fixed_to_metal_(fixed_to_metal), metal_to_mobile_(metal_to_mobile), metal_site_(metal_site)
{
	Mover::type( "ZincHeterodimerMover" );
	generate_scorefunctions();
	generate_factory();
}

ZincHeterodimerMover::~ZincHeterodimerMover(){}

}//metal_interface
}//protocols
