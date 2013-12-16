// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/analysis/InterfaceAnalyzerMover.cc
/// @brief InterfaceAnalyzerMover implementation - class for in-depth interface quality analysis
/// @author Steven Lewis, Bryan Der, Ben Stranges

// Unit Headers
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <protocols/analysis/InterfaceAnalyzerMoverCreator.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
//#include <core/scoring/Energies.hh>
//#include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh> //for centroid switch

// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>

#include <protocols/rigid/RigidBodyMover.hh>

#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MutateResidue.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/InterGroupNeighborsCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>

#include <core/scoring/packstat/compute_sasa.hh>
// AUTO-REMOVED #include <protocols/analysis/PackStatMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <basic/options/keys/run.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/packstat.OptionKeys.gen.hh>

// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/tag/Tag.hh>


// Utility Headers
#include <ObjexxFCL/FArray1D.hh> //necessary for fold tree tricks
#include <ObjexxFCL/FArray2D.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
//#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <sstream>
#include <set>
#include <string>

#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.analysis.InterfaceAnalyzerMover");
static basic::Tracer TRinterface("protocols.analysis.InterfaceAnalyzerMover.interface_selection");
static basic::Tracer TRhbonds("protocols.analysis.InterfaceAnalyzerMover.missing_hbonds");



///stupid helper function needed because ternary operator does not allow variable return types
std::ostream & which_ostream( std::ostream & ost, std::ostream & oss, bool const tracer){
  if(tracer) return ost;
  return oss;
}

namespace protocols{
namespace analysis{

using namespace protocols::moves;

std::string
InterfaceAnalyzerMoverCreator::keyname() const
{
	return InterfaceAnalyzerMoverCreator::mover_name();
}

protocols::moves::MoverOP
InterfaceAnalyzerMoverCreator::create_mover() const {
	return new InterfaceAnalyzerMover();
}

std::string
InterfaceAnalyzerMoverCreator::mover_name()
{
	return "InterfaceAnalyzerMover";
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////InterfaceAnalyzerMover///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
InterfaceAnalyzerMover::InterfaceAnalyzerMover(
	core::Size interface_jump,
	bool const tracer,
	core::scoring::ScoreFunctionCOP sf,
	bool compute_packstat,
	bool pack_input,
	bool pack_separated,
	bool use_jobname
) :
  Mover(),
  interface_jump_(interface_jump),
  sf_(sf),
  chain1_(0),
  chain2_(0),
  tracer_(tracer),
  calcs_ready_(false),
  compute_packstat_(compute_packstat),
	compute_interface_sc_(true), //would prefer a default of false
  compute_separated_sasa_(true),
	compute_interface_energy_(true),
	calc_hbond_sasaE_(true),
	compute_interface_delta_hbond_unsat_(true),
	skip_reporting_(false),
	multichain_constructor_(false),
  pack_input_(pack_input),
  pack_separated_(pack_separated),
  use_jobname_(use_jobname),
  use_resfile_(false),
  use_centroid_(false),
  ddG_(0),
  total_sasa_(0),
  interface_delta_sasa_(0),
  n_interface_res_(0),
  per_residue_energy_(0),
  complex_energy_(0),
  separated_interface_energy_(0),
  crossterm_interface_energy_(0),
  separated_interface_energy_ratio_(0),
  crossterm_interface_energy_ratio_(0),
  interface_packstat_(0),
  gly_dG_(0),
  centroid_dG_(0),
  delta_unsat_hbond_counter_(0),
  //pymol_sel_interface_,
  //pymol_sel_hbond_unsat_,
  //pymol_sel_packing_,
  side1_score_(0),
  side2_score_(0),
  nres_(0),
  side1_nres_(0),
  side2_nres_(0),
  //hbond_exposure_ratio_(0),
  //total_hb_sasa_(0),
  total_hb_E_(0),
	sc_value_(0.0),
  task_(0)//,
{
	protocols::moves::Mover::type( "InterfaceAnalyzer" );
	upstream_chains_.insert(0);
	downstream_chains_.insert(0);
}

//Alternate constructor for multichain poses
//Takes a set of ints for the chain nums
InterfaceAnalyzerMover::InterfaceAnalyzerMover(
	std::set<int> fixed_chains,
	bool const tracer,
	core::scoring::ScoreFunctionCOP sf,
	bool compute_packstat,
	bool pack_input,
	bool pack_separated,
	bool use_jobname
) :
  Mover(),
  interface_jump_(0),
  fixed_chains_(fixed_chains),
  sf_(sf),
  //posename_,
  //posename_base_,
  chain1_(0),
  chain2_(0),
  //chain1_char_;
  //chain2_char_;
  tracer_(tracer),
  calcs_ready_(false),
  compute_packstat_(compute_packstat),
	compute_interface_sc_(true), //would prefer a default of false
  compute_separated_sasa_(true),
	compute_interface_energy_(true),
	calc_hbond_sasaE_(true),
	compute_interface_delta_hbond_unsat_(true),
	skip_reporting_(false),
	multichain_constructor_(true),
  pack_input_(pack_input),
  pack_separated_(pack_separated),
  use_jobname_(use_jobname),
  use_resfile_(false),
  use_centroid_(false),
  ddG_(0),
  total_sasa_(0),
  interface_delta_sasa_(0),
  n_interface_res_(0),
  per_residue_energy_(0),
  complex_energy_(0),
  separated_interface_energy_(0),
  crossterm_interface_energy_(0),
  separated_interface_energy_ratio_(0),
  crossterm_interface_energy_ratio_(0),
  interface_packstat_(0),
  gly_dG_(0),
  centroid_dG_(0),
  delta_unsat_hbond_counter_(0),
  //pymol_sel_interface_,
  //pymol_sel_hbond_unsat_,
  //pymol_sel_packing_,
  side1_score_(0),
  side2_score_(0),
  nres_(0),
  side1_nres_(0),
  side2_nres_(0),
  //hbond_exposure_ratio_(0),
  //total_hb_sasa_(0),
  total_hb_E_(0),
	sc_value_(0.0),
  //interface_set_,
  //chain_groups_,
  task_(0)//,
  //Sasa_,
  //InterfaceNeighborDefinition_,
  //InterfaceSasaDefinition_,
  //InterfaceDeltaEnergetics_,
  //NumberHBonds_,
  //BuriedUnsatisfiedPolars_,
  //InterGroupNeighborsCalculator_,
{
	protocols::moves::Mover::type( "InterfaceAnalyzer" );
	upstream_chains_.insert(0);
	downstream_chains_.insert(0);
}


InterfaceAnalyzerMover::~InterfaceAnalyzerMover() {}

core::Real InterfaceAnalyzerMover::get_interface_ddG() const {return ddG_; } //previous functionality: redundant with get_separated_interface_energy, but supports other protocols

///@details getters
core::Real InterfaceAnalyzerMover::get_total_sasa() { return total_sasa_; }
core::Real InterfaceAnalyzerMover::get_interface_delta_sasa( ) { return interface_delta_sasa_; }
core::Real InterfaceAnalyzerMover::get_separated_interface_energy( ) { return separated_interface_energy_; }
core::Size InterfaceAnalyzerMover::get_num_interface_residues( ) {return n_interface_res_;}
core::Real InterfaceAnalyzerMover::get_complex_energy( ) { return complex_energy_; }
core::Real InterfaceAnalyzerMover::get_per_residue_energy( ) { return per_residue_energy_;}
core::Real InterfaceAnalyzerMover::get_crossterm_interface_energy( ) { return crossterm_interface_energy_; }
core::Real InterfaceAnalyzerMover::get_separated_interface_energy_ratio( ) { return separated_interface_energy_ratio_; }
core::Real InterfaceAnalyzerMover::get_crossterm_interface_energy_ratio( ) { return crossterm_interface_energy_ratio_; }
core::Real InterfaceAnalyzerMover::get_interface_packstat() { return interface_packstat_; }
core::Size InterfaceAnalyzerMover::get_interface_delta_hbond_unsat() { return delta_unsat_hbond_counter_; }
std::string InterfaceAnalyzerMover::get_pymol_sel_interface() { return pymol_sel_interface_; }
std::string InterfaceAnalyzerMover::get_pymol_sel_hbond_unsat() { return pymol_sel_hbond_unsat_; }
std::string InterfaceAnalyzerMover::get_pymol_sel_packing() { return pymol_sel_packing_; }
bool InterfaceAnalyzerMover::get_multichain_constructor() {return multichain_constructor_; }
std::set<int> InterfaceAnalyzerMover::get_fixed_chains() {return fixed_chains_; }
std::set<core::Size> InterfaceAnalyzerMover::get_interface_set() {return interface_set_;}
InterfaceAnalyzerMover::group_set InterfaceAnalyzerMover::get_chain_groups () { return chain_groups_; }
bool InterfaceAnalyzerMover::get_pack_input(){return pack_input_; }
core::pack::task::PackerTaskOP InterfaceAnalyzerMover::get_packer_task(){return task_ ;}
core::Real InterfaceAnalyzerMover::get_gly_interface_energy(){ return gly_dG_ ; }
core::Real InterfaceAnalyzerMover::get_centroid_dG(){ return centroid_dG_ ; }
//core::Real InterfaceAnalyzerMover::get_interface_Hbond_sasa() {return total_hb_sasa_; }
//core::Real InterfaceAnalyzerMover::get_Hbond_exposure_ratio() {return hbond_exposure_ratio_; }
core::Real InterfaceAnalyzerMover::get_total_Hbond_E(){return total_hb_E_;}

bool InterfaceAnalyzerMover::get_use_resfile() const {return use_resfile_;}
bool InterfaceAnalyzerMover::get_use_centroid_dG() const {return use_centroid_;}

///@details setters
void InterfaceAnalyzerMover::set_use_resfile(bool const use_resfile) {use_resfile_ = use_resfile;}
void InterfaceAnalyzerMover::set_use_centroid_dG(bool const use_centroid) {use_centroid_ = use_centroid;}
void InterfaceAnalyzerMover::set_compute_packstat(bool const compute_packstat) {compute_packstat_ = compute_packstat;}
void InterfaceAnalyzerMover::set_pack_input(bool const pack_input) {pack_input_ = pack_input;}
void InterfaceAnalyzerMover::set_interface_jump(core::Size const interface_jump) {interface_jump_ = interface_jump;}
void InterfaceAnalyzerMover::set_tracer(bool const tracer) {tracer_ = tracer;}
//void InterfaceAnalyzerMover::set_calcs_ready(bool const calcs_ready) {calcs_ready_ = calcs_ready;}
void InterfaceAnalyzerMover::set_use_jobname(bool const use_jobname) {use_jobname_ = use_jobname;}
void InterfaceAnalyzerMover::set_pack_separated(bool const pack_separated) {pack_separated_ = pack_separated;}


///@details InterfaceAnalyzerMover computes various interface statistics and makes them available through getters
void InterfaceAnalyzerMover::apply( core::pose::Pose & pose )
{
	using namespace core;
	if( pose.conformation().num_chains() < 2 ){
		utility_exit_with_message_status( "InterfaceAnalyzerMover: pose has only one chain, aborting analysis \n", 1 );
		//remaining code works for actual interfaces
	}
	
	//If we've specified a ligand chainn, then every other chain is the fixed chain
	//This isn't a very efficient way of figuring this out but this mover is full of
	//fairly slow stuff so another trip through the residues isn't going to kill anyone
	if(ligand_chain_.size() != 0)
	{
		char this_chain (ligand_chain_[0]);
		for (core::Size i = 1; i<=pose.total_residue(); ++i){
			if (pose.pdb_info()->chain( i ) != this_chain){
				fixed_chains_.insert( pose.chain(i) );
			}
		}
		
	}
	//check for multichain poses
	if(pose.conformation().num_chains() > 2){
		if (multichain_constructor_){
			

			//fix the foldtree to reflect the fixed chains we want
			core::Size newjump;
			std::set< int > fixedchains( get_fixed_chains() );
			newjump = reorder_foldtree_find_jump( pose, fixedchains );
			interface_jump_ = newjump ;
			set_pose_info( pose );
			//need calculators to register here if not already around
			//make_multichain_interface_set(pose, fixedchains);
			if(!calcs_ready_){
				register_calculators();
				calcs_ready_ = true;
			}
		}
		//if multi chains with wrong constructor, work but print a warming
		else{
			TR<< "WARNING: more than 2 chains present w/o using the right constructor!  Values might be over the wrong jump." << std::endl;
			//sets up the pose for calculations
			set_pose_info( pose );
			//register calculators here if need be
			if(!calcs_ready_){
				register_calculators();
				calcs_ready_ = true;
			}
		}
	}//end if greater than 2 chains

	//if only 2 chains
	else {
		//sets up the pose for calculations
		//If the user has specified "fixedchains" and provided a 2 chain pose then the interface jump hasn't been set properly.
		//The interface jump ought to be 1, but we'll just do it properly and get it from the chain id
		if(interface_jump_ == 0)
		{
			std::set< int > fixedchains( get_fixed_chains() );
			if(fixedchains.size() != 1)
			{
				utility_exit_with_message_status( "Pose only has 2 chains, but 'fixedchains' option in Interface analyzer specified more than one chain \n", 1 );
			}
			
			core::Size newjump = reorder_foldtree_find_jump( pose, fixedchains );
			interface_jump_ = newjump ;
		}

		
		set_pose_info( pose );
		//register calculators here if need be
		if(!calcs_ready_){
			register_calculators();
			calcs_ready_ = true;
		}
	}// end if only two chains

	//setup an interface set
	if(interface_set_.empty())  //should always be empty at this point
		make_interface_set(pose);
	//set up the packer task for later
	setup_task(pose);
	//init compexed and separated pose objects
	core::pose::Pose complexed_pose( pose );
	core::pose::Pose separated_pose( make_separated_pose( complexed_pose ) );
	//actual computation here
	if(compute_separated_sasa_) compute_separated_sasa( complexed_pose, separated_pose );
	if(compute_interface_energy_) compute_interface_energy( complexed_pose, separated_pose );
	//mut_to_gly(complexed_pose, separated_pose );  this didn't help
	if(use_centroid_) calc_centroid_dG( complexed_pose, separated_pose );
	if(calc_hbond_sasaE_) calc_hbond_sasaE( complexed_pose ); //get hbond E at interface
	if( compute_packstat_) compute_interface_packstat( complexed_pose );
	//find the change in interface hbonds
	if(compute_interface_delta_hbond_unsat_) compute_interface_delta_hbond_unsat( complexed_pose, separated_pose );
	//find the shape complementarity stats for the interface
	if( compute_interface_sc_ ) compute_interface_sc(interface_jump_, complexed_pose);

	if( !skip_reporting_ ){
		//always will fill a selection option to get the pymol selection
		print_pymol_selection_of_interface_residues( complexed_pose, interface_set_);

		//report to tracer or job all this cool stuff we calculated
		report_data();
	}

	//why on earth are we deleting this set?  Nobody can use it if we delete it... SML 11/2/12
	// //get the interface set again setup some output and clear it.
	// if( ! (interface_set_.empty()) ){
	// 	//TR << "Total interface residues for " << posename_base_ << ":  " << interface_set_.size() << std::endl;
	// 	//clear any old values before return
	// 	interface_set_.clear();
	// 	//set_interface_set( interface_set_ );
	// }


	return;
}//end apply

std::string
InterfaceAnalyzerMover::get_name() const {
	return "InterfaceAnalyzerMover";
}

///@details sets up the pose information such as the name and chain ids
void InterfaceAnalyzerMover::set_pose_info( core::pose::Pose const & pose ) {
	if( use_jobname_ ){
		protocols::jd2::JobOP current_job(protocols::jd2::JobDistributor::get_instance()->current_job());
		posename_base_ = protocols::jd2::JobDistributor::get_instance()->job_outputter()->output_name( current_job );
	}
	else {
		posename_ = pose.pdb_info()->name();
		posename_base_ = posename_.base();
	}

	chain1_ = pose.residue(pose.fold_tree().upstream_jump_residue(interface_jump_)).chain();
	chain2_ = pose.residue(pose.fold_tree().downstream_jump_residue(interface_jump_)).chain();
	chain1_char_ = pose.pdb_info()->chain(pose.conformation().chain_end(chain1_));
	chain2_char_ = pose.pdb_info()->chain(pose.conformation().chain_end(chain2_));
	//set here_ will be later replaced by multichain constructor if need be
	upstream_chains_.clear();
	downstream_chains_.clear();
	upstream_chains_.insert(chain1_);
	downstream_chains_.insert( chain2_);
}

///@details makes the interface sets for either constructor
//always will make a new interface set, make sure this is what you want
void InterfaceAnalyzerMover::make_interface_set( core::pose::Pose & pose ){
	if(multichain_constructor_){
		make_multichain_interface_set(pose, fixed_chains_);
	}
	//std::set< core::Size > interface_set ( get_interface_set() );
	else {
		TR << "Making an interface set with default calculator." << std::endl;
		basic::MetricValue< std::set< core::Size > > mv_interface_set;
		pose.metric(InterfaceNeighborDefinition_, "interface_residues", mv_interface_set);
		interface_set_ = mv_interface_set.value();
		//set_interface_set( interface_set_ );
	}
	if ( interface_set_.empty() )
		TR << "NO INTERFACE FOR: " << posename_base_ << std::endl;
}

///@details reorder the fold tree to allow multichain interfaces to be evaluated returns the new chain for the jump
core::Size InterfaceAnalyzerMover::reorder_foldtree_find_jump( core::pose::Pose & pose ,
	std::set<int> & fixed_chains){
	using namespace core;
	using core::kinematics::Edge;
	core::Size mobile_jump(1);  //setup, will be changed later
	if(fixed_chains.empty())
		utility_exit_with_message_status( "Can't find fixed chains.  Exiting...\n", 1 );
	//now get info about chains
	Size numchains ( pose.conformation().num_chains() );
	utility::vector1<Size> chain_starts;
	utility::vector1<Size> chain_ends;
	for(Size ii = 1; ii<= numchains; ++ii){
		chain_starts.push_back( pose.conformation().chain_begin( ii ) );
		chain_ends.push_back( pose.conformation().chain_end( ii ) );
	}

	//find a non mobile chain
	Size anchor_chain (1); //switch later if needed
	for(Size ii = 1; ii<= numchains; ++ii){
		if ( !fixed_chains.count( ii ) ){
			anchor_chain = ii;
			break;
		}
	}

	//now build our fold tree
	//this will move the mobile chains together
	kinematics::FoldTree foldtree (pose.fold_tree());
	Size this_jump(1);
	foldtree.clear();
	bool previous_mobile (false);
	for( Size ii=1; ii<= numchains; ++ii){
		foldtree.add_edge( Edge(chain_starts[ii], chain_ends[ii], kinematics::Edge::PEPTIDE) );
		if (ii != anchor_chain) {
			if ( fixed_chains.count(ii) && previous_mobile  ){ //not in the first mobile chain
				foldtree.add_edge( Edge(chain_starts[ii], chain_ends[*(fixed_chains.begin())], this_jump) );
			}
			else if ( fixed_chains.count(ii) ){ //ie in the first mobile chain
				foldtree.add_edge( Edge(chain_starts[ii], chain_ends[anchor_chain], this_jump) );
				mobile_jump = this_jump; //sets this jump to be the new mobile one
				previous_mobile=true;
			}
			else foldtree.add_edge( Edge(chain_starts[ii], chain_ends[anchor_chain], this_jump) );
			++this_jump;
		}
	}
	//now set all the new values!
	foldtree.reorder(chain_starts[anchor_chain]);
	pose.fold_tree(foldtree);
	//TR << "Remade foldtree:\n"<< foldtree << std::endl;

	return mobile_jump;

}//end reorder foldtree

// sets the interface set for a multichain pose
// uses everything in the fixed chains as one side, everything else is the other side
void InterfaceAnalyzerMover::make_multichain_interface_set( core::pose::Pose & pose,
	std::set<int> & fixed_chains){
	using namespace core;
	using namespace utility;
	std::set<Size> fixed_side_residues, other_side_residues;
	std::set<core::Size> fixed_chain_nums ;
	std::set<core::Size> other_chain_nums;
	//itterate over all residues determine what part of the interface they are
	//also select what chain(s) are upstream and downstream of the interface
	for( Size ii = 1; ii<= pose.total_residue(); ++ii){
		if( fixed_chains.count( pose.chain( ii ) ) ){
			fixed_side_residues.insert( ii );
			fixed_chain_nums.insert( pose.chain( ii ) );
		}
		else{
			other_side_residues.insert( ii );
			other_chain_nums.insert( pose.chain( ii ) );
		}
	}
	//now assign the correct chains
	upstream_chains_ = fixed_chain_nums;
	downstream_chains_= other_chain_nums;
	//debugging
	//TR << "Fixed residues: " << fixed_side_residues.size() << "   Other side residues: " << other_side_residues.size() << std::endl;

	//prep a vector of a pair of these residue sets for Steven's calculator
	std::pair< std::set<Size>, std::set<Size> > side_pairs;
	side_pairs.first = fixed_side_residues;
	side_pairs.second = other_side_residues;
	InterfaceAnalyzerMover::group_set interface_deffinition_vector;
	interface_deffinition_vector.push_back( side_pairs );
	chain_groups_ = interface_deffinition_vector ;

	//get calculator for multichain poses
	register_intergroup_calculator();

	//std::set<Size> multichain_interface;
	basic::MetricValue< std::set<Size> > mv_interface_set;
	pose.metric(  InterGroupNeighborsCalculator_, "neighbors", mv_interface_set);
	//set_interface_set( mv_interface_set.value() );
	interface_set_ =  mv_interface_set.value();
	//debugging
	//TR << "Interface set residues total: " << interface_set_.size() << std::endl;

}//end make_multichain_interface_set


///@details  makes the complexed and separated poses for the mover
core::pose::Pose InterfaceAnalyzerMover::make_separated_pose( core::pose::Pose & pose ){
	using namespace core;

	if(!sf_) {
		TR << "NULL scorefunction.  Initialize from cmd line." << std::endl;
		using namespace core::scoring;
		sf_ = getScoreFunction();
	}
	(*sf_)(pose); //shits, giggles, and segfault prevention
	//complexed_pose_ = pose;
	pose::Pose separated_pose( pose );
	protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( separated_pose, interface_jump_ ) );
	translate->step_size( 1000.0 );
	translate->apply( separated_pose );
	(*sf_)(separated_pose);
	//debugging step to make sure the jump is right
	//separated_pose_.dump_pdb("IAM_test.pdb");
	return separated_pose;
}


///@details computes the SASA by finding difference between complex and separated SASA
/// also does the same thing for hydrophobic/polar SASA
void InterfaceAnalyzerMover::compute_separated_sasa(core::pose::Pose & complexed_pose,
	core::pose::Pose & separated_pose){
	using namespace core;
	basic::MetricValue< core::Real > mv_complex_sasa;
	basic::MetricValue< core::Real > mv_separated_sasa;

	complexed_pose.metric(Sasa_, "total_sasa", mv_complex_sasa);
	separated_pose.metric(Sasa_, "total_sasa", mv_separated_sasa);
	Real const total_sasa_value (mv_complex_sasa.value());
	Real const separated_sasa_value ( mv_separated_sasa.value() );
	total_sasa_ = total_sasa_value ;
	interface_delta_sasa_ = separated_sasa_value - total_sasa_value ;

	//need to get the hydrophobic sasa at the interface, only need to subtract them
	utility::vector1< core::Real > complexed_residue_sasa( complexed_pose.total_residue(), 0.0 );
	utility::vector1< core::Real > separated_residue_sasa( separated_pose.total_residue(), 0.0 );
	utility::vector1< core::Real > complexed_residue_hsasa( complexed_pose.total_residue(), 0.0 ); // hydrophobic SASA only
	utility::vector1< core::Real > separated_residue_hsasa( separated_pose.total_residue(), 0.0 ); // hydrophobic SASA only
	core::Real probe_radius = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius];
	core::Real complexed_hSASA = core::scoring::calc_per_res_hydrophobic_sasa( complexed_pose, complexed_residue_sasa,
		complexed_residue_hsasa, probe_radius,
		false /*no n_acccess_radii*/  );
	core::Real separated_hSASA = core::scoring::calc_per_res_hydrophobic_sasa( separated_pose, separated_residue_sasa,
		separated_residue_hsasa, probe_radius,
		false /*no n_acccess_radii*/ );

	//Now assume anything that isn't hydrophobic is polar
	interface_hsasa_ = separated_hSASA - complexed_hSASA;
	interface_polar_sasa_ = interface_delta_sasa_ - interface_hsasa_;
	return;
}


///@details computes the interface energy of the interface
void InterfaceAnalyzerMover::compute_interface_energy( core::pose::Pose & complexed_pose,
	core::pose::Pose & separated_pose )
{
	//separated interface energy and ratio
	score_separated_chains(complexed_pose, separated_pose);  //sets separated_interface_energy_
	separated_interface_energy_ratio_ = separated_interface_energy_ / interface_delta_sasa_ ;
	calc_per_residue_energy( complexed_pose ); //find avg residue E at interface
	//crossterm interface energy and ratio
	basic::MetricValue< core::Real > mv_delta_total;
	complexed_pose.metric(InterfaceDeltaEnergetics_, "weighted_total", mv_delta_total);
	crossterm_interface_energy_ = mv_delta_total.value() ;
	crossterm_interface_energy_ratio_ = crossterm_interface_energy_ / interface_delta_sasa_ ;

	// Get the first and last resnum of each chain, using name_of_InterfaceNeighborDefinitionCalculator_
	basic::MetricValue<Size> mv_size;
	separated_pose.metric(InterfaceNeighborDefinition_,"first_chain_first_resnum",mv_size);
	core::Size ch1_begin_num = mv_size.value();
	separated_pose.metric(InterfaceNeighborDefinition_,"first_chain_last_resnum",mv_size);
	core::Size ch1_end_num = mv_size.value();
	separated_pose.metric(InterfaceNeighborDefinition_,"second_chain_first_resnum",mv_size);
	core::Size ch2_begin_num = mv_size.value();
	separated_pose.metric(InterfaceNeighborDefinition_,"second_chain_last_resnum",mv_size);
	core::Size ch2_end_num = mv_size.value();

	side1_score_ = 0;
	side1_nres_ = ch1_end_num - ch1_begin_num+1;
	for( core::Size i(ch1_begin_num); i<= ch1_end_num; ++i) {
		side1_score_ += separated_pose.energies().residue_total_energy(i);
	}

	side2_score_ = 0;
	side2_nres_ = ch2_end_num - ch2_begin_num+1;
	for( core::Size i(ch2_begin_num); i<= ch2_end_num; ++i) {
		side2_score_ += separated_pose.energies().residue_total_energy(i);
	}

	nres_ = complexed_pose.total_residue();

	return;
}
///@details actual function to separate the chains based on the chosen jump and score
void InterfaceAnalyzerMover::score_separated_chains( core::pose::Pose & complexed_pose,
	core::pose::Pose & separated_pose ){
	//using namespace core::pack::task;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::moves;

	//core::pose::Pose copy(pose);
	//PackerTaskOP const task( get_packer_task() );
	//we know there is a score function from above
	core::Real complex_score( (*sf_)(complexed_pose) );
	core::Real separated_score( (*sf_)(separated_pose) );
	complex_energy_ = complex_score ;

	// TR<< "Initial complex_score: " << complex_score << "  sep_score: "
	// 		<< separated_score << std::endl;

	//setup mover to pack interface
	//core::Size ndruns = option[packing::ndruns];
	protocols::simple_moves::PackRotamersMoverOP repacker = new protocols::simple_moves::PackRotamersMover(sf_, task_ );

	//do we pack the complex?
	if ( pack_input_ ){
		repacker->apply(complexed_pose);
		complex_score = (*sf_)(complexed_pose);
		complex_energy_ = complex_score ;
	}
	//do we pack the separated pose
	if (pack_separated_){
		repacker->apply( separated_pose );
		separated_score = (*sf_)(separated_pose) ;
	}
	//debugging step to make sure the jump is right
	//copy.dump_pdb("IAM_test.pdb");
	// 	TR<< "Post Packing complex_score: " << complex_score << "  sep_score: "
	// 		<< separated_score << std::endl;
	ddG_ = complex_score - separated_score;
	separated_interface_energy_ = ddG_ ; //these are the same thing
}

///@details reports all the cool stuff we calculate to tracer output OR puts it into the job object.
void InterfaceAnalyzerMover::report_data(){
	//make output
	protocols::jd2::JobOP current_job(protocols::jd2::JobDistributor::get_instance()->current_job());
	//std::string const posename(protocols::jd2::JobDistributor::get_instance()->job_outputter()->output_name( current_job ) );
	//std::ostringstream results_oss;
	//std::ostream & results = which_ostream(T(posename_base_), results_oss, tracer_); //easy swap between tracer/job output

	//std::ostringstream interface_sele, missingHbond;

	//set up what data is reported in a string stream
	//report to job
	if (tracer_){
		//TR<<"Debugging print interface info:" << std::endl;
		T(posename_base_) << "TOTAL SASA: " << total_sasa_ << std::endl;
		T(posename_base_) << "NUMBER OF RESIDUES: " << n_interface_res_ << std::endl;
		T(posename_base_) << "AVG RESIDUE ENERGY: " << per_residue_energy_ << std::endl;
		T(posename_base_) << "INTERFACE DELTA SASA: " << interface_delta_sasa_ << std::endl;
		T(posename_base_) << "INTERFACE HYDROPHOBIC SASA: " << interface_hsasa_ << std::endl;
		T(posename_base_) << "INTERFACE POLAR SASA: " << interface_polar_sasa_ << std::endl;
		T(posename_base_) << "CROSS-INTERFACE ENERGY SUMS: " << crossterm_interface_energy_ << std::endl;
		T(posename_base_) << "SEPARATED INTERFACE ENERGY DIFFERENCE: " << separated_interface_energy_ << std::endl;
		T(posename_base_) << "INTERFACE DELTA SASA/CROSS-INTERFACE ENERGY: " << crossterm_interface_energy_ratio_ << std::endl;
		T(posename_base_) << "INTERFACE DELTA SASA/SEPARATED INTERFACE ENERGY: " << separated_interface_energy_ratio_ << std::endl;
		T(posename_base_) << "DELTA UNSTAT HBONDS: " << delta_unsat_hbond_counter_ << std::endl;
		//T(posename_base_) << "ALL Gly INTERFACE ENERGY:  " << gly_dG_ << std::endl; does not help
		if(use_centroid_)
			T(posename_base_) << "CENTORID dG: " << centroid_dG_ << std::endl;
		//T(posename_base_) << "AVG HBOND EXPOSURE RATIO: " << hbond_exposure_ratio_ << std::endl; does not help
		//T(posename_base_) << "HBOND SASA / INTERFACE dSASA: " << total_hb_sasa_ / interface_delta_sasa_ << std::endl;
		T(posename_base_) << "HBOND ENERGY: " << total_hb_E_ << std::endl;
		T(posename_base_) << "HBOND ENERGY/ SEPARATED INTERFACE ENERGY: " << total_hb_E_ / separated_interface_energy_ << std::endl;
		if (compute_packstat_)
			T(posename_base_) << "INTERFACE PACK STAT: " << interface_packstat_ << std::endl;
		if (compute_interface_sc_)
			T(posename_base_) << "SHAPE COMPLEMENTARITY VALUE: " << sc_value_ << std::endl;
		//results <<  pymol_sel_interface_;
		//	results <<  pymol_sel_hbond_unsat_;
		//if(compute_packstat_)
		//	results <<  pymol_sel_packing_;

	}
	//or report to job
	else{
		current_job->add_string_real_pair("dSASA_int", interface_delta_sasa_);
		current_job->add_string_real_pair("dSASA_polar", interface_polar_sasa_);
		current_job->add_string_real_pair("dSASA_hphobic", interface_hsasa_);
		current_job->add_string_real_pair("dG_separated", separated_interface_energy_);
		current_job->add_string_real_pair("dG_separated/dSASAx100", separated_interface_energy_ratio_*100.0 );//purpose of 1000 is to move decimal point to make more significant digits appear in fixed-precision scorefile (0.022 is only 2 digits, 2.234 is more useful)
		current_job->add_string_real_pair("delta_unsatHbonds", delta_unsat_hbond_counter_ );
		current_job->add_string_real_pair("packstat", interface_packstat_ );
		current_job->add_string_real_pair("dG_cross", crossterm_interface_energy_ );
		current_job->add_string_real_pair("dG_cross/dSASAx100", crossterm_interface_energy_ratio_*100.0);//as above
		//current_job->add_string_real_pair("AllGly_dG", gly_dG_ ); does not help
		if(use_centroid_)
			current_job->add_string_real_pair("cen_dG", centroid_dG_ );
		current_job->add_string_real_pair("nres_int", n_interface_res_ );
		current_job->add_string_real_pair("per_residue_energy_int", per_residue_energy_ );
		current_job->add_string_real_pair("side1_score", side1_score_ );
		current_job->add_string_real_pair("side2_score", side2_score_ );
		current_job->add_string_real_pair("nres_all", nres_ );
		current_job->add_string_real_pair("side1_normalized", (side1_score_/(core::Real(side1_nres_))) );
		current_job->add_string_real_pair("side2_normalized", (side2_score_/(core::Real(side2_nres_))) );
		current_job->add_string_real_pair("complex_normalized", complex_energy_/(core::Real(nres_)) );
		current_job->add_string_real_pair("hbond_E_fraction", total_hb_E_/separated_interface_energy_);
		current_job->add_string_real_pair("sc_value", sc_value_);

	}

}

///@details Averaged packstat of interface residues only
void InterfaceAnalyzerMover::compute_interface_packstat( core::pose::Pose & pose )
{

	TR << "Computing interface packstats." << std::endl;
	//check if there is already an interface set, if so use it, if not, make own
	//will make own using the basic calculator
	//std::set< core::Size > interface_set ( get_interface_set() );
	if( interface_set_.empty() )
		make_interface_set( pose );

	//calculates the packstat scores for the interface set
	utility::vector1< core::Real > interface_pack_scores;
	core::Real interface_pack_score_sum (0.0);
	core::Size interface_num_res (0);
	interface_pack_scores =  core::scoring::packstat::compute_residue_packing_scores( pose, basic::options::option[ basic::options::OptionKeys::packstat::oversample ]() );

	for( std::set< Size >::const_iterator it( interface_set_.begin() ), end( interface_set_.end());
			 it != end; ++it){
		interface_pack_score_sum += interface_pack_scores[*it];
		++interface_num_res;
	}
	//fills the selection for pymol output, doesn't print anything here
	print_pymol_selection_of_packing( pose, interface_pack_scores );
	interface_packstat_ = interface_pack_score_sum / interface_num_res ;
}

///@details If a polar atom at the interface is also "buried unsat" in the monomer, we don't count this one
void InterfaceAnalyzerMover::compute_interface_delta_hbond_unsat( core::pose::Pose & complexed_pose,
	core::pose::Pose & separated_pose)
{

	//already should be an interface set but check anyway
	if(interface_set_.empty() )
		make_interface_set(complexed_pose);
	//do the calculations
	TR << "Computing delta unsat polar residues..." << std::endl;


	//core::pose::Pose complexed_pose_(pose);
	//core::pose::Pose separated_pose_(pose);
	//protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( separated_pose_, interface_jump_ ) );
	//translate->step_size( 1000.0 );
	//translate->apply( separated_pose_ );

	basic::MetricValue< utility::vector1< core::Size > > mv_complexed_unsat_res;
	basic::MetricValue< utility::vector1< core::Size > > mv_separated_unsat_res;
	basic::MetricValue< core::id::AtomID_Map< bool > > mv_complexed_unsat_map;
	basic::MetricValue< core::id::AtomID_Map< bool > > mv_separated_unsat_map;

	complexed_pose.metric(BuriedUnsatisfiedPolars_, "residue_bur_unsat_polars", mv_complexed_unsat_res);
	separated_pose.metric(BuriedUnsatisfiedPolars_, "residue_bur_unsat_polars", mv_separated_unsat_res);
	utility::vector1< core::Size > const complexed_unsat_res(mv_complexed_unsat_res.value());
	utility::vector1< core::Size > const separated_unsat_res(mv_separated_unsat_res.value());

	complexed_pose.metric(BuriedUnsatisfiedPolars_, "atom_bur_unsat", mv_complexed_unsat_map);
	separated_pose.metric(BuriedUnsatisfiedPolars_, "atom_bur_unsat", mv_separated_unsat_map);
	core::id::AtomID_Map< bool > const complexed_unsat_map(mv_complexed_unsat_map.value());
	core::id::AtomID_Map< bool > const separated_unsat_map(mv_separated_unsat_map.value());

	//loop over the interface set and figure out what's burried/unsat
	core::Size delta_unsat_hbond_counter( 0 );
	utility::vector1< core::id::AtomID > delta_unsat_hbond_atid_vector;
	for( std::set< core::Size >::const_iterator it(interface_set_.begin()), end(interface_set_.end());
			 it != end; ++it){
		//TR << "UnsatHbond res " << *it << std::endl;
		//iterate over all its atoms, check if they're in the map of missing, and print their string name
		core::chemical::ResidueType const & res(complexed_pose.residue_type(*it));
		for( core::Size i=1 ; i <= res.natoms(); ++i ){
			core::id::AtomID const atid(i, *it);

			//TR << "UnsatHbond atom " << i << std::endl;

			//if atom is buried unsat in complexed pose but not separated pose, count it as delta unsat
			bool unsat_complex(complexed_unsat_map.has(atid) && complexed_unsat_map(atid) );
			bool unsat_separated(separated_unsat_map.has(atid) && separated_unsat_map(atid) );

			if( unsat_complex && !unsat_separated ) {
				//results << " " << res.atom_name(i) ;
				delta_unsat_hbond_atid_vector.push_back(atid);
				++delta_unsat_hbond_counter;
			}
		}//end for each atom
	}//end loop of interface_set_

	delta_unsat_hbond_counter_ = delta_unsat_hbond_counter ; //sets the total
	//print_pymol_selection_of_interface_residues( pose, interface_set_);
	//calculated here but need to call get the selection to output it
	print_pymol_selection_of_hbond_unsat( complexed_pose, delta_unsat_hbond_atid_vector );

	return;
}


///@details prints tracer output of pymol selction of interface residues, also builds a pymol selection that can be used from a file.
void InterfaceAnalyzerMover::print_pymol_selection_of_interface_residues( core::pose::Pose const & pose, std::set< core::Size > const interface_set )
{

	//make output
	protocols::jd2::JobOP current_job(protocols::jd2::JobDistributor::get_instance()->current_job());
	//for tracer or job output
	std::ostringstream  interface_oss;
	std::ostream & pymol_interface = which_ostream(TRinterface, interface_oss, tracer_);

	//setup naming of pymol objects
	std::ostringstream interface_sele;
	interface_sele << std::endl;
	std::string pymol_obj_for_interface_sel;
	if( compute_packstat_ ) {
		pymol_obj_for_interface_sel = posename_base_ + "_fullpose_pack";
	}
	else {
		pymol_obj_for_interface_sel = posename_base_;
	}
	//setup the tracer output
	pymol_interface << "pymol-style selection for interface res \n"
	<< "select " << posename_base_ << "_interface, ";
	//itterate through the interface set and build the selection syntaxt
	bool first_sel_complete( false );
	core::Size resnum;
	char chain_char, chain_char_last('z');
	for( std::set< core::Size >::const_iterator it(interface_set.begin()), end(interface_set.end());
			 it != end; ++it){
		//sets the current values
		resnum = pose.pdb_info()->number(*it);
		chain_char = pose.pdb_info()->chain(*it);
		//special print if the first time through
		if( !first_sel_complete ) {
			interface_sele << "cmd.select(\"/" << pymol_obj_for_interface_sel << "//" << chain_char << "/"
			<< resnum << "\")" << std::endl;
			pymol_interface << "/" << posename_base_ << "//" << chain_char << "/" << resnum << "+" ;
			first_sel_complete = true;
		}
		else if (chain_char != chain_char_last){
			interface_sele << "cmd.select(\"sele + /" << pymol_obj_for_interface_sel << "//" << chain_char << "/"
			<< resnum << "\")" << std::endl;
			pymol_interface <<" + "<< "/" << posename_base_ << "//" << chain_char << "/" << resnum << "+";
		}
		else {
			interface_sele << "cmd.select(\"sele + /" << pymol_obj_for_interface_sel << "//" << chain_char << "/"
			<< resnum << "\")" << std::endl;
			pymol_interface << resnum << "+";
		}
		chain_char_last = chain_char;
	} //end itterate over interface
    //finish up
	pymol_interface << std::endl;
	interface_sele << "cmd.create(\"" << posename_base_ << "_interface_sel\", \"sele\")" << std::endl;
	pymol_sel_interface_ =  interface_sele.str() ;
	//job output if wanted
	if(!tracer_)
		current_job->add_string( interface_oss.str() );
	return;

}//end


///@details This function reports a few things: a pymol sytle selection of the unstat atoms and reports to the tracer or job what these atoms are.  The app InterfaceAnalyzer gets the multi-line string to write a file or print the selection.  Unsat hbonds to be shown as Spheres
void InterfaceAnalyzerMover::print_pymol_selection_of_hbond_unsat( core::pose::Pose & pose, utility::vector1< core::id::AtomID > delta_unsat_hbond_atid_vector )
{
	//make output
	protocols::jd2::JobOP current_job(protocols::jd2::JobDistributor::get_instance()->current_job());
	//for tracer or job output
	std::ostringstream results_oss, unsathbond_oss;
	std::ostream & results = which_ostream(T(posename_base_), results_oss, tracer_); //easy swap between tracer/job output
	std::ostream & unsathbond = which_ostream(TRhbonds, unsathbond_oss, tracer_);
	results << "Residues missing H-bonds:" << std::endl;
	results << "Residue \t Chain \t Atom " << std::endl;
	bool first_sel_complete( false );
	std::ostringstream missingHbond; // for pymol selection
	missingHbond << std::endl;

	//setup the tracer output
	unsathbond << "pymol-style selection for unstat hbond res \n"
	<< "select " << posename_base_ << "_unsat, ";
	//setup for looping over all unstat hbonds
	core::Size resnum;
	char chain_char, chain_char_last ('z');
	std::string atomname ;
	for( core::Size i(1); i <= delta_unsat_hbond_atid_vector.size(); i++ ) {
		core::id::AtomID const id ( delta_unsat_hbond_atid_vector[ i ] );
		resnum = pose.pdb_info()->number(id.rsd());
		chain_char = pose.pdb_info()->chain( id.rsd() );
		atomname = pose.residue(id.rsd()).atom_name(id.atomno());
		//get rid of whitespace in the atomname
		std::string temp;
		for (unsigned int j = 0; j < atomname.length(); j++) {
			if (atomname[j] != ' ') { temp += atomname[j]; }
		}
		atomname = temp;
		//do the tracer/job output
		results << resnum << " \t " << chain_char << " \t "<< atomname << std::endl;
		//now setup pymol output
		if( !first_sel_complete ) {
			missingHbond << "cmd.select(\"/" << posename_base_ << "//" << chain_char << "/" << resnum << "/" << atomname << "\")"<< std::endl;
			unsathbond << "/" << posename_base_ << "//" << chain_char << "/" << resnum << "+" ;
			first_sel_complete = true;
		}
		else  if (chain_char != chain_char_last){
			missingHbond << "cmd.select(\"sele + /" << posename_base_ << "//" << chain_char << "/" << resnum << "/" << atomname << "\")"<< std::endl;
			unsathbond <<" + "<< "/" << posename_base_ << "//" << chain_char << "/" << resnum << "+";
		}
		else {
			missingHbond << "cmd.select(\"sele + /" << posename_base_ << "//" << chain_char << "/" << resnum << "/" << atomname << "\")"<< std::endl;
			unsathbond << resnum << "+";
		}
		chain_char_last = chain_char;
	} //end itterate over all unsat AtomIDs
	unsathbond << std::endl;
	//finalize output
	missingHbond << "cmd.create(\"" << posename_base_ << "_unsat_hbond\", \"sele\")" << std::endl;
	missingHbond << "cmd.show(\"spheres\", \"" << posename_base_ << "_unsat_hbond\")" << std::endl;
	pymol_sel_hbond_unsat_ = missingHbond.str() ;

	//if we're not doing tracer output report this stuff to the job
	if(!tracer_){
		current_job->add_string( results_oss.str() );
		current_job->add_string( unsathbond_oss.str() );
	}

	return;
}


///@details This function doesn't do the printing itself.  The app InterfaceAnalyzer gets the multi-line string to write a file or print the selection
///@details From best packing to worse packing, colors go as Blue, Purple, Pink, Red
void InterfaceAnalyzerMover::print_pymol_selection_of_packing( core::pose::Pose const & pose,	utility::vector1< core::Real > interface_pack_scores) {
	std::ostringstream pymol_packing;
	pymol_packing << std::endl;

	std::string pymol_object_fullpose_pack = posename_base_ + "_fullpose_pack";
	pymol_packing << "cmd.create(\"" << pymol_object_fullpose_pack << "\", \"" << posename_base_ << "\")" << std::endl;

	///////////////////////////////////////////// WHOLE POSE COLOR BY PACKING /////////////////////
	for( core::Size i(1); i <= pose.total_residue(); i++ ) {
		core::Size resnum = pose.pdb_info()->number(i);
		char chain_char = pose.pdb_info()->chain(i);
		core::Size color;
		if      ( interface_pack_scores[i] >= 0.75 ) { color = 2; }  //blue
		else if ( interface_pack_scores[i] >= 0.50 ) { color = 16; } //purple
		else if ( interface_pack_scores[i] >= 0.25 ) { color = 12; } //pink
		else if ( interface_pack_scores[i] >  0.00 ) { color = 4; }  //red
		else { color = 24; } //gray, something went wrong

		pymol_packing << "cmd.select(\"/" << pymol_object_fullpose_pack << "//" << chain_char << "/" << resnum << "\")" << std::endl;
		pymol_packing << "cmd.color(" << color << ", \"sele\")" << std::endl;
	}

	pymol_sel_packing_ = pymol_packing.str() ;

	return;
}


///@detail Only want to register the calculators once, thus the 'if' statement in apply
void InterfaceAnalyzerMover::register_calculators()
{
	using namespace core::pose::metrics;
	using namespace protocols::toolbox::pose_metric_calculators;
	//determine name
	std::ostringstream interface_jump_cast;
	interface_jump_cast << interface_jump_;
	std::string const ijump(interface_jump_cast.str());


	//TR << "Running register_calculators" << std::endl;

	//this sucks but I can't figure out a way to iterate...
	Sasa_ = "Sasa_" + ijump;
	if( CalculatorFactory::Instance().check_calculator_exists( Sasa_ ) ){
		Warning() << "In InterfaceAnalyzerMover, calculator " << Sasa_
		<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( Sasa_, new core::pose::metrics::simple_calculators::SasaCalculator);
	}

	InterfaceNeighborDefinition_ = "InterfaceNeighborDefinition_" + ijump;
	if( CalculatorFactory::Instance().check_calculator_exists( InterfaceNeighborDefinition_ ) ){
		Warning() << "In InterfaceAnalyzerMover, calculator " << InterfaceNeighborDefinition_
		<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( InterfaceNeighborDefinition_, new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator(chain1_, chain2_));
	}

	InterfaceSasaDefinition_ = "InterfaceSasaDefinition_" + ijump;
	if( CalculatorFactory::Instance().check_calculator_exists( InterfaceSasaDefinition_ ) ){
		Warning() << "In InterfaceAnalyzerMover, calculator " << InterfaceSasaDefinition_
		<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( InterfaceSasaDefinition_, new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator(chain1_, chain2_));
	}

	InterfaceDeltaEnergetics_ = "InterfaceDeltaEnergetics_" + ijump;
	if( CalculatorFactory::Instance().check_calculator_exists( InterfaceDeltaEnergetics_ ) ){
		Warning() << "In InterfaceAnalyzerMover, calculator " << InterfaceDeltaEnergetics_
		<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( InterfaceDeltaEnergetics_, new core::pose::metrics::simple_calculators::InterfaceDeltaEnergeticsCalculator(InterfaceNeighborDefinition_));
	}

	NumberHBonds_ = "NumberHBonds_" + ijump;
	if( CalculatorFactory::Instance().check_calculator_exists( NumberHBonds_ ) ){
		Warning() << "In InterfaceAnalyzerMover, calculator " << NumberHBonds_
		<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( NumberHBonds_, new NumberHBondsCalculator);
	}

	BuriedUnsatisfiedPolars_ = "BuriedUnsatisfiedPolars_" + ijump;
	if( CalculatorFactory::Instance().check_calculator_exists( BuriedUnsatisfiedPolars_ ) ){
		Warning() << "In InterfaceAnalyzerMover, calculator " << BuriedUnsatisfiedPolars_
		<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator(  BuriedUnsatisfiedPolars_, new BuriedUnsatisfiedPolarsCalculator(Sasa_, NumberHBonds_));
	}



	return;
}//register_calculators

///@detail register calculator for multichain poses
void InterfaceAnalyzerMover::register_intergroup_calculator(){
	using namespace core::pose::metrics;
	using namespace protocols::toolbox::pose_metric_calculators;
	InterGroupNeighborsCalculator_ = "InterGroupNeighborsCalculator_" + posename_base_ ;
	if( CalculatorFactory::Instance().check_calculator_exists( InterGroupNeighborsCalculator_ ) ){
		Warning() << "In InterfaceAnalyzerMover, calculator " << InterGroupNeighborsCalculator_
		<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator(  InterGroupNeighborsCalculator_, new InterGroupNeighborsCalculator(chain_groups_, basic::options::option[basic::options::OptionKeys::pose_metrics::inter_group_neighbors_cutoff] ) );
	}
}

///@details calculate the average energy per residue in the interface
void InterfaceAnalyzerMover::calc_per_residue_energy( core::pose::Pose pose ){
	using namespace core::scoring;
	using namespace core;
	//using namespace core::scoring::Energies;
	//pose.update_residue_neighbors();
	(*sf_) (pose); //shits, giggles, and segfault prevention
	//itterate over all residues calc total energy for each then take avg
	core::Real total= 0.0;
	for( std::set< core::Size >::const_iterator it(interface_set_.begin()), end(interface_set_.end());
			 it != end; ++it){
		total+=pose.energies().residue_total_energies( *it )[ ScoreType( total_score ) ];
	}
	per_residue_energy_ = ( total / interface_set_.size() );
	n_interface_res_ = interface_set_.size();
}
///@details calculate the hbond energy and dampen it by exposure
void InterfaceAnalyzerMover::calc_hbond_sasaE( core::pose::Pose pose ){
	using namespace core::scoring::hbonds;
	using namespace core;
	using namespace core::chemical;
	(*sf_) (pose); //shits, giggles, and segfault prevention
	//calculate the per atom sasa
	// an atomID map is needed for the calc_per_atom_sasa method; it stores the actual calculated sasa for every atom
	//core::id::AtomID_Map< core::Real > atom_sasa;
	//core::pose::initialize_atomid_map( atom_sasa, pose, (core::Real)0.0 ); // initialize to 0.0 for "not computed"

	//utility::vector1< Real > rsd_sasa( pose.total_residue(), 0.0 );
	//core::Real probe_radius = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius];

	// create an atom_subset mask such that only the heavy atoms will have their sasa computed (ignore H's to make it faster)
	core::id::AtomID_Map< bool > atom_subset;
	atom_subset.clear();
	atom_subset.resize( pose.n_residue() );
	for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {
		atom_subset.resize( ii, pose.residue_type(ii).natoms(), false );
		for ( Size jj = 1; jj <= pose.residue_type(ii).nheavyatoms(); ++jj ) {
			atom_subset[ ii ][ jj ] = true;
		}
	}
	//now actually calculate the SASA for heavy atoms only
	//   core::Real total_sasa = 0.0;
	//   total_sasa = core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false /* no big polar H */, atom_subset );
	// 	//some tracer output to make sure things work
	// 	for ( core::Size ii=1; ii <= pose.total_residue(); ++ii ) {
	// 		core::conformation::Residue const & rsd = pose.residue( ii );
	// 		TR << "residue " << rsd.name3() << pose.pdb_info()->number(ii) << " atom_sasas: [ ";
	// 		for ( Size at=1; at <= rsd.nheavyatoms(); ++at ) {
	// 			core::id::AtomID atid( at, ii );
	// 			TR << utility::trim( rsd.atom_name( at ) ) << ":" << atom_sasa[ atid ] << ", ";
	// 		}
	// 		TR << "], total residue SASA: " << rsd_sasa[ ii ] << std::endl;
	// 	}

	//setup a vector of the sasa radii
	//Real const four_pi = 4.0f * Real( numeric::constants::d::pi );
	//AtomTypeSet const & atom_type_set = * ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
	//utility::vector1< Real > radii( atom_type_set.n_atomtypes(), 0.0 );
	//core::Size SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index( "SASA_RADIUS" );

	// TR << "radii: [ ";
	// for ( core::Size ii=1; ii <= atom_type_set.n_atomtypes(); ++ii ) {
	// 	radii[ ii ] = atom_type_set[ ii ].extra_parameter( SASA_RADIUS_INDEX );
	// 	TR << radii[ ii ] << ", ";
	// }
	// TR << "]" << std::endl;

	//EM options for bb-bb hbond output
	core::scoring::ScoreFunctionOP new_sf = scoring::getScoreFunction();
	scoring::methods::EnergyMethodOptions energymethodoptions( new_sf->energy_method_options() );
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
	new_sf->set_energy_method_options( energymethodoptions );

	//make an HbondSet
	//figure out energy statistics
	core::scoring::hbonds::HBondSet hbond_set;
	pose.update_residue_neighbors();
	//fill_hbond_set( pose,
	//						 false /*calc_deriv*/,
	//						 hbond_set,
	//						 false /*bb only*/ );
	hbond_set.setup_for_residue_pair_energies(pose, false, false);

	//itterate through all hbonds and figure out which ones are bb-bb betas
	Real n_crosschain_hbonds( 0 );
	//Real total_ratios( 0 );
	//total_hb_sasa_ = 0 ;
	total_hb_E_ = 0;
	for(Size ii=1; ii <= hbond_set.nhbonds(); ++ii){
		core::scoring::hbonds::HBond hbond ( hbond_set.hbond(ii) );
		Size don_resnum = hbond.don_res(); Size acc_resnum = hbond.acc_res();
		//need to take special consideration for multichain constructor
		if (multichain_constructor_){
			//TR<< "Do multichain hbonds eval..." << std::endl;
			//if they are different chains
			if( pose.chain( don_resnum ) != pose.chain( acc_resnum )){
				//if the acceptor or donor is in a fixed chain
				if ( fixed_chains_.count( pose.chain( don_resnum ) ) || fixed_chains_.count( pose.chain( acc_resnum ) )  ){
					//if the acceptor and donor are not both fixed chains
					if ( fixed_chains_.count( pose.chain( don_resnum ) ) != fixed_chains_.count( pose.chain( acc_resnum ) ) ){
						TR << "Found Hbond between chains: "<< pose.chain( don_resnum ) << " and " << pose.chain( acc_resnum ) <<std::endl;
						//now copy the same stuff as for the normal case:
						n_crosschain_hbonds += 1;
						//TR << "Hbond number: "<< ii << " is between chains." << std::endl;
						//now look at atoms involved
						// Size hatm = hbond.don_hatm() ;
						// conformation::Residue don_rsd = pose.residue(  don_resnum );
						// conformation::Residue acc_rsd = pose.residue(  acc_resnum );
						// Size don_atm = don_rsd.atom_base( hatm );
						// //Size don_atm = hbond.don_hatm() ;
						// Size acc_atm = hbond.acc_atm();

						// //Make AtomIDs
						// core::id::AtomID don_atid( don_atm, don_resnum );
						// core::id::AtomID acc_atid( acc_atm, acc_resnum );
						// Real acc_SASA( atom_sasa[ acc_atid ] );
						// Real don_SASA( atom_sasa[ don_atid ] );
						// Real hbond_SASA ( acc_SASA + don_SASA );

						// //now find possible sasa
						// Real const don_rad = radii[ don_rsd.atom(don_atm).type() ] + probe_radius;
						// Real const acc_rad = radii[ acc_rsd.atom(acc_atm).type() ] + probe_radius;
						// Real const possible_sasa( four_pi * (don_rad*don_rad + acc_rad*acc_rad ) );

						// //Adujust Hbond energy
						// Real adjust_weight( 2 );
						// Real sasa_ratio( hbond_SASA / possible_sasa );
						// Real const adjustment = (1 - adjust_weight * sasa_ratio);
						// //Real const new_Hbond_E = adjustment * hbond.energy();

						// //Print some stuff to make sure
						// const std::string & don_atomname = don_rsd.atom_name( don_atm );
						// const std::string & acc_atomname = acc_rsd.atom_name( acc_atm );
						// // TR << "Donor Res: " <<  don_rsd.name3() << pose.pdb_info()->number( don_resnum ) << " Donor Atm: " << don_atid << " (" << don_atomname << ")"
						// // 	 << " Acc Res: " << acc_rsd.name3() << pose.pdb_info()->number( acc_resnum ) << " Acc Atm: " << acc_atid << " (" << acc_atomname << ")"
						// // 	 << " Hbond_SASA: " << hbond_SASA << " Possible SASA: " << possible_sasa << " SASA ratio: " << sasa_ratio
						// // 	 << " Initial HBond E: " << hbond.energy() << " Adjusted: "<< new_Hbond_E
						// // 	 << std::endl ;
						// //add up the ratios and total hbond sasa
						// total_ratios += sasa_ratio;
						// total_hb_sasa_ += hbond_SASA;
						total_hb_E_ += hbond.energy();
					}
				}
			}//end if different chains
		}
		else { //not multichain constructor
			if ( pose.chain( hbond.don_res() ) != pose.chain( hbond.acc_res() ) ) {
				n_crosschain_hbonds += 1;
				// //TR << "Hbond number: "<< ii << " is between chains." << std::endl;
				// //now look at atoms involved
				// Size hatm = hbond.don_hatm() ;
				// conformation::Residue don_rsd = pose.residue(  don_resnum );
				// conformation::Residue acc_rsd = pose.residue(  acc_resnum );
				// Size don_atm = don_rsd.atom_base( hatm );
				// //Size don_atm = hbond.don_hatm() ;
				// Size acc_atm = hbond.acc_atm();

				// //Make AtomIDs
				// core::id::AtomID don_atid( don_atm, don_resnum );
				// core::id::AtomID acc_atid( acc_atm, acc_resnum );
				// Real acc_SASA( atom_sasa[ acc_atid ] );
				// Real don_SASA( atom_sasa[ don_atid ] );
				// Real hbond_SASA ( acc_SASA + don_SASA );

				// //now find possible sasa
				// Real const don_rad = radii[ don_rsd.atom(don_atm).type() ] + probe_radius;
				// Real const acc_rad = radii[ acc_rsd.atom(acc_atm).type() ] + probe_radius;
				// Real const possible_sasa( four_pi * (don_rad*don_rad + acc_rad*acc_rad ) );

				// //Adujust Hbond energy
				// Real adjust_weight( 2 );
				// Real sasa_ratio( hbond_SASA / possible_sasa );
				// Real const adjustment = (1 - adjust_weight * sasa_ratio);
				// Real const new_Hbond_E = adjustment * hbond.energy();

				// // //Print some stuff to make sure
				// // const std::string & don_atomname = don_rsd.atom_name( don_atm );
				// // const std::string & acc_atomname = acc_rsd.atom_name( acc_atm );
				// // // TR << "Donor Res: " <<  don_rsd.name3() << pose.pdb_info()->number( don_resnum ) << " Donor Atm: " << don_atid << " (" << don_atomname << ")"
				// // // 	 << " Acc Res: " << acc_rsd.name3() << pose.pdb_info()->number( acc_resnum ) << " Acc Atm: " << acc_atid << " (" << acc_atomname << ")"
				// // // 	 << " Hbond_SASA: " << hbond_SASA << " Possible SASA: " << possible_sasa << " SASA ratio: " << sasa_ratio
				// // // 	 << " Initial HBond E: " << hbond.energy() << " Adjusted: "<< new_Hbond_E
				// // // 	 << std::endl ;
				// // //add up the ratios and total hbond sasa
				// total_ratios += sasa_ratio;
				// total_hb_sasa_ += hbond_SASA;
				total_hb_E_ += hbond.energy();
			} //end if chains not equal
		} //end if not multichain
	}//end loop over all hbonds
	//get the avg exposure of the hbonds
	//hbond_exposure_ratio_ = total_ratios / n_crosschain_hbonds;
}//end function def

// //helper function for above to calculate hbond stats based on an hbond input
// void hbond_info_calculate( core::scoring::hbonds::HBond hbond ){
// 	//does nothing for now, figure out later
// }

void
InterfaceAnalyzerMover::compute_interface_sc( core::Size &, core::pose::Pose const & complexed_pose){

	core::scoring::sc::ShapeComplementarityCalculator sc_calc;
	// Split PDB into two surfaces
	for(core::Size i = 1; i <= complexed_pose.n_residue(); i++) {
		if(upstream_chains_.count( complexed_pose.chain( i ) ) )
			sc_calc.AddResidue(0, complexed_pose.residue(i));
		else if(downstream_chains_.count( complexed_pose.chain( i ) ) )
			sc_calc.AddResidue(1, complexed_pose.residue(i));
		else
			continue;
	}
	//now calculate and print results
	TR << "Computing Shape Complementarity Score..." << std::endl;
	TR << "Upstream chain(s) numbers: ";
	for( std::set< core::Size >::const_iterator it(upstream_chains_.begin()), end(upstream_chains_.end());
			 it != end; ++it){
		TR << *it << ", ";
	}
	TR << std::endl;

	TR << "Downstream chain(s) numbers: ";
	for( std::set< core::Size >::const_iterator it(downstream_chains_.begin()), end(downstream_chains_.end());
			 it != end; ++it){
		TR << *it << ", ";
	}
	TR << std::endl;

	//actual calculate function
	sc_calc.Calc();
	core::scoring::sc::RESULTS const results = sc_calc.GetResults();
	sc_value_ = results.sc;

}//end compute_interface_sc


///@details  Mutate all residues to GlY rescore complex energy and separated energy
void InterfaceAnalyzerMover::mut_to_gly( core::pose::Pose complex_pose, core::pose::Pose separated_pose ) {
	using namespace core;
	//need a copy of the pose to avoid screwing up the good one
	pose::Pose copy_complex( complex_pose );
	pose::Pose copy_separate( separated_pose );
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;

	//setupt task info
	core::pack::task::PackerTaskOP task(core::pack::task::TaskFactory::create_packer_task((copy_complex)));
	utility::vector1_bool packable(copy_complex.total_residue(), false); //false = nobody is packable
	utility::vector1< bool > allowed_aa( chemical::num_canonical_aas, false ); //no allowed residues
	allowed_aa[ core::chemical::aa_from_oneletter_code( 'G' ) ] = true; //allow gly only
	//allow all interface residues to be mutated to Gly
	for( std::set< core::Size >::const_iterator it(interface_set_.begin()), end(interface_set_.end());
			 it != end; ++it){
		task->nonconst_residue_task( *it ).restrict_absent_canonical_aas(allowed_aa);
		packable[ *it ] = true;
	}
	task->restrict_to_residues(packable);  //prevents non interface res from changing

#ifndef NDEBUG
	TR<< "GLY Packer Task: " << *(task) << std::endl;
#endif

	//apply mutations
	protocols::simple_moves::PackRotamersMoverOP packrot_mover( new protocols::simple_moves::PackRotamersMover( sf_ , task) );
	packrot_mover->apply ( copy_complex );
	packrot_mover->apply( copy_separate );
	gly_dG_ =  (*sf_) ( copy_complex ) - (*sf_) ( copy_separate )  ;
}

///@details
void InterfaceAnalyzerMover::calc_centroid_dG ( core::pose::Pose complex_pose, core::pose::Pose separated_pose ){
	if(!use_centroid_){
		centroid_dG_ = 0;
		return;
	}
	core::pose::Pose copy_complex( complex_pose );
	core::pose::Pose copy_separated( separated_pose );
	core::util::switch_to_residue_type_set( copy_complex , core::chemical::CENTROID );
	core::util::switch_to_residue_type_set( copy_separated , core::chemical::CENTROID );
	// use score3 but turn of RG
	core::scoring::ScoreFunctionOP scorefxn  = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	scorefxn->set_weight( core::scoring::rg, 0.0 );
	//Debugging:
	TR << "Centroid score of complex: " << (*scorefxn) ( copy_complex ) << std::endl;
	TR << "Centroid score of separated: " << (*scorefxn) ( copy_separated ) << std::endl;

	centroid_dG_ =  (*scorefxn) ( copy_complex ) - (*scorefxn ) ( copy_separated )  ;
}


///@details  sets up the packer task for the interface
void InterfaceAnalyzerMover::setup_task(core::pose::Pose & pose){
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;
	//set up the task to match this calculation
	//set up a packer task
	TaskFactoryOP tf = new TaskFactory();
	tf->push_back( new InitializeFromCommandline() );
	//force include current to prevent wonky results
	tf->push_back( new IncludeCurrent() );
	tf->push_back( new RestrictToRepacking() );
	if( use_resfile_ ) tf->push_back( new ReadResfile() );
	//use the same logic for calculators to restrict to the interface
	utility::vector1< std::pair< std::string, std::string> > calculators_used;
	if(multichain_constructor_){  //different calculator for different constructor
		std::pair< std::string, std::string> multichain_strings ( InterGroupNeighborsCalculator_, "neighbors" );
		calculators_used.push_back( multichain_strings );
	}
	else{ //if not multichain constructor
		std::pair< std::string, std::string> chain_strings ( InterfaceNeighborDefinition_, "interface_residues" );
		calculators_used.push_back( chain_strings );
	}
	tf->push_back( new RestrictByCalculatorsOperation( calculators_used ) );
	core::pack::task::PackerTaskOP task( tf->create_task_and_apply_taskoperations( pose ) );
	task_ = task ;
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
InterfaceAnalyzerMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const &,
	Movers_map const &,
	Pose const & pose
)
{
	if ( tag->getName() != "InterfaceAnalyzerMover" ) {
		//TR(t_warning) << " received incompatible Tag " << tag << std::endl;
		// assert(false);
		//   return;
	}
	sf_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );

	set_pack_separated(tag->getOption<bool>("pack_separated", false));
	set_use_resfile(tag->getOption<bool>("resfile", false ) );
	set_compute_packstat(tag->getOption<bool> ("packstat", false));
	set_compute_interface_sc(tag->getOption<bool> ("interface_sc", false));
	set_pack_input(tag->getOption<bool> ("pack_input", false));
	set_tracer(tag->getOption("tracer", false));
	set_use_jobname(tag->getOption("use_jobname", false));

	if (tag->hasOption("jump") && tag->hasOption("fixedchains"))
    {
    	throw utility::excn::EXCN_RosettaScriptsOption("Jump and fixedchains are mutually exclusive. Use either jump or fixedchains");
    }
	
	if ( (tag->hasOption("jump") || tag->hasOption("fixedchains") ) && tag->hasOption("ligandchain") )
	{
		throw utility::excn::EXCN_RosettaScriptsOption("if you specify Jump or fixedchains you cannot also specify ligandchain");
	}
	
	if (tag->hasOption("fixedchains"))
    {
    	set_interface_jump(0);
    	std::string chains_string = tag->getOption<std::string>("fixedchains");
			utility::vector1<std::string> fixed_chains_string = utility::string_split(chains_string,',');
			//parse the fixed chains to figure out pose chain nums
			//std::set< int > fixed_chains_ ; //This is a set of the CHAIN IDs, not residue ids
			TR << "Fixed chains are: " ;
			for(core::Size j = 1; j <= fixed_chains_string.size(); ++j){
				char this_chain (fixed_chains_string[ j ][0]);
				for (core::Size i = 1; i<=pose.total_residue(); ++i){
					if (pose.pdb_info()->chain( i ) == this_chain){
						fixed_chains_.insert( pose.chain(i) );
						break; //once we know something about the chain we can skip - we just need the chain id
					}
				}
				TR << this_chain << ", ";
			}
			TR << "these will be moved together." << std::endl;

			multichain_constructor_ = true;
			//fixed_chains_(fixed_chains)
    }else if(tag->hasOption("ligandchain"))
	{
		ligand_chain_ = tag->getOption<std::string>("ligandchain");
		multichain_constructor_ = true;
		multichain_constructor_ = true;
	}
	else
    {
    	set_interface_jump(tag->getOption("jump", 1));
    }
	//      tracer_(false), //output to tracer
	//      calcs_ready_(false), //calculators are not ready
	//      use_jobname_(false), //use the pose name
}

///@brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
InterfaceAnalyzerMover::fresh_instance() const
{
	return new InterfaceAnalyzerMover;
}

///@brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
InterfaceAnalyzerMover::clone() const
{
	return new InterfaceAnalyzerMover( *this );
}


}//analysis
}//protocols
