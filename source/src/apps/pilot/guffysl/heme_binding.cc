// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/guffysl/heme_binding.cc
/// @brief Application to redesign an enzyme to bind an alternative ligand
/// @author Sharon Guffy


// Headers
#include <protocols/moves/Mover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/database/open.hh>
//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static basic::Tracer TR("apps.pilot.guffysl.heme_binding");


//Define local options
namespace local {
basic::options::IntegerOptionKey const num_iterations ( "num_iterations" );
basic::options::RealOptionKey const prob_br_move ( "prob_br_move" );
basic::options::BooleanOptionKey const add_ligand_constraints ( "add_ligand_constraints" );
basic::options::StringVectorOptionKey const ligand_cst_atoms ( "ligand_cst_atoms" ); //Vector of atoms to constrain
basic::options::RealOptionKey const dist_cst_weight ( "dist_cst_weight" );
basic::options::RealOptionKey const dihedral_cst_weight ( "dihedral_cst_weight" );
basic::options::RealOptionKey const mc_temperature ( "mc_temperature" );
}







//Header file information
class HemeBindingMover: public protocols::moves::Mover{
public:
	//Default constructor
	HemeBindingMover();
	//Constructor with options
	HemeBindingMover(core::scoring::ScoreFunctionOP const & sfxn, protocols::backrub::BackrubMoverOP & brm,
		protocols::simple_moves::PackRotamersMoverOP & prm, core::Size n_it = 1000, core::Real temperature = 0.6);
	//Constructor with different options
	HemeBindingMover(core::scoring::ScoreFunctionOP const & sfxn, core::pack::task::TaskFactoryOP const & tf, core::Size n_it, core::Real temperature = 0.6);
	//Copy constructor
	HemeBindingMover (HemeBindingMover const & other);
	//Destructor
	virtual ~HemeBindingMover();
	//Set private data
	void score_function(core::scoring::ScoreFunctionOP other_score_function);
	void backrub_mover(protocols::backrub::BackrubMoverOP br_mover);
	void pack_rotamers_mover(protocols::simple_moves::PackRotamersMoverOP pr_mover);
	void num_iterations(core::Size n_it);
	void temperature(core::Real temp);
	void task_factory(core::pack::task::TaskFactoryOP tf);
	//Get private data
	core::scoring::ScoreFunctionCOP score_function() const;
	protocols::backrub::BackrubMoverCOP backrub_mover() const;
	protocols::simple_moves::PackRotamersMoverCOP pack_rotamers_mover() const;
	core::Size num_iterations() const;
	core::Real temperature() const;
	core::pack::task::TaskFactoryCOP task_factory() const;
	//Virtual methods
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual void apply(core::pose::Pose & pose);
	virtual std::string get_name() const;
	//For RosettaScripts
	//Any other methods
	void add_initial_constraints( core::pose::Pose & pose, utility::vector1< core::Size > designable_residues, utility::vector1< core::Size > packable_residues );
	void initialize_from_options();
private:
	//private data
	core::scoring::ScoreFunctionOP score_function_;
	core::scoring::ScoreFunctionOP design_score_function_;
	core::scoring::ScoreFunctionOP repack_score_function_;
	//Store OPs to the other movers so we can set their properties directly
	protocols::backrub::BackrubMoverOP backrub_mover_;
	protocols::simple_moves::PackRotamersMoverOP pack_rotamers_mover_;
	core::Size num_iterations_;
	core::Real temperature_;
	core::Real prob_br_;
	core::pack::task::TaskFactoryOP tf_;
	protocols::simple_moves::PackRotamersMoverOP design_mover_;
	utility::vector1< std::string > atoms_to_constrain_;
	core::scoring::constraints::BoundFuncOP dist_func_;
	core::scoring::func::CircularHarmonicFuncOP dihedral_func_;
	core::scoring::func::CircularHarmonicFuncOP dihedral_func_2_;
	core::scoring::constraints::ConstraintCOPs current_constraints_;
};




//Owning pointers
typedef utility::pointer::shared_ptr< HemeBindingMover > HemeBindingMoverOP;
typedef utility::pointer::shared_ptr< HemeBindingMover const > HemeBindingMoverCOP;


//local options
//namespace local{
// basic::options::FileOptionKey const num_it( "num_it");
/* basic::options::FileOptionKey const scaffold_protein( "scaffold_protein" );
basic::options::IntegerOptionKey const hbondres1("hbondres1");
basic::options::IntegerOptionKey const hbondres2("hbondres2");*/
//}//local

//Define methods
HemeBindingMover::HemeBindingMover() :
	protocols::moves::Mover("HemeBindingMover"),
	score_function_(core::scoring::get_score_function()),
	backrub_mover_ (protocols::backrub::BackrubMoverOP(new protocols::backrub::BackrubMover) ),
	pack_rotamers_mover_ (protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover) ),
	num_iterations_(1000),
	temperature_(0.6),
	//Set up a TaskFactory
	tf_ (core::pack::task::TaskFactoryOP(new core::pack::task::TaskFactory) ),
	design_mover_ (protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover))
{
	this->initialize_from_options();
	TR << "HemeBindingMover initialization complete" << std::endl;
	design_score_function_ = score_function_->clone();
	repack_score_function_ = score_function_->clone();
	design_score_function_->add_weights_from_file(basic::database::full_name( "scoring/weights/soft_rep_design.wts"));
	repack_score_function_->add_weights_from_file(basic::database::full_name( "scoring/weights/soft_rep.wts"));

}
HemeBindingMover::HemeBindingMover(core::scoring::ScoreFunctionOP const & sfxn, protocols::backrub::BackrubMoverOP & brm, protocols::simple_moves::PackRotamersMoverOP & prm, core::Size n_it, core::Real temperature):
	protocols::moves::Mover("HemeBindingMover")
{
	score_function_ = sfxn->clone();
	backrub_mover_ = protocols::backrub::BackrubMoverOP(new protocols::backrub::BackrubMover(*brm) );
	pack_rotamers_mover_ = protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover(*prm) );
	num_iterations_ = n_it;
	temperature_ = temperature;
	tf_ = prm->task_factory()->clone();
	design_mover_ = protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover );
	this->initialize_from_options();
	design_score_function_ = score_function_->clone();
	repack_score_function_ = score_function_->clone();
	design_score_function_->add_weights_from_file(basic::database::full_name( "scoring/weights/soft_rep_design.wts"));
	repack_score_function_->add_weights_from_file(basic::database::full_name( "scoring/weights/soft_rep.wts"));
}
HemeBindingMover::HemeBindingMover(core::scoring::ScoreFunctionOP const & sfxn, core::pack::task::TaskFactoryOP const & tf, core::Size n_it, core::Real temperature):
	protocols::moves::Mover("HemeBindingMover")
{
	score_function_ = sfxn->clone();
	tf_ = tf->clone();
	temperature_ = temperature;
	num_iterations_ = n_it;
	backrub_mover_ = protocols::backrub::BackrubMoverOP(new protocols::backrub::BackrubMover );
	pack_rotamers_mover_ = protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover);
	design_mover_ = protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover);
	this->initialize_from_options();
	design_score_function_ = score_function_->clone();
	repack_score_function_ = score_function_->clone();
	design_score_function_->add_weights_from_file(basic::database::full_name( "scoring/weights/soft_rep_design.wts"));
	repack_score_function_->add_weights_from_file(basic::database::full_name( "scoring/weights/soft_rep.wts"));
}

HemeBindingMover::HemeBindingMover(HemeBindingMover const & other):
	protocols::moves::Mover(other)
{
	score_function_ = other.score_function()->clone();
	backrub_mover_ = protocols::backrub::BackrubMoverOP(new protocols::backrub::BackrubMover(*other.backrub_mover()) );
	pack_rotamers_mover_ = protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover(*other.pack_rotamers_mover()));
	num_iterations_ = other.num_iterations();
	temperature_ = other.temperature();
	tf_ = other.task_factory()->clone();
	design_mover_ = protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover);
	this->initialize_from_options();
	design_score_function_ = score_function_->clone();
	repack_score_function_ = score_function_->clone();
	design_score_function_->add_weights_from_file(basic::database::full_name( "scoring/weights/soft_rep_design.wts"));
	repack_score_function_->add_weights_from_file(basic::database::full_name( "scoring/weights/soft_rep.wts"));
}
HemeBindingMover::~HemeBindingMover(){}
void HemeBindingMover::score_function(core::scoring::ScoreFunctionOP other_score_function){
	score_function_ = other_score_function->clone();
}
void HemeBindingMover::backrub_mover(protocols::backrub::BackrubMoverOP br_mover){
	backrub_mover_ = protocols::backrub::BackrubMoverOP(new protocols::backrub::BackrubMover(*br_mover) );
}
void HemeBindingMover::pack_rotamers_mover(protocols::simple_moves::PackRotamersMoverOP pr_mover){
	pack_rotamers_mover_ = protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover(*pr_mover) );
}
void HemeBindingMover::num_iterations(core::Size n_it){
	num_iterations_ = n_it;
}
void HemeBindingMover::temperature(core::Real temp){
	temperature_ = temp;
}
void HemeBindingMover::task_factory(core::pack::task::TaskFactoryOP tf){
	tf_ = tf->clone();
}
core::scoring::ScoreFunctionCOP HemeBindingMover::score_function() const{
	return score_function_;
}
protocols::backrub::BackrubMoverCOP HemeBindingMover::backrub_mover() const{
	return backrub_mover_;
}
protocols::simple_moves::PackRotamersMoverCOP HemeBindingMover::pack_rotamers_mover() const{
	return pack_rotamers_mover_;
}
core::Size HemeBindingMover::num_iterations() const{
	return num_iterations_;
}
core::Real HemeBindingMover::temperature() const{
	return temperature_;
}
core::pack::task::TaskFactoryCOP HemeBindingMover::task_factory() const{
	return tf_;
}
//Virtual methods
protocols::moves::MoverOP HemeBindingMover::clone() const{
	return new HemeBindingMover(*this);
}
protocols::moves::MoverOP HemeBindingMover::fresh_instance() const{
	return new HemeBindingMover;
}
std::string HemeBindingMover::get_name() const{
	return "HemeBindingMover";
}
void HemeBindingMover::initialize_from_options(){
	/*****************************************************************************************/
	//Set up from command line options


	num_iterations_ = basic::options::option[ local::num_iterations ].value();



	//Check options to see if they gave a probability for backrub moves

	//Only set it if the value is less than one since it's a probability
	if ( basic::options::option[ local::prob_br_move ].value() <= 1 ) {
		prob_br_ = basic::options::option[ local::prob_br_move ].value();
	} else {
		TR << "Cannot have a backrub move probability greater than 1" << std::endl;
		prob_br_ = 1.0;
	}
	//Set temperature
	temperature_ = basic::options::option[ local::mc_temperature ].value();

	//Set up functions in case of constraints
	//First = bounded function for atompair constraints
	core::Real lower_distance_bound = 2.8;
	core::Real upper_distance_bound = 4.7;
	core::Real distance_sd = 0.5; //This is somewhat arbitrary--look for an experimental value in Seeliger 2007 paper
	core::Real dihedral_1 = 0;
	core::Real dihedral_2 = numeric::conversions::radians(180.);
	core::Real dihedral_sd = numeric::conversions::radians(20.);

	dist_func_ = core::scoring::constraints::BoundFuncOP(new core::scoring::constraints::BoundFunc ( lower_distance_bound, upper_distance_bound, distance_sd, 0.5, "atom_pair_constraint" ) );
	//Second = circular harmonic function for dihedral constraints
	//We want to allow it to be 0 or 3.14 (anything planar)
	//Therefore we actually want two CircularHarmonicFuncs (ambiguous constraint)
	dihedral_func_ = core::scoring::func::CircularHarmonicFuncOP(new core::scoring::func::CircularHarmonicFunc ( dihedral_1, dihedral_sd ) );
	dihedral_func_2_ = core::scoring::func::CircularHarmonicFuncOP(new core::scoring::func::CircularHarmonicFunc ( dihedral_2, dihedral_sd ) );



}
void HemeBindingMover::add_initial_constraints( core::pose::Pose & pose, utility::vector1< core::Size > designable_residues, utility::vector1< core::Size > packable_residues ){
	//We know the ligand atoms that need to be constrained (stored in atoms_to_constrain_ ) and the residues the constraints will go to (in packable_residues)
	//For each atom to be constrained, we will have an AmbiguousConstraint filled with MultiConstraints, each with an AtomPairConstraint and a DihedralConstraint
	//We will use the residue object's atom_base( int atom# ) method to make the Dihedral constraints
	utility::vector1 < core::scoring::constraints::AmbiguousConstraintOP > constraints_to_add;


	//Start of a HUGE nest of loops--beware!
	for ( core::Size i = 1; i <= atoms_to_constrain_.size(); ++i ) {
		//Hope that the ligand is the last residue, but double check
		core::id::AtomID ligand_atom;
		if ( pose.residue( pose.total_residue() ).is_ligand() ) {
			ligand_atom = core::id::AtomID(pose.residue( pose.total_residue() ).atom_index(atoms_to_constrain_[i] ) , pose.total_residue()) ;
		} else {
			for ( core::Size resnum = 1; resnum <= pose.total_residue(); ++resnum ) {
				if ( pose.residue( resnum ).is_ligand() ) {
					ligand_atom = core::id::AtomID(pose.residue( resnum ).atom_index(atoms_to_constrain_[i] ) , resnum) ;
					break;
				}
			}
		}
		//For each atom in the ligand, we will make an ambiguous constraint
		core::scoring::constraints::AmbiguousConstraintOP possible_constraints (new core::scoring::constraints::AmbiguousConstraint );
		//We will now loop through all of the packable residues
		runtime_assert ( designable_residues.size() <= packable_residues.size() ); //I'm assuming designable_residues is a subset, so better make sure that's true
		for ( core::Size j = 1; j <= packable_residues.size(); ++j ) {
			//We now need to loop through all of the heavy atoms in this residue
			core::conformation::Residue current_res = pose.residue( j );
			for ( core::Size atomno = 1; atomno <= current_res.nheavyatoms(); ++atomno ) {
				//Make an AtomPairConstraint between ligand_atom and (atomno, j)
				core::scoring::constraints::AtomPairConstraintOP dist_cst (new core::scoring::constraints::AtomPairConstraint(ligand_atom, core::id::AtomID(atomno, j ), dist_func_ ) );

				//Now make a DihedralConstraint between base atom for ligand_atom, ligand_atom, (atomno, j), and base atom for atomno on current_res using dihedral_func_ or dihedral_func_2_
				core::id::AtomID ligand_base(pose.residue( ligand_atom.rsd() ).atom_base( ligand_atom.atomno() ) ,  ligand_atom.rsd() );
				core::id::AtomID res_base(current_res.atom_base( atomno ), j );
				core::scoring::constraints::DihedralConstraintOP dihedral_cst (new core::scoring::constraints::DihedralConstraint ( ligand_base, ligand_atom, core::id::AtomID( atomno, j ), res_base , dihedral_func_ ) );
				core::scoring::constraints::DihedralConstraintOP dihedral_cst_2 (new core::scoring::constraints::DihedralConstraint ( ligand_base, ligand_atom, core::id::AtomID( atomno, j ), res_base , dihedral_func_2_ ) );

				//Now make a MultiConstraint that contains both of these constraints

				core::scoring::constraints::MultiConstraintOP this_atom_csts (new core::scoring::constraints::MultiConstraint );
				this_atom_csts->add_individual_constraint(dist_cst);
				this_atom_csts->add_individual_constraint(dihedral_cst);


				core::scoring::constraints::MultiConstraintOP this_atom_csts_2 (new core::scoring::constraints::MultiConstraint );
				this_atom_csts->add_individual_constraint(dist_cst);
				this_atom_csts->add_individual_constraint(dihedral_cst_2);

				//Now add this MultiConstraint to our AmbiguousConstraint
				possible_constraints->add_individual_constraint( this_atom_csts );
				possible_constraints->add_individual_constraint( this_atom_csts_2 );
			}

		}
		//Only add this AmbiguousConstraint if it is not empty
		if ( possible_constraints->size() != 0 ) {
			constraints_to_add.push_back( possible_constraints );
		}
	}

	TR << "Adding initial constraints to pose: " << std::endl;
	current_constraints_ = pose.add_constraints( constraints_to_add );
}



void HemeBindingMover::apply(core::pose::Pose & pose){
	/*****************************************************************************************/
	//Set up TaskOperations for the TaskFactory
	tf_->push_back( core::pack::task::operation::InitializeFromCommandlineOP (new core::pack::task::operation::InitializeFromCommandline ) );
	core::pack::task::operation::ReadResfileOP resfile_op (new core::pack::task::operation::ReadResfile );
	//Get name of resfile to use in setting up backrub
	//There must be some step required in here to actually read the resfile into the ReadResfile object
	resfile_op->default_filename();
	std::string resfile_name = resfile_op->filename();
	TR << "Resfile in " << resfile_name << std::endl;
	tf_->push_back(resfile_op);
	/*****************************************************************************************/
	//Find which residues are being designed and repacked
	//Make a PackerTask
	core::pack::task::PackerTaskOP pt = tf_->create_task_and_apply_taskoperations(pose);
	//Find out which residues are being designed/repacked
	utility::vector1<bool> being_designed = pt->designing_residues();
	utility::vector1<bool> being_repacked = pt->repacking_residues();
	//Create vectors to contain list of residues being designed and residues being repacked
	utility::vector1< core::Size > residues_to_design;
	utility::vector1< core::Size > residues_to_repack;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( being_designed[i] ) {
			residues_to_design.push_back(i);
		}
		//We also want to include designed residues in repackable residues
		if ( being_repacked[i] ) {
			residues_to_repack.push_back(i);
		}
	}//End for loop, i out of scope
	/*****************************************************************************************/

	//Constraints will not work for this--remove this code
	/*
	//Set any constraints for the score function here
	//If the user asked for constraints, first add them to the pose
	if ( basic::options::option [ local::add_ligand_constraints ].value() ) {

	//Find which atoms to constrain
	atoms_to_constrain_ = basic::options::option [ local::ligand_cst_atoms ].value();
	if ( atoms_to_constrain_.size() == 0 )  {
	TR << "Error: Did not give any atoms to constrain! No constraints will be added." << std::endl;
	}
	else {
	add_initial_constraints( pose, residues_to_design, residues_to_repack );

	//Then add the weights for those constraints to the score function
	score_function_->set_weight(core::scoring::atom_pair_constraint, basic::options::option [ local::dist_cst_weight ].value() );
	score_function_->set_weight(core::scoring::dihedral_constraint, basic::options::option [ local::dihedral_cst_weight ].value() );
	TR << "Done setting up constraints" << std::endl;
	}
	}
	else{
	TR << "No constraints were requested." << std::endl;
	}
	//Might need to add mm_bend term to score function--check if weight is zero, if so, set to 1
	if ( score_function_->get_weight(core::scoring::mm_bend) == 0 ){
	TR << "Setting weight for mm_bend term to 1.0 (default)" << std::endl;
	score_function_->set_weight( core::scoring::mm_bend, 1.0);
	}
	*/
	//Set the scorefunction for the pack_rotamers_mover_ and design_mover_
	pack_rotamers_mover_->score_function(repack_score_function_);
	design_mover_->score_function(design_score_function_);

	/*****************************************************************************************/
	//Set up the BackrubMover
	//Set the maximum angle for our backrub moves to 20 degrees
	core::Real max_rotation_angle = numeric::conversions::radians(20.);
	core::Size segment_added;
	//Make a PoseOP that will hold a copy of our pose to use for setting input pose for backrub_mover_
	core::pose::PoseOP br_pose = pose.clone();
	backrub_mover_->set_input_pose(br_pose);
	//Create the set of possible segments for backrub moves (Include only redesigned residues)
	for ( core::Size i = 1; i <= residues_to_design.size(); ++ i ) {
		//Get the current residue number
		core::Size current_resnum = residues_to_design[i];
		segment_added = 0;
		//If this is not the first residue and is not the last residue
		if ( current_resnum != 1 && current_resnum != pose.total_residue() ) {
			core::id::AtomID prev_id(pose.residue(current_resnum - 1).atom_index(" CA "), current_resnum - 1);
			core::id::AtomID next_id(pose.residue(current_resnum + 1).atom_index(" CA "), current_resnum + 1);
			segment_added = backrub_mover_->add_segment(prev_id, next_id, max_rotation_angle);

		} else {
			//If this is the first residue or the last residue, it will go from N to C of this residue
			//Get AtomID for N
			core::id::AtomID n_id(pose.residue(current_resnum).atom_index(" N  "), current_resnum);
			//Get AtomID for C
			core::id::AtomID c_id(pose.residue(current_resnum).atom_index(" C  "), current_resnum);
			backrub_mover_->add_segment( n_id, c_id, max_rotation_angle );
		}
		if ( segment_added == 0 ) {
			TR << "ERROR: Segment not added for residue " << current_resnum << std::endl;
		} else {
			TR << "Segment added for residue " << current_resnum << std::endl;
		}
	} //End for loop, i out of scope
	TR << "Backrub segments have now been added" << std::endl;
	/*****************************************************************************************/
	//Set up minimization
	//First create movemap to minimize backbone and sidechains
	core::kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );
	//Then create the MinimizerOptions object
	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true); //Tolerance = 0.01, use_nb_list_in true
	//Now create the AtomTreeMinimizer
	core::optimization::AtomTreeMinimizer atm;

	//Set up MonteCarlo object
	protocols::moves::MonteCarloOP mc (new protocols::moves::MonteCarlo(pose, *score_function_, temperature_) );


	//****************************************************************************************
	//Finding neighbors of designed residues
	std::map< core::Size, std::set <core::Size > > neighbors;//Maps each designable residue
	//to a set of its neighbors
	//Create a point graph that stores coordinate information from the pose
	core::conformation::PointGraphOP coords (new core::conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation(pose.conformation(), *coords);
	//How far away do we want these neighbors to be?  Let's say 10 angstroms
	//Set our distance cutoff to 10 angstroms
	core::Real neighbor_cutoff = 10.;
	core::conformation::find_neighbors< core::conformation::PointGraphVertexData, core::conformation::PointGraphEdgeData >(coords, neighbor_cutoff);
	//This makes edges between residues separated by no more than our cutoff distance
	//Iterate through the edges for each repackable residue
	core::Size designable_res;
	for ( core::Size ii = 1; ii < residues_to_design.size(); ++ii ) {
		designable_res = residues_to_design[ii];
		for ( core::conformation::PointGraph::UpperEdgeListConstIter iter = coords->get_vertex( designable_res ).const_upper_edge_list_begin();
				iter != coords->get_vertex( designable_res ).const_upper_edge_list_end();
				++iter ) {
			neighbors[designable_res].insert(iter->upper_vertex());
		}
	}
	/*****************************************************************************************/
	//Initialize pose copy
	core::pose::PoseOP pose_copy;
	//Begin main loop
	core::Size accepted = 0;
	for ( core::Size i = 1; i <= num_iterations_; ++i ) {
		//****************************************************************************************
		//************** BACKRUB STAGE ************************************************************
		//****************************************************************************************
		pose_copy = pose.clone();

		//Find a random number (uniform)
		core::Real prob_check = numeric::random::uniform();
		//Check if the random number is less than prob_br_

		//If so, execute this code
		//Make sure the mover has the correct input pose for foldtree (otherwise would reset segments)
		//delete br_pose;
		core::pose::PoseOP br_pose = pose.clone();
		backrub_mover_->set_input_pose(br_pose);
		core::Size res1;
		core::Size res2;
		if ( prob_check <= prob_br_ ) {
			backrub_mover_->apply(pose);
			//For debugging purposes, print which residues were just moved
			res1 = backrub_mover_->segment(backrub_mover_->last_segment_id()).start_atomid().rsd();
			res2 = backrub_mover_->segment(backrub_mover_->last_segment_id()).end_atomid().rsd();
			TR << "Backrub move applied to residues " << res1 << " to " << res2 << " angle " << backrub_mover_->last_angle() << std::endl;
		} else {
			//Otherwise we still need to pick a residue to redesign
			//Find how many segments there are total
			core::Size total_segments = backrub_mover_->num_segments();
			//Pick a random segment
			core::Size segment_number = numeric::random::random_range( 1, total_segments );
			//Get res1 and res2 from that segment
			res1 = backrub_mover_->segment( segment_number ).start_atomid().rsd();
			res2 = backrub_mover_->segment( segment_number ).end_atomid().rsd();
			TR << "No backrub move this iteration. Segment " << res1 << " to " << res2 << " will redesign." << std::endl;
		}
		//We will only redesign (if applicable) res1 and res2
		//****************************************************************************************
		//************** DESIGN STAGE ************************************************************
		//****************************************************************************************
		//Make a new PackerTask from the TaskFactory
		core::pack::task::PackerTaskOP design_task = tf_->create_task_and_apply_taskoperations(pose);
		//Turn OFF design and repacking for all residues except res1 and res2
		//Make a boolean vector which is true only at the residue where we made the backrub move
		//If res1 == res2, then this is the residue that was repacked
		//Otherwise, the residue was (res1 + res2) /2
		core::Size repacked_res;
		if ( res1 == res2 ) {
			repacked_res = res1;
		} else {
			repacked_res = (res1 + res2 ) / 2;
		}
		utility::vector1 <bool> is_backrub_res;
		for ( core::Size j = 1; j <= pose.total_residue(); ++j ) {
			if ( j == repacked_res ) {
				is_backrub_res.push_back(true);
			} else {
				is_backrub_res.push_back(false);
			}
		}
		design_task->restrict_to_residues(is_backrub_res);
		//Initialize a PackRotamersMover from this PackerTask and our score_function_
		design_mover_->task(design_task);
		//Call apply on the PackRotamersMover
		design_mover_->apply(pose);
		//**********************************************************************************
		//Check if the residue was mutated & constraints are being used

		//Constraints don't work--remove this code


		/*
		bool mutation_occurred = ( pose.residue( repacked_res ).name3() == pose_copy->residue( repacked_res ).name3() );
		if ( mutation_occurred &&  basic::options::option [ local::add_ligand_constraints ].value() ) {
		TR << "Replacing constraints for mutated residue " << repacked_res << std::endl;
		//In this case, we need to remove all of the old constraints that apply to repacked_res
		//Make a non-const copy of current_constraints_
		core::scoring::constraints::ConstraintOPs new_constraints;
		core::scoring::constraints::AmbiguousConstraintOP new_amb_constraint = new core::scoring::constraints::AmbiguousConstraint;
		//Loop through each AmbiguousConstraint
		for (core::Size amb_cst_no(1); amb_cst_no <= current_constraints_.size(); ++amb_cst_no ){
		//We will be making a new AmbiguousConstraint to replace it
		new_amb_constraint = new core::scoring::constraints::AmbiguousConstraint;
		core::scoring::constraints::AmbiguousConstraintCOP current_amb_constraint = dynamic_cast< core::scoring::constraints::AmbiguousConstraint const * > (current_constraints_[ amb_cst_no ]() );
		//Loop through each MultiConstraint in the AmbiguousConstraint
		core::scoring::constraints::ConstraintCOPs member_csts =  current_amb_constraint->member_constraints();
		for ( core::Size multi_cst_no(1); multi_cst_no <= member_csts.size(); ++multi_cst_no ){
		//Call residues() for each of the MultiConstraints within--gives vector1 of residue positions that are included
		core::scoring::constraints::MultiConstraintCOP current_multi_cst = dynamic_cast< core::scoring::constraints::MultiConstraint const * > (member_csts[ multi_cst_no ]() ) ;
		utility::vector1< core::Size > cst_resnums = current_multi_cst->residues();
		//If residues do not include repacked_res, add that MultiConstraint to the new AmbiguousConstraint
		if ( cst_resnums.has_value( repacked_res ) ) {
		new_amb_constraint->add_individual_constraint( member_csts[ multi_cst_no ] );
		}
		else {
		TR << "Removing old constraint on " << repacked_res << std::endl;
		}
		}
		//Now add in new MultiConstraints for all heavy atoms of the new residue (use similar code to add_initial_constraints() block )
		//First find the AtomID for the ligand atom from one of the existing AtomPairConstraints--this might be a little tricky
		//Look at the first MultiConstraint in member_csts & loop through *its* member constraints
		core::scoring::constraints::MultiConstraintCOP first_multi_cst =  dynamic_cast< core::scoring::constraints::MultiConstraint const * > (member_csts[ 1 ]() );
		core::scoring::constraints::ConstraintCOPs csts_to_check = first_multi_cst->member_constraints();
		core::id::AtomID ligand_atom_id;
		for (core::Size ii = 1; ii <= csts_to_check.size(); ++ii ) {
		if (csts_to_check[ ii ]->type() == "AtomPair" ){
		//There should be just two AtomIDs associated with this constraint.  Get the one that is for a ligand atom
		if (pose.residue( csts_to_check[ ii ]->atom( 1 ).rsd() ).is_ligand() ){
		ligand_atom_id = csts_to_check[ ii ]->atom( 1 );
		break;
		}
		else if (pose.residue( csts_to_check[ ii ]->atom( 2 ).rsd() ).is_ligand() ) {
		ligand_atom_id = csts_to_check[ ii ]->atom( 2 );
		break;
		}
		else{
		TR << "Did not find ligand atom in the AtomPairConstraint.  Check add_initial_constraints and lines 556-600." << std::endl;
		}
		}
		}
		//Now that we have the id for the ligand atom, we can create the new MultiConstraints for all heavy atoms in residue repacked_res
		core::conformation::Residue repacked_res_object = pose.residue( repacked_res );
		for ( core::Size atomno = 1; atomno <= repacked_res_object.nheavyatoms(); ++atomno ) {
		//Make an AtomPairConstraint between ligand_atom and (atomno, j)
		core::scoring::constraints::AtomPairConstraintOP dist_cst = new core::scoring::constraints::AtomPairConstraint(ligand_atom_id, core::id::AtomID(atomno, repacked_res ), dist_func_ );

		//Now make a DihedralConstraint between base atom for ligand_atom, ligand_atom, (atomno, j), and base atom for atomno on current_res using dihedral_func_ or dihedral_func_2_
		core::id::AtomID ligand_base(pose.residue( ligand_atom_id.rsd() ).atom_base( ligand_atom_id.atomno() ) ,  ligand_atom_id.rsd() );
		core::id::AtomID res_base(repacked_res_object.atom_base( atomno ), repacked_res );
		core::scoring::constraints::DihedralConstraintOP dihedral_cst = new core::scoring::constraints::DihedralConstraint ( ligand_base, ligand_atom_id, core::id::AtomID( atomno, repacked_res ), res_base , dihedral_func_ );
		core::scoring::constraints::DihedralConstraintOP dihedral_cst_2 = new core::scoring::constraints::DihedralConstraint ( ligand_base, ligand_atom_id, core::id::AtomID( atomno, repacked_res ), res_base , dihedral_func_2_ );

		//Now make a MultiConstraint that contains both of these constraints

		core::scoring::constraints::MultiConstraintOP this_atom_csts = new core::scoring::constraints::MultiConstraint;
		this_atom_csts->add_individual_constraint(dist_cst);
		this_atom_csts->add_individual_constraint(dihedral_cst);

		core::scoring::constraints::MultiConstraintOP this_atom_csts_2 = new core::scoring::constraints::MultiConstraint;
		this_atom_csts->add_individual_constraint(dist_cst);
		this_atom_csts->add_individual_constraint(dihedral_cst_2);

		//Now add this MultiConstraint to our AmbiguousConstraint
		new_amb_constraint->add_individual_constraint( this_atom_csts );
		new_amb_constraint->add_individual_constraint( this_atom_csts_2 );
		}

		//Finally add the new AmbiguousConstraint to new_constraints
		new_constraints.push_back( new_amb_constraint );
		}
		//Next we need to add the constraints back (replace the pose's old constraint set with pose.add_constraints()
		pose.remove_constraints();
		//Store the returned constant constraint set in current_constraints_
		current_constraints_ = pose.add_constraints( new_constraints );

		}
		*/



		//****************************************************************************************
		//************** REPACKING STAGE ************************************************************
		//****************************************************************************************
		//Find neighbors for the two residues that were used in the move
		//utility::vector1 <core::Size> neighbors;
		//We will only repack residues that are in neighbors
		//Use neighbors to make a vector1 of booleans stating if that residue will be repacked
		utility::vector1 <bool> is_neighbor;

		for ( core::Size j = 1; j <= pose.total_residue(); ++j ) {
			if ( neighbors[ repacked_res ].count( j ) != 0 && being_repacked[ j ] ) { //If this residue is a neighbor and is repackable
				is_neighbor.push_back(true);
			} else {
				is_neighbor.push_back(false);
			}
		}
		//We will make a PackerTask that can only repack and that only repacks the residues
		core::pack::task::PackerTaskOP neighbor_packer_task = tf_->create_task_and_apply_taskoperations(pose);
		neighbor_packer_task->restrict_to_repacking();
		neighbor_packer_task->restrict_to_residues(is_neighbor);
		//Prepare the PackRotamersMover
		pack_rotamers_mover_->task(neighbor_packer_task);
		//Repack neighbors
		pack_rotamers_mover_->apply(pose);
		//****************************************************************************************
		//*****************************Minimization***********************************************
		//****************************************************************************************
		atm.run(pose, mm,  *score_function_, min_opts);




		//****************************************************************************************
		//******************************Monte Carlo check***********************************************
		//****************************************************************************************

		TR << "Design round " << i << " complete" << std::endl;
		//MonteCarlo check
		if ( mc->boltzmann(pose) ) {
			TR << "Move accepted " << std::endl;
			++accepted;
		}
	}
	//****************************************************************************************
	//*******************************End main loop******************************************
	//Score the pose for output score file
	//Recover the lowest scoring pose from the MC object
	//Consider changing score function first (remove any constraints)
	mc->recover_low( pose );
	core::scoring::ScoreFunctionOP basic_score_function = core::scoring::get_score_function();
	( *basic_score_function )(pose);
	TR << "Final acceptance rate: " << (1.0*accepted)/num_iterations_ << std::endl;
}
int main( int argc, char* argv[] ){
	try{
		basic::options::option.add( local::num_iterations, "number of iterations").def(1000);
		basic::options::option.add( local::prob_br_move, "probability of making a backrub move before each design run" ).def(1);
		basic::options::option.add( local::add_ligand_constraints, "the specified ligand atoms are constrained to interact with protein" ).def(false);
		basic::options::option.add( local::ligand_cst_atoms, "ligand atoms for which we will constrain distances/angles with the protein" ); //No default
		basic::options::option.add( local::dist_cst_weight, "weight for atom pair constraints" ).def(1.0);
		basic::options::option.add( local::dihedral_cst_weight, "weight for dihedral constraints" ).def(1.0);
		basic::options::option.add( local::mc_temperature, "Temperature to use in MonteCarlo object" ).def(0.6);
		/* basic::options::BooleanOptionKey const add_ligand_constraints ( "add_ligand_constraints" );
		basic::options::StringVectorOptionKey const ligand_cst_atoms ( "ligand_cst_atoms" ); //Vector of atoms to constrain
		basic::options::RealOptionKey const dist_cst_weight ( "dist_cst_weight" );
		basic::options::RealOptionKey const dihedral_cst_weight ( "dihedral_cst_weight" );*/
		// option.add( local::scaffold_protein, "scaffold protein filename").def("scaffold.pdb");
		//    option.add( local::hbondres2, "second h bond residue").def(2);
		//Add any local options
		devel::init(argc, argv);
		TR << "Initialization complete" << std::endl;
		HemeBindingMoverOP hb_mover (new HemeBindingMover );
		//Apply the HemeBindingMover using JD2
		TR << "Preparing to run: " << std::endl;
		protocols::jd2::JobDistributor::get_instance()->go(hb_mover);
		return 0;
	}
catch ( utility::excn::Exception const &e ){
	std::cout << "Caught exception " << e.msg() << std::endl;
	return -1;
}
}

