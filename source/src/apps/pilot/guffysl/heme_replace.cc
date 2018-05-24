// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/guffysl/heme_binding.cc
/// @brief Mover to replace heme alternative ligand with heme of the appropriate variant and rescore the pose
/// @author Sharon Guffy

//Headers
#include <protocols/moves/Mover.hh>
#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

//Define tracer
using basic::Error;
using basic::Warning;
using basic::T;
static basic::Tracer TR("apps.pilot.guffysl.heme_replace");

//Should take the following options:
//Original pdb file (through in:file:s)
//PDB file with coordinates for the appropriate form of heme
//Params files for heme and for the new alternative ligand

//Define local options
namespace local {
basic::options::FileOptionKey const heme_conformation( "heme_conformation" );
basic::options::BooleanOptionKey const no_relax( "no_relax" );
basic::options::IntegerOptionKey const num_trials( "num_trials" );
basic::options::RealOptionKey const repacking_distance( "repacking_distance" );
}//local


/// @brief Replaces original porphyrin ligand with heme and rescores the pose. Prints energy difference between two ligands.
class HemeReplace : public protocols::moves::Mover {
public:
	HemeReplace(){
		//Set up the score function
		score_function_ = core::scoring::get_score_function();
		//Set up FastRelax
		fast_relax_ = protocols::relax::FastRelaxOP( new protocols::relax::FastRelax(score_function_) ); //Will set up other options from command line and defaults
		//Set native pose for fast_relax_
		if ( basic::options::option [ basic::options::OptionKeys::in::file::native ] .user() ) {
			//If a native pose was given, set it
			core::pose::PoseOP native_pose = core::import_pose::pose_from_pdb( basic::options::option[ basic::options::OptionKeys::in::file::native ].value() );
			fast_relax_->set_native_pose( native_pose );

		}
		//Load in the coordinates for heme
		heme_coords_ = core::import_pose::pose_from_pdb( basic::options::option[ local::heme_conformation ].value() );

		//Set up the mover for initial repacking
		core::pack::task::TaskFactoryOP initial_repack_tf (new core::pack::task::TaskFactory );
		initial_repack_tf->push_back(new core::pack::task::operation::InitializeFromCommandline);
		initial_repack_mover_ = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover(score_function_) );
		initial_repack_mover_->task_factory(initial_repack_tf);
		restricted_repack_mover_ = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover(score_function_) );//Its tf will be set later
	}
	virtual ~HemeReplace(){
	}
	//Getters
	//Setters
	//Pure virtual methods
	virtual std::string get_name() const {
		return "HemeReplace";
	}
	virtual void apply(core::pose::Pose & pose){
		//Get name of the structure from PDBInfo
		std::string structure_name = pose.pdb_info()->name();

		find_ligand(pose);
		setup_task_factory(pose);
		runtime_assert ( porphyrin_position_ != 0); //It found the porphyrin residue in the pose

		//Overall structure for apply:
		//First repack everything (according to the command line) with fixed backbone
		initial_repack_mover_->apply(pose);

		//Initialize vectors to hold scores
		utility::vector1< core::Real > porphyrin_scores;
		utility::vector1< core::Real > heme_scores;
		//Save the starting pose for the loop
		core::pose::PoseCOP loop_starting_pose = pose.clone();

		//Declare OPs to hold lowest energy poses with each ligand
		core::pose::PoseOP lowest_por_pose;
		core::pose::PoseOP lowest_heme_pose;
		//Now start the loop
		for ( int iteration = 1; iteration <= basic::options::option[ local::num_trials ].value(); ++iteration ) {
			//Restore pose to the starting conformation
			pose = *loop_starting_pose;
			//Relax original pose to ensure lowest score
			if ( !basic::options::option[ local::no_relax ].value() ) { //If we *are* relaxing
				fast_relax_->apply( pose );
			}
			//Score the pose
			porphyrin_scores.push_back( ( *score_function_ )( pose ) );

			//************************
			//If this is the first iteration or if this is the lowest porphyrin_score so far, set lowest_por_pose to a copy of pose
			if ( iteration == 1 || utility::min( porphyrin_scores ) == porphyrin_scores[ iteration ] ) {
				lowest_por_pose = pose.clone();

			}
			//************************

			//Reset the pose
			pose = *loop_starting_pose;

			//Replace the ligand
			replace_ligand( pose );
			//Relax the structure using FastRelax
			if ( !basic::options::option[ local::no_relax ].value() ) { //If we *are* relaxing
				fast_relax_->apply( pose );
			} else {
				restricted_repack_mover_->apply( pose );
			}
			//Rescore the pose
			heme_scores.push_back( ( *score_function_ )( pose ) );

			//************************
			//If this is the first iteration or if this is the lowest porphyrin_score so far, set lowest_por_pose to a copy of pose
			if ( iteration == 1 || utility::min( heme_scores ) == heme_scores[ iteration ] ) {
				lowest_heme_pose = pose.clone();

			}
			//************************



			//End for loop here
		}

		//************************
		//Dump the lowest-scoring poses with each ligand
		std::string por_name = structure_name + "_por.pdb";
		lowest_por_pose->dump_pdb(por_name);
		std::string heme_name = structure_name + "_heme.pdb";
		lowest_heme_pose->dump_pdb(heme_name);
		//************************


		//************************
		//Print the avg. & std. dev of porphyrin, heme, & differences

		TR << "Summary for " << structure_name << ": " << std::endl;
		core::Real heme_total = 0;
		core::Real por_total = 0;
		core::Real diff_total = 0;
		for ( core::Size ii = 1; ii <= heme_scores.size(); ++ii ) {
			heme_total += heme_scores[ ii ];
			por_total += porphyrin_scores[ ii ];
			diff_total += ( porphyrin_scores[ ii ] - heme_scores[ ii ] );
		}
		core::Real heme_average = heme_total / heme_scores.size();
		core::Real por_average = por_total / porphyrin_scores.size();
		core::Real diff_average = diff_total / heme_scores.size();

		//Calculate standard deviations
		core::Real heme_sd = 0;

		core::Real por_sd = 0;

		core::Real diff_sd = 0;
		for ( core::Size ii = 1; ii <= heme_scores.size(); ++ii ) {
			por_sd += std::pow( ( porphyrin_scores[ii] - por_average ), 2 );
			heme_sd += std::pow( ( heme_scores[ii] - por_average ), 2 );
			diff_sd += std::pow( (porphyrin_scores[ ii ] - heme_scores[ ii ] - diff_average ), 2 );
		}
		por_sd = por_sd / (heme_scores.size() - 1 );
		heme_sd = heme_sd / (heme_scores.size() - 1 );
		diff_sd = diff_sd / (heme_scores.size() - 1 );


		TR << "Porphyrin average " << por_average << " Std dev: " << por_sd << std::endl;
		TR << "Heme average " << heme_average << " Std dev: " << heme_sd << std::endl;
		TR << "Difference average " << diff_average << " Std dev: " << diff_sd << std::endl;
		TR<< "*******************************************" << std::endl;

		//************************
		//  TR << "Start score: " << original_score << "  Heme score: " << heme_score << std::endl;
	}
private:
	//Private method to make apply more modular
	void replace_ligand( core::pose::Pose & pose){
		//pose contains structure with alternate ligand bound
		//heme_coords_ contains coordinates for the heme to replace it with
		//First we should align the two ligands
		//Alternate porphyrin will appear in pdb file as POR
		core::conformation::Residue hem_residue = heme_coords_->residue(1);//Should be the only residue in this pose
		runtime_assert(hem_residue.name3() == "HEM");
		utility::vector1<std::pair < std::string, std::string > > atom_pairs;
		//Define pairs of atoms for alignment--make sure they will be consistent for all of our heme derivatives
		//Use 3 of our coordinating nitrogen atoms
		atom_pairs.push_back(std::make_pair( " NA ", " NA ") );
		atom_pairs.push_back(std::make_pair( " NB ", " NB ") );
		atom_pairs.push_back(std::make_pair( " NC ", " NC ") );
		//Find the porphyrin residue in the original pose
		runtime_assert ( porphyrin_position_ != 0); //It found the porphyrin residue in the pose
		hem_residue.orient_onto_residue( pose.residue(porphyrin_position_), atom_pairs );
		//The residues are now aligned, so we can just call pose.replace_residue()
		pose.replace_residue(porphyrin_position_, hem_residue, atom_pairs );
		return;
	}
	void find_ligand(core::pose::Pose const & pose ){
		//Find the sequence position of the porphyrin
		porphyrin_position_ = 0;
		for ( core::Size i(1); i <= pose.total_residue(); ++i ) {
			if ( pose.residue(i).name3() == "POR" ) {
				porphyrin_position_ = i;
				break;
			}
		}
		runtime_assert ( porphyrin_position_ != 0); //It found the porphyrin residue in the pose
	}
	void setup_task_factory(core::pose::Pose const & pose ){
		core::pack::task::TaskFactoryOP task_factory (new core::pack::task::TaskFactory );
		task_factory->push_back(new core::pack::task::operation::InitializeFromCommandline ); //Will now take resfile, etc.
		//NOTE:: This must be moved to apply because it requires the pose
		if ( basic::options::option[ local::repacking_distance ].value() != 0 ) { //If we are going to relax the structure and restrict repacking
			//Set up the TaskFactory for fast_relax so it will only repack within repacking_distance of heme
			core::Real heme_neighbor_distance = basic::options::option[ local::repacking_distance ].value();

			//Set the operation to restrict repacking distance
			//Create a point graph that stores coordinate information from the pose
			core::conformation::PointGraphOP coords (new core::conformation::PointGraph );
			core::conformation::residue_point_graph_from_conformation(pose.conformation(), *coords);
			//Now find neighbors
			core::conformation::find_neighbors< core::conformation::PointGraphVertexData, core::conformation::PointGraphEdgeData >(coords, heme_neighbor_distance);


			//Create a vector of booleans with the same length as the pose
			utility::vector1 < bool > repack_residue;
			//Initially make all values false
			for ( core::Size resnum = 1; resnum <= pose.total_residue(); ++resnum ) {
				repack_residue[ resnum ] = false;
			}
			//Iterate through all neighbors of ligand
			runtime_assert ( porphyrin_position_ != 0); //It should know the ligand position; if not, something is wrong
			for ( core::conformation::PointGraph::UpperEdgeListConstIter iter = coords->get_vertex( porphyrin_position_ ).const_upper_edge_list_begin(); iter != coords->get_vertex( porphyrin_position_ ).const_upper_edge_list_end(); ++iter ) {
				repack_residue[ iter->upper_vertex() ] = true;

			}
			repack_residue[ porphyrin_position_ ] = true; //Make sure the ligand always repacks no matter what

			//Will need to create a PreventRepacking task operation and then for all non-neighbor residues use include_residue
			core::pack::task::operation::PreventRepacking prev_repack_t_o;
			for ( core::Size resnum = 1; resnum <= pose.total_residue(); ++resnum ) {
				if ( !repack_residue[ resnum ] ) { //If this residue is not allowed to repack
					prev_repack_t_o.include_residue( resnum );
				}
			}
			task_factory->push_back( prev_repack_t_o.clone() );



		}
		//Set this as the tf for fast_relax and the restricted_repack_mover
		fast_relax_->set_task_factory(task_factory);
		restricted_repack_mover_->task_factory( task_factory );

	}
	//Score function
	core::scoring::ScoreFunctionOP score_function_;
	//Heme coordinates
	core::pose::PoseOP heme_coords_;
	//FastRelax mover
	protocols::relax::FastRelaxOP fast_relax_;
	//PackRotamers mover for initial repacking
	protocols::simple_moves::PackRotamersMoverOP initial_repack_mover_;
	//PackRotamers mover for fixed bb ligand replacement
	protocols::simple_moves::PackRotamersMoverOP restricted_repack_mover_;
	//Porphyrin sequence position
	core::Size porphyrin_position_;
};

//Define owning pointers
typedef utility::pointer::shared_ptr< HemeReplace > HemeReplaceOP;

int main( int argc, char* argv[] ){
	basic::options::option.add( local::heme_conformation, "Input pdb file containing only the heme conformation to replace ligand").def("heme.pdb");
	basic::options::option.add( local::no_relax, "Perform fixed bb ligand replacement without relaxing structure").def(false);
	basic::options::option.add( local::num_trials, "Number of trials over which to average ligand replacement scores").def(5);
	basic::options::option.add( local::repacking_distance, "Distance from heme (in Angstroms) at which repacking will be cut off").def(12.0);
	try{
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go( new HemeReplace );

		TR << "***************************d**o**n**e***********************************" << std::endl;
		return 0;
	}
catch (utility::excn::Exception const & e ){
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
