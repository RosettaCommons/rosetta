// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/movers/AssemblyMover.cc
/// @brief an interface for making Movers that deal with Assemblies
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <protocols/sewing/movers/AssemblyMover.hh>
#include <protocols/sewing/movers/AssemblyMoverCreator.hh>

//Package headers
#include <protocols/sewing/scoring/AssemblyScorerFactory.hh>
#include <protocols/sewing/scoring/AssemblyScorerCreator.hh>
#include <protocols/sewing/scoring/AssemblyScorer.hh>
#include <protocols/sewing/requirements/AssemblyRequirement.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <protocols/sewing/requirements/AssemblyRequirementCreator.hh>
#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <protocols/sewing/hashing/EdgeMapGenerator.hh>

#include <protocols/moves/mover_schemas.hh>
// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
static basic::Tracer TR( "protocols.sewing.movers.AssemblyMover" );

namespace protocols {
namespace sewing {
namespace movers {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
AssemblyMover::AssemblyMover():
	protocols::moves::Mover( AssemblyMover::class_name() )
{
	//Don't initialize from options here (RosettaScripts)
	//basis_map_generator_ = hashing::BasisMapGeneratorOP ( new hashing::BasisMapGenerator );
	min_cycles_ = 1000;
	max_cycles_ = 1000;
	start_temperature_ = 1.0;
	end_temperature_ = 0.01;
	add_probability_ = 0.05;
	delete_probability_ = 0.005;
	conformer_switch_probability_ = 0;
	output_pose_per_move_ = false;
	output_pose_count_ = 0;
}
std::string
AssemblyMover::class_name(){
	return "AssemblyMover";
}

///////////////////////////////////////////////////////////////////////////// ///@brief Copy constructor
AssemblyMover::AssemblyMover( AssemblyMover const & src ):
	protocols::moves::Mover( src )
{
	edge_map_generator_ = src.get_edge_map_generator();
	basis_map_generator_ = src.get_basis_map_generator();
	model_file_name_ = src.get_model_file_name();
	edge_file_name_ = src.get_edge_file_name();
	assembly_scorers_ = src.get_assembly_scorers();
	requirements_ = src.get_requirements();
	max_num_segments_ = src.get_max_num_segments();
	hashed_ = src.get_hashed();
	max_segment_length_ = src.get_max_segment_length();
	min_cycles_ = src.get_min_cycles();
	max_cycles_ = src.get_max_cycles();
	add_probability_ = src.get_add_probability();
	delete_probability_ = src.get_delete_probability();
	conformer_switch_probability_ = src.get_conformer_switch_probability();
	start_temperature_ = src.get_start_temperature();
	end_temperature_ = src.get_end_temperature();
	segment_vector_ = src.get_segment_vector();
	output_pose_per_move_ = src.get_output_pose_per_move();
	output_pose_count_ = src.get_output_pose_count();
	window_width_ = src.get_window_width();
	recover_lowest_assembly_ = src.get_recover_lowest_assembly();

}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////



void
AssemblyMover::set_up_hashing(){
	if ( hashed_ ) {
		edge_map_generator_ = hashing::EdgeMapGeneratorOP(new hashing::EdgeMapGenerator(edge_file_name_));
		basis_map_generator_ = hashing::BasisMapGeneratorOP( new hashing::BasisMapGenerator );
		basis_map_generator_->segment_vector( segment_vector_ );
		basis_map_generator_->set_edge_file_reader( edge_map_generator_ );

	}
}
/// @brief Apply the mover
void
AssemblyMover::apply( core::pose::Pose& pose){
	//Try to generate an assembly
	core::Size starttime = time(NULL);

	data_storage::SmartAssemblyOP assembly = set_up_assembly( pose );
	if ( max_cycles_ > 0 ) {
		generate_assembly(assembly,pose);
	}
	//If we failed, set to FAIL_RETRY
	if ( assembly == nullptr ) {
		TR << "Failed to generate an Assembly" << std::endl;
		set_last_move_status(protocols::moves::FAIL_RETRY);
		return;
	}

	//Now convert the assembly to a pose
	//(not actually sure how this function will look)

	core::Size endttime = time(NULL);
	TR << "Assembly successfully generated in " << endttime - starttime << " seconds" << std::endl;

	print_statistics( assembly );

	pose = assembly->to_pose("fa_standard");

	add_motif_scorers_to_score_file( pose, assembly );


	//At some point in here we'll need to output the side chain coordinates
	//Again, no idea what this function will look like
	assembly->dump_side_chains();

	//
}

data_storage::SmartAssemblyOP
AssemblyMover::set_up_assembly( core::pose::Pose &){

	data_storage::SmartAssemblyOP assembly;
	if ( hashed_ ) {
		data_storage::HashedSmartAssemblyOP hashed_assembly( new data_storage::HashedSmartAssembly( segment_vector_ ) );
		hashed_assembly->set_basis_map_generator( basis_map_generator_->clone() );
		assembly = hashed_assembly;
	} else {
		assembly = data_storage::SmartAssemblyOP( new data_storage::SmartAssembly( segment_vector_, window_width_ ) );
	}
	assembly->pick_random_starting_segment();
	return assembly;

}

void
AssemblyMover::print_statistics( data_storage::SmartAssemblyOP assembly) const{

	TR <<"N-terminal segment is: " << assembly->get_n_terminal_segment()->get_segment_id() << std::endl;
	TR <<"C-terminal segment is: " << assembly->get_c_terminal_segment()->get_segment_id() << std::endl;
	TR << "Current Forward Assembly is: " << assembly->get_comprehensive_forward_assembly() << std::endl;
	TR << "Current Reverse Assembly is: " << assembly->get_comprehensive_reverse_assembly() << std::endl;
	TR << "Assembly size is " << assembly->get_size() << " segments." << std::endl;
	TR << "Assembly length is " << assembly->get_length() << " residues." << std::endl;

}
void
AssemblyMover::add_motif_scorers_to_score_file( core::pose::Pose & pose, data_storage::SmartAssemblyOP & assembly){
	/*
	for( std::pair< std::string, std::pair< core::Real, scoring::AssemblyScorerOP > > scorer: assembly_scorers_ ){
	//temp_score_file << scorer.second.second->get_last_score() << " ";
	core::pose::setPoseExtraScore( pose, scorer.first, scorer.second.second->get_last_score() );
	}
	*/
	core::Size counter = 1;
	core::Real latest_score = score( assembly ); //Makes sure all the scores are up to date
	for ( scoring::AssemblyScorerOP scorer : assembly_scorers_ ) {
		//core::pose::setPoseExtraScore( pose, scorer->get_name(), scorer->get_last_score() );
		core::pose::setPoseExtraScore( pose, scorer->get_name() + "_" + utility::to_string( counter ), scorer->get_last_score() );
		++counter;
	}
	core::pose::setPoseExtraScore( pose, "nres", pose.total_residue() );
	//core::Real sew_tot_score = score(asembly);
	core::pose::setPoseExtraScore( pose, "sew_total", latest_score );
}



void
AssemblyMover::generate_assembly(data_storage::SmartAssemblyOP assembly, core::pose::Pose& pose){
	score_ = 1000;
	lowest_score_ = 1000;
	bool revertable_move = true;
	bool revert = false;
	bool continue_assembly = true;
	bool scored_this_cycle = false;
	core::Size cycle = 0;
	output_pose_count_ = 0;
	//assembly->set_basis_map_generator(basis_map_generator_);
	while ( continue_assembly ) {
		scored_this_cycle = false;
		cycle = cycle+1;
		old_score_ = score_;
		TR << "Cycle: " << cycle << std::endl;
		print_statistics( assembly );
		// short-circuiting logic happens here
		if ( !assembly->can_delete() ) {
			TR.Debug << "N terminus " << assembly->get_n_terminal_segment()->get_segment_id() << "returns " << assembly->get_n_terminal_segment()->is_vital() << std::endl;
			TR.Debug << "C terminus " << assembly->get_c_terminal_segment()->get_segment_id() << "returns " << assembly->get_c_terminal_segment()->is_vital() << std::endl;
			TR.Debug << "All termini vital. Adding." << std::endl;
			current_add_probability_ = 1.0 - conformer_switch_probability_;
		} else if ( assembly->get_size() >= max_num_segments_ ) {
			TR.Debug << "Assembly at max size. Cannot add." << std::endl;
			current_add_probability_ = 0.0;
		} else {
			current_add_probability_ = add_probability_;
		}
		random_action_ = numeric::random::rg().uniform();
		//TR << "random action " << random_action_ << std::endl;
		// and now we actually move
		if ( random_action_ < current_add_probability_ ) {
			TR.Debug << "Adding a segment" << std::endl;
			if ( not(assembly->add_segment(assembly->get_modifiable_terminus('A'))) ) { // add happens here
				revertable_move = false;
				cycle--;
			} else {
				revertable_move = true;
			}
		} else if ( random_action_ < (current_add_probability_+delete_probability_) ) {
			TR.Debug << "Deleting a segment" << std::endl;
			if ( not(assembly->delete_segment(assembly->get_modifiable_terminus('D'))) ) { //delete happens here
				revertable_move = false;
				cycle--;
			} else {
				revertable_move = true;
			}
		} else if ( random_action_ < (current_add_probability_ + delete_probability_ + conformer_switch_probability_ ) ) {
			TR.Debug << "Switching ligand conformers" << std::endl;
			if ( !( assembly->sample_ligand() ) ) {
				revertable_move = false;
				cycle--;
			} else {
				revertable_move = true;
			}
		} else {
			TR.Debug << "Switching a segment" << std::endl;
			if ( not(assembly->switch_segment(assembly->get_modifiable_terminus('S'))) ) {
				revertable_move = false;
				cycle--;
			} else {
				revertable_move = true;
			}
		}
		if ( output_pose_per_move_ ) {
			this->output_pose( assembly , pose);//pose will only be output if output_pose_per_move is set
		}
		revert = false;
		if ( revertable_move ) {
			for ( core::Size i=1; i<=requirements_.size(); i++ ) {
				if ( !requirements_[i]->test(assembly).first ) {
					TR.Debug << "Assembly failed requirement " << requirements_[ i ]->get_name() << std::endl;
					revert = true;
					break;
				}
			}
			if ( !revert ) {
				score_ = score(assembly);
				scored_this_cycle = true;
				TR.Debug << "Old assembly score: " << old_score_ << std::endl;
				TR.Debug << "New assembly score: " << score_ << std::endl;
				core::Real current_temperature = std::max(0.001, ((end_temperature_ - start_temperature_)/min_cycles_)*cycle + start_temperature_ );
				// TR << "temperature " << current_temperature << "   cycle  "  << cycle << std::endl;
				//if(score_ > old_score_ && numeric::random::rg().uniform() > std::exp( std::min (40.0, std::max(-40.0,((old_score_ - score_)/temperature_))))){
				//  revert = true;
				//}
				if ( score_ > old_score_ && numeric::random::rg().uniform() > std::exp( std::min (40.0, std::max(-40.0,((old_score_ - score_)/current_temperature)))) ) {
					revert = true;
				} else if ( score_ < lowest_score_ && this->assembly_meets_requirements( assembly ) ) {

					// lowest_scoring_assembly_ = assembly->get_abbreviated_assembly();
					// lowest_score_ = score_;
					lowest_score_ = score_;
					best_assembly_ = data_storage::SmartAssemblyOP( new data_storage::SmartAssembly(*assembly));
				}
			}
		}
		//print_statistics( assembly );
		if ( revert ) {
			TR.Debug << "Reverting change." << std::endl;
			//print_statistics( assembly );
			if ( scored_this_cycle ) {
				/*    for ( std::pair< std::string, std::pair< core::Real, scoring::AssemblyScorerOP > > entry : assembly_scorers_ ){
				entry.second.second->set_last_score( entry.second.second->get_old_last_score() );
				}
				*/    revert_score();
			}
			assembly->revert();
			for ( core::Size i=1; i<=requirements_.size(); i++ ) {
				if ( !requirements_[i]->test(assembly).first ) {
					TR << "WARNING: Assembly failed requirement on revert!" << requirements_[ i ]->get_name() << std::endl;
					// revert = true;
					// break;
				}
			}
			score_ = old_score_;
			if ( output_pose_per_move_ ) {
				this->output_pose( assembly , pose);//pose will only be output if output_pose_per_move is set
			}
			//print_statistics( assembly );
		}
		//TR << "Does the assembly have chain breaks?: " << !assembly->is_continuous() << std::endl;
		//TR << "Cycle: " << cycle << std::endl;
		// TR << "END CYCLE" << std::endl;
		//now see if we should run another cycle
		if ( cycle > max_cycles_ ) {
			continue_assembly = false;
		} else if ( cycle > min_cycles_ ) {
			continue_assembly = !this->assembly_meets_requirements( assembly );
			//   continue_assembly = false;
			//   for(core::Size i=1;i<=requirements_.size();i++){
			//    if( !requirements_[i]->test(assembly).second ){
			//     continue_assembly = true;
			//     break;
			//    }
			//   }
		}
	}
	// The following commented out lines implement how to recover the lowest scoring assembly. However, in append assembly runs with single helices it causes isolated chimaeras during theadd phase of the recovery. I wasn't able to track the bug, but my assumption is it occurs during the deletion phase.
	if ( score_ > lowest_score_ && recover_lowest_assembly_ ) {
		TR << "Current score ( " << score_ << " ) is greater than lowest score ( " << lowest_score_ << " )" << std::endl;
		// print_statistics( assembly );
		TR << "Recovering lowest scoring assembly" << std::endl;
		assembly = best_assembly_;
		//assembly->reconstitute_assembly_from_string(lowest_scoring_assembly_string_);
		//assembly->revert_to_abbreviated_assembly( lowest_scoring_assembly_ );
		old_score_ = score_;
		score_ = lowest_score_;
	} else {
		TR << "Final assembly is lowest scoring assembly" << std::endl;
		//assembly->reconstitute_assembly_from_string(lowest_scoring_assembly_string_);
	}

	TR << "Final score: " << score_ << "\t";
	TR << "Old score: " << old_score_ << std::endl;
	TR << "Final Stats: " << std::endl;
	//print_statistics( assembly );
}
////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
AssemblyMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

/// @brief Get the name of the Mover
std::string
AssemblyMover::get_name() const{
	return "AssemblyMover";
}


////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
AssemblyMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & filtermap,
	protocols::moves::Movers_map const & movermap,
	core::pose::Pose const & pose)
{
	//get the model file name
	min_cycles_ = tag->getOption<core::Size>("minimum_cycles",10000);
	max_cycles_ = tag->getOption<core::Size>("maximum_cycles",100000);
	start_temperature_ = tag->getOption<core::Real>("start_temperature",0.6);
	end_temperature_ = tag->getOption<core::Real>("end_temperature",0.6);
	add_probability_ = tag->getOption<core::Real>("add_probability",0.05);
	delete_probability_ = tag->getOption<core::Real>("delete_probability",0.005);
	conformer_switch_probability_ = tag->getOption< core::Real >( "conformer_switch_probability", 0 );
	output_pose_per_move_ = tag->getOption<bool>("output_pose_per_move", false);
	recover_lowest_assembly_ = tag->getOption<bool>("recover_lowest_assembly", true);
	if ( tag->hasOption( "model_file_name" ) ) {
		model_file_name_ = tag->getOption<std::string>("model_file_name");
	} else {
		utility_exit_with_message( "Please provide a model file name in your AssemblyMover tag! (i.e. model_file_name=\"my_models.models\")" );
	}
	if ( tag->hasOption( "edge_file_name" ) ) {
		edge_file_name_ = tag->getOption< std::string>( "edge_file_name", "" );
	}
	if ( tag->hasOption( "hashed" ) ) {
		hashed_ = tag->getOption< bool >( "hashed", false );
		if ( tag->hasOption( "window_width" ) ) {
			TR << "Warning: The specified window width will not be used. AssemblyMover will use the hash score cutoff specified in the edge file." << std::endl;
		}
	}

	max_num_segments_ = tag->getOption< core::Size >( "max_segments", 100 );
	//Are we in any way checking/enforcing this during assembly? A: No, we are not. We will, though, with a DsspSpecificLengthRequirement.
	max_segment_length_ = tag->getOption< core::Size >( "max_segment_length", 100 );
	window_width_ = tag->getOption< core::Size >( "window_width", 4 ); //This will only be used in hashless sewing
	// core::Size recursive_depth = tag->getOption< core::Size>( "max_recursion", 1 );
	// utility::vector1< core::Size > dummy_vector;
	// dummy_vector.clear();
	// basis_map_generator_->set_alignment_settings( dummy_vector, dummy_vector, dummy_vector, recursive_depth );
	//Check for AssemblyScorers subtag and Requirements subtag
	//If we have an AssemblyScorers subtag, parse scorers
	if ( tag->hasTag( "AssemblyScorers" ) ) {
		parse_assembly_scorers( tag->getTag( "AssemblyScorers" ), datamap, filtermap, movermap, pose );
	} else {
		//Otherwise set to default score function
		set_default_assembly_scorers();
	}
	//If we have a Requirements tag, parse requirements
	if ( tag->hasTag( "AssemblyRequirements" ) ) {
		parse_requirements( tag->getTag( "AssemblyRequirements" ), datamap, filtermap, movermap, pose ) ;
	} else {
		//Otherwise set to default requirements
		set_default_requirements();
	}
	segment_vector_ = hashing::ModelFileReader::read_model_file( model_file_name_ );
	set_up_hashing();


}

void
AssemblyMover::parse_requirements(
	utility::tag::TagCOP requirements_tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & filtermap,
	protocols::moves::Movers_map const & movermap,
	core::pose::Pose const & pose)
{
	//Each requirement will need to have its own parse_tag function that returns an instance of that requirement with the appropriate options set
	requirements::AssemblyRequirementOP current_requirement;
	for ( utility::tag::TagCOP requirement: requirements_tag->getTags() ) {
		current_requirement = requirements::AssemblyRequirementFactory::get_instance()->get_requirement( requirement->getName() );
		current_requirement->set_options_from_tag( requirement, datamap, filtermap, movermap, pose );
		add_requirement( current_requirement );
	}
}



void
AssemblyMover::parse_assembly_scorers(
	utility::tag::TagCOP scorers_tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & filtermap,
	protocols::moves::Movers_map const & movermap,
	core::pose::Pose const & pose)
{
	//Loop through all of the subtags
	scoring::AssemblyScorerOP current_scorer;
	for ( utility::tag::TagCOP scorer: scorers_tag->getTags() ) {
		//Create the scorer
		//  core::Real weight = scorer->getOption<core::Real>( "weight" );
		current_scorer = scoring::AssemblyScorerFactory::get_instance()->get_assembly_scorer( scorer->getName() );
		current_scorer->set_options_from_tag( scorer, datamap, filtermap, movermap, pose );
		// *****Later scorers may have additional functionality besides just a weight--give them functions to parse themselves to accommodate this!!******
		/*  //Get the weight--if none is provided, will get exception
		if( !scorer->hasOption( "weight" )){
		TR << "No weight provided for scorer " << scorer->getName() << std::endl;
		utility_exit_with_message( "Weights must be provided for all AssemblyScorers!" );
		}
		//Add the scorer to the map
		add_scorer( current_scorer, weight );
		*/
		add_scorer( current_scorer );
	}

}
void
AssemblyMover::set_default_assembly_scorers()
{
	//Current default only has MotifScorer and InterModelMotifScorer
	add_scorer( scoring::AssemblyScorerFactory::get_instance()->get_assembly_scorer( "MotifScorer" ) );
	add_scorer( scoring::AssemblyScorerFactory::get_instance()->get_assembly_scorer( "InterModelMotifScorer" ) );
	// add_scorer( scoring::AssemblyScorerFactory::get_instance()->get_assembly_scorer( "MotifScorer" ), 1.0 );
	// add_scorer( scoring::AssemblyScorerFactory::get_instance()->get_assembly_scorer( "InterModelMotifScorer" ), 10.0 );
}

bool
AssemblyMover::assembly_meets_requirements( data_storage::SmartAssemblyOP assembly ) {
	for ( core::Size i=1; i<=requirements_.size(); i++ ) {
		if ( !requirements_[i]->test(assembly).second ) {
			return false;
		}
	}
	return true;
}



void
AssemblyMover::set_default_requirements()
{
	//By default, the only requirement we should have is a ClashRequirement
	//Ideally, we should create this with the factory as before
	add_requirement( requirements::AssemblyRequirementFactory::get_instance()->get_requirement( "ClashRequirement" ) );
	TR << "No requirements specified. Adding default requirement" << std::endl;
	TR << "ClashRequirement clash_radius=5.0 maximum_clashes_allowed=0" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AssemblyMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new AssemblyMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AssemblyMover::clone() const
{
	return protocols::moves::MoverOP( new AssemblyMover( *this ) );
}

///@brief Add requirements to the mover from the outside (not sure if we'll need this?
void
AssemblyMover::add_requirement( requirements::AssemblyRequirementOP requirement )
{
	requirements_.push_back( requirement );
}



std::ostream &
operator<<( std::ostream & os, AssemblyMover const & mover )
{
	mover.show(os);
	return os;
}


//Getters


//Setters



















//XML Schema





////////////////////////////////////////////////////////////////////////////////
/// Creator ///
///////////////

/////////////// Creator ///////////////


protocols::moves::MoverOP
AssemblyMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AssemblyMover );
}

std::string
AssemblyMoverCreator::keyname() const
{
	return AssemblyMover::mover_name();
}


//not really sure why this function exists
std::string
AssemblyMoverCreator::class_name()
{
	return AssemblyMover::mover_name();
}
std::string
AssemblyMover::mover_name()
{
	return "AssemblyMover";
}


//Getters
core::Size
AssemblyMover::get_window_width() const{
	return window_width_;
}

hashing::SegmentVectorCOP
AssemblyMover::get_segment_vector() const{
	return segment_vector_;
}

hashing::BasisMapGeneratorOP
AssemblyMover::get_basis_map_generator() const{
	return basis_map_generator_;
}

hashing::EdgeMapGeneratorOP
AssemblyMover::get_edge_map_generator() const{
	return edge_map_generator_;
}

std::string
AssemblyMover::get_model_file_name() const{
	return model_file_name_;
}

std::string
AssemblyMover::get_edge_file_name() const{
	return edge_file_name_;
}

// ScorerMap
utility::vector1< scoring::AssemblyScorerOP >
AssemblyMover::get_assembly_scorers() const{
	return assembly_scorers_;
}

utility::vector1< requirements::AssemblyRequirementOP >
AssemblyMover::get_requirements() const{
	return requirements_;
}

bool
AssemblyMover::get_hashed() const{
	return hashed_;
}

core::Real
AssemblyMover::get_start_temperature() const{
	return start_temperature_;
}

core::Real
AssemblyMover::get_end_temperature() const{
	return end_temperature_;
}

core::Size
AssemblyMover::get_max_num_segments() const{
	return max_num_segments_;
}

core::Size
AssemblyMover::get_max_segment_length() const{
	return max_segment_length_;
}

core::Size
AssemblyMover::get_min_cycles() const{
	return min_cycles_;
}

core::Size
AssemblyMover::get_max_cycles() const{
	return max_cycles_;
}

core::Real
AssemblyMover::get_add_probability() const{
	return add_probability_;
}

core::Real
AssemblyMover::get_delete_probability() const{
	return delete_probability_;
}

core::Real
AssemblyMover::get_conformer_switch_probability() const{
	return conformer_switch_probability_;
}

bool
AssemblyMover::get_output_pose_per_move() const{
	return output_pose_per_move_;
}

core::Size
AssemblyMover::get_output_pose_count() const{
	return output_pose_count_;
}

core::Real
AssemblyMover::get_lowest_score() const{
	return lowest_score_;
}

//Setters
void
AssemblyMover::set_window_width( core::Size width ){
	window_width_ = width;
}

void
AssemblyMover::set_segment_vector( hashing::SegmentVectorCOP segvec ){
	segment_vector_ = segvec;
}

void
AssemblyMover::set_basis_map_generator( hashing::BasisMapGeneratorOP bmg ){
	basis_map_generator_ = bmg;
}

void
AssemblyMover::set_edge_map_generator( hashing::EdgeMapGeneratorOP emg ){
	edge_map_generator_ = emg;
}

void
AssemblyMover::set_model_file_name( std::string models ){
	model_file_name_ = models;
}

void
AssemblyMover::set_edge_file_name( std::string edges ){
	edge_file_name_ = edges;
}

void
AssemblyMover::set_assembly_scorers( utility::vector1< scoring::AssemblyScorerOP > assembly_scorers ){
	assembly_scorers_ = assembly_scorers;
}

void
AssemblyMover::set_requirements( utility::vector1 <requirements::AssemblyRequirementOP > requirements ){
	requirements_ = requirements;
}

void
AssemblyMover::set_hashed( bool hashed ){
	hashed_ = hashed;
}

void
AssemblyMover::set_start_temperature( core::Real temp ){
	start_temperature_ = temp;
}

void
AssemblyMover::set_end_temperature( core::Real temp ){
	end_temperature_ = temp;
}

void
AssemblyMover::set_max_num_segments( core::Size max_segs ){
	max_num_segments_ = max_segs;
}

void
AssemblyMover::set_max_segment_length( core::Size max_seg_length ){
	max_segment_length_ = max_seg_length;
}

void
AssemblyMover::set_min_cycles( core::Size min_cycles ){
	min_cycles_ = min_cycles;
}

void
AssemblyMover::set_max_cycles( core::Size max_cycles ){
	max_cycles_ = max_cycles;
}

void
AssemblyMover::set_add_probability( core::Real add_probability ){
	add_probability_ = add_probability;
}

void
AssemblyMover::set_delete_probability( core::Real delete_probability ){
	delete_probability_ = delete_probability;
}

void
AssemblyMover::set_conformer_switch_probability( core::Real prob ){
	conformer_switch_probability_ = prob;
}

void
AssemblyMover::set_output_pose_per_move( bool output_pose_per_move ){
	output_pose_per_move_ = output_pose_per_move;
}

void
AssemblyMover::output_pose( data_storage::SmartAssemblyOP assembly, core::pose::Pose& pose){
	pose = assembly->to_pose( "fa_standard" );
}

void
AssemblyMover::set_lowest_score( core::Real lowest_score ){
	lowest_score_ = lowest_score;
}
void
AssemblyMover::set_recover_lowest_assembly(bool new_setting){
	recover_lowest_assembly_ = new_setting;
}
bool
AssemblyMover::get_recover_lowest_assembly() const{
	return recover_lowest_assembly_;
}






//Scoring methods

//Main method for scoring an assembly
core::Real
AssemblyMover::score( data_storage::SmartAssemblyOP assembly ){
	core::Real total_score = 0;
	//Iterate over all score terms
	/*
	for ( std::pair< std::string, std::pair< core::Real, scoring::AssemblyScorerOP > > entry : assembly_scorers_ ){
	//Add the weighted score
	total_score = total_score +  ( entry.second.first * entry.second.second->score( assembly ) );
	}
	return total_score;
	*/
	for ( scoring::AssemblyScorerOP scorer : assembly_scorers_ ) {
		total_score = total_score + ( scorer->score( assembly ) * scorer->get_weight() );
	}
	return total_score;
}

void
AssemblyMover::revert_score(){
	for ( scoring::AssemblyScorerOP scorer : assembly_scorers_ ) {
		scorer->set_last_score( scorer->get_old_last_score() );
	}
}


///@details  We give add_scorer the scorer and not the name to accommodate possible future cases in which the scorer has additional options that we may want to set and incase there are multiple copies. Returns the index of the added scorer.
core::Size
//AssemblyMover::add_scorer( scoring::AssemblyScorerOP scorer, core::Real weight)
AssemblyMover::add_scorer( scoring::AssemblyScorerOP scorer )
{
	// std::string new_scorer_name = scorer->get_name();
	// //Add the scorer to the ScorerMap
	// assembly_scorers_[ new_scorer_name ] = std::make_pair( weight, scorer );
	assembly_scorers_.push_back( scorer );
	return assembly_scorers_.size();
}

core::Real
//AssemblyMover::get_unweighted_score_term( data_storage::SmartAssemblyOP assembly, std::string scoretype )
AssemblyMover::get_unweighted_score_term( data_storage::SmartAssemblyOP assembly, core::Size scorer_index )
{
	//return assembly_scorers_[ scoretype ].second->score( assembly );
	return assembly_scorers_[ scorer_index ]->score( assembly );
}


void
AssemblyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	//Define xsd for our AssemblyScorers using utility function
	scoring::AssemblyScorerFactory::get_instance()->define_assembly_scorer_subtag( xsd );
	//Define xsd for our Requirements using utility function
	requirements::AssemblyRequirementFactory::get_instance()->define_assembly_requirement_subtag( xsd );

	//Define the subelement list for the mover (consists of AssemblyScorers and AssemblyRequirements )
	XMLSchemaSimpleSubelementList mover_subelements;
	mover_subelements
		.add_already_defined_subelement( "AssemblyScorers", & assembly_mover_subtag_ct_namer )
		.add_already_defined_subelement( "AssemblyRequirements", & assembly_mover_subtag_ct_namer );

	AttributeList attributes;
	define_generic_assembly_mover_attributes( attributes );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.add_attributes( attributes )
		.add_optional_name_attribute( "Name to identify this mover" )
		.set_subelements_single_appearance_optional( mover_subelements )
		.description( "Basic SEWING mover for generating assemblies using random substructures. Uses Monte Carlo sampling and scores based on motif score measuring potential packing interactions." )
		.element_name( mover_name() )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}

std::string
AssemblyMover::assembly_mover_subtag_ct_namer( std::string tag_name ){
	return "assembly_mover_" + tag_name + "_complex_type";
}



void
AssemblyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	AssemblyMover::provide_xml_schema( xsd );
}

void
AssemblyMover::define_generic_assembly_mover_attributes( utility::tag::AttributeList & attributes ){
	using namespace utility::tag;
	//Define non negative real type

	/*XMLSchemaRestriction non_negative_float;
	non_negative_float.name( "non_negative_float" );
	non_negative_float.base_type( xsct_real );
	non_negative_float.add_restriction( xsr_minInclusive, "0" );
	xsd.add_top_level_element( non_negative_float );*/

	//Add attributes

	attributes
		+ XMLSchemaAttribute::attribute_w_default( "start_temperature", xsct_real, "Temperature at start of simulated annealing", "0.6")
		+ XMLSchemaAttribute::attribute_w_default( "end_temperature", xsct_real, "Temperature at end of simulated annealing", "0.6")
		+ XMLSchemaAttribute::attribute_w_default( "add_probability", xsct_real, "Probability of adding a triplet of segments at any given step during assembly", "0.05" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_probability", xsct_real,  "Probability of deleting a terminal triplet of segments at any given step during assembly", "0.005" )
		+ XMLSchemaAttribute::attribute_w_default( "conformer_switch_probability", xsct_real,  "Probability of switching ligand conformers during assembly. This should only be used if a ligand is present AND if you have provided conformers for that ligand.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "window_width", xsct_positive_integer, "Required number of overlapping residues for two segments to be considered a match. Used in hashless SEWING only (for hashed SEWING, this is determined by the hasher settings used when generating the edge file).", "4" )
		//  + XMLSchemaAttribute::attribute_w_default( "cycles", xsct_non_negative_integer,  "Number of Monte Carlo cycles for assembly", "1000" )
		+ XMLSchemaAttribute::attribute_w_default( "minimum_cycles", xsct_non_negative_integer,  "Minimum number of Monte Carlo cycles for assembly before completion requirements are checked.", "10000" )
		+ XMLSchemaAttribute::attribute_w_default( "maximum_cycles", xsct_non_negative_integer,  "Maximum number of Monte Carlo cycles for assembly before forced termination.", "100000" )
		+ XMLSchemaAttribute::required_attribute( "model_file_name", xs_string, "Path to file defining segments to use during assembly" )
		+ XMLSchemaAttribute::attribute_w_default( "hashed", xsct_rosetta_bool,  "Use the hasher during assembly to check overlap of all atoms? Requires an input edge file.", "false" )
		+ XMLSchemaAttribute( "edge_file_name", xs_string, "Path to edge file to use during assembly (only used if hashed is set to true)" )
		+ XMLSchemaAttribute::attribute_w_default( "max_segments", xsct_non_negative_integer,  "Maximum number of segments to include in the final assembly", "100" )
		+ XMLSchemaAttribute::attribute_w_default( "max_segment_length", xsct_non_negative_integer,  "Maximum number of residues to include in a segment", "100")
		+ XMLSchemaAttribute::attribute_w_default( "output_pose_per_move", xsct_rosetta_bool,  "Setting to true will output a pose after each move/revert.", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "recover_lowest_assembly", xsct_rosetta_bool,  "Setting to true will output the lowest assembly in the final pose", "true" );

}









////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

} //protocols
} //sewing
} //movers

