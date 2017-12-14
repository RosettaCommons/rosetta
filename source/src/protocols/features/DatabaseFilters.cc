// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/DatabaseFilters.hh
/// @brief  report atom-atom pair geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/DatabaseFilters.hh>
#include <protocols/features/DatabaseStatements.hh>

//External

// Project Headers
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh> // for total_nr_jobs
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>


// C++ Headers
#include <string>


//Auto Headers
#include <utility/excn/Exceptions.hh>
namespace protocols {
namespace features {

static basic::Tracer TR( "protocols.features.DatabaseFilters" );


core::Real get_current_model_score(core::pose::Pose const & pose, core::Size score_type_id){
	//Set up to pull scores out of the energy map
	core::scoring::Energies const & energies(pose.energies());
	core::scoring::EnergyMap emap;
	core::scoring::ScoreFunctionOP score_function(core::scoring::get_score_function());
	auto score_type = static_cast<core::scoring::ScoreType>(score_type_id);
	utility::vector1< bool > relevant_residues(pose.size(), true);
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose(pose, score_function);
	score_function->get_sub_score(pose,relevant_residues,emap);
	return energies.weights()[score_type] * emap[score_type];
}

DatabaseFilterOP get_DB_filter_ptr(){
	if ( ! basic::options::option[basic::options::OptionKeys::out::database_filter].user() ) {
		return nullptr;// just leave the database_filter pointer null
	}
	utility::vector1<std::string> filter_option=
		basic::options::option[basic::options::OptionKeys::out::database_filter];
	auto begin = filter_option.begin();
	std::string type= *begin;
	++begin;
	utility::vector1<std::string> arguments(begin, filter_option.end());

	if ( type == "TopPercentOfEachInput" ) return DatabaseFilterOP( new TopPercentOfEachInput(arguments) );
	if ( type == "TopPercentOfAllInputs" ) return DatabaseFilterOP( new TopPercentOfAllInputs(arguments) );
	if ( type == "TopCountOfEachInput" ) return DatabaseFilterOP( new TopCountOfEachInput(arguments) );
	if ( type == "TopCountOfAllInputs" ) return DatabaseFilterOP( new TopCountOfAllInputs(arguments) );

	utility_exit_with_message(type+" is not a valid Database Filter name");
	return nullptr; // To keep the compiler happy
}

WriteDeletePair get_write_delete_pair(
	core::pose::Pose const & pose,
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id,
	core::Size const & count,
	std::string const & score_term,
	std::string current_input=""
){
	core::Size n_models = get_current_structure_count(db_session,protocol_id,current_input);

	//store all the structures until you have at least percentile_count worth of models
	if ( n_models < count ) {
		return WriteDeletePair(true, utility::vector1<StructureID>());
	}
	// else we need to delete some structure(s) from the database
	StructureID struct_id;
	core::Real cutoff_score = 0.0;
	core::Real current_model_score = 0.0;
	core::Size score_type_id = 0;
	try
{
		score_type_id = get_score_type_id_from_score_term(db_session,protocol_id,score_term);
	}catch (utility::excn::Exception &)
{
		TR << "no score type term, looking in the job data map" <<std::endl;
	}
	//Some applications (most ligand docking) store score terms in the job data.
	//If this is the case score_type_id will return 0
	//otherwise, we know that the score term in question is a normal scoring function term

	StructureID struct_id_to_remove;

	if ( score_type_id != 0 ) {
		current_model_score = get_current_model_score(pose, score_type_id);

		//Get the structure ID and associated score from the database
		struct_id = get_struct_id_with_nth_lowest_score_from_score_data(db_session,score_type_id,count, protocol_id, current_input);
		cutoff_score = get_score_for_struct_id_and_score_term_from_score_data(db_session,struct_id,score_type_id);
		struct_id_to_remove = get_struct_id_with_highest_score_from_score_data(db_session,score_type_id,protocol_id, current_input);

	} else {
		//get the current score out of the job data map
		std::map< std::string, core::Real > const & srpairmap( protocols::jd2::get_string_real_pairs_from_current_job() );
		if ( srpairmap.count( score_term ) ) {
			current_model_score = srpairmap.at( score_term );
		}

		struct_id = get_struct_id_with_nth_lowest_score_from_job_data(db_session, score_term, count,protocol_id, current_input);
		cutoff_score = get_score_for_struct_id_and_score_term_from_job_data(db_session,struct_id,score_term);

		//See if our current model is better
		struct_id_to_remove = get_struct_id_with_highest_score_from_job_data(db_session,score_term,protocol_id,current_input);
	}
	//See if our current model is better
	if ( current_model_score < cutoff_score ) {
		TR <<"Current score is " <<current_model_score << ". Best score is " << cutoff_score <<
			". Deleting struct_id " <<struct_id_to_remove <<std::endl;
		utility::vector1<StructureID> struct_ids_to_delete(1,struct_id_to_remove);
		return WriteDeletePair(true, struct_ids_to_delete);
	} else {
		TR <<"Current score is " <<current_model_score << ". Best score is " << cutoff_score <<
			". Discarding current pose" <<std::endl;
		return WriteDeletePair(false, utility::vector1<StructureID>());
	}

}

TopPercentOfEachInput::TopPercentOfEachInput(
	utility::vector1<std::string> arguments
):
	DatabaseFilter(),
	top_count_of_each_input_()
{
	if ( arguments.size() != 2 ) {
		utility_exit_with_message("TopPercentOfEachInput option takes 2 arguments");
	}
	top_count_of_each_input_.score_term_ = arguments[1];

	core::Real percent = utility::from_string(arguments[2], core::Real());
	core::Size n_structs = basic::options::option[basic::options::OptionKeys::out::nstruct];
	auto percentile_count = static_cast<core::Size>(floor(percent*n_structs));
	top_count_of_each_input_.count_ = percentile_count;
}

WriteDeletePair
TopPercentOfEachInput::operator()(
	core::pose::Pose const & pose,
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id
){
	return top_count_of_each_input_(pose, db_session, protocol_id);
}

TopPercentOfAllInputs::TopPercentOfAllInputs(utility::vector1<std::string> arguments):
	DatabaseFilter()
{
	if ( arguments.size() != 2 ) {
		utility_exit_with_message("TopPercentOfAllInputs option takes 2 arguments");
	}
	top_count_of_all_inputs_.score_term_ = arguments[1];

	core::Real percent = utility::from_string(arguments[2], core::Real());
	if ( protocols::jd2::jd2_used() ) {
		core::Size total_nr_jobs = protocols::jd2::JobDistributor::get_instance()->total_nr_jobs();
		core::Size percentile_count( floor(percent*total_nr_jobs) );
		top_count_of_all_inputs_.count_ = percentile_count;
	} else {
		top_count_of_all_inputs_.count_ = 0;
	}
}

WriteDeletePair
TopPercentOfAllInputs::operator()(
	core::pose::Pose const & pose,
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id
){
	return top_count_of_all_inputs_(pose, db_session, protocol_id);
}

TopCountOfEachInput::TopCountOfEachInput():
	DatabaseFilter()
{}

TopCountOfEachInput::TopCountOfEachInput(utility::vector1<std::string> arguments):
	DatabaseFilter(),
	count_(0)
{
	if ( arguments.size() != 2 ) {
		utility_exit_with_message("TopCountOfEachInput option takes 2 arguments");
	}
	score_term_ = arguments[1];
	count_ = utility::from_string(arguments[2], core::Size());
}

WriteDeletePair
TopCountOfEachInput::operator()(
	core::pose::Pose const & pose,
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id
){
	std::string current_input( protocols::jd2::current_input_tag() );
	return get_write_delete_pair(pose, db_session, protocol_id, count_, score_term_, current_input);
}

TopCountOfAllInputs::TopCountOfAllInputs():
	DatabaseFilter()
{}

TopCountOfAllInputs::TopCountOfAllInputs(utility::vector1<std::string> arguments):
	DatabaseFilter(),
	count_(0)
{
	if ( arguments.size() != 2 ) {
		utility_exit_with_message("TopCountOfAllInputs option takes 2 arguments");
	}
	score_term_ = arguments[1];
	count_ = utility::from_string(arguments[2], core::Size());
}

WriteDeletePair
TopCountOfAllInputs::operator()(
	core::pose::Pose const & pose,
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id
){
	return get_write_delete_pair(pose, db_session, protocol_id, count_, score_term_);
}

} // features
} // protocols

