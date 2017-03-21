// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/qsarOptFunc.cc
/// @author Sam DeLuca


//This is a bit of a mess, probably refactor it some time

#include <protocols/qsar/qsarOptFunc.hh>
#include <basic/database/sql_utils.hh>
#include <numeric/roc_curve.hh>
#include <utility>
#include <utility/exit.hh>

//Auto Headers
#include <utility/vector1.hh>

namespace protocols {
namespace qsar {

qsarOptFunc::qsarOptFunc(
	utility::sql_database::sessionOP db_session,
	core::optimization::Multivec const & initial_values,
	std::map<std::string,core::Size>  grid_indices) : core::optimization::Multifunc(), initial_values_(initial_values), grid_indices_(std::move(grid_indices)), cutoff_(0.0)
{
	std::string value_string =
		"SELECT job_string_real_data.data_value\n"
		"\tFROM job_string_real_data\n"
		"\tWHERE\n"
		"\t\tjob_string_real_data.data_key = ?\n"
		"\tAND\n"
		"\t\tjob_string_real_data.struct_id = ?;";

	score_selection_ = basic::database::safely_prepare_statement(value_string,db_session);

	std::string struct_id_string = "SELECT structures.struct_id FROM structures;";

	struct_id_selection_ = basic::database::safely_prepare_statement(struct_id_string,db_session);

	std::string tag_activity_string =
		"SELECT structures.tag, structure_activity.activity\n"
		"\tFROM structures\n"
		"\tINNER JOIN structure_activity ON structures.input_tag = structure_activity.input_tag\n"
		"\tWHERE structures.struct_id = ?;";

	tag_activity_selection_ =  basic::database::safely_prepare_statement(tag_activity_string,db_session);

}

void qsarOptFunc::setup_data_map()
{
	data_map_.clear();

	cppdb::result struct_id_result(basic::database::safely_read_from_database(struct_id_selection_));
	while ( struct_id_result.next() )
			{
		core::Size struct_id;
		struct_id_result >> struct_id;
		data_map_.push_back(get_struct_data(struct_id));
	}

}

void qsarOptFunc::set_initial_values(core::optimization::Multivec const & initial_values)
{
	initial_values_ = initial_values;
}

core::Real qsarOptFunc::operator() (core::optimization::Multivec const & vars) const
{
	debug_assert(vars.size() == grid_indices_.size());

	numeric::RocCurve roc_curve;

	for ( auto const & data_it : data_map_ ) {

		core::Real total_score = 0.0;

		for ( const auto & score_it : data_it.score_map ) {
			core::Size vec_index = grid_indices_.find(score_it.first)->second;
			core::Real initial_weight = initial_values_[vec_index];
			core::Real current_weight = vars[vec_index];
			core::Real component_score = score_it.second;

			total_score += (component_score/initial_weight)*current_weight;
		}

		bool predicted(total_score < cutoff_);
		roc_curve.insert_point(predicted,data_it.activity,data_it.tag,total_score);
	}

	roc_curve.generate_roc_curve();
	return roc_curve.calculate_auc();

}

void qsarOptFunc::dfunc( core::optimization::Multivec const &, core::optimization::Multivec & ) const
{
	utility_exit_with_message("haven't implemented dfunc sorry bye");
}

void qsarOptFunc::dump( core::optimization::Multivec const &, core::optimization::Multivec const & ) const
{
	utility_exit_with_message("haven't implemented dump sorry bye");
}

qsarOptData qsarOptFunc::get_struct_data(core::Size const & struct_id )
{

	tag_activity_selection_.bind(1,struct_id);
	std::string tag;
	core::Size activity;
	cppdb::result tag_activity_result(basic::database::safely_read_from_database(tag_activity_selection_));
	if ( tag_activity_result.next() ) {
		tag_activity_result >> tag >> activity;
	} else {
		utility_exit_with_message("Unable to read tag and activity from database.");
	}

	score_selection_.bind(2,struct_id);
	std::map<std::string,core::Real> score_map;

	for ( auto & grid_indice : grid_indices_ ) {
		score_selection_.bind(1,grid_indice.first+"_score_X");
		cppdb::result score_result(basic::database::safely_read_from_database(score_selection_));
		if ( score_result.next() ) {
			core::Real component_score;
			score_result >> component_score;
			score_map.insert(std::make_pair(grid_indice.first,component_score));

		}
	}


	qsarOptData new_point;
	new_point.activity = static_cast<bool>(activity);
	new_point.tag = tag;
	new_point.score_map = score_map;

	return new_point;
}

}
}
