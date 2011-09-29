// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/features/ProteinSilentReport_util.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_features_ProteinSilentReprot_util_HH_
#define INCLUDED_protocols_features_ProteinSilentReprot_util_HH_

#include <core/types.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <string>


namespace protocols {
namespace features {

core::Size get_current_structure_count_by_input_tag(
	utility::sql_database::sessionOP db_session,
	std::string const & input_tag);

core::Size get_score_type_id_from_score_term(
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id,
	std::string const & score_term
);

core::Size get_struct_id_with_lowest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	std::string const & input_tag );

core::Size get_struct_id_with_lowest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	std::string const & input_tag );

core::Size get_struct_id_with_highest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	std::string const & input_tag );

core::Size get_struct_id_with_highest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	std::string const & input_tag );

core::Size get_struct_id_with_nth_lowest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	std::string const & input_tag,
	core::Size const & cutoff_index);

core::Size get_struct_id_with_nth_lowest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	std::string const & input_tag,
	core::Size const & cutoff_index);

core::Real get_score_for_struct_id_and_score_term_from_job_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & struct_id,
	std::string const & score_term);

core::Real get_score_for_struct_id_and_score_term_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & struct_id,
	core::Size const & score_type_id);

}
}


#endif /* PROTEINSILENTREPORT_UTIL_HH_ */
