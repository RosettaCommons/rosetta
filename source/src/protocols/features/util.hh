// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_util_hh
#define INCLUDED_protocols_features_util_hh

#include <utility/sql_database/DatabaseSessionManager.hh>

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>

#include <string>
#include <map>

#include <protocols/features/FeaturesReporter.fwd.hh>

namespace protocols{
namespace features{

/// @brief write the given protocol and batch ids to the database. The protocol and batches
///features reporters will check for an existing entry with the same key, and write if one
///does not exist. Not recommended for parallel use as it is subject to race conditions (due
///to the nature of 'insert or ignore' type database writing)
void
set_protocol_and_batch_id(
	core::Size protocol_id,
	core::Size batch_id,
	std::string const & batch_name,
	std::string const & batch_description,
	utility::vector1<FeaturesReporterOP> features_reporters,
	utility::sql_database::sessionOP db_session
);

/// @brief Get the protocol and batch ids or create them if they don't
///yet exist. For MPI protocols, only allow the head node to create
///protocol or batch ids and have the other nodes ask the head node
///for the info.
std::pair<core::Size, core::Size>
get_protocol_and_batch_id(
	std::string const & batch_name,
	std::string const & batch_description,
	utility::vector1<FeaturesReporterOP> features_reporters,
	utility::sql_database::sessionOP db_session);

void
write_features_reporters_table(
	utility::vector1<FeaturesReporterOP> features_reporters,
	utility::sql_database::sessionOP db_session
);

/// @brief write the linking table between features reporters
///and batches. This happens here so that the protocol/batch
///id framework can be used to prevent duplicate key entries.
///This function gets called when the batch id is written.
void
write_batch_reports_table(
	utility::vector1<FeaturesReporterOP> features_reporters,
	core::Size batch_id,
	utility::sql_database::sessionOP db_session
);

std::pair<core::Size, core::Size> deserialize_db_listener_data(std::string data);

std::string serialize_ids(int protocol_id, std::string batch_name, core::Size batch_id);

core::Size
get_batch_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session
);

utility::vector1<StructureID>
struct_ids_from_tag(
	utility::sql_database::sessionOP db_session,
	std::string const & tag);

std::string serialize_residue_xyz_coords(core::conformation::Residue const & residue);

utility::vector1< numeric::xyzVector<core::Real> > deserialize_xyz_coords(std::string const & data, core::Size natoms);

/// @brief Returns (?,?,?) With question marks of length n to help create database query.
std::string
get_question_mark_string(core::Size const n);

} //features
} //protocols

#endif
