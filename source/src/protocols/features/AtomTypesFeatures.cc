// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/AtomTypesFeatures.cc
/// @brief  report AtomType properties to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/AtomTypesFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomTypeDatabaseIO.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/database/sql_utils.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>

// C++ Headers
#include <sstream>

namespace protocols {
namespace features {

using std::string;
using std::stringstream;
using std::endl;
using std::pair;
using std::set;
using core::Size;
using core::Real;
using core::pose::Pose;
using utility::sql_database::sessionOP;
using core::chemical::AtomType;
using core::chemical::AtomTypeSet;
using utility::vector1;
using basic::Tracer;

static Tracer TR("protocols.features.AtomTypesFeatures");

AtomTypesFeatures::AtomTypesFeatures() {}

AtomTypesFeatures::~AtomTypesFeatures() {}

string
AtomTypesFeatures::type_name() const { return "AtomTypesFeatures"; }

void
AtomTypesFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	// NOTE: To support building feature databases in parallel, the
	// ResidueTypeSet and ResidueType objects must be identified by
	// their names rather then assigning them a unique id.
	return atom_type_dbio_.write_schema_to_db(db_session);
}

vector1<string>
AtomTypesFeatures::features_reporter_dependencies() const {
	vector1<string> dependencies;
	return dependencies;
}

Size
AtomTypesFeatures::report_features(
	Pose const & pose,
	vector1< bool > const &,
	StructureID const,
	sessionOP db_session
){
	if ( pose.total_residue() == 0 ) return 0;

	// All residues in a pose must have the same atom type set. So look
	// up the atom type set of the first residue if it exists
	AtomTypeSet const & atom_type_set(pose.residue(1).atom_type_set());

	atom_type_dbio_.write_atom_type_set_to_database(atom_type_set, db_session);
	return 0;
}


} // namesapce
} // namespace
