// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueTypesFeatures.cc
/// @brief  report ResidueTypes to features Statistics Scientific Benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/ResidueTypesFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueDatabaseIO.hh>
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
#define foreach BOOST_FOREACH

// C++ Headers
#include <set>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using std::pair;
using std::set;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::BondName;
using core::chemical::AtomIndices;
using core::chemical::ResidueType;
using core::chemical::ResidueTypeCOP;
using utility::sql_database::sessionOP;
using utility::vector1;
using basic::Tracer;

static Tracer TR("protocols.features.ResidueTypesFeatures");

ResidueTypesFeatures::ResidueTypesFeatures()
{
	version_ = residue_dbio_.get_version();
}

ResidueTypesFeatures::~ResidueTypesFeatures() {}

string
ResidueTypesFeatures::type_name() const { return "ResidueTypesFeatures"; }

void
ResidueTypesFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	// NOTE: To support building feature databases in parallel, the
	// ResidueTypeSet and ResidueType objects must be identified by
	// their names rather then assigning them a unique id.
	return residue_dbio_.write_schema_to_db(db_session);
}

utility::vector1<std::string>
ResidueTypesFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	return dependencies;
}


Size
ResidueTypesFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const,
	sessionOP db_session
){

	// Get a set of the unique residue types that are used in this pose
	set< ResidueTypeCOP > res_types;
	for(Size i=1; i <= pose.n_residue(); ++i){
		if(!check_relevant_residues(relevant_residues, i)) continue;
		res_types.insert(&pose.residue_type(i));
	}

	foreach( ResidueTypeCOP res_type, res_types){
		string const & residue_type_set_name(res_type->residue_type_set().name());

		// Is this residue type already in the database?

		residue_dbio_.write_residuetype_to_database(residue_type_set_name,*res_type,db_session);
	}
	return 0;
}


} // namesapce
} // namespace
