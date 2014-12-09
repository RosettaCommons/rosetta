// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ConcurrencyTest.hh
///
/// @brief
/// @author Tim Jacobs

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <protocols/features/helixAssembly/ConcurrencyTest.hh>

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

//Devel
#include <protocols/features/helixAssembly/HelixBundleFeatures.hh>
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility and basic
#include <numeric/random/random.hh>
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//C++
#include <string>
#include <math.h>

//External Headers
#include <cppdb/frontend.h>

//Basic
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
//#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>

namespace protocols {
namespace features {
namespace helixAssembly {

	void
	ConcurrencyTest::write_schema_to_db(utility::sql_database::sessionOP db_session) const{

		using namespace basic::database::schema_generator;

		PrimaryKey id(Column("id", DbDataTypeOP( new DbBigInt() ), false));
		Column random_number(Column("description", DbDataTypeOP( new DbInteger() )));

		Schema concurrency_test("concurrency_test", id);
		concurrency_test.add_column(random_number);

		concurrency_test.write(db_session);

	}

	///@brief collect all the feature data for the pose
	core::Size
	ConcurrencyTest::report_features(
		core::pose::Pose const &,
		utility::vector1<bool> const &,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	){


		std::string test_insert =  "INSERT INTO concurrency_test (id, random_num) VALUES (?);";
		for(int i=1; i<=100000; i++){
			cppdb::statement test_stmt(basic::database::safely_prepare_statement(test_insert,db_session));
			test_stmt.bind(1,struct_id);
			test_stmt.bind(2,numeric::random::random_range(0,INT_MAX));
			basic::database::safely_write_to_database(test_stmt);
		}
		return 0;
	}

}
}
}
