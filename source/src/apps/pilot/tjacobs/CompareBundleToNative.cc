// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DatabaseEntryWorkUnit.cc
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

#include <protocols/sewing/BundlePairRmsdWorkUnit.hh>

#include <basic/gpu/GPU.hh>
#include <basic/gpu/Timer.hh>

#include <core/pose/util.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitManager.hh>
#include <protocols/wum/MPI_WorkUnitManager_Slave.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

#include <cstdio>

#include <devel/init.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>

// ObjexxFCL libraries
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//TEST
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

//Protocol headers
#include <protocols/features/ProteinSilentReport.hh>

static basic::Tracer TR("BundleDumper");

static int BLOCK_SIZE = 512;

namespace BundleDumper {
	basic::options::IntegerOptionKey const num_bundles( "num_bundles" ); // fraction to run
}

typedef utility::vector1< std::pair<protocols::features::StructureID, core::Size> > BundleKeys;

BundleKeys
get_random_bundle_ids(
	utility::sql_database::sessionOP db_session,
	core::Size n
){
	std::string struct_id_select_string =
		"SELECT struct_id, bundle_id\n"
		"FROM helix_bundles\n"
		"ORDER BY random()\n"
		"LIMIT ?\n";

	cppdb::statement struct_id_stmt=basic::database::safely_prepare_statement(struct_id_select_string, db_session);
	struct_id_stmt.bind(1, n);

	cppdb::result struct_id_res=basic::database::safely_read_from_database(struct_id_stmt);

	protocols::features::StructureID struct_id;
	core::Size bundle_id;
	BundleKeys ids;
	while(struct_id_res.next())
	{
		struct_id_res >> struct_id >> bundle_id;
		ids.push_back(std::make_pair(struct_id, bundle_id));
	}

	return ids;
}

void
dump_bundle_pose(
	utility::sql_database::sessionOP db_session,
	protocols::features::StructureID struct_id,
	core::Size bundle_id
){
	std::string select_resnums =
		"SELECT residue_begin, residue_end\n"
		"FROM bundle_helices r\n"
		"WHERE\n"
		"	struct_id = ? AND"
		"	bundle_id = ?;";

	cppdb::statement select_resnums_stmt =
		basic::database::safely_prepare_statement(select_resnums, db_session);
	
	protocols::features::ProteinSilentReportOP protein_silent_report = new protocols::features::ProteinSilentReport();
	
	select_resnums_stmt.bind(1, struct_id);
	select_resnums_stmt.bind(2, bundle_id);
	cppdb::result res = basic::database::safely_read_from_database(select_resnums_stmt);
	std::set<core::Size> resnums;
	while(res.next()){
		core::Size residue_begin, residue_end;
		res >> residue_begin >> residue_end;
		for(core::Size i=std::min(residue_begin, residue_end); i<=std::max(residue_begin, residue_end); ++i){
			resnums.insert(i);
		}
	}

	core::pose::Pose pose;
	protein_silent_report->load_pose(db_session, struct_id, resnums, pose);

	std::stringstream pdb_name;
	pdb_name << struct_id << "_" << bundle_id << ".pdb";
	pose.dump_pdb(pdb_name.str());
}

int
main( int argc, char * argv [] )
{
	using cppdb::statement;
	using cppdb::result;
	using namespace basic::database::schema_generator;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::wum;

	using namespace core;

	using ObjexxFCL::FArray2D;
	using ObjexxFCL::FArray1D;

	option.add( BundleDumper::num_bundles, "The number of bundles to output");

	// initialize core
	devel::init(argc, argv);

	core::Size num_bundles = option[BundleDumper::num_bundles].def(10);

	// Initialize DB
	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );
	
	BundleKeys keys = get_random_bundle_ids(db_session, num_bundles);
	for(core::Size i=1; i<=keys.size(); ++i){
		dump_bundle_pose(db_session, keys[i].first, keys[i].second);
	}
}
