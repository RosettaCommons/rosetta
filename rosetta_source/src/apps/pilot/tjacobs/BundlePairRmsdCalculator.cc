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

#include <protocols/wum/DatabaseQueryWorkUnitManager.hh>
#include <protocols/wum/DatabaseEntryWorkUnit.hh>

#include <devel/helixAssembly/BundlePairRmsdWorkUnit.hh>


#include <core/pose/util.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitManager.hh>
#include <protocols/wum/MPI_WorkUnitManager_Slave.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <cstdio>

#include <core/init.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//TEST
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

int
main( int argc, char * argv [] )
{
    using cppdb::statement;
    using cppdb::result;
    
	// initialize core
	core::init(argc, argv);
    
    using namespace basic::options;
    using namespace protocols::wum;
        
    // Initialize DB - hardcoded currently
    utility::sql_database::sessionOP db_session(basic::database::get_db_session(option[OptionKeys::inout::database_filename], option[OptionKeys::inout::database_mode], false, false));
    
    
    WorkUnitList wulist;
    
    // Create a workunit for waiting - why is this necessary?
    protocols::wum::WorkUnit_WaitOP wait_wu = new protocols::wum::WorkUnit_Wait( 2 );
    wulist.register_work_unit( "waitwu", wait_wu );
    
    //Create a workunit for processing db results
    std::string db_wu_name = "bundle_pair_rmsd_wu";
    DatabaseEntryWorkUnitOP bundle_pair_rmsd_wu = new BundlePairRmsdWorkUnit( db_session );
    wulist.register_work_unit(db_wu_name, bundle_pair_rmsd_wu );
    
    //Create the query string whose results will be distributed across nodes
    std::string select_string = 
    "SELECT * FROM helix_bundles";
    
    std::string bundle_pairs_creation_string =
    
    "CREATE TABLE IF NOT EXISTS bundle_pairs (\n"
    "	pair_id INTEGER,\n"
    "	bundle_id INTEGER,\n"
    "	helix_1 INTEGER,\n"
    "	helix_2 INTEGER,\n"
    "	FOREIGN KEY (bundle_id)\n"
    "		REFERENCES helix_bundles (bundle_id)\n"
    "		DEFERRABLE INITIALLY DEFERRED,\n"
    "	CONSTRAINT dist_is_nonnegative CHECK (count >= 0),\n"
    "	PRIMARY KEY (struct_id, resType1, atmType1, resType2, atmType2, distBin));\n"
   
    "CREATE TABLE IF NOT EXISTS bundle_pairs AS \n" 
    "SELECT\n"
    "   bundles.struct_id,\n"
    "	bundles.bundle_id,\n"
    "   helices1.helix_id,\n"
    "	helices1.residue_begin AS helix_1_begin,\n"
    "	helices1.residue_end AS helix_1_end,\n"
    "   helices2.helix_id,\n"
    "	helices2.residue_begin AS helix_2_begin,\n"
    "	helices2.residue_end AS helix_2_end\n"
    "FROM\n"
    "	helix_bundles AS bundles\n"
    "JOIN bundle_helices AS helices1 ON\n"
    "   bundles.bundle_id = helices1.bundle_id\n"
    "JOIN bundle_helices AS helices2 ON\n"
    "   bundles.bundle_id = helices2.bundle_id\n"
    "WHERE\n"
    "   helices1.bundle_id = helices2.bundle_id AND\n"
    "   helices1.helix_id <> helices2.helix_id;";

    
    // Now define the structure of the algorithm, i.e. assign roles to individual nodes
    WorkUnitManagerOP wu_manager;
    core::Size master_id = 0; /* only 1 master for this protocol */
    if ( mpi_rank() == 0 ){
        wu_manager = new DatabaseQueryWorkUnitManager<BundlePairRmsdWorkUnit>(master_id, db_session, select_string, db_wu_name);
//        wu_manager = new DatabaseQueryWorkUnitManager(master_id, db_session, select_string, "dbwu");
    }
    else
    {
        wu_manager = new MPI_WorkUnitManager_Slave(master_id);
    }
    
    // make sure all the necessary work unit have been properly initiated.
    wu_manager->register_work_units( wulist );
    
    // now launch this node
    wu_manager->go();
    
#ifdef USEMPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
#endif
	return 0;
}




