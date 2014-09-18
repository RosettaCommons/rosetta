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
/// @brief A helper application that converts data in the format expected by the GPU RmsdCalculator

/// @author Tim Jacobs

//Basic
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>

#include <basic/Tracer.hh>

//Utility
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>


//Core
#include <core/init/init.hh>
#include <core/types.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>

// External libraries
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

//C++
#include <cstdio>

static thread_local basic::Tracer TR( "NodePairRmsdCalculator" );

void
generate_schema(
	utility::sql_database::sessionOP db_session
){
	using namespace basic::database::schema_generator;
	using utility::vector1;

	Column struct_id("struct_id", new DbUUID(), false);
	
	//******substructures******//
	Column substructure_id_pkey("substructure_id", new DbInteger(), false, true);
	
//	utility::vector1<Column> substructure_pkey_cols;
//	substructure_pkey_cols.push_back(struct_id);
//	substructure_pkey_cols.push_back(substructure_id);
	
	utility::vector1<Column> substructure_fkey_1_cols;
	substructure_fkey_1_cols.push_back(struct_id);
	
	vector1< std::string > substructure_reference_1_cols;
	substructure_reference_1_cols.push_back("struct_id");
	
	ForeignKey substructure_fkey_1(substructure_fkey_1_cols, "structures", substructure_reference_1_cols, false);
	
	Schema substructures("substructures", PrimaryKey(substructure_id_pkey));
	substructures.add_column(struct_id);
	substructures.add_foreign_key(substructure_fkey_1);
	
	substructures.write(db_session);
	
	//******nodes******//
	Column node_id_pkey("node_id", new DbInteger(), false, true);
	Column substructure_id_fkey("substructure_id", new DbInteger(), false, false);
	
//	utility::vector1<Column> node_pkey_cols;
//	node_pkey_cols.push_back(struct_id);
//	node_pkey_cols.push_back(node_id);
	
	utility::vector1<Column> node_fkey_1_cols;
	node_fkey_1_cols.push_back(substructure_id_fkey);
//	node_fkey_1_cols.push_back(struct_id);
	
	vector1< std::string > node_reference_1_cols;
	node_reference_1_cols.push_back("substructure_id");
//	node_reference_1_cols.push_back("struct_id");
	
	ForeignKey node_fkey_1(node_fkey_1_cols, "substructures", node_reference_1_cols, false);
	
	Schema substructure_nodes("substructure_nodes", PrimaryKey(node_id_pkey));
	substructure_nodes.add_column(struct_id);
	substructure_nodes.add_foreign_key(node_fkey_1);
	
	substructure_nodes.write(db_session);
	
	//*****node_elements*****//
	Column element_id("element_id", new DbInteger(), false, false/*don't auto increment, it's really a fk*/);
	Column index_in_node("index_in_node", new DbInteger(), false);
	Column node_id_fkey("node_id", new DbInteger(), false, false);
	
	utility::vector1<Column> node_element_pkey_cols;
	node_element_pkey_cols.push_back(struct_id);
	node_element_pkey_cols.push_back(element_id);
	node_element_pkey_cols.push_back(node_id_fkey);
	
	utility::vector1<Column> node_element_fkey_1_cols;
	node_element_fkey_1_cols.push_back(node_id_fkey);
//	node_element_fkey_1_cols.push_back(struct_id);
	
	vector1< std::string > node_element_reference_1_cols;
	node_element_reference_1_cols.push_back("node_id");
//	node_element_reference_1_cols.push_back("struct_id");
	
	ForeignKey node_element_fkey_1(node_element_fkey_1_cols, "substructure_nodes", node_element_reference_1_cols, false);
	
	Schema node_elements("node_elements", PrimaryKey(node_element_pkey_cols));
	node_elements.add_column(struct_id);
	node_elements.add_foreign_key(node_element_fkey_1);
	node_elements.add_column(index_in_node);
	
	node_elements.write(db_session);
}
	
void
prepare_smotifs(
	utility::sql_database::sessionOP db_session
){
	using cppdb::statement;
	using cppdb::result;
	
	///****Insert statements****///
	std::string substructure_insert =
		"INSERT INTO substructures (struct_id) VALUES (?)";
	statement substructure_insert_stmt =
		basic::database::safely_prepare_statement(substructure_insert, db_session);
		
	std::string node_insert =
		"INSERT INTO substructure_nodes (struct_id, substructure_id) VALUES (?,?)";
	statement node_insert_stmt =
		basic::database::safely_prepare_statement(node_insert, db_session);
		
	std::string node_element_insert =
		"INSERT INTO node_elements (element_id, struct_id, node_id, index_in_node) VALUES (?,?,?,?)";
	statement node_element_insert_stmt =
		basic::database::safely_prepare_statement(node_element_insert, db_session);

	///****Select statements****///
	std::string select_smotifs =
		"SELECT struct_id, smotif_id, secondary_struct_segment_id_1, secondary_struct_segment_id_2\n"
		"FROM smotifs\n"
		"ORDER BY struct_id, smotif_id";
	
	statement select_smotifs_stmt =
		basic::database::safely_prepare_statement(select_smotifs, db_session);
		
	///****Go****///
	result smotifs_res =
		basic::database::safely_read_from_database( select_smotifs_stmt );
	
//	if ( ! smotifs_res.next() ) {
//		utility_exit_with_message( "Failed to retrieve smotifs" );
//	}
	
	boost::uuids::uuid struct_id;
	core::Size smotif_id;
	core::Size secondary_struct_segment_id_1, secondary_struct_segment_id_2;
	db_session->begin();
	while(smotifs_res.next())
	{
		smotifs_res >> struct_id >> smotif_id >> secondary_struct_segment_id_1 >> secondary_struct_segment_id_2;
		
		std::string struct_id_string = boost::uuids::to_string(struct_id);
		
		//each smotif corresponds to one substructure
		substructure_insert_stmt.bind(1, struct_id);
		basic::database::safely_write_to_database( substructure_insert_stmt );

		int substructure_id = substructure_insert_stmt.sequence_last("substructures_substructure_id_seq");
		
		//add the first node
		node_insert_stmt.bind(1, struct_id);
		node_insert_stmt.bind(2, substructure_id);
		basic::database::safely_write_to_database( node_insert_stmt );
		
		int node_id_1 = node_insert_stmt.sequence_last("substructure_nodes_node_id_seq");

		//two comparison elements for each smotifs (the secondary structure segments)
		node_element_insert_stmt.bind(1, secondary_struct_segment_id_1);
		node_element_insert_stmt.bind(2, struct_id);
		node_element_insert_stmt.bind(3, node_id_1);
		node_element_insert_stmt.bind(4, 1);
		basic::database::safely_write_to_database( node_element_insert_stmt );
		
		node_element_insert_stmt.bind(1, secondary_struct_segment_id_2);
		node_element_insert_stmt.bind(2, struct_id);
		node_element_insert_stmt.bind(3, node_id_1);
		node_element_insert_stmt.bind(4, 2);
		basic::database::safely_write_to_database( node_element_insert_stmt );
		
		//increment node id and add second node
		node_insert_stmt.bind(1, struct_id);
		node_insert_stmt.bind(2, substructure_id);
		basic::database::safely_write_to_database( node_insert_stmt );
		
		int node_id_2 = node_insert_stmt.sequence_last("substructure_nodes_node_id_seq");
		
		//two comparison elements for each smotifs (the secondary structure segments)
		node_element_insert_stmt.bind(1, secondary_struct_segment_id_2);
		node_element_insert_stmt.bind(2, struct_id);
		node_element_insert_stmt.bind(3, node_id_2);
		node_element_insert_stmt.bind(4, 1);
		basic::database::safely_write_to_database( node_element_insert_stmt );
		
		node_element_insert_stmt.bind(1, secondary_struct_segment_id_1);
		node_element_insert_stmt.bind(2, struct_id);
		node_element_insert_stmt.bind(3, node_id_2);
		node_element_insert_stmt.bind(4, 2);
		basic::database::safely_write_to_database( node_element_insert_stmt );
	}
	db_session->commit();
}

void
prepare_helix_bundles(
	utility::sql_database::sessionOP db_session
){
	using cppdb::statement;
	using cppdb::result;

	///****Insert statements****///
	std::string substructure_insert =
		"INSERT INTO substructures (struct_id) VALUES (?)";
	statement substructure_insert_stmt =
		basic::database::safely_prepare_statement(substructure_insert, db_session);
		
	std::string node_insert =
		"INSERT INTO substructure_nodes (struct_id, substructure_id) VALUES (?,?)";
	statement node_insert_stmt =
		basic::database::safely_prepare_statement(node_insert, db_session);
		
	std::string node_element_insert =
		"INSERT INTO node_elements (element_id, struct_id, node_id, index_in_node) VALUES (?,?,?,?)";
	statement node_element_insert_stmt =
		basic::database::safely_prepare_statement(node_element_insert, db_session);

	///****Select statements****///
	std::string select_helices =
		"SELECT struct_id, bundle_id, helix_id\n"
		"FROM bundle_helices\n"
		"ORDER BY struct_id, bundle_id";
	
	statement select_helices_stmt =
		basic::database::safely_prepare_statement(select_helices, db_session);
		
	///****Go****///
	result helices_res =
		basic::database::safely_read_from_database( select_helices_stmt );
	
	boost::uuids::uuid struct_id;
	core::Size bundle_id, helix_id;
	core::Size prev_bundle_id=0;
	utility::vector1<core::Size> bundle_helices;
	while(helices_res.next())
	{
		helices_res >> struct_id >> bundle_id >> helix_id;
		
		std::string struct_id_string = boost::uuids::to_string(struct_id);
		TR << "Working with struct_id: " << struct_id_string << std::endl;
		
		//new bundle
		if(bundle_id!=prev_bundle_id)
		{
			substructure_insert_stmt.bind(1, struct_id);
			basic::database::safely_write_to_database(substructure_insert_stmt);
			
			//write to substructures
			//write to nodes
			
			bundle_helices.clear();
		}
		bundle_helices.push_back(helix_id);
		prev_bundle_id=bundle_id;
	}
}

prepare_beta_sandwiches(
	utility::sql_database::sessionOP db_session
){
}
	
int
main( int argc, char * argv [] )
{
	using utility::vector1;
	using namespace ObjexxFCL;
	
	// initialize core
	core::init::init(argc, argv);
	
	// Initialize DB
	utility::sql_database::sessionOP db_session(
		basic::database::get_db_session());
		
	generate_schema(db_session);
	prepare_smotifs(db_session);
//	prepare_helix_bundles(db_session);
		
	//db_session->commit();

	return 0;
}
