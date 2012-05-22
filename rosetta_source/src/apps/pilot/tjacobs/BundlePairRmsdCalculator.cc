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

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

#include <cstdio>

#include <core/init.hh>

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

static basic::Tracer TRACER("BundlePairRmsdCalculator");

namespace BundlePairRmsdCalculator {
	basic::options::IntegerOptionKey const total_fractions( "total_fractions" ); // fraction to run
	basic::options::IntegerOptionKey const fraction( "fraction" ); // fraction to run
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
	
	option.add( BundlePairRmsdCalculator::total_fractions, "The number of fractions to split the selection into");
	option.add( BundlePairRmsdCalculator::fraction, "The fraction of results to run RMSD comparisons on");
	
	// initialize core
	core::init(argc, argv);
	
	if(!option[OptionKeys::inout::database_mode].user()){
		utility_exit_with_message("ERROR: Must provide database mode using inout::database_mode. Allowed modes are sqlite3, mysql, and postgres");
	}
	std::string db_mode=option[OptionKeys::inout::database_mode].value();
	std::cout << "database mode: " << db_mode << std::endl;
	
	
	if(!option[OptionKeys::inout::database_filename].user()){
		utility_exit_with_message("ERROR: Must provide database file name using inout::database_filename");
	}
	std::string db_name=option[OptionKeys::inout::database_filename].value();
	
	std::cout << "database name " << db_name << std::endl;
	
	// Initialize DB
	utility::sql_database::sessionOP db_session(basic::database::get_db_session(db_name, db_mode, false, false));
	
//	WorkUnitList wulist;
//	
//	// Create a workunit for waiting - why is this necessary?
//	protocols::wum::WorkUnit_WaitOP wait_wu = new protocols::wum::WorkUnit_Wait( 2 );
//	wulist.register_work_unit( "waitwu", wait_wu );
//	
//	//Create a workunit for processing db results
//	std::string db_wu_name = "bundle_pair_rmsd_wu";
//	DatabaseEntryWorkUnitOP bundle_pair_rmsd_wu = new BundlePairRmsdWorkUnit( db_session );
//	wulist.register_work_unit(db_wu_name, bundle_pair_rmsd_wu );

//	std::string bundle_pairs_generation =
//	"CREATE TABLE bundle_pairs AS \n" 
//	"SELECT\n"
//	"   bundles.struct_id,\n"
//	"	bundles.bundle_id,\n"
//	"   helices1.helix_id AS helix_id_1,\n"
//	"	helices1.residue_begin AS helix_1_begin,\n"
//	"	helices1.residue_end AS helix_1_end,\n"
//	"   helices2.helix_id AS helix_id_2,\n"
//	"	helices2.residue_begin AS helix_2_begin,\n"
//	"	helices2.residue_end AS helix_2_end\n"
//	"FROM\n"
//	"	helix_bundles AS bundles\n"
//	"JOIN bundle_helices AS helices1 ON\n"
//	"   bundles.bundle_id = helices1.bundle_id\n"
//	"JOIN bundle_helices AS helices2 ON\n"
//	"   bundles.bundle_id = helices2.bundle_id\n"
//	"WHERE\n"
//	"   helices1.bundle_id = helices2.bundle_id AND\n"
//	"   helices1.helix_id < helices2.helix_id AND\n"
//	"   helices1.helix_id <> helices2.helix_id;";
//	
//	TRACER << "Bundle selection string " << bundle_pairs_generation << std::endl;
//	
//	cppdb::statement bundle_pairs_statement(basic::database::safely_prepare_statement(bundle_pairs_generation,db_session));
//	basic::database::safely_write_to_database(bundle_pairs_statement);
	
	//database query to distribute
	
	int fraction = option[BundlePairRmsdCalculator::fraction].def(1);
	int total_fractions = option[BundlePairRmsdCalculator::total_fractions].def(1000);
	
	std::string select_string = 
"SELECT * FROM(\n"
	"SELECT\n"
	"		pair1.struct_id AS struct_id_1,\n"
	"		pair1.pair_id AS pair_id_1,\n"
	"		pair1.helix_1_flipped AS pair_1_helix_1_flipped,\n"
	"		pair1.helix_1_begin AS pair_1_helix_1_begin,\n"
	"		pair1.helix_1_end AS pair_1_helix_1_end,\n"
	"		pair1.helix_2_flipped AS pair_1_helix_2_flipped,\n"
	"		pair1.helix_2_begin AS pair_1_helix_2_begin,\n"
	"		pair1.helix_2_end AS pair_1_helix_2_end,\n"
	"		pair1.helix_3_flipped AS pair_1_helix_3_flipped,\n"
	"		pair1.helix_3_begin AS pair_1_helix_3_begin,\n"
	"		pair1.helix_3_end AS pair_1_helix_3_end,\n"
	"\n"
	"		pair2.struct_id AS struct_id_2,\n"
	"		pair2.pair_id AS pair_id_2,\n"
	"		pair2.helix_1_flipped AS pair_2_helix_1_flipped,\n"
	"		pair2.helix_1_begin AS pair_2_helix_1_begin,\n"
	"		pair2.helix_1_end AS pair_2_helix_1_end,\n"
	"		pair2.helix_2_flipped AS pair_2_helix_2_flipped,\n"
	"		pair2.helix_2_begin AS pair_2_helix_2_begin,\n"
	"		pair2.helix_2_end AS pair_2_helix_2_end,\n"
	"		pair2.helix_3_flipped AS pair_2_helix_3_flipped,\n"
	"		pair2.helix_3_begin AS pair_2_helix_3_begin,\n"
	"		pair2.helix_3_end AS pair_2_helix_3_end\n"
", ntile(" + utility::to_string(total_fractions) + ") over (order by pair1.struct_id) as fraction\n"
	"FROM bundle_pairs pair1\n"
	"JOIN bundle_pairs pair2 ON\n"
	"		pair1.struct_id <> pair2.struct_id AND\n"
	"		pair1.pair_id < pair2.pair_id\n"
") as temp WHERE fraction = " + utility::to_string(fraction);
	;
	
	//begin a transaction for all this nonsense
	db_session->begin();
	
	cppdb::statement select_stmt=basic::database::safely_prepare_statement(select_string, db_session);
	cppdb::result res=basic::database::safely_read_from_database(select_stmt);
	
	boost::uuids::uuid struct_id_1;
	boost::uuids::uuid struct_id_2;
	
	core::Size pair_id_1;
	
	core::Size pair_1_helix_1_flipped;
	core::Size pair_1_helix_1_begin;
	core::Size pair_1_helix_1_end;
	
	core::Size pair_1_helix_2_flipped;
	core::Size pair_1_helix_2_begin;
	core::Size pair_1_helix_2_end;
	
	core::Size pair_1_helix_3_flipped;
	core::Size pair_1_helix_3_begin;
	core::Size pair_1_helix_3_end;
	
	core::Size pair_id_2;
	
	core::Size pair_2_helix_1_flipped;
	core::Size pair_2_helix_1_begin; 
	core::Size pair_2_helix_1_end; 
	
	core::Size pair_2_helix_2_flipped;
	core::Size pair_2_helix_2_begin; 
	core::Size pair_2_helix_2_end;
	
	core::Size pair_2_helix_3_flipped;
	core::Size pair_2_helix_3_begin; 
	core::Size pair_2_helix_3_end;
	
	std::string helix_coords_select =
	"SELECT\n"
	"	res.x,\n"
	"	res.y,\n"
	"	res.z\n"
	"FROM residue_atom_coords res\n"
	"WHERE\n"
	"	res.struct_id = ? AND\n"
	"	res.atomno = 2 AND\n" //atomno 2 is the CA. It is bad practice to use this directly, but, I do what I want
	"	res.seqpos BETWEEN ? AND ?\n"
	"ORDER BY seqpos ASC;";
	
	statement select_coords_stmt(basic::database::safely_prepare_statement(helix_coords_select,db_session));
	
	std::string comparison_insert = "INSERT INTO pair_comparisons(pair_id_1, pair_id_2,rmsd,third_helix_rmsd) VALUES(?,?,?,?);";
	cppdb::statement comparison_insert_stmt(basic::database::safely_prepare_statement(comparison_insert,db_session));
	while(res.next()){
	
		res >> struct_id_1 >> pair_id_1 >> pair_1_helix_1_flipped >> pair_1_helix_1_begin >> pair_1_helix_1_end 
			>> pair_1_helix_2_flipped >> pair_1_helix_2_begin >> pair_1_helix_2_end
			>> pair_1_helix_3_flipped >> pair_1_helix_3_begin >> pair_1_helix_3_end
			>> struct_id_2 >> pair_id_2 >> pair_2_helix_1_flipped >> pair_2_helix_1_begin >> pair_2_helix_1_end 
			>> pair_2_helix_2_flipped >> pair_2_helix_2_begin >> pair_2_helix_2_end
			>> pair_2_helix_3_flipped >> pair_2_helix_3_begin >> pair_2_helix_3_end;
		
		core::Size helix_size = pair_1_helix_1_end-pair_1_helix_1_begin;//all helices are the same size
				
		select_coords_stmt.bind(1,struct_id_1);
		select_coords_stmt.bind(2, pair_1_helix_1_begin);
		select_coords_stmt.bind(3, pair_1_helix_1_end);
		result pair_1_helix_1_res(basic::database::safely_read_from_database(select_coords_stmt));
		
		select_coords_stmt.bind(1,struct_id_2);
		select_coords_stmt.bind(2, pair_2_helix_1_begin);
		select_coords_stmt.bind(3, pair_2_helix_1_end);
		result pair_2_helix_1_res(basic::database::safely_read_from_database(select_coords_stmt));
		
		core::Size counter = 0;
		utility::vector1< numeric::xyzVector<core::Real> > pair_1_helix_1_coords;
		utility::vector1< numeric::xyzVector<core::Real> > pair_2_helix_1_coords;
		while(pair_1_helix_1_res.next()){
			pair_2_helix_1_res.next();
			
			core::Real x1,y1,z1,x2,y2,z2;
			pair_1_helix_1_res >> x1 >> y1 >> z1;
			pair_2_helix_1_res >> x2 >> y2 >> z2;
			
			//		TR << "Coords for pos " << counter << " of pair 1: " << x1 << " " << y1 << " " << z1 << std::endl;
			//		TR << "Coords for pos " << counter << " of pair 2: " << x2 << " " << y2 << " " << z2 << std::endl;
			pair_1_helix_1_coords.push_back(numeric::xyzVector< core::Real >(x1,y1,z1));
			pair_2_helix_1_coords.push_back(numeric::xyzVector< core::Real >(x2,y2,z2));
			++counter;
		}
		if(pair_1_helix_1_flipped){reverse(pair_1_helix_1_coords.begin(), pair_1_helix_1_coords.end());}
		if(pair_2_helix_1_flipped){reverse(pair_2_helix_1_coords.begin(), pair_2_helix_1_coords.end());}
		
		
		select_coords_stmt.bind(1,struct_id_1);
		select_coords_stmt.bind(2, pair_1_helix_2_begin);
		select_coords_stmt.bind(3, pair_1_helix_2_end);
		result pair_1_helix_2_res(basic::database::safely_read_from_database(select_coords_stmt));
		
		select_coords_stmt.bind(1,struct_id_2);
		select_coords_stmt.bind(2, pair_2_helix_2_begin);
		select_coords_stmt.bind(3, pair_2_helix_2_end);
		result pair_2_helix_2_res(basic::database::safely_read_from_database(select_coords_stmt));
		
		utility::vector1< numeric::xyzVector<core::Real> > pair_1_helix_2_coords;
		utility::vector1< numeric::xyzVector<core::Real> > pair_2_helix_2_coords;
		while(pair_1_helix_2_res.next()){
			pair_2_helix_2_res.next();
			
			core::Real x1,y1,z1,x2,y2,z2;
			pair_1_helix_2_res >> x1 >> y1 >> z1;
			pair_2_helix_2_res >> x2 >> y2 >> z2;
			
			pair_1_helix_2_coords.push_back(numeric::xyzVector< core::Real >(x1,y1,z1));
			pair_2_helix_2_coords.push_back(numeric::xyzVector< core::Real >(x2,y2,z2));
			++counter;
		}
		if(pair_1_helix_2_flipped){reverse(pair_1_helix_2_coords.begin(), pair_1_helix_2_coords.end());}
		if(pair_2_helix_2_flipped){reverse(pair_2_helix_2_coords.begin(), pair_2_helix_2_coords.end());}
		
		
		utility::vector1< numeric::xyzVector<core::Real> > pair_1_coords;
		pair_1_coords.insert( pair_1_coords.end(), pair_1_helix_1_coords.begin(), pair_1_helix_1_coords.end() );
		pair_1_coords.insert( pair_1_coords.end(), pair_1_helix_2_coords.begin(), pair_1_helix_2_coords.end() );
		
		utility::vector1< numeric::xyzVector<core::Real> > pair_2_coords;
		pair_2_coords.insert( pair_2_coords.end(), pair_2_helix_1_coords.begin(), pair_2_helix_1_coords.end() );
		pair_2_coords.insert( pair_2_coords.end(), pair_2_helix_2_coords.begin(), pair_2_helix_2_coords.end() );
		
//		core::Real bundle_pair_rmsd = numeric::model_quality::calc_rms(pair_1_coords,pair_2_coords);
		
		
		/* NOW COMES THE TRICKEIR PART*/
		
		//Convert to FArrays
		runtime_assert(pair_1_coords.size() == pair_2_coords.size());
		
		//Save coords to FArrays
		FArray2D< numeric::Real > p1_coords( 3, pair_1_coords.size() );
		FArray2D< numeric::Real > p2_coords( 3, pair_2_coords.size() );
		for ( Size i = 1; i <= pair_1_coords.size(); ++i ) {
			for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
				p1_coords(k,i) = pair_1_coords[i](k);
				p2_coords(k,i) = pair_2_coords[i](k);
			}
		}
		
		//uu is rotational matrix that transforms p2coords onto p1coords in a way that minimizes rms. Do this for ONLY the first two helices
		FArray1D< numeric::Real > ww( pair_1_coords.size(), 1.0 );//weight matrix, all 1 for my purposes
		FArray2D< numeric::Real > uu( 3, 3, 0.0 );//transformation matrix
		numeric::Real ctx;
		numeric::model_quality::findUU( p1_coords, p2_coords, ww, pair_1_coords.size(), uu, ctx );
		
		//get coords for third helices of both pairs
		select_coords_stmt.bind(1,struct_id_1);
		select_coords_stmt.bind(2, pair_1_helix_3_begin);
		select_coords_stmt.bind(3, pair_1_helix_3_end);
		result pair_1_helix_3_res(basic::database::safely_read_from_database(select_coords_stmt));
		
		select_coords_stmt.bind(1,struct_id_2);
		select_coords_stmt.bind(2, pair_2_helix_3_begin);
		select_coords_stmt.bind(3, pair_2_helix_3_end);
		result pair_2_helix_3_res(basic::database::safely_read_from_database(select_coords_stmt));
		
		counter = 0;
		utility::vector1< numeric::xyzVector<core::Real> > pair_1_helix_3_coords;
		utility::vector1< numeric::xyzVector<core::Real> > pair_2_helix_3_coords;
		while(pair_1_helix_3_res.next()){
			pair_2_helix_3_res.next();
			
			core::Real x1,y1,z1,x2,y2,z2;
			pair_1_helix_3_res >> x1 >> y1 >> z1;
			pair_2_helix_3_res >> x2 >> y2 >> z2;
			
			pair_1_helix_3_coords.push_back(numeric::xyzVector< core::Real >(x1,y1,z1));
			pair_2_helix_3_coords.push_back(numeric::xyzVector< core::Real >(x2,y2,z2));
			++counter;
		}
		if(pair_1_helix_3_flipped){reverse(pair_1_helix_3_coords.begin(), pair_1_helix_3_coords.end());}
		if(pair_2_helix_3_flipped){reverse(pair_2_helix_3_coords.begin(), pair_2_helix_3_coords.end());}
		
		//Add third helix to pair vectors, these are now essentiall bundle-vectors
		pair_1_coords.insert( pair_1_coords.end(), pair_1_helix_3_coords.begin(), pair_1_helix_3_coords.end() );
		pair_2_coords.insert( pair_2_coords.end(), pair_2_helix_3_coords.begin(), pair_2_helix_3_coords.end() );
		
		//Save coords, now with third helix, to new FArrays
		FArray2D< numeric::Real > b1_coords( 3, pair_1_coords.size() );
		FArray2D< numeric::Real > b2_coords( 3, pair_2_coords.size() );
		for ( Size i = 1; i <= pair_1_coords.size(); ++i ) {
			for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
				b1_coords(k,i) = pair_1_coords[i](k);
				b2_coords(k,i) = pair_2_coords[i](k);
			}
		}		
		
		//move bundle to the origin using the center of mass for only the first two helices
		for ( int k = 1; k <= 3; ++k ) {
			numeric::Real bundle_1_offset = 0.0;
			numeric::Real bundle_2_offset = 0.0;
			
			for ( int j = 1; j <= pair_2_coords.size()-pair_2_helix_3_coords.size(); ++j ) {
				bundle_1_offset += b1_coords(k,j);
				bundle_2_offset += b2_coords(k,j);
			}
			bundle_1_offset /= (pair_2_coords.size()-pair_2_helix_3_coords.size());
			bundle_2_offset /= (pair_2_coords.size()-pair_2_helix_3_coords.size());
			
			for ( int j = 1; j <= pair_1_coords.size(); ++j ) {
				b1_coords(k,j) -= bundle_1_offset;
				b2_coords(k,j) -= bundle_2_offset;
			}
		}
		
		//transform the coords of the entire bundle using rotational matrix produced by findUU for the first two helices
		FArray2D< numeric::Real > b2_coords_transformed(3, pair_2_coords.size());
		for ( int i = 1; i <= pair_2_coords.size(); ++i ) {//iterate through all points
			b2_coords_transformed(1,i) = ( uu(1,1)*b2_coords(1,i) )+( uu(1,2)*b2_coords(2,i) ) +( uu(1,3)*b2_coords(3,i) );
			b2_coords_transformed(2,i) = ( uu(2,1)*b2_coords(1,i) )+( uu(2,2)*b2_coords(2,i) ) +( uu(2,3)*b2_coords(3,i) );
			b2_coords_transformed(3,i) = ( uu(3,1)*b2_coords(1,i) )+( uu(3,2)*b2_coords(2,i) ) +( uu(3,3)*b2_coords(3,i) );
		}
				
		//Calculate RMSD for the first two helices
		numeric::Real tot = 0;
		for ( int i = 1; i <= pair_2_coords.size()-pair_2_helix_3_coords.size(); ++i ) {
			for ( int j = 1; j <= 3; ++j ) {
				tot += std::pow( b1_coords(j,i) - b2_coords_transformed(j,i), 2 );
			}
		}
		core::Real bundle_pair_rmsd = std::sqrt(tot/(pair_2_coords.size()-pair_2_helix_3_coords.size()));
		
		//calc rms between transformed bundle_2_helix_3 pts and bundle_1_helix_3 pts 
		tot = 0;
		for ( int i = pair_1_helix_1_coords.size()+pair_1_helix_2_coords.size()+1; i <= pair_2_coords.size(); ++i ) {
			for ( int j = 1; j <= 3; ++j ) {
				tot += std::pow( b1_coords(j,i) - b2_coords_transformed(j,i), 2 );
			}
		}
		core::Real third_helix_rmsd = std::sqrt(tot/pair_2_helix_3_coords.size());

//		TRACER << "ALL THE RMSDS: " << bundle_pair_rmsd << " " << test_rmsd << " " << third_helix_rmsd << std::endl;

		comparison_insert_stmt.bind(1, pair_id_1);
		comparison_insert_stmt.bind(2, pair_id_2);
		comparison_insert_stmt.bind(3, bundle_pair_rmsd);
		comparison_insert_stmt.bind(4, third_helix_rmsd);
		basic::database::safely_write_to_database(comparison_insert_stmt);	
		
	}
	db_session->commit();
	
	//Currently, fuck WUM, message sending is way to slow compared to the RMSD comparisons
	// Now define the structure of the algorithm, i.e. assign roles to individual nodes
//	WorkUnitManagerOP wu_manager;
//	core::Size master_id = 0; /* only 1 master for this protocol */
//	if ( mpi_rank() == 0 ){
//		wu_manager = new DatabaseQueryWorkUnitManager<BundlePairRmsdWorkUnit>(master_id, db_session, select_string, db_wu_name);
////		wu_manager = new DatabaseQueryWorkUnitManager(master_id, db_session, select_string, "dbwu");
//	}
//	else
//	{
//		wu_manager = new MPI_WorkUnitManager_Slave(master_id);
//	}
//	
//	// make sure all the necessary work unit have been properly initiated.
//	wu_manager->register_work_units( wulist );
//	
//	// now launch this node
//	wu_manager->go();
//	
//#ifdef USEMPI
//	MPI_Barrier( MPI_COMM_WORLD );
//	MPI_Finalize();
//#endif
	return 0;
}




