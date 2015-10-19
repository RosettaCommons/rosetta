// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BundlePairRmsdWorkUnit.cc
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

//Unit
#include <devel/sewing/BundlePairRmsdWorkUnit.hh>

//Basic
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>

//Utility
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>

//C++
#include <string>
#include <map>

namespace devel {
namespace sewing {
	
static THREAD_LOCAL basic::Tracer TR( "BundlePairRmsdWorkUnit" );

using namespace std;

BundlePairRmsdWorkUnit::BundlePairRmsdWorkUnit(utility::sql_database::sessionOP db_session):
DatabaseEntryWorkUnit(db_session)
{}

BundlePairRmsdWorkUnit::BundlePairRmsdWorkUnit( std::map<std::string,std::string> row_map ):
DatabaseEntryWorkUnit(row_map)
{}

/// @brief Calculate the pair RMSD using data from the results_map_
void
BundlePairRmsdWorkUnit::run(){
	using namespace cppdb;
	using namespace utility;

	TR << "BundlePairRMSD work unit recieved following key/value pairs:" << endl;

	for( map<string,string>::const_iterator it = row_map_.begin();
		it != row_map_.end(); ++it){
		TR << it->first << ", " << it->second << endl;
	}

	
  core::Size struct_id_1 = string2int(row_map_["struct_id_1"]);
	core::Size struct_id_2 = string2int(row_map_["struct_id_2"]);

	core::Size pair_id_1 = string2int(row_map_["pair_id_1"]);

	bool pair_1_helix_1_flipped = string2int(row_map_["pair_1_helix_1_flipped"]);
	core::Size pair_1_helix_1_begin = string2int(row_map_["pair_1_helix_1_begin"]);
	core::Size pair_1_helix_1_end = string2int(row_map_["pair_1_helix_1_end"]);

	bool pair_1_helix_2_flipped = string2int(row_map_["pair_1_helix_2_flipped"]);
	core::Size pair_1_helix_2_begin = string2int(row_map_["pair_1_helix_2_begin"]);
	core::Size pair_1_helix_2_end = string2int(row_map_["pair_1_helix_2_end"]);

	core::Size pair_id_2 = string2int(row_map_["pair_id_2"]);

	bool pair_2_helix_1_flipped = string2int(row_map_["pair_2_helix_1_flipped"]);
	core::Size pair_2_helix_1_begin = string2int(row_map_["pair_2_helix_1_begin"]);
	core::Size pair_2_helix_1_end = string2int(row_map_["pair_2_helix_1_end"]);

	bool pair_2_helix_2_flipped = string2int(row_map_["pair_2_helix_2_flipped"]);
	core::Size pair_2_helix_2_begin = string2int(row_map_["pair_2_helix_2_begin"]);
	core::Size pair_2_helix_2_end = string2int(row_map_["pair_2_helix_2_end"]);

	//Need to take into account the flipped status of a helix

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

	statement select_coords_stmt(basic::database::safely_prepare_statement(helix_coords_select,db_session_));

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

//		TR << "Coords for pos " << counter << " of pair 1: " << x1 << " " << y1 << " " << z1 << std::endl;
//		TR << "Coords for pos " << counter << " of pair 2: " << x2 << " " << y2 << " " << z2 << std::endl;
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

	core::Real bundle_pair_rmsd = numeric::model_quality::calc_rms(pair_1_coords,pair_2_coords);

	std::string comparison_insert = "INSERT INTO pair_comparisons(pair_id_1, pair_id_2,rmsd) VALUES(?,?,?);";
	cppdb::statement comparison_insert_stmt(basic::database::safely_prepare_statement(comparison_insert,db_session_));
	comparison_insert_stmt.bind(1, pair_id_1);
	comparison_insert_stmt.bind(2, pair_id_2);
	comparison_insert_stmt.bind(3, bundle_pair_rmsd);
	basic::database::safely_write_to_database(comparison_insert_stmt);

	TR << "Bundle pair " << pair_id_1 << " and pair " << pair_id_2 << " RMSD: " << bundle_pair_rmsd << std::endl;
}

} //sewing namespace
} //devel namespace
