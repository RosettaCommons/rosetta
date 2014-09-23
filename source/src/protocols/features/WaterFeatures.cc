// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/WaterFeatures.cc
/// @brief
/// @author Kevin Houlihan

// Unit Headers
#include <protocols/features/WaterFeatures.hh>

//External
#include <cppdb/frontend.h>

// Platform Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <string>
#include <sstream>

#include <cppdb/frontend.h>

#include <basic/Tracer.hh>
using basic::Tracer;
static thread_local basic::Tracer TR( "protocols.features.WaterFeatures" );

namespace protocols{
namespace features{

//using std::string;
//using std::stringstream;
//using std::endl;
//using core::Size;
//using core::Real;
//using core::Distance;
using core::Vector;
using core::Length;
using core::Angle;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;
using utility::tag::TagCOP;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;

using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

WaterFeatures::WaterFeatures() :
	FeaturesReporter(),
	acc_dist_cutoff_(3),
	don_dist_cutoff_(3)
{}

WaterFeatures::WaterFeatures(
	core::Length acc_dist_cutoff, core::Length don_dist_cutoff) :
	FeaturesReporter(),
	acc_dist_cutoff_(acc_dist_cutoff),
	don_dist_cutoff_(don_dist_cutoff)
{}

WaterFeatures::WaterFeatures(
	WaterFeatures const & src ) :
	FeaturesReporter(),
	acc_dist_cutoff_(src.acc_dist_cutoff_),
	don_dist_cutoff_(src.don_dist_cutoff_)
{}

std::string
WaterFeatures::type_name() const { return "WaterFeatures"; }

void
WaterFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_water_hbond_geom_table_schema(db_session);
}

void
WaterFeatures::write_water_hbond_geom_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column unrecognized_atom_res("unrecognized_atom_res", DbDataTypeOP( new DbInteger() ));
	Column unrecognized_atom_name("unrecognized_atom_name", DbDataTypeOP( new DbText(255) ));
	Column partner_site_id("partner_site_id", DbDataTypeOP( new DbInteger() ));
	// donor only
	Column whdist("WHdist", DbDataTypeOP( new DbReal() ));
	Column whd("WHD", DbDataTypeOP( new DbReal() ));
	// acceptors only
	Column awdist("AWdist", DbDataTypeOP( new DbReal() ));
	Column baw("BAW", DbDataTypeOP( new DbReal() ));
	Column chi("chi", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(unrecognized_atom_res);
	primary_key_columns.push_back(unrecognized_atom_name);
	primary_key_columns.push_back(partner_site_id);
	PrimaryKey primary_key1(primary_key_columns);
	PrimaryKey primary_key2(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(unrecognized_atom_res);
	foreign_key_columns1.push_back(unrecognized_atom_name);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("residue_number");
	reference_columns1.push_back("atom_name");
	ForeignKey foreign_key1(foreign_key_columns1, "unrecognized_atoms", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(partner_site_id);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("site_id");
	ForeignKey foreign_key2(foreign_key_columns2, "hbond_sites", reference_columns2, true);

	// TODO Add back foreign keys for both tables
	Schema table1("water_hbond_acceptors", primary_key1);
	//table1.add_foreign_key(foreign_key1);
	//table1.add_foreign_key(foreign_key2);
	table1.add_column(whdist);
	table1.add_column(whd);

	table1.write(db_session);

	Schema table2("water_hbond_donors", primary_key2);
	//table2.add_foreign_key(foreign_key1);
	//table2.add_foreign_key(foreign_key2);
	table2.add_column(awdist);
	table2.add_column(baw);
	table2.add_column(chi);

	table2.write(db_session);
}


utility::vector1<std::string>
WaterFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("UnrecognizedAtomFeatures");
	return dependencies;
}

void
WaterFeatures::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {

	// TODO add logic to parse targets=name3:atom_name,... tag option
  	if ( tag->hasOption("targets") ) {
		//stringstream error_msg;
		//error_msg
		//	<< "The " << type_name() << " reporter requires a 'targets' tag:" << endl
		//	<< endl
		//	<< "    <feature name=" << type_name() <<" targets=comma_separated_resName:atmName_list />" << endl;
		//throw utility::excn::EXCN_RosettaScriptsOption(error_msg.str());
		std::string water_names_list = tag->getOption< std::string >("targets");
		utility::vector0< std::string > const water_names( utility::string_split( water_names_list, ',' ) );
  		for ( utility::vector0< std::string>::const_iterator water_name( water_names.begin() ),
				end( water_names.end() ); water_name != end; ++water_name ) {
			utility::vector0< std::string > const resn_atomn( utility::string_split( *water_name, ':' ) );
			std::string resname = resn_atomn[0];
			std::string atomname = resn_atomn[1];
			names_for_water_.push_back(std::make_pair(resname, atomname));
		}
	}
	else {
		// resn:atomn, resn:atomn ...
		// HOH:O, DOD:O, WAT:O
		names_for_water_.push_back(std::make_pair("HOH", "O"));
		names_for_water_.push_back(std::make_pair("DOD", "O"));
		names_for_water_.push_back(std::make_pair("WAT", "O"));
	}
	for ( utility::vector1< std::pair<std::string, std::string> >::iterator
			ur_ua_it = names_for_water_.begin(), ur_ua_preend = names_for_water_.end();
			ur_ua_it != ur_ua_preend; ++ur_ua_it) {
		TR << "Name for water: " << "resName: " << ur_ua_it->first << ", atomName: " << ur_ua_it->second << std::endl;
	}

	acc_dist_cutoff_ = tag->getOption<core::Length>("acc_dist_cutoff", 3.0);
	don_dist_cutoff_ = tag->getOption<core::Length>("don_dist_cutoff", 3.0);
	ahd_cutoff_ = tag->getOption<core::Length>("ahd_cutoff", 90.0);
}


core::Size
WaterFeatures::report_features(
	Pose const & pose,
	vector1< bool > const &,
	StructureID struct_id,
	sessionOP db_session
){
	core::pose::PDBInfoCOP pdb_info(pose.pdb_info());
	if(!pdb_info) return 0;

	InsertGenerator igen_wat_accepts("water_hbond_acceptors");
	igen_wat_accepts.add_column("struct_id");
	igen_wat_accepts.add_column("unrecognized_atom_res");
	igen_wat_accepts.add_column("unrecognized_atom_name");
	igen_wat_accepts.add_column("partner_site_id");
	igen_wat_accepts.add_column("WHdist");
	igen_wat_accepts.add_column("WHD");

	InsertGenerator igen_wat_donates("water_hbond_donors");
	igen_wat_donates.add_column("struct_id");
	igen_wat_donates.add_column("unrecognized_atom_res");
	igen_wat_donates.add_column("unrecognized_atom_name");
	igen_wat_donates.add_column("partner_site_id");
	igen_wat_donates.add_column("AWdist");
	igen_wat_donates.add_column("BAW");
	igen_wat_donates.add_column("chi");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id", struct_id) );

	// locate candidate polar sites from the hbond_sites table
	std::stringstream water_ss;

	water_ss <<
		"SELECT\n"
		"	ua.residue_number,\n"
		"	ua.atom_name,\n"
		"   ua.coord_x,\n"
		"   ua.coord_y,\n"
		"   ua.coord_z \n"
		"FROM\n"
		"	unrecognized_atoms AS ua\n"
		"JOIN\n"
		"   unrecognized_residues AS ur\n"
		"ON\n"
		"   ua.struct_id = ur.struct_id AND\n"
		"   ua.residue_number = ur.residue_number\n"
		"WHERE\n"
		"	ua.struct_id = ? AND\n"
		"	(\n";

	for ( utility::vector1< std::pair<std::string, std::string> >::iterator
			ur_ua_it = names_for_water_.begin(), ur_ua_preend = names_for_water_.end() - 1;
			ur_ua_it != ur_ua_preend; ++ur_ua_it) {
		water_ss << "   (ur.name3 = ? AND ua.atom_name = ?) OR\n";
	}
	water_ss << "   (ur.name3 = ? and ua.atom_name = ?) );\n";
	std::string water_string(water_ss.str());
	statement water_stmt(basic::database::safely_prepare_statement(water_string,db_session));
//	water_stmt.bind(1, struct_id);
//	for (Size i = 1; i <= names_for_water_.size(); ++i) {
//		water_stmt.bind(2*i,names_for_water_[i].first);
//		water_stmt.bind(2*i+1,names_for_water_[i].second);
//	}

	TR << "Water query template:" << water_string << std::endl;

	water_stmt << struct_id;
	for (Size i = 1; i <= names_for_water_.size(); ++i) {
		water_stmt << names_for_water_[i].first;
		water_stmt << names_for_water_[i].second;
	}


	std::string hbond_site_string =
		"SELECT\n"
		"	hs.site_id,\n"
		"	hs.is_donor,\n"
		"	hsa.atm_x,\n"
		"	hsa.atm_y,\n"
		"	hsa.atm_z,\n"
		"	hsa.base_x,\n"
		"	hsa.base_y,\n"
		"	hsa.base_z,\n"
		"	hsa.base2_x,\n"
		"	hsa.base2_y,\n"
		"	hsa.base2_z \n"
		"FROM\n"
		"	hbond_sites AS hs\n"
		"JOIN\n"
		"	hbond_site_atoms AS hsa\n"
		"ON\n"
		"	hs.struct_id = hsa.struct_id AND\n"
		"	hs.site_id = hsa.site_id\n"
		"WHERE\n"
		"	hs.struct_id = ?;\n";

	statement hbond_site_stmt(basic::database::safely_prepare_statement(hbond_site_string,db_session));
	hbond_site_stmt.bind(1,struct_id);

	result wat_result(basic::database::safely_read_from_database(water_stmt));

	//core::Real dist, cosBAH, cosAHD, chi;


	TR << "Checking water hbond site combinations for struct_id" << struct_id << std::endl;
	while(wat_result.next()){
		TR << "Checking a water site for nearby hbond-sites." << std::endl;
		Size unrecognized_atom_res;
		std::string unrecognized_atom_name;
		wat_result >> unrecognized_atom_res;
		wat_result >> unrecognized_atom_name;
		TR << "Looking at water with resNum " << unrecognized_atom_res <<
			" and name " << unrecognized_atom_name << std::endl;
		Length wat_x, wat_y, wat_z;
		wat_result >> wat_x;
		wat_result >> wat_y;
		wat_result >> wat_z;
		Vector wat_xyz(wat_x, wat_y, wat_z);

		result hbond_site_result(basic::database::safely_read_from_database(hbond_site_stmt));
		while(hbond_site_result.next()){
			//TR << "Testing water-hbond site for potential interaction." << std::endl;
			Size partner_site_id;
			hbond_site_result >> partner_site_id;
			TR << "Looking at hbond site with site_id " << partner_site_id << std::endl;
			Size hbond_site_is_donor_int;
			hbond_site_result >> hbond_site_is_donor_int;
			bool hbond_site_is_donor = static_cast<bool>(hbond_site_is_donor_int);
			Length hbond_site_x, hbond_site_y, hbond_site_z;
			hbond_site_result >> hbond_site_x;
			hbond_site_result >> hbond_site_y;
			hbond_site_result >> hbond_site_z;
			Vector hbond_site_xyz(hbond_site_x, hbond_site_y, hbond_site_z);

			Length hbond_site_base_x, hbond_site_base_y, hbond_site_base_z;
			hbond_site_result >> hbond_site_base_x;
			hbond_site_result >> hbond_site_base_y;
			hbond_site_result >> hbond_site_base_z;
			Vector hbond_site_base_xyz(hbond_site_base_x, hbond_site_base_y, hbond_site_base_z);

			Length hbond_site_base2_x, hbond_site_base2_y, hbond_site_base2_z;
			hbond_site_result >> hbond_site_base2_x;
			hbond_site_result >> hbond_site_base2_y;
			hbond_site_result >> hbond_site_base2_z;
			Vector hbond_site_base2_xyz(hbond_site_base2_x, hbond_site_base2_y, hbond_site_base2_z);

			if (hbond_site_is_donor) {
				core::Length whdist = hbond_site_xyz.distance(wat_xyz);
				TR << "Hbond site is a donor." << std::endl;
				if (whdist > don_dist_cutoff_) continue;
				TR << "Hbond site is within interaction range." << std::endl;
				Angle whd = numeric::angle_degrees(hbond_site_base_xyz, hbond_site_xyz,
						                              wat_xyz);

				TR << "Attempting write to water_hbonds_acceptors with struct_id="
					<< struct_id << ", unrecognized_atom_res="
					<< unrecognized_atom_res << ", unrecognized_atom_name="
					<< unrecognized_atom_name << ", partner_site_id="
					<< partner_site_id << std::endl;
				igen_wat_accepts.add_row(
						utility::tools::make_vector(
							struct_id_data,
							RowDataBaseOP( new RowData<Size>(
									"unrecognized_atom_res",
									unrecognized_atom_res) ),
							RowDataBaseOP( new RowData<std::string>(
									"unrecognized_atom_name",
									unrecognized_atom_name) ),
							RowDataBaseOP( new RowData<core::Real>("partner_site_id", partner_site_id) ),
							RowDataBaseOP( new RowData<core::Real>("WHdist", whdist) ),
							RowDataBaseOP( new RowData<core::Real>("WHD", whd) )));
			}
			else {
				core::Length awdist = hbond_site_xyz.distance(wat_xyz);
				TR << "Hbond site is an acceptor." << std::endl;
				if (awdist > acc_dist_cutoff_) continue;
				TR << "Hbond site is within interaction range." << std::endl;
				Angle baw = numeric::angle_degrees(hbond_site_base_xyz, hbond_site_xyz,
						                              wat_xyz);
				Angle chi = numeric::angle_degrees(hbond_site_base2_xyz,
						                              hbond_site_base_xyz,
						                              hbond_site_xyz, wat_xyz);
				TR << "Attempting write to water_hbonds_donors with struct_id="
					<< struct_id << ", unrecognized_atom_res="
					<< unrecognized_atom_res << ", unrecognized_atom_name="
					<< unrecognized_atom_name << ", partner_site_id="
					<< partner_site_id << std::endl;
				igen_wat_donates.add_row(
						utility::tools::make_vector(
							struct_id_data,
							RowDataBaseOP( new RowData<Size>(
									"unrecognized_atom_res",
									unrecognized_atom_res) ),
							RowDataBaseOP( new RowData<std::string>(
									"unrecognized_atom_name",
									unrecognized_atom_name) ),
							RowDataBaseOP( new RowData<Size>(
									"partner_site_id",
									partner_site_id) ),
							RowDataBaseOP( new RowData<core::Real>("AWdist", awdist) ),
							RowDataBaseOP( new RowData<core::Real>("BAW", baw) ),
							RowDataBaseOP( new RowData<core::Real>("chi", chi) )));
			}
		} //while(hbond_site_result.next()){
	} //while(wat_result.next()){

	igen_wat_accepts.write_to_database(db_session);
	igen_wat_donates.write_to_database(db_session);

	return 0;
}

} //namesapce
} //namespace
