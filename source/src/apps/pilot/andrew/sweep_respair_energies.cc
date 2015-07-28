// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/andrew/sweep_respair_energies
/// @brief  Output a database with the sets of residue pair energies swept over parameters given in a match constraint file
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <devel/init.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/ResidueKinWriter.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/import_pose/import_pose.hh>

// Protocol headers
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.hh>
#include <protocols/toolbox/match_enzdes_util/LigandConformer.hh>

// Numeric includes
#include <numeric/random/random.hh>
#include <numeric/constants.hh>

// Basic includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

// Utility headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/LexicographicalIterator.hh>

// C++ headers
#include <fstream>
#include <sstream>

// External Headers
#include <cppdb/frontend.h>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <ObjexxFCL/format.hh>

OPT_1GRP_KEY( Boolean, sweep_respair_energies, output_bestpose_only )
OPT_1GRP_KEY( Boolean, sweep_respair_energies, find_local_minima )

class hbgeom_to_measure {
public:
	core::Size r1atind; // the acceptor
	core::Size r2atind; // the donor
	utility::vector1< std::string > hbgeoms_;
};

class local_minimum {
public:
	core::conformation::ResidueOP res1_conformation;
	core::Real score_;
};

typedef utility::vector1< hbgeom_to_measure > hbgeoms_to_measure;
typedef numeric::HomogeneousTransform< core::Real > HTReal;


void
create_schema(
	utility::sql_database::sessionOP db_session
)
{
	using namespace basic::database;
	using namespace basic::database::schema_generator;

	Column respair_id( "respair_id", DbDataTypeOP( new DbUUID ) );
	Column respair_name( "respair_name", DbDataTypeOP( new DbText ) );
	Column res1name( "res1name", DbDataTypeOP( new DbText ) );
	Column res2name( "res2name", DbDataTypeOP( new DbText ) );
	Column focused_geom_param( "focused_geom_param", DbDataTypeOP( new DbText ) );
	PrimaryKey respair_id_pk( respair_id );
	Schema table_residue_pairs( "residue_pairs", respair_id_pk );
	table_residue_pairs.add_column( respair_name );
	table_residue_pairs.add_column( res1name );
	table_residue_pairs.add_column( res2name );
	table_residue_pairs.add_column( focused_geom_param );
	table_residue_pairs.write( db_session );

	Column conf_id( "conf_id", DbDataTypeOP( new DbUUID ) );
	PrimaryKey conf_id_pk( conf_id );
	ForeignKey respair_id_fk( respair_id, "residue_pairs", "respair_id" );
	Schema table_conformations( "conformations", conf_id_pk );
	table_conformations.add_foreign_key( respair_id_fk );
	table_conformations.write( db_session );

	Column geom_name( "geom_name", DbDataTypeOP( new DbText ) );
	PrimaryKey geom_name_pk( geom_name );
	Column geom_desc( "geom_desc", DbDataTypeOP( new DbText ) );
	Schema table_geom_types( "geom_types", geom_name_pk );
	table_geom_types.add_column( geom_desc );
	table_geom_types.write( db_session );

	Columns respair_id_and_geom_name;
	respair_id_and_geom_name.push_back( respair_id );
	respair_id_and_geom_name.push_back( geom_name );
	PrimaryKey respair_id_and_geom_name_pk( respair_id_and_geom_name );
	ForeignKey geom_name_fk( geom_name, "geom_types", "geom_name" );
	Column at1( "at1", DbDataTypeOP( new DbText ) );
	Column res1( "res1", DbDataTypeOP( new DbInteger ) );
	Column at2( "at2", DbDataTypeOP( new DbText ) );
	Column res2( "res2", DbDataTypeOP( new DbInteger ) );
	Column at3( "at3", DbDataTypeOP( new DbText ) );
	Column res3( "res3", DbDataTypeOP( new DbInteger ) );
	Column at4( "at4", DbDataTypeOP( new DbText ) );
	Column res4( "res4", DbDataTypeOP( new DbInteger ) );
	Schema table_geom_definitions( "geom_definitions", respair_id_and_geom_name_pk );
	table_geom_definitions.add_foreign_key( geom_name_fk );
	table_geom_definitions.add_foreign_key( respair_id_fk );
	table_geom_definitions.add_column( at1 );
	table_geom_definitions.add_column( res1 );
	table_geom_definitions.add_column( at2 );
	table_geom_definitions.add_column( res2 );
	table_geom_definitions.add_column( at3 );
	table_geom_definitions.add_column( res3 );
	table_geom_definitions.add_column( at4 );
	table_geom_definitions.add_column( res4 );
	table_geom_definitions.write( db_session );

	Columns conf_id_and_geom_name;
	conf_id_and_geom_name.push_back( conf_id );
	conf_id_and_geom_name.push_back( geom_name );
	PrimaryKey conf_id_and_geom_name_pk( conf_id_and_geom_name );
	ForeignKey conf_id_fk( conf_id, "conformations", "conf_id" );
	Column geom_value( "geom_value", DbDataTypeOP( new DbReal ) );
	Schema table_geometries( "geometries", conf_id_and_geom_name_pk );
	table_geometries.add_foreign_key( conf_id_fk );
	table_geometries.add_foreign_key( geom_name_fk );
	table_geometries.add_column( geom_value );
	table_geometries.write( db_session );

	Column score_name( "score_name", DbDataTypeOP( new DbText ) );
	PrimaryKey score_name_pk( score_name );
	Schema table_score_types( "score_types", score_name_pk );
	table_score_types.write( db_session );

	Columns conf_id_and_score_name;
	conf_id_and_score_name.push_back( conf_id );
	conf_id_and_score_name.push_back( score_name );
	PrimaryKey conf_id_and_score_name_pk( conf_id_and_score_name );
	ForeignKey score_name_fk( score_name, "score_types", "score_name" );
	Column score_value( "score_value", DbDataTypeOP( new DbReal ) );
	Schema table_conf_scores( "conf_scores", conf_id_and_score_name_pk );
	table_conf_scores.add_foreign_key( conf_id_fk );
	table_conf_scores.add_foreign_key( score_name_fk );
	table_conf_scores.add_column( score_value );
	table_conf_scores.write( db_session );

}

/// @details write to database all the score types possible; there may be overlap if this is called
/// multiple times on a single database, (e.g. multiple invocations) so I'm using the "INSERT OR IGNORE"
/// SQL statement.  This code is basically a complete rip off of ScoreTypeFeatures::insert_score_type_rows
void
write_score_types_to_database(
	utility::sql_database::sessionOP db_session
)
{
	using cppdb::statement;
	using namespace core::scoring;

	std::string statement_string;

	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3:
		statement_string = "INSERT OR IGNORE INTO score_types (score_name) VALUES (?);";
		break;
	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres:
		statement_string = "INSERT IGNORE INTO score_tyeps (score_name) VALUES (?);";
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for ( core::Size score_type_id=1; score_type_id <= n_score_types; ++score_type_id ) {
		ScoreType type(static_cast<ScoreType>(score_type_id));
		std::string const score_type( ScoreTypeManager::name_from_score_type(type) );
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);
	}

}

void
write_hbgeom_types_to_database(
	utility::sql_database::sessionOP db_session

)
{
	using cppdb::statement;

	std::string statement_string( "INSERT OR IGNORE INTO geom_types (geom_name,geom_desc) VALUES (?,?);" );
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	std::string geom_name, geom_desc;

	geom_name = "AHD";
	geom_desc = "Angle in a hydrogen bond between the acceptor, the hydrogen, and the donor heavy atom, in radians";
	stmt.bind(1,geom_name);
	stmt.bind(2,geom_desc);
	basic::database::safely_write_to_database(stmt);

	geom_name = "BAH";
	geom_desc = "Angle in a hydrogen bond between the acceptor-base, the acceptor, and the hydrogen in radians";
	stmt.bind(1,geom_name);
	stmt.bind(2,geom_desc);
	basic::database::safely_write_to_database(stmt);

	geom_name = "AHdist";
	geom_desc = "Distance in a hydrogen bond between the acceptor and the hydrogen in Angstroms";
	stmt.bind(1,geom_name);
	stmt.bind(2,geom_desc);
	basic::database::safely_write_to_database(stmt);

	geom_name = "chi";
	geom_desc = "Dihedral in a hydrogen bond defined by the abase2, the acceptor-base, the acceptor and the hydrogen in radians";
	stmt.bind(1,geom_name);
	stmt.bind(2,geom_desc);
	basic::database::safely_write_to_database(stmt);

}

/// @brief return the new respair_id uuid created by this insertion
boost::uuids::uuid
write_new_respair_to_database(
	utility::sql_database::sessionOP db_session,
	protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo const & mcfi
)
{
	using cppdb::statement;
	using namespace protocols::toolbox::match_enzdes_util;

	std::string acc3let = mcfi.allowed_res_name3s(1)[1];
	std::string don3let = mcfi.allowed_res_name3s(2)[1];

	/// 1. Figure out what parameters are being swept
	ExternalGeomSamplerOP exgeom = mcfi.create_exgs();
	utility::vector1< core::Size > nexgeoms( 6, 0 );
	nexgeoms[1] = exgeom->n_tor_U3D1_samples();
	nexgeoms[2] = exgeom->n_ang_U2D1_samples(); //
	nexgeoms[3] = exgeom->n_dis_U1D1_samples();
	nexgeoms[4] = exgeom->n_ang_U1D2_samples();
	nexgeoms[5] = exgeom->n_tor_U2D2_samples();
	nexgeoms[6] = exgeom->n_tor_U1D3_samples();

	std::string sweep_string( "sweep" );
	if ( nexgeoms[1] > 1 ) { sweep_string += "_toru3d1"; }
	if ( nexgeoms[2] > 1 ) { sweep_string += "_AHD"; }
	if ( nexgeoms[3] > 1 ) { sweep_string += "_AHdist"; }
	if ( nexgeoms[4] > 1 ) { sweep_string += "_BAH"; }
	if ( nexgeoms[5] > 1 ) { sweep_string += "_toru2d2"; }
	if ( nexgeoms[6] > 1 ) { sweep_string += "_chi"; }

	boost::uuids::basic_random_generator<numeric::random::RandomGenerator> uuids_rng(numeric::random::rg());

	/// 2. insert into residue_pairs
	boost::uuids::uuid respair_id = uuids_rng();
	std::string statement_string( "INSERT INTO residue_pairs (respair_id,respair_name,res1name,res2name,focused_geom_param) VALUES (?,?,?,?,?);" );
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	std::string respair_name = acc3let + "_" + don3let + "_" + sweep_string;
	stmt.bind(1, respair_id );
	stmt.bind(2, respair_name);
	stmt.bind(3, acc3let );
	stmt.bind(4, don3let );
	stmt.bind(5, sweep_string );
	basic::database::safely_write_to_database(stmt);

	return respair_id;
}

void
write_hbgeom_to_db(
	utility::sql_database::sessionOP db_session,
	boost::uuids::uuid const conf_id,
	std::string const & geom_name,
	core::Real geometry
)
{
	using cppdb::statement;

	std::string statement_string( "INSERT INTO geometries (conf_id,geom_name,geom_value) VALUES (?,?,?);" );
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,conf_id);
	stmt.bind(2,geom_name);
	stmt.bind(3,geometry);
	basic::database::safely_write_to_database(stmt);
}

void
write_ahdist_to_db(
	utility::sql_database::sessionOP db_session,
	boost::uuids::uuid const conf_id,
	core::conformation::Residue const & acc_res,
	core::conformation::Residue const & don_res,
	core::Size acc_ind,
	core::Size hdon_ind
){
	core::Real dist = acc_res.xyz( acc_ind ).distance( don_res.xyz( hdon_ind ) );
	std::string geom_name( "AHdist" );
	write_hbgeom_to_db( db_session, conf_id, geom_name, dist );
}

void
write_AHD_to_db(
	utility::sql_database::sessionOP db_session,
	boost::uuids::uuid const conf_id,
	core::conformation::Residue const & acc_res,
	core::conformation::Residue const & don_res,
	core::Size acc_ind,
	core::Size hdon_ind,
	core::Size don_ind
){
	core::Real AHDang = numeric::angle_radians( acc_res.xyz( acc_ind ), don_res.xyz( hdon_ind ), don_res.xyz( don_ind ));
	std::string geom_name( "AHD" );

	write_hbgeom_to_db( db_session, conf_id, geom_name, AHDang );
}

void
write_BAH_to_db(
	utility::sql_database::sessionOP db_session,
	boost::uuids::uuid const conf_id,
	core::conformation::Residue const & acc_res,
	core::conformation::Residue const & don_res,
	core::Size acc_base,
	core::Size acc_ind,
	core::Size hdon_ind
){
	core::Real BAHang = numeric::angle_radians( acc_res.xyz( acc_base ), acc_res.xyz( acc_ind ), don_res.xyz( hdon_ind ));
	std::string geom_name( "BAH" );
	write_hbgeom_to_db( db_session, conf_id, geom_name, BAHang );
}

void
write_chi_to_db(
	utility::sql_database::sessionOP db_session,
	boost::uuids::uuid const conf_id,
	core::conformation::Residue const & acc_res,
	core::conformation::Residue const & don_res,
	core::Size acc_base2,
	core::Size acc_base,
	core::Size acc_ind,
	core::Size hdon_ind
){
	core::Real chi = numeric::dihedral_radians(
		acc_res.xyz( acc_base2 ),
		acc_res.xyz( acc_base ),
		acc_res.xyz( acc_ind ),
		don_res.xyz( hdon_ind ));
	std::string geom_name( "chi" );
	write_hbgeom_to_db( db_session, conf_id, geom_name, chi );
}


void
write_hbgeoms_to_database(
	utility::sql_database::sessionOP db_session,
	boost::uuids::uuid const conf_id,
	core::pose::Pose const & pose,
	hbgeoms_to_measure const & geoms,
	core::Size res1ind,
	core::Size res2ind
)
{
	for ( core::Size ii = 1; ii <= geoms.size(); ++ii ) {
		core::Size acc_ind  = geoms[ii].r1atind;
		core::Size accbase  = pose.residue_type(res1ind).atom_base( acc_ind );
		core::Size accbase2 = pose.residue_type(res1ind).abase2( acc_ind );
		core::Size hdon_ind = geoms[ii].r2atind;
		core::Size don_ind  = pose.residue_type(res2ind).atom_base(hdon_ind);
		core::conformation::Residue const & r1 = pose.residue( res1ind );
		core::conformation::Residue const & r2 = pose.residue( res2ind );
		for ( core::Size jj = 1; jj <= geoms[ii].hbgeoms_.size(); ++jj ) {
			std::string const & jjgeom = geoms[ii].hbgeoms_[jj];
			if ( jjgeom == "AHdist" ) {
				write_ahdist_to_db( db_session, conf_id, r1, r2, acc_ind, hdon_ind );
			}
			if ( jjgeom == "AHD" ) {
				write_AHD_to_db( db_session, conf_id, r1, r2, acc_ind, hdon_ind, don_ind );
			}
			if ( jjgeom == "BAH" ) {
				write_BAH_to_db( db_session, conf_id, r1, r2, accbase, acc_ind, hdon_ind );
			}
			if ( jjgeom == "chi" ) {
				write_chi_to_db( db_session, conf_id, r1, r2, accbase2, accbase, acc_ind, hdon_ind );
			}
		}
	}
}

void
write_respair_scores_to_database(
	utility::sql_database::sessionOP db_session,
	boost::uuids::uuid const conf_id,
	core::scoring::ScoreFunction const & sfxn,
	core::pose::Pose & pose,
	core::Size res1ind,
	core::Size res2ind
)
{
	using cppdb::statement;
	using namespace core::scoring;
	sfxn( pose );
	EnergyMap emap;
	sfxn.eval_ci_2b( pose.residue( res1ind ), pose.residue( res2ind ), pose, emap );
	sfxn.eval_cd_2b( pose.residue( res1ind ), pose.residue( res2ind ), pose, emap );
	// MISSING LONG RANGE ENERGIES

	emap *= sfxn.weights();

	// Now write out the terms to the database
	std::string statement_string( "INSERT INTO conf_scores (conf_id,score_name,score_value) VALUES (?,?,?);" );
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	std::string score_name; core::Real score_val(0), total(0);
	stmt.bind(1,conf_id);

	//std::cout << "ahdist: " << pose.residue( res1ind ).atom( "OG" ).xyz().distance(
	//	pose.residue( res2ind ).atom( "1HD2" ).xyz() ) << " ";

	for ( core::Size ii = 1; ii <= n_score_types; ++ii ) {
		ScoreType iist = ScoreType(ii);
		if ( (sfxn.weights())[ iist ] == 0 ) continue;
		score_name = ScoreTypeManager::name_from_score_type( iist );
		score_val = emap[ iist ];
		total += score_val;
		//std::cout << iist << " " << score_val << " ";
		stmt.bind(2,score_name);
		stmt.bind(3,score_val);
		basic::database::safely_write_to_database(stmt);
	}
	//std::cout << "total: " << total << std::endl;
	score_name = "total";
	score_val = total;
	stmt.bind(2,score_name);
	stmt.bind(3,score_val);
	basic::database::safely_write_to_database(stmt);

}

void
write_scores_and_hbond_geoms_to_database(
	utility::sql_database::sessionOP db_session,
	boost::uuids::uuid const respair_id,
	core::scoring::ScoreFunction const & sfxn,
	core::pose::Pose & pose,
	hbgeoms_to_measure const & geoms,
	core::Size res1ind,
	core::Size res2ind
)
{
	using cppdb::statement;

	boost::uuids::basic_random_generator<numeric::random::RandomGenerator> uuids_rng(numeric::random::rg());

	// create a new conformation id for this residue pair
	// and insert this conformation into the "conformations" table
	boost::uuids::uuid conf_id = uuids_rng();
	{
		std::string statement_string( "INSERT INTO conformations (conf_id, respair_id) VALUES (?,?);" );
		statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
		stmt.bind(1,conf_id);
		stmt.bind(2,respair_id);
		basic::database::safely_write_to_database(stmt);
	}

	/// now insert the geometries for this residue pair into the geometries table
	write_hbgeoms_to_database( db_session, conf_id, pose, geoms, res1ind, res2ind );
	// now insert the residue pair scores into the conf_scores table
	write_respair_scores_to_database( db_session, conf_id, sfxn, pose, res1ind, res2ind );

}


/// A function to populate the database with fake data
/// in order to write the R scripts that will be necessary
/// to visualize the data that will be generated by this
/// program.
void
populate_database_with_dummy_data(
	utility::sql_database::sessionOP db_session
)
{

	using namespace basic::database;
	using cppdb::statement;
	//using namespace basic::database::schema_generator;

	// define score types, geom types
	// 1. score types
	{
		std::string statement_string( "INSERT INTO score_types (score_name) VALUES (?);" );
		statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
		std::string score_type;

		score_type = "fa_atr";
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);

		score_type = "fa_rep";
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);

		score_type = "fa_sol";
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);

		score_type = "fa_elec";
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);

		score_type = "hb_sr_bb";
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);

		score_type = "hb_lr_bb";
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);

		score_type = "hb_bb_sc";
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);

		score_type = "hb_sc";
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);

		score_type = "total";
		stmt.bind(1,score_type);
		basic::database::safely_write_to_database(stmt);

	}

	// geom types
	{
		std::string statement_string( "INSERT INTO geom_types (geom_name,geom_desc) VALUES (?,?);" );
		statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
		std::string geom_name, geom_desc;

		geom_name = "AHD";
		geom_desc = "Angle in a hydrogen bond between the acceptor, the hydrogen, and the donor heavy atom, in radians";
		stmt.bind(1,geom_name);
		stmt.bind(2,geom_desc);
		basic::database::safely_write_to_database(stmt);

		geom_name = "BAH";
		geom_desc = "Angle in a hydrogen bond between the acceptor-base, the acceptor, and the hydrogen in radians";
		stmt.bind(1,geom_name);
		stmt.bind(2,geom_desc);
		basic::database::safely_write_to_database(stmt);

		geom_name = "AHdist";
		geom_desc = "Distance in a hydrogen bond between the acceptor and the hydrogen in Angstroms";
		stmt.bind(1,geom_name);
		stmt.bind(2,geom_desc);
		basic::database::safely_write_to_database(stmt);
	}

	/// ok, now let's insert data for a fake hbond distance sweeps.
	/// 1. Iterate over pairs of residues
	///   2. describe the respair (insert into residue_pairs)
	///   3. -- skip -- give the atoms that define the interaction (insert into respair_id_and_geom_name)
	///   4. iterate:
	///     5. add a new confromation
	///     6. describe the geometries for this conformation
	///     7. describe the scores for this conformation

	/// residue pairs list:
	utility::vector1< std::string > acc_names;
	utility::vector1< std::string > acc_hbtypes;
	acc_names.push_back( "ASP" ); acc_names.push_back("ASN" ); acc_names.push_back( "SER" ); acc_names.push_back( "HIS" ); acc_names.push_back( "GLY" );
	acc_hbtypes.push_back( "CXL" ); acc_hbtypes.push_back( "CXA" ); acc_hbtypes.push_back( "HXL" );  acc_hbtypes.push_back( "IME" ); acc_hbtypes.push_back( "PBA" );
	utility::vector1< std::string > don_names;
	utility::vector1< std::string > don_hbtypes;
	don_names.push_back( "SER" ); don_names.push_back( "TYR" ); don_names.push_back( "HIS" ); don_names.push_back( "ARG" ); don_names.push_back( "LYS" ); don_names.push_back( "GLY" );
	don_hbtypes.push_back( "HXL" ); don_hbtypes.push_back( "AHX" ); don_hbtypes.push_back( "IME" );don_hbtypes.push_back( "GDH" );don_hbtypes.push_back( "AMO" );don_hbtypes.push_back( "PBA" );

	boost::uuids::basic_random_generator<numeric::random::RandomGenerator> uuids_rng(numeric::random::rg());


	for ( core::Size ii = 1; ii <= acc_names.size(); ++ii ) {
		std::string const & ii_acc_name = acc_names[ii];
		std::string const ii_acc_hbtype = acc_hbtypes[ii];
		for ( core::Size jj = 1; jj <= don_names.size(); ++jj ) {
			std::string const & jj_don_name = don_names[jj];
			std::string const jj_don_hbtype = don_hbtypes[jj];

			/// 2. insert into residue_pairs
			boost::uuids::uuid respair_id = uuids_rng();
			{
				std::string statement_string( "INSERT INTO residue_pairs (respair_id,respair_name,res1name,res2name,focused_geom_param) VALUES (?,?,?,?,?);" );
				statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
				std::string respair_name = ii_acc_hbtype + "_" + jj_don_hbtype + "_sweep_AHdist";
				std::string ahdist( "AHdist" );
				stmt.bind(1, respair_id );
				stmt.bind(2, respair_name);
				stmt.bind(3, ii_acc_name );
				stmt.bind(4, jj_don_name );
				stmt.bind(5, ahdist );
				basic::database::safely_write_to_database(stmt);
			}
			///// 2. describe the atoms in the geometries
			//{
			//	std::string statement_string( "INSERT INTO geom_definitions (respair_id,geom_name,at1,res1,at2,res2,at3,res3,at4,res4) VALUES (?,?,?,?,?,?,?,?,?,?);" );
			//	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
			//	std::string geom_name,at1,at2,at3,at4;
			//	core::Size res1(0),res2(0),res3(0),res4(0);
			//
			//	stmt.bind(1,respair_id);
			//
			//	geom_name = "AHD";   stmt.bind(2,geom_name);
			//	at1 = "OD1"; res1=1; stmt.bind( 3,at1); stmt.bind( 4,res1);
			//	at2 = "HG";  res2=2; stmt.bind( 5,at2); stmt.bind( 6,res2);
			//	at3 = "OG";  res3=2; stmt.bind( 7,at3); stmt.bind( 8,res3);
			//	at4 = "";    res4=0; stmt.bind( 9,at4); stmt.bind(10,res4);
			//	basic::database::safely_write_to_database(stmt);
			//
			//	geom_name = "BAH";   stmt.bind(2,geom_name);
			//	at1 = "CG";  res1=1; stmt.bind( 3,at1); stmt.bind( 4,res1);
			//	at2 = "OD1"; res2=1; stmt.bind( 5,at2); stmt.bind( 6,res2);
			//	at3 = "HG";  res3=2; stmt.bind( 7,at3); stmt.bind( 8,res3);
			//	at4 = "";    res4=0; stmt.bind( 9,at4); stmt.bind(10,res4);
			//	basic::database::safely_write_to_database(stmt);
			//
			//	geom_name = "AHdist";   stmt.bind(2,geom_name);
			//	at1 = "OD1"; res1=1; stmt.bind( 3,at1); stmt.bind( 4,res1);
			//	at2 = "HG";  res2=2; stmt.bind( 5,at2); stmt.bind( 6,res2);
			//	at3 = "";    res3=0; stmt.bind( 7,at3); stmt.bind( 8,res3);
			//	at4 = "";    res4=0; stmt.bind( 9,at4); stmt.bind(10,res4);
			//	basic::database::safely_write_to_database(stmt);
			//}

			/// 3. Let's iterate across a distance span
			for ( core::Real dist = 1.5; dist < 4.0; dist += 0.1 ) {
				/// Insert a new conformation into the conformations table
				boost::uuids::uuid conf_id = uuids_rng();
				{
					std::string statement_string( "INSERT INTO conformations (conf_id, respair_id) VALUES (?,?);" );
					statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
					stmt.bind(1,conf_id);
					stmt.bind(2,respair_id);
					basic::database::safely_write_to_database(stmt);
				}
				/// Insert the distance geometry into the geometries table
				{
					std::string statement_string( "INSERT INTO geometries (conf_id,geom_name,geom_value) VALUES (?,?,?);" );
					statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
					std::string geom_name( "AHdist" );
					stmt.bind(1,conf_id);
					stmt.bind(2,geom_name);
					stmt.bind(3,dist);
					basic::database::safely_write_to_database(stmt);
				}
				/// Now insert scores for all the terms
				{
					std::string statement_string( "INSERT INTO conf_scores (conf_id,score_name,score_value) VALUES (?,?,?);" );
					statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
					std::string score_name; core::Real score_val;

					score_name = "fa_rep";
					score_val = std::max( 0., 2-dist);
					stmt.bind(1,conf_id);
					stmt.bind(2,score_name);
					stmt.bind(3,score_val);
					basic::database::safely_write_to_database(stmt);

					score_name = "fa_atr";
					score_val = dist > 2.0 ? -4.0 / std::pow( dist, 6 ) : -4.0 / std::pow( 2.0, 6 );
					stmt.bind(1,conf_id);
					stmt.bind(2,score_name);
					stmt.bind(3,score_val);
					basic::database::safely_write_to_database(stmt);

					score_name = "fa_elec";
					score_val = dist > 2.0 ? 4.0 / dist : 4.0 / 2.0;
					stmt.bind(1,conf_id);
					stmt.bind(2,score_name);
					stmt.bind(3,score_val);
					basic::database::safely_write_to_database(stmt);

					score_name = ii_acc_hbtype == "PBA" ? ( jj_don_hbtype == "PBA" ? "hb_lr_bb" : "hb_bb_sc" ) : ( jj_don_hbtype == "PBA" ? "hb_bb_sc" : "hb_sc_sc" );
					score_val = std::min( 0., -1*std::sin( (dist-1.5) / 1.5 * numeric::constants::d::pi ) );
					stmt.bind(1,conf_id);
					stmt.bind(2,score_name);
					stmt.bind(3,score_val);
					basic::database::safely_write_to_database(stmt);
				}
			}
		}
	}
}

/// @brief assumes that the geometry is specified between the donor and the acceptor;
/// this is not necessarily true!  FIX THIS so that the hbond geometry can be defined
/// independently of the mat
hbgeoms_to_measure
create_fullcoverage_hbgeoms_to_measure(
	core::pose::Pose const & pose,
	protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo const & mcfi,
	core::Size acc_id,
	core::Size don_id
)
{
	hbgeoms_to_measure geoms;
	hbgeom_to_measure geom;
	geom.r1atind = mcfi.template_atom_inds(1,1,pose.residue_type(acc_id))[1];
	geom.r2atind = mcfi.template_atom_inds(2,1,pose.residue_type(don_id))[1];
	geom.hbgeoms_.push_back( "AHdist" );
	geom.hbgeoms_.push_back( "AHD" );
	geom.hbgeoms_.push_back( "BAH" );
	geom.hbgeoms_.push_back( "chi" );
	geoms.push_back( geom );
	return geoms;
}

core::Real
score_residue_pair(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::Size const res1_index,
	core::Size const res2_index
)
{
	core::scoring::EnergyMap emap;
	sfxn.eval_ci_2b( pose.residue( res1_index ), pose.residue( res2_index ), pose, emap );
	sfxn.eval_cd_2b( pose.residue( res1_index ), pose.residue( res2_index ), pose, emap );
	core::Real this_score = emap.dot( sfxn.weights() );
	return this_score;
}

core::Real
build_and_score_res1_against_res2(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	HTReal const & htend,
	core::Size res1_index,
	core::Size res2_index,
	protocols::toolbox::match_enzdes_util::LigandConformer const & res1_geom
)
{
	for ( core::Size ii=1, iiend = pose.residue(2).natoms(); ii <= iiend; ++ii ) {
		pose.set_xyz( core::id::AtomID( ii, res1_index ), res1_geom.coordinate_in_global_frame( ii, htend ) );
	}
	return score_residue_pair( pose, sfxn, res1_index, res2_index );
}

bool
determine_if_at_local_minimum_numeric(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	protocols::toolbox::match_enzdes_util::ExternalGeomSamplerCOP exgeom,
	HTReal const & launch,
	core::Size res1_index,
	core::Size res2_index,
	utility::LexicographicalIterator const & lex,
	protocols::toolbox::match_enzdes_util::LigandConformer const & res1_geom,
	core::Real & score
)
{
	/// ok: for each of the 6 DOFs, construct the position to both "below" and "above"
	/// of the given sample configuration and score those conformation.
	/// If the energies of these conformations are higher
	/// than the energy of the sample conformation, then count this structure as being at
	/// a local minimum.
	using namespace protocols::toolbox::match_enzdes_util;

	HTReal tu3d1 = exgeom->transform( HT_tor_U3D1, lex[1] );
	HTReal au2d1 = exgeom->transform( HT_ang_U2D1, lex[2] );
	HTReal du1d1;  du1d1.walk_along_z( exgeom->dis_U1D1_samples()[lex[3]] );
	HTReal tu2d2 = exgeom->transform( HT_tor_U2D2, lex[4] );
	HTReal au1d2 = exgeom->transform( HT_ang_U1D2, lex[5] );
	HTReal tu1d3 = exgeom->transform( HT_tor_U1D3, lex[6] );

	HTReal htend = launch*tu3d1*au2d1*du1d1*tu2d2*au1d2*tu1d3;
	score = build_and_score_res1_against_res2( pose, sfxn, htend, res1_index, res2_index, res1_geom );

	/// 1. tu3d1 -- periodic
	if ( lex.dimsize(1) != 1 ) {
		HTReal left_ht = launch;
		HTReal right_ht = au2d1*du1d1*tu2d2*au1d2*tu1d3;
		// the conformation below
		HTReal below_ht = left_ht * exgeom->transform( HT_tor_U3D1, lex[1] == 1 ? lex.dimsize(1) : lex[1]-1 ) * right_ht;
		core::Real below_score = build_and_score_res1_against_res2( pose, sfxn, below_ht, res1_index, res2_index, res1_geom );
		//pose.dump_pdb( "test_minima_tu3d1_below_" + utility::to_string( lex.index() ) + ".pdb" );
		if ( below_score < score ) return false;

		HTReal above_ht = left_ht * exgeom->transform( HT_tor_U3D1, lex[1] == lex.dimsize(1) ? 1 : lex[1]+1 ) * right_ht;
		core::Real above_score = build_and_score_res1_against_res2( pose, sfxn, above_ht, res1_index, res2_index, res1_geom );
		//pose.dump_pdb( "test_minima_tu3d1_above_" + utility::to_string( lex.index() ) + ".pdb" );
		if ( above_score < score ) return false;
	}

	// 2. au2d1 -- non-periodic
	if ( lex.dimsize(2) != 1 ) {
		HTReal left_ht = launch*tu3d1;
		HTReal right_ht = du1d1*tu2d2*au1d2*tu1d3;
		if ( lex[2] != 1 ) {
			HTReal below_ht = left_ht * exgeom->transform( HT_ang_U2D1, lex[2]-1 ) * right_ht;
			core::Real below_score = build_and_score_res1_against_res2( pose, sfxn, below_ht, res1_index, res2_index, res1_geom );
			//pose.dump_pdb( "test_minima_au2d1_below_" + utility::to_string( lex.index() ) + ".pdb" );
			if ( below_score < score ) return false;
		}
		if ( lex[2] != lex.dimsize(2) ) {
			HTReal above_ht = left_ht * exgeom->transform( HT_ang_U2D1, lex[2]+1 ) * right_ht;
			core::Real above_score = build_and_score_res1_against_res2( pose, sfxn, above_ht, res1_index, res2_index, res1_geom );
			//pose.dump_pdb( "test_minima_au2d1_above_" + utility::to_string( lex.index() ) + ".pdb" );
			if ( above_score < score ) return false;
		}
	}

	// 3. du1d1 -- non-periodic
	if ( lex.dimsize(3) != 1 ) {
		HTReal left_ht = launch*tu3d1*au2d1;
		HTReal right_ht = tu2d2*au1d2*tu1d3;
		if ( lex[3] != 1 ) {
			HTReal below_du1d1; below_du1d1.walk_along_z( exgeom->dis_U1D1_samples()[ lex[3]-1 ] );
			HTReal below_ht = left_ht * below_du1d1 * right_ht;
			core::Real below_score = build_and_score_res1_against_res2( pose, sfxn, below_ht, res1_index, res2_index, res1_geom );
			//pose.dump_pdb( "test_minima_du1d1_below_" + utility::to_string( lex.index() ) + ".pdb" );
			if ( below_score < score ) return false;
		}
		if ( lex[3] != lex.dimsize(3) ) {
			HTReal above_du1d1; above_du1d1.walk_along_z( exgeom->dis_U1D1_samples()[ lex[3]+1 ] );
			HTReal above_ht = left_ht * above_du1d1 * right_ht;
			core::Real above_score = build_and_score_res1_against_res2( pose, sfxn, above_ht, res1_index, res2_index, res1_geom );
			//pose.dump_pdb( "test_minima_du1d1_above_" + utility::to_string( lex.index() ) + ".pdb" );
			if ( above_score < score ) return false;
		}
	}

	/// 4. tu2d2 -- periodic
	if ( lex.dimsize(4) != 1 ) {
		HTReal left_ht = launch*tu3d1*au2d1*du1d1;
		HTReal right_ht = au1d2*tu1d3;
		// the conformation below
		HTReal below_ht = left_ht * exgeom->transform( HT_tor_U2D2, lex[4] == 1 ? lex.dimsize(4) : lex[4]-1 ) * right_ht;
		core::Real below_score = build_and_score_res1_against_res2( pose, sfxn, below_ht, res1_index, res2_index, res1_geom );
		//pose.dump_pdb( "test_minima_tu2d2_below_" + utility::to_string( lex.index() ) + ".pdb" );
		if ( below_score < score ) return false;

		HTReal above_ht = left_ht * exgeom->transform( HT_tor_U2D2, lex[4] == lex.dimsize(4) ? 1 : lex[4]+1 ) * right_ht;
		core::Real above_score = build_and_score_res1_against_res2( pose, sfxn, above_ht, res1_index, res2_index, res1_geom );
		//pose.dump_pdb( "test_minima_tu2d2_above_" + utility::to_string( lex.index() ) + ".pdb" );
		if ( above_score < score ) return false;
	}

	// 5. au1d2 -- non-periodic
	if ( lex.dimsize(5) != 1 ) {
		HTReal left_ht = launch*tu3d1*au2d1*du1d1*tu2d2;
		HTReal right_ht = tu1d3;
		if ( lex[5] != 1 ) {
			HTReal below_ht = left_ht * exgeom->transform( HT_ang_U1D2, lex[5]-1 ) * right_ht;
			core::Real below_score = build_and_score_res1_against_res2( pose, sfxn, below_ht, res1_index, res2_index, res1_geom );
			//pose.dump_pdb( "test_minima_au1d2_below_" + utility::to_string( lex.index() ) + ".pdb" );
			if ( below_score < score ) return false;
		}
		if ( lex[5] != lex.dimsize(5) ) {
			HTReal above_ht = left_ht * exgeom->transform( HT_ang_U1D2, lex[5]+1 ) * right_ht;
			core::Real above_score = build_and_score_res1_against_res2( pose, sfxn, above_ht, res1_index, res2_index, res1_geom );
			//pose.dump_pdb( "test_minima_au1d2_above_" + utility::to_string( lex.index() ) + ".pdb" );
			if ( above_score < score ) return false;
		}
	}

	/// 6. tu1d3 -- periodic
	if ( lex.dimsize(6) != 1 ) {
		HTReal left_ht = launch*tu3d1*au2d1*du1d1*tu2d2*au1d2;
		// the conformation below
		HTReal below_ht = left_ht * exgeom->transform( HT_tor_U1D3, lex[6] == 1 ? lex.dimsize(6) : lex[6]-1 );
		core::Real below_score = build_and_score_res1_against_res2( pose, sfxn, below_ht, res1_index, res2_index, res1_geom );
		//pose.dump_pdb( "test_minima_tu1d3_below_" + utility::to_string( lex.index() ) + ".pdb" );
		if ( below_score < score ) return false;

		HTReal above_ht = left_ht * exgeom->transform( HT_tor_U1D3, lex[6] == lex.dimsize(6) ? 1 : lex[6]+1 );
		core::Real above_score = build_and_score_res1_against_res2( pose, sfxn, above_ht, res1_index, res2_index, res1_geom );
		//pose.dump_pdb( "test_minima_tu1d3_above_" + utility::to_string( lex.index() ) + ".pdb" );
		if ( above_score < score ) return false;
	}

	/// OK -- we have a winner -- put this residue back onto the pose (and score it 'cause that's easier than not scoring it)
	build_and_score_res1_against_res2( pose, sfxn, htend, res1_index, res2_index, res1_geom );

	return true;
}

void
create_mingraph_for_focused_residue_pair(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::Size res1_index,
	core::Size res2_index,
	core::optimization::MinimizerMapOP & minimizer_map,
	core::scoring::MinimizationGraphOP & mingraph
)
{
	using namespace core;

	kinematics::MoveMap mm;
	mm.set_bb( res1_index, true );
	mm.set_bb( res2_index, true );
	minimizer_map = core::optimization::MinimizerMapOP( new optimization::MinimizerMap );
	minimizer_map->setup( pose, mm );

	mingraph = core::scoring::MinimizationGraphOP( new scoring::MinimizationGraph( pose.total_residue() ) );

	mingraph->add_edge( res1_index, res2_index );
	scoring::EnergyMap dummy;

	sfxn.setup_for_minimizing_for_node(
		* mingraph->get_minimization_node( res1_index ), pose.residue( res1_index ),
		* minimizer_map, pose, false, dummy );
	sfxn.setup_for_minimizing_for_node(
		* mingraph->get_minimization_node( res2_index ), pose.residue( res2_index ),
		* minimizer_map, pose, false, dummy );

	sfxn.setup_for_minimizing_sr2b_enmeths_for_minedge(
		pose.residue( res1_index ), pose.residue( res2_index ),
		* mingraph->find_minimization_edge( res1_index, res2_index ),
		* minimizer_map, pose, true, false,
		static_cast< scoring::EnergyEdge const * > ( 0 ), dummy );

}

core::Real
gradient_magnitude_for_conformation(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::Size res1_index,
	core::Size res2_index,
	core::optimization::MinimizerMapOP & minimizer_map,
	core::scoring::MinimizationGraphOP & mingraph
)
{
	minimizer_map->zero_torsion_vectors();

	core::scoring::MinimizationNode & minnode1( * mingraph->get_minimization_node( res1_index ) );
	core::scoring::MinimizationNode & minnode2( * mingraph->get_minimization_node( res2_index ) );

	minnode1.setup_for_minimizing( pose.residue( res1_index ), pose, sfxn, *minimizer_map );
	minnode2.setup_for_minimizing( pose.residue( res2_index ), pose, sfxn, *minimizer_map );

	minnode1.setup_for_derivatives( pose.residue( res1_index ), pose, sfxn );
	minnode2.setup_for_derivatives( pose.residue( res2_index ), pose, sfxn );

	core::scoring::MinimizationEdge & minedge = static_cast< core::scoring::MinimizationEdge & >
		( *mingraph->find_minimization_edge( res1_index, res2_index ) );

	minedge.setup_for_minimizing(
		pose.residue( res1_index ),
		pose.residue( res2_index ),
		pose, sfxn, *minimizer_map );

	minedge.setup_for_derivatives(
		pose.residue( res1_index ),
		pose.residue( res2_index ),
		pose, sfxn );

	core::scoring::eval_atom_derivatives_for_minedge(
		minedge, pose.residue( res1_index ), pose.residue( res2_index ),
		minnode1.res_min_data(), minnode2.res_min_data(),
		pose, sfxn.weights(),
		minimizer_map->atom_derivatives( res1_index ),
		minimizer_map->atom_derivatives( res2_index ));

	core::Vector deriv_sum( 0.0 );
	utility::vector1< core::scoring::DerivVectorPair > const & r1derivs = minimizer_map->atom_derivatives( res1_index );
	for ( core::Size ii = 1; ii <= r1derivs.size(); ++ii ) {
		deriv_sum += r1derivs[ii].f2();
	}

	return deriv_sum.length();

}

bool
determine_if_at_local_minimum_analytic(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	protocols::toolbox::match_enzdes_util::ExternalGeomSamplerCOP exgeom,
	HTReal const & launch,
	core::Size res1_index,
	core::Size res2_index,
	utility::LexicographicalIterator const & lex,
	protocols::toolbox::match_enzdes_util::LigandConformer const & res1_geom,
	core::optimization::MinimizerMapOP minimizer_map,
	core::scoring::MinimizationGraphOP mingraph,
	core::Real & score
)
{
	core::Real const gradient_magnitude_cutoff( 0.5 ); // I have no idea what to use here
	core::Real const score_function_magnitude_cutoff( -0.8 );

	using namespace protocols::toolbox::match_enzdes_util;

	HTReal tu3d1 = exgeom->transform( HT_tor_U3D1, lex[1] );
	HTReal au2d1 = exgeom->transform( HT_ang_U2D1, lex[2] );
	HTReal du1d1;  du1d1.walk_along_z( exgeom->dis_U1D1_samples()[lex[3]] );
	HTReal tu2d2 = exgeom->transform( HT_tor_U2D2, lex[4] );
	HTReal au1d2 = exgeom->transform( HT_ang_U1D2, lex[5] );
	HTReal tu1d3 = exgeom->transform( HT_tor_U1D3, lex[6] );

	HTReal htend = launch*tu3d1*au2d1*du1d1*tu2d2*au1d2*tu1d3;
	score = build_and_score_res1_against_res2( pose, sfxn, htend, res1_index, res2_index, res1_geom );

	if ( score > score_function_magnitude_cutoff ) return false;

	core::Real gradient_magnitude = gradient_magnitude_for_conformation(
		pose, sfxn, res1_index, res2_index, minimizer_map, mingraph );

	return gradient_magnitude < gradient_magnitude_cutoff;
}


struct
sort_minima_by_scores{
	bool operator() ( local_minimum const & a, local_minimum const & b ) { return a.score_ < b.score_; }
};


bool
same_coordinates(
	core::conformation::Residue const & r1,
	core::conformation::Residue const & r2,
	core::Real tolerance
)
{
	core::Real d2sum = 0;
	for ( core::Size ii = 1, iiend = r1.natoms(); ii <= iiend; ++ii ) {
		d2sum += r1.xyz( ii ).distance_squared( r2.xyz( ii ) );
	}
	return d2sum < tolerance;
}

std::list< local_minimum >
trim_redundant_local_minima(
	std::list< local_minimum > const & input_minima,
	core::Real tolerance
)
{
	utility::vector1< local_minimum > minima( input_minima.size() );
	std::copy( input_minima.begin(), input_minima.end(), minima.begin() );
	utility::vector1< int > good( minima.size(), 1 );
	std::list< local_minimum > output_minima;
	for ( core::Size ii = 1; ii <= minima.size(); ++ii ) {
		if ( !good[ii] ) continue;
		output_minima.push_back( minima[ii] );
		for ( core::Size jj = ii+1; jj <= minima.size(); ++jj ) {
			if ( ! good[jj] ) continue;
			if ( same_coordinates( *minima[ii].res1_conformation, *minima[jj].res1_conformation, tolerance ) ) {
				good[ jj ] = 0;
			}
		}
	}
	return output_minima;
}

void
output_local_minima(
	core::pose::Pose & pose,
	std::list< local_minimum > & local_minima,
	std::string const & acc3let,
	std::string const & don3let,
	core::Size res1_index
)
{
	core::id::AtomID_Map< core::Real > bfactors;

	sort_minima_by_scores sorter;
	local_minima.sort( sorter );

	local_minima = trim_redundant_local_minima( local_minima, 0.01 );

	core::pose::initialize_atomid_map( bfactors, pose, core::Real(0.0) );

	core::Size max_digits( core::Size( std::ceil( std::log10( local_minima.size() ) ) ) );
	core::Size count( 0 );
	for ( std::list< local_minimum >::const_iterator
			iter = local_minima.begin(), iter_end = local_minima.end();
			iter != iter_end; ++iter ) {
		core::Real bfactor = iter->score_ * -10;
		if ( bfactor > 100 ) bfactor = 100;
		for ( core::Size ii = 1; ii <= pose.residue( res1_index ).natoms(); ++ii ) {
			bfactors[ core::id::AtomID( ii, res1_index ) ] = bfactor;
		}
		pose.replace_residue( res1_index, *(iter->res1_conformation ), false );

		std::string out_fname = acc3let + "_" + don3let + "_minima_" + ObjexxFCL::format::I( max_digits, max_digits, ++count ) + ".pdb";
		std::ofstream out( out_fname.c_str() );
		core::io::pdb::dump_bfactor_pdb( pose, bfactors, out, "model_" + utility::to_string(count));
	}
}


void
sweep_params_from_match_constraint_file(
	utility::sql_database::sessionOP db_session,
	core::scoring::ScoreFunction const & sfxn,
	std::string const & cst_fname
)
{
	using namespace core::chemical;
	using namespace protocols::toolbox::match_enzdes_util;

	bool const output_bestpose_only( basic::options::option[ basic::options::OptionKeys::sweep_respair_energies::output_bestpose_only ] );
	bool const find_local_minima( basic::options::option[ basic::options::OptionKeys::sweep_respair_energies::find_local_minima ] );
	bool const no_database_output(
		basic::options::option[ basic::options::OptionKeys::sweep_respair_energies::output_bestpose_only ] ||
		basic::options::option[ basic::options::OptionKeys::sweep_respair_energies::find_local_minima ] );

	ResidueTypeSetCOP fa_restypeset( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	/// Read in the match constraint file info (list)
	MatchConstraintFileInfoList mcfil( fa_restypeset );
	utility::io::izstream infile( cst_fname );
	mcfil.read_data( infile );

	/// Construct the pose that represents the amino acids described in the match constraint file
	char res1n1 = oneletter_code_from_aa( aa_from_name( mcfil.mcfi(1)->allowed_res_name3s(1)[1] ));
	char res2n1 = oneletter_code_from_aa( aa_from_name( mcfil.mcfi(1)->allowed_res_name3s(2)[1] ));
	core::pose::Pose pose;
	// 2=res1, 8=res2
	core::Size const res1_index(2), res2_index(8);
	std::string pose_seq = std::string("G") + res1n1 + "GGGGG" + res2n1 + "G";
	std::cout << "creating pose from sequence: " << pose_seq << std::endl;
	core::pose::make_pose_from_sequence( pose, pose_seq, *fa_restypeset );
	for ( core::Size ii = 1; ii <= 9; ++ii ) {
		pose.set_phi( ii, -121 );
		pose.set_psi( ii,  136 );
		pose.set_omega(ii, 180 );
	}

	/// UGLY HACK -- make acceptor histadines HIS_Ds (to accept at NE2) and
	/// donor histadines HISs (to donate at HE2)
	{
		if ( res1n1 == 'H' ) {
			std::cout << "putting HIS_D as acceptor" << std::endl;
			core::chemical::ResidueType const & his_d = fa_restypeset->name_map( "HIS_D" );
			core::conformation::ResidueOP acc_his = core::conformation::ResidueFactory::create_residue(
				his_d, pose.residue( res1_index ), pose.conformation(), false );
			pose.replace_residue( res1_index, *acc_his, false );
		}
		if ( res2n1 == 'H' ) {
			std::cout << "putting HIS(_E) as donor" << std::endl;
			core::chemical::ResidueType const & his_e = fa_restypeset->name_map( "HIS" );
			core::conformation::ResidueOP don_his = core::conformation::ResidueFactory::create_residue(
				his_e, pose.residue( res2_index ), pose.conformation(), false );
			pose.replace_residue( res2_index, *don_his, false );
		}
	}

	boost::uuids::uuid respair_id;
	hbgeoms_to_measure geoms;
	if ( ! no_database_output ) {
	  /// Insert a new entry into the respairs table of the database
	  /// representing the parameters that have been specified in the input match
	  respair_id = write_new_respair_to_database( db_session, *mcfil.mcfi(1) );
	  geoms = create_fullcoverage_hbgeoms_to_measure( pose, *mcfil.mcfi(1), res1_index, res2_index );
	}

	/// Set DOFs
	std::map< std::string, utility::vector1< std::string > > const &
		alg_info( mcfil.mcfi(1)->algorithm_inputs());
	if ( alg_info.find( "match" ) != alg_info.end() ) {
		utility::vector1< std::string > const & chi_to_set( alg_info.find( "match" )->second );
		for ( core::Size ii = 1; ii <= chi_to_set.size(); ++ii ) {
			std::istringstream iiline( chi_to_set[ii] );
			std::string res, chi;
			core::Size which_res, which_chi;
			core::Real chival;
			iiline >> res;
			iiline >> which_res;
			iiline >> chi;
			iiline >> which_chi;
			iiline >> chival;
			pose.set_chi( which_chi, which_res == 1 ? res1_index : res2_index, chival );
			//std::cout << "setting chi: " << (which_res == 1 ? res1_index : res2_index) << " " << which_chi << " " << chival << std::endl;
		}
	}

	//pose.dump_pdb( "pose_after_setting_chi_from_mcfi.pdb" );

	utility::vector1< core::Size > r1ats( 3, 0 );
	r1ats[1] = mcfil.mcfi(1)->template_atom_inds(1,1,pose.residue_type(res1_index) )[1];
	r1ats[2] = mcfil.mcfi(1)->template_atom_inds(1,2,pose.residue_type(res1_index) )[1];
	r1ats[3] = mcfil.mcfi(1)->template_atom_inds(1,3,pose.residue_type(res1_index) )[1];

	utility::vector1< core::Size > r2ats( 3, 0 );
	r2ats[1] = mcfil.mcfi(1)->template_atom_inds(2,1,pose.residue_type(res2_index) )[1];
	r2ats[2] = mcfil.mcfi(1)->template_atom_inds(2,2,pose.residue_type(res2_index) )[1];
	r2ats[3] = mcfil.mcfi(1)->template_atom_inds(2,3,pose.residue_type(res2_index) )[1];

	LigandConformer res1_geom;
	res1_geom.initialize_from_residue(
		r1ats[1], r1ats[2], r1ats[3],
		r1ats[1], r1ats[2], r1ats[3],
		pose.residue(res1_index) );

	core::conformation::Residue const & res2( pose.residue(res2_index) );
	HTReal launch( res2.xyz( r2ats[3] ), res2.xyz( r2ats[2] ), res2.xyz( r2ats[1] ) );

	// prepare to enumerate all combinations of external geometries
	ExternalGeomSamplerOP exgeom = mcfil.mcfi(1)->create_exgs();
	exgeom->set_dis_D1D2( res1_geom.atom1_atom2_distance() );
	exgeom->set_dis_D2D3( res1_geom.atom2_atom3_distance() );
	exgeom->set_ang_D1D2D3( res1_geom.atom1_atom2_atom3_angle() );
	exgeom->precompute_transforms();

	utility::vector1< core::Size > nexgeoms( 6, 0 );
	nexgeoms[1] = exgeom->n_tor_U3D1_samples();
	nexgeoms[2] = exgeom->n_ang_U2D1_samples();
	nexgeoms[3] = exgeom->n_dis_U1D1_samples();
	nexgeoms[4] = exgeom->n_tor_U2D2_samples();
	nexgeoms[5] = exgeom->n_ang_U1D2_samples();
	nexgeoms[6] = exgeom->n_tor_U1D3_samples();
	utility::LexicographicalIterator lex( nexgeoms );

	sfxn( pose );
	core::pose::Pose best_pose = pose;
	core::Real best_score( 1234 );


	core::optimization::MinimizerMapOP minimizer_map;
	core::scoring::MinimizationGraphOP mingraph;
	if ( find_local_minima ) {
		create_mingraph_for_focused_residue_pair(
			pose, sfxn, res1_index, res2_index,
			minimizer_map,	mingraph );
	}


	core::Size count( 0 );

	std::string acc3let = mcfil.mcfi(1)->allowed_res_name3s(1)[1];
	std::string don3let = mcfil.mcfi(1)->allowed_res_name3s(2)[1];

	std::list< local_minimum > local_minima;

	for( lex.begin(); ! lex.at_end(); ++lex ) {
		++count;

		if ( find_local_minima ) {
			core::Real score(0);
			//bool at_local_minimum = determine_if_at_local_minimum(
			//pose, sfxn, exgeom, launch, res1_index, res2_index, lex, res1_geom, score );

			bool at_local_minimum = determine_if_at_local_minimum_analytic(
				pose, sfxn, exgeom, launch, res1_index, res2_index, lex, res1_geom,
				minimizer_map, mingraph, score );

			if ( at_local_minimum ) {
				// record the conformation and the score
				local_minimum lm;
				lm.res1_conformation = pose.residue( res1_index ).clone();
				lm.score_ = score;
				local_minima.push_back( lm );
			}
		} else {

			HTReal tu3d1 = exgeom->transform( HT_tor_U3D1, lex[1] );
			HTReal au2d1 = exgeom->transform( HT_ang_U2D1, lex[2] );
			HTReal du1d1;  du1d1.walk_along_z( exgeom->dis_U1D1_samples()[lex[3]] );
			HTReal tu2d2 = exgeom->transform( HT_tor_U2D2, lex[4] );
			HTReal au1d2 = exgeom->transform( HT_ang_U1D2, lex[5] );
			HTReal tu1d3 = exgeom->transform( HT_tor_U1D3, lex[6] );

			HTReal htend = launch*tu3d1*au2d1*du1d1*tu2d2*au1d2*tu1d3;
			for ( core::Size ii=1, iiend = pose.residue(2).natoms(); ii <= iiend; ++ii ) {
				pose.set_xyz( core::id::AtomID( ii, res1_index ), res1_geom.coordinate_in_global_frame( ii, htend ) );
			}

			if ( no_database_output ) {
				core::scoring::EnergyMap emap;
				sfxn.eval_ci_2b( pose.residue( res1_index ), pose.residue( res2_index ), pose, emap );
				sfxn.eval_cd_2b( pose.residue( res1_index ), pose.residue( res2_index ), pose, emap );
				core::Real this_score = emap.dot( sfxn.weights() );

				if ( count == 0 || this_score < best_score ) {
					best_score = this_score;
					best_pose = pose;
				}
			} else {
				write_scores_and_hbond_geoms_to_database(
					db_session, respair_id, sfxn, pose,
					geoms, res1_index, res2_index );
			}

			if ( count == lex.num_states_total() / 2 ) {
			  pose.dump_pdb( acc3let + "_" + don3let + "_" + utility::to_string(count) + ".pdb" );
			}
		}
	}

	if ( output_bestpose_only ) {
		best_pose.dump_pdb( acc3let + "_" + don3let + "_lowest_energy.pdb" );
	}
	if ( find_local_minima ) {
		output_local_minima( pose, local_minima, acc3let, don3let, res1_index );
	}

	//{
	//	utility::vector1< core::Size > const & template_atom_inds_1_1( mcfil.mcfi(1)->template_atom_inds(1,1,pose.residue_type(2) ));
	//	for ( core::Size ii = 1; ii <= template_atom_inds_1_1.size(); ++ii ) {
	//		std::cout << "Template atom1 for ASP: " << template_atom_inds_1_1[ ii ] << " " << pose.residue_type(2).atom_name( template_atom_inds_1_1[ii] ) << std::endl;
	//	}
	//}
	//{
	//	utility::vector1< core::Size > const & template_atom_inds_1_2( mcfil.mcfi(1)->template_atom_inds(1,2,pose.residue_type(2) ));
	//	for ( core::Size ii = 1; ii <= template_atom_inds_1_2.size(); ++ii ) {
	//		std::cout << "Template atom2 for ASP: " << template_atom_inds_1_2[ ii ] << " " << pose.residue_type(2).atom_name( template_atom_inds_1_2[ii] ) << std::endl;
	//	}
	//}
	//{
	//	utility::vector1< core::Size > const & template_atom_inds_1_3( mcfil.mcfi(1)->template_atom_inds(1,3,pose.residue_type(2) ));
	//	for ( core::Size ii = 1; ii <= template_atom_inds_1_3.size(); ++ii ) {
	//		std::cout << "Template atom3 for ASP: " << template_atom_inds_1_3[ ii ] << " " << pose.residue_type(2).atom_name( template_atom_inds_1_3[ii] ) << std::endl;
	//	}
	//}
	//
	//{
	//	utility::vector1< core::Size > const & template_atom_inds_2_1( mcfil.mcfi(1)->template_atom_inds(2,1,pose.residue_type(8) ));
	//	for ( core::Size ii = 1; ii <= template_atom_inds_2_1.size(); ++ii ) {
	//		std::cout << "Template atom1 for SER: " << template_atom_inds_2_1[ ii ] << " " << pose.residue_type(8).atom_name( template_atom_inds_2_1[ii] ) << std::endl;
	//	}
	//}
	//{
	//	utility::vector1< core::Size > const & template_atom_inds_2_2( mcfil.mcfi(1)->template_atom_inds(2,2,pose.residue_type(8) ));
	//	for ( core::Size ii = 1; ii <= template_atom_inds_2_2.size(); ++ii ) {
	//		std::cout << "Template atom2 for SER: " << template_atom_inds_2_2[ ii ] << " " << pose.residue_type(8).atom_name( template_atom_inds_2_2[ii] ) << std::endl;
	//	}
	//}
	//
	//{
	//	utility::vector1< core::Size > const & template_atom_inds_2_3( mcfil.mcfi(1)->template_atom_inds(2,3,pose.residue_type(8) ));
	//	for ( core::Size ii = 1; ii <= template_atom_inds_2_3.size(); ++ii ) {
	//		std::cout << "Template atom3 for SER: " << template_atom_inds_2_3[ ii ] << " " << pose.residue_type(8).atom_name( template_atom_inds_2_3[ii] ) << std::endl;
	//	}
	//}

}


int main( int argc, char * argv [] )
{
	try{

	using namespace core;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::id;
	using namespace core::io::pdb;
	using namespace core::graph;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::database;

	NEW_OPT( sweep_respair_energies::output_bestpose_only, "write out the best geometries ", false );
	NEW_OPT( sweep_respair_energies::find_local_minima, "step to the left and right from each conformation along all 6 DOFs and output multi-model PDBs for the structures at local minima", false );

	devel::init( argc, argv );

	bool const no_database_output(
		basic::options::option[ basic::options::OptionKeys::sweep_respair_energies::output_bestpose_only ] ||
		basic::options::option[ basic::options::OptionKeys::sweep_respair_energies::find_local_minima ] );

	utility::sql_database::sessionOP db_session;
	if ( ! no_database_output ) {
		db_session = get_db_session( "test.db3" );

		basic::database::set_cache_size( db_session, 1024 );

		create_schema( db_session );
		//populate_database_with_dummy_data( db_session );

		db_session->begin();
		write_score_types_to_database( db_session );
		write_hbgeom_types_to_database( db_session );
	}

	/// Read score fumction from the command line.
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	using namespace core;
	scoring::methods::EnergyMethodOptionsOP emopts( new scoring::methods::EnergyMethodOptions( sfxn->energy_method_options() ) );
	emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
	emopts->hbond_options().use_hb_env_dep( false );
	sfxn->set_energy_method_options( *emopts );

	utility::vector1< utility::file::FileName > match_files( option[ in::file::l ] );

	utility::vector1< std::string > fnames;
	for ( core::Size ii = 1; ii <= match_files.size(); ++ii ) {
		std::string iimatchfilelist = match_files[ii].name();
		std::ifstream list_file( iimatchfilelist.c_str() );
		while ( list_file ) {
			std::string fname;
			list_file >> fname;
			if ( ! list_file.bad() && fname[0] != '#' && fname.size() != 0 ) {
				fnames.push_back( fname );
			}
		}
	}
	for ( core::Size ii = 1; ii <= fnames.size(); ++ii ) {
		sweep_params_from_match_constraint_file( db_session, *sfxn, fnames[ ii ] );
	}

	if ( ! no_database_output ) {
		db_session->commit();
	}

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
}

