// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/RotamerFeatures.cc
/// @brief  report rotamer library features to the database
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/RotamerFeatures.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
// Project Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>


//Basic Headers
#include <basic/basic.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/Tracer.hh>

// External Headers
#include <cppdb/frontend.h>

namespace protocols {
namespace features {

using std::string;
using std::endl;
using std::copy;
using core::Size;
using core::Real;
using core::Vector;
using core::chemical::AA;
using core::chemical::FA_STANDARD;
using core::chemical::num_canonical_aas;
using core::conformation::Residue;
using core::pack::dunbrack::RotamerLibrary;
using core::pack::dunbrack::RotamerLibraryScratchSpace;
using core::pack::dunbrack::PackedDunbrackRotamer;
using core::pack::dunbrack::ChiVector;
using core::pack::dunbrack::subtract_chi_angles;
using core::pack::dunbrack::ONE;
using core::pack::dunbrack::TWO;
using core::pack::dunbrack::THREE;
using core::pack::dunbrack::FOUR;
using core::pose::Pose;
using basic::periodic_range;
using basic::database::safely_prepare_statement;
using basic::database::safely_write_to_database;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using cppdb::result;

static THREAD_LOCAL basic::Tracer TR( "protocols.features.RotamerFeatures" );

string
RotamerFeatures::type_name() const {
	return "RotamerFeatures";
}

void
RotamerFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_residue_rotamers_table_schema(db_session);
}

void
RotamerFeatures::write_residue_rotamers_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column residue_number("residue_number", DbDataTypeOP( new DbInteger() ), false);
	Column rotamer_bin("rotamer_bin", DbDataTypeOP( new DbInteger() ), false);
	Column nchi("nchi", DbDataTypeOP( new DbInteger() ), false);
	Column semi_rotameric("semi_rotameric", DbDataTypeOP( new DbInteger() ), false);
	Column chi1_mean("chi1_mean", DbDataTypeOP( new DbReal() ), true);
	Column chi2_mean("chi2_mean", DbDataTypeOP( new DbReal() ), true);
	Column chi3_mean("chi3_mean", DbDataTypeOP( new DbReal() ), true);
	Column chi4_mean("chi4_mean", DbDataTypeOP( new DbReal() ), true);
	Column chi1_standard_deviation("chi1_standard_deviation", DbDataTypeOP( new DbReal() ), true);
	Column chi2_standard_deviation("chi2_standard_deviation", DbDataTypeOP( new DbReal() ), true);
	Column chi3_standard_deviation("chi3_standard_deviation", DbDataTypeOP( new DbReal() ), true);
	Column chi4_standard_deviation("chi4_standard_deviation", DbDataTypeOP( new DbReal() ), true);
	Column chi1_deviation("chi1_deviation", DbDataTypeOP( new DbReal() ), true);
	Column chi2_deviation("chi2_deviation", DbDataTypeOP( new DbReal() ), true);
	Column chi3_deviation("chi3_deviation", DbDataTypeOP( new DbReal() ), true);
	Column chi4_deviation("chi4_deviation", DbDataTypeOP( new DbReal() ), true);
	Column rotamer_bin_probability("rotamer_bin_probability", DbDataTypeOP( new DbReal() ), false);

	vector1<Column> primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(residue_number);
	PrimaryKey primary_key(primary_key_columns);

	vector1<Column> residues_foreign_key_columns;
	residues_foreign_key_columns.push_back(struct_id);
	residues_foreign_key_columns.push_back(residue_number);

	vector1<string> residues_foreign_key_reference_columns;
	residues_foreign_key_reference_columns.push_back("struct_id");
	residues_foreign_key_reference_columns.push_back("resNum");
	ForeignKey residues_foreign_key(
		residues_foreign_key_columns,
		"residues",
		residues_foreign_key_reference_columns,
		true);


	Schema residue_rotamers("residue_rotamers", primary_key);
	residue_rotamers.add_column(struct_id);
	residue_rotamers.add_column(residue_number);
	residue_rotamers.add_column(rotamer_bin);
	residue_rotamers.add_column(nchi);
	residue_rotamers.add_column(semi_rotameric);
	residue_rotamers.add_column(chi1_mean);
	residue_rotamers.add_column(chi2_mean);
	residue_rotamers.add_column(chi3_mean);
	residue_rotamers.add_column(chi4_mean);
	residue_rotamers.add_column(chi1_standard_deviation);
	residue_rotamers.add_column(chi2_standard_deviation);
	residue_rotamers.add_column(chi3_standard_deviation);
	residue_rotamers.add_column(chi4_standard_deviation);
	residue_rotamers.add_column(chi1_deviation);
	residue_rotamers.add_column(chi2_deviation);
	residue_rotamers.add_column(chi3_deviation);
	residue_rotamers.add_column(chi4_deviation);
	residue_rotamers.add_column(rotamer_bin_probability);
	residue_rotamers.add_foreign_key(residues_foreign_key);
	residue_rotamers.write(db_session);
}


vector1<string>
RotamerFeatures::features_reporter_dependencies() const {
	vector1<string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
RotamerFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	string library_name;
	vector1< AA > rotameric_amino_acids;
	vector1< Size > rotameric_n_chi;
	vector1< Size > rotameric_n_bb;

	vector1< AA > sraa;
	vector1< Size > srnchi;
	vector1< Size > srnbb;
	vector1< bool > scind;
	vector1< bool > sampind;
	vector1< bool > sym;
	vector1< Real > astr;

	if ( option[ corrections::score::dun10 ] ) {
		library_name = "dun10";
		RotamerLibrary::initialize_dun10_aa_parameters(
			rotameric_amino_acids, rotameric_n_chi, rotameric_n_bb,
			sraa, srnchi, srnbb, scind, sampind, sym, astr );

	} else {
		library_name = "dun02";
		RotamerLibrary::initialize_dun02_aa_parameters(
			rotameric_amino_acids, rotameric_n_chi, rotameric_n_bb);
	}

	string insert_sql(
		"INSERT INTO residue_rotamers ("
		"\tstruct_id,"
		"\tresidue_number,"
		"\trotamer_bin,"
		"\tnchi,"
		"\tsemi_rotameric,"
		"\tchi1_mean,"
		"\tchi2_mean,"
		"\tchi3_mean,"
		"\tchi4_mean,"
		"\tchi1_standard_deviation,"
		"\tchi2_standard_deviation,"
		"\tchi3_standard_deviation,"
		"\tchi4_standard_deviation,"
		"\tchi1_deviation,"
		"\tchi2_deviation,"
		"\tchi3_deviation,"
		"\tchi4_deviation,"
		"\trotamer_bin_probability) VALUES"
		"\t(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);");
	statement insert_stmt(safely_prepare_statement(insert_sql, db_session));

	RotamerLibraryScratchSpace scratch;

	for (
			Size residue_number=1;
			residue_number <= pose.size();
			++residue_number ) {

		// only report features for the specified subset of residues
		if ( !check_relevant_residues(relevant_residues, residue_number) ) continue;

		Residue const & residue(pose.residue(residue_number));
		if ( residue.type().residue_type_set()->category() != core::chemical::FULL_ATOM_t ) {
			utility_exit_with_message("Attempting to report protein rotamer features for a pose that is not full-atom.");
		}

		if ( ! residue.is_protein() || residue.type().aa() > num_canonical_aas ) {
			TR.Warning << "Unable to report protein rotamer features for residue " << residue_number << " because it is not a canonical amino acid: '" << residue.type().name() << "'" << endl;
			continue;
		}

		if ( residue.type().rotamer_library_specification() && residue.type().rotamer_library_specification()->keyname() == "NCAA" ) {
			TR.Warning << "Currently the RotamerFeatures only supports canonical amino acids, but this residue type is '" << residue.type().name() << "'" << endl;
			continue;
		}


		Size rotamer_bin(0);
		bool recognized_residue_type;

		ChiVector const & chis(residue.chi());

		bool semi_rotameric(true);
		Size nchi = 0;
		Size n_bb = 0;
		for ( Size ii=1; ii <= rotameric_amino_acids.size(); ++ii ) {
			if ( rotameric_amino_acids[ii] == residue.aa() ) {
				semi_rotameric = false;
				nchi = rotameric_n_chi[ ii ];
				n_bb = rotameric_n_bb[ ii ];
				break;
			}
		}

		if ( semi_rotameric ) {
			for ( Size ii=1; ii <= sraa.size(); ++ii ) {
				if ( sraa[ii] == residue.aa() ) {
					nchi = srnchi[ii];
					n_bb = srnbb[ ii ];
					break;
				}
			}
		}


		if ( nchi == 0 ) {
			continue;
		} else if ( nchi == ONE ) {
			if ( n_bb == ONE ) {
				recognized_residue_type = RotamerInitializer<ONE, ONE>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == TWO ) {
				recognized_residue_type = RotamerInitializer<ONE, TWO>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == THREE ) {
				recognized_residue_type = RotamerInitializer<ONE, THREE>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == FOUR ) {
				recognized_residue_type = RotamerInitializer<ONE, FOUR>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else {
				continue;
			}
		} else if ( nchi == TWO ) {
			if ( n_bb == ONE ) {
				recognized_residue_type = RotamerInitializer<TWO, ONE>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == TWO ) {
				recognized_residue_type = RotamerInitializer<TWO, TWO>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == THREE ) {
				recognized_residue_type = RotamerInitializer<TWO, THREE>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == FOUR ) {
				recognized_residue_type = RotamerInitializer<TWO, FOUR>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else {
				continue;
			}
		} else if ( nchi == THREE ) {
			if ( n_bb == ONE ) {
				recognized_residue_type = RotamerInitializer<THREE, ONE>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == TWO ) {
				recognized_residue_type = RotamerInitializer<THREE, TWO>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == THREE ) {
				recognized_residue_type = RotamerInitializer<THREE, THREE>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == FOUR ) {
				recognized_residue_type = RotamerInitializer<THREE, FOUR>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else {
				continue;
			}
		} else if ( nchi == FOUR ) {
			if ( n_bb == ONE ) {
				recognized_residue_type = RotamerInitializer<FOUR, ONE>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == TWO ) {
				recognized_residue_type = RotamerInitializer<FOUR, TWO>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == THREE ) {
				recognized_residue_type = RotamerInitializer<FOUR, THREE>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else if ( n_bb == FOUR ) {
				recognized_residue_type = RotamerInitializer<FOUR, FOUR>::initialize_rotamer( residue, scratch, rotamer_bin);
			} else {
				continue;
			}
		} else {
			continue;
		}

		if ( !recognized_residue_type ) {
			TR.Warning << "Unable to report protein rotamer features for residue " << residue_number << " because there is no rotamer library defined for this residue type: '" << residue.type().name() << "'" << endl;
			continue;
		}


		vector1< Real > chi_deviations(0, 4);

		for ( Size chi = 1; chi <= nchi; ++chi ) {
			if ( library_name.compare("dun02") == 0 ) {
				chi_deviations.push_back(
					subtract_chi_angles(
					chis[chi], scratch.chimean()[chi], residue.aa(), chi));
			} else {
				chi_deviations.push_back(
					periodic_range(chis[chi] - scratch.chimean()[chi], 360));
			}
		}

		insert_stmt.bind(1, struct_id);
		insert_stmt.bind(2, residue_number);
		insert_stmt.bind(3, rotamer_bin);
		insert_stmt.bind(4, nchi);
		insert_stmt.bind(5, semi_rotameric);

		for ( Size chi = 1; chi <= 4; ++chi ) {
			if ( chi <= nchi ) {
				insert_stmt.bind(5+chi, scratch.chimean()[chi]);
				insert_stmt.bind(9+chi, scratch.chisd()[chi]);
				insert_stmt.bind(13+chi, chi_deviations[chi]);
			} else {
				insert_stmt.bind_null(5+chi);
				insert_stmt.bind_null(9+chi);
				insert_stmt.bind_null(13+chi);
			}
		}

		insert_stmt.bind(18, scratch.rotprob());
		safely_write_to_database(insert_stmt);
	}
	return 0;
}

void
RotamerFeatures::delete_record(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session
){

	statement conf_stmt(basic::database::safely_prepare_statement("DELETE FROM residue_rotamers WHERE struct_id = ?;\n",db_session));
	conf_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(conf_stmt);
}

} // namespace
} // namespace
