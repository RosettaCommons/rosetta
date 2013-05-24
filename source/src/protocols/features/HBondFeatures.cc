// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/HBondFeatures.cc
/// @brief  report HBond geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/HBondFeatures.hh>

//External

// Package Headers
#include <core/scoring/hbonds/HBondSet.hh>

// Project Headers
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondTypeManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/types.hh>
#include <protocols/moves/DataMap.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>


// Utility Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/assert.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector0.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <cmath>
#include <utility/excn/Exceptions.hh>
#include <algorithm>


//Auto Headers


namespace protocols{
namespace features{

using std::string;
using std::endl;
using std::sort;
using std::sqrt;
using std::stringstream;
using core::Size;
using core::Real;
using core::Vector;
using core::chemical::AtomIndices;
using core::chemical::Hybridization;
using core::conformation::Residue;
using core::id::AtomID_Map;
using core::id::AtomID;
using core::pose::Pose;
using core::pose::initialize_atomid_map;
using core::graph::EdgeListConstIterator;
using core::graph::Graph;
using core::scoring::calc_per_atom_sasa;
using core::scoring::EnergyMap;
using core::scoring::getScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunction;
using core::scoring::ScoringManager;
using core::scoring::ScoreTypes;
using core::scoring::TenANeighborGraph;
using core::scoring::etable::AnalyticEtableEnergy;
using core::scoring::etable::TableLookupEtableEnergy;
using core::scoring::hbonds::HBDonChemType;
using core::scoring::hbonds::HBAccChemType;
using core::scoring::hbonds::HBondOptions;
using core::scoring::hbonds::HBondTypeManager;
using core::scoring::hbonds::get_hb_don_chem_type;
using core::scoring::hbonds::get_hb_acc_chem_type;
using core::scoring::hbonds::make_hbBasetoAcc_unitvector;
using core::scoring::hbonds::hb_eval_type_weight;
using core::scoring::hbonds::HBond;
using core::scoring::hbonds::HBondCOP;
using core::scoring::hbonds::HBondSet;
using core::scoring::fa_atr;
using core::scoring::fa_rep;
using core::scoring::fa_sol;
using numeric::xyzVector;
using numeric::dihedral_radians;
using protocols::filters::Filters_map;
using protocols::moves::DataMap;
using protocols::moves::Movers_map;
using utility::sql_database::sessionOP;
using utility::tag::TagPtr;
using cppdb::statement;
using utility::vector1;
using basic::Tracer;
using basic::database::insert_or_ignore;

static Tracer TR("protocols.features.HBondFeatures");

HBondFeatures::HBondFeatures() :
	scfxn_(getScoreFunction())
{}

HBondFeatures::HBondFeatures(
	ScoreFunctionOP scfxn) :
	scfxn_(scfxn),
	definition_type_(hbdef_ENERGY),
	definition_threshold_(0)
{}

HBondFeatures::HBondFeatures(HBondFeatures const & src) :
	FeaturesReporter(),
	scfxn_(src.scfxn_),
	definition_type_(src.definition_type_),
	definition_threshold_(src.definition_threshold_)
{}

HBondFeatures::~HBondFeatures() {}

string
HBondFeatures::type_name() const { return "HBondFeatures"; }

void
HBondFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_hbond_chem_types_table_schema(db_session);
	write_hbond_sites_table_schema(db_session);
	write_hbond_sites_pdb_table_schema(db_session);
	write_hbond_site_environment_table_schema(db_session);
	write_hbond_site_atoms_table_schema(db_session);
	write_hbonds_table_schema(db_session);
	write_hbond_lennard_jones_table_schema(db_session);
	write_hbond_geom_coords_table_schema(db_session);
	write_hbond_dehydrons_table_schema(db_session);
}

void
HBondFeatures::write_hbond_chem_types_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;
	using boost::assign::list_of;

	Column chem_type("chem_type", new DbText(255));
	Column label("label", new DbText(255));

	Columns primary_key_columns;
	primary_key_columns.push_back(chem_type);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("hbond_chem_types", primary_key);
	table.add_column(label);
	table.write(db_session);

	//insert static values
	string const t("hbond_chem_types");
	std::vector< string > c;
	c.push_back("chem_type");
	c.push_back("label");
 	insert_or_ignore(t, c, list_of("'hbacc_NONE'")("'aNONE'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_PBA'")("'aPBA: bb'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_CXA'")("'aCXA: n,q'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_CXL'")("'aCXL: d,e'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_IMD'")("'aIMD: h'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_IME'")("'aIME: h'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_AHX'")("'aAHX: t'"), db_session);
 	insert_or_ignore(t, c, list_of("'hbacc_HXL'")("'aHXL: s,t'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_PCA_DNA'")("'aPCA_DNA: O{1,2}P'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_PES_DNA'")("'aPES_DNA: O{3,5}*'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_RRI_DNA'")("'aRRI_DNA: O4*'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_PCA_RNA'")("'aPCA_RNA: O{1,2}P'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_PES_RNA'")("'aPES_RNA: O{3,5}*'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_RRI_RNA'")("'aRRI_RNA: O4*'"), db_session);
 	insert_or_ignore(t, c, list_of("'hbacc_H2O'")("'aH2O'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_GENERIC_SP2BB'")("'aGEN: sp2 bb'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_GENERIC_SP2SC'")("'aGEN: sp2 sc'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_GENERIC_SP3BB'")("'aGEN: sp3 bb'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_GENERIC_SP3SC'")("'aGEN: sp3 sc'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_GENERIC_RINGBB'")("'aGEN: ring bb'"), db_session);
	insert_or_ignore(t, c, list_of("'hbacc_GENERIC_RINGSC'")("'aGEN: ring sc'"), db_session);

	insert_or_ignore(t, c, list_of("'hbdon_NONE'")("'dNONE'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_PBA'")("'dPBA: bb'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_CXA'")("'dCXA: n,q'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_IMD'")("'dIMD: h'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_IME'")("'dIME: h'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_IND'")("'dIND: w'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_AMO'")("'dAMO: k'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_GDE'")("'dGDE: r'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_GDH'")("'dGDH: r'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_AHX'")("'dAHX: s,t'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_HXL'")("'dHXL: y'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_H2O'")("'dH2O'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_GENERIC_BB'")("'dGEN: bb'"), db_session);
	insert_or_ignore(t, c, list_of("'hbdon_GENERIC_SC'")("'dGEN: sc'"), db_session);
}

void
HBondFeatures::write_hbond_sites_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column site_id("site_id", new DbInteger());
	Column resNum("resNum", new DbInteger());
	Column atmNum("atmNum", new DbInteger());
	Column is_donor("is_donor", new DbInteger());
	Column chain("chain", new DbInteger());
	Column resType("resType", new DbText());
	Column atmType("atmType", new DbText());
	Column HBChemType("HBChemType", new DbText(255));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(site_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(resNum);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(HBChemType);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("chem_type");
	ForeignKey foreign_key2(foreign_key_columns2, "hbond_chem_types", reference_columns2, true);


	Schema table("hbond_sites", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);

	table.add_column(atmNum);
	table.add_column(is_donor);
	table.add_column(chain);
	table.add_column(resType);
	table.add_column(atmType);


	table.write(db_session);
}

void
HBondFeatures::write_hbond_sites_pdb_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column site_id("site_id", new DbInteger());
	Column chain("chain", new DbText(1));
	Column resNum("resNum", new DbInteger());
	Column iCode("iCode", new DbText(1));
	Column heavy_atom_temperature("heavy_atom_temperature", new DbReal());
	Column heavy_atom_occupancy("heavy_atom_occupancy", new DbReal());


	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(site_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(site_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("site_id");
	ForeignKey foreign_key(foreign_key_columns, "hbond_sites", reference_columns, true);

	Schema table("hbond_sites_pdb", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(chain);
	table.add_column(resNum);
	table.add_column(iCode);
	table.add_column(heavy_atom_temperature);
	table.add_column(heavy_atom_occupancy);

	table.write(db_session);
}

void
HBondFeatures::write_hbond_site_environment_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column site_id("site_id", new DbInteger());
	Column sasa_r100("sasa_r100", new DbReal());
	Column sasa_r140("sasa_r140", new DbReal());
	Column sasa_r200("sasa_r200", new DbReal());
	Column hbond_energy("hbond_energy", new DbReal());
	Column num_hbonds("num_hbonds", new DbInteger());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(site_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(site_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("site_id");
	ForeignKey foreign_key(foreign_key_columns, "hbond_sites", reference_columns, true);

	Schema table("hbond_site_environment", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(sasa_r100);
	table.add_column(sasa_r140);
	table.add_column(sasa_r200);
	table.add_column(hbond_energy);
	table.add_column(num_hbonds);

	table.write(db_session);
}

void
HBondFeatures::write_hbond_site_atoms_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column site_id("site_id", new DbInteger());
	Column atm_x("atm_x", new DbReal());
	Column atm_y("atm_y", new DbReal());
	Column atm_z("atm_z", new DbReal());
	Column base_x("base_x", new DbReal());
	Column base_y("base_y", new DbReal());
	Column base_z("base_z", new DbReal());
	Column bbase_x("bbase_x", new DbReal());
	Column bbase_y("bbase_y", new DbReal());
	Column bbase_z("bbase_z", new DbReal());
	Column base2_x("base2_x", new DbReal());
	Column base2_y("base2_y", new DbReal());
	Column base2_z("base2_z", new DbReal());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(site_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(site_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("site_id");
	ForeignKey foreign_key(foreign_key_columns, "hbond_sites", reference_columns, true);

	Schema table("hbond_site_atoms", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(atm_x);
	table.add_column(atm_y);
	table.add_column(atm_z);
	table.add_column(base_x);
	table.add_column(base_y);
	table.add_column(base_z);
	table.add_column(bbase_x);
	table.add_column(bbase_y);
	table.add_column(bbase_z);
	table.add_column(base2_x);
	table.add_column(base2_y);
	table.add_column(base2_z);

	table.write(db_session);
}

void
HBondFeatures::write_hbonds_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column hbond_id("hbond_id", new DbInteger());
	Column don_id("don_id", new DbInteger());
	Column acc_id("acc_id", new DbInteger());
	Column HBEvalType("HBEvalType", new DbInteger());
	Column energy("energy", new DbReal());
	Column envWeight("envWeight", new DbReal());
	Column score_weight("score_weight", new DbReal());
	Column donRank("donRank", new DbInteger());
	Column accRank("accRank", new DbInteger());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(hbond_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(don_id);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("site_id");
	ForeignKey foreign_key1(foreign_key_columns1, "hbond_sites", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(acc_id);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("site_id");
	ForeignKey foreign_key2(foreign_key_columns2, "hbond_sites", reference_columns2, true);

	Schema table("hbonds", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_column(don_id);
	table.add_column(acc_id);
	table.add_column(HBEvalType);
	table.add_column(energy);
	table.add_column(envWeight);
	table.add_column(score_weight);
	table.add_column(donRank);
	table.add_column(accRank);

	table.write(db_session);
}

void
HBondFeatures::write_hbond_lennard_jones_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column hbond_id("hbond_id", new DbInteger());
	Column don_acc_atrE("don_acc_atrE", new DbReal());
	Column don_acc_repE("don_acc_repE", new DbReal());
	Column don_acc_solv("don_acc_solv", new DbReal());
	Column don_acc_base_atrE("don_acc_base_atrE", new DbReal());
	Column don_acc_base_repE("don_acc_base_repE", new DbReal());
	Column don_acc_base_solv("don_acc_base_solv", new DbReal());
	Column h_acc_atrE("h_acc_atrE", new DbReal());
	Column h_acc_repE("h_acc_repE", new DbReal());
	Column h_acc_solv("h_acc_solv", new DbReal());
	Column h_acc_base_atrE("h_acc_base_atrE", new DbReal());
	Column h_acc_base_repE("h_acc_base_repE", new DbReal());
	Column h_acc_base_solv("h_acc_base_solv", new DbReal());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(hbond_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(hbond_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("hbond_id");
	ForeignKey foreign_key(foreign_key_columns, "hbonds", reference_columns, true);

	Schema table("hbond_lennard_jones", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(don_acc_atrE);
	table.add_column(don_acc_repE);
	table.add_column(don_acc_solv);
	table.add_column(don_acc_base_atrE);
	table.add_column(don_acc_base_repE);
	table.add_column(don_acc_base_solv);
	table.add_column(h_acc_atrE);
	table.add_column(h_acc_repE);
	table.add_column(h_acc_solv);
	table.add_column(h_acc_base_atrE);
	table.add_column(h_acc_base_repE);
	table.add_column(h_acc_base_solv);

	table.write(db_session);
}

void
HBondFeatures::write_hbond_geom_coords_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column hbond_id("hbond_id", new DbInteger());
	Column AHdist("AHdist", new DbReal());
	Column cosBAH("cosBAH", new DbReal());
	Column cosAHD("cosAHD", new DbReal());
	Column chi("chi", new DbReal());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(hbond_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(hbond_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("hbond_id");
	ForeignKey foreign_key(foreign_key_columns, "hbonds", reference_columns, true);

	Schema table("hbond_geom_coords", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(AHdist);
	table.add_column(cosBAH);
	table.add_column(cosAHD);
	table.add_column(chi);

	table.write(db_session);
}

void
HBondFeatures::write_hbond_dehydrons_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column hbond_id("hbond_id", new DbInteger());
	Column wrapping_count("wrapping_count", new DbInteger());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(hbond_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(hbond_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("hbond_id");
	ForeignKey foreign_key(foreign_key_columns, "hbonds", reference_columns, true);

	Schema table("hbond_dehydrons", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(wrapping_count);

	table.write(db_session);
}

utility::vector1<std::string>
HBondFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

void
HBondFeatures::parse_my_tag(
	TagPtr const tag,
	DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	if(tag->hasOption("scorefxn")){
		string const scorefxn_name(tag->getOption<string>("scorefxn"));
		scfxn_ = data.get<ScoreFunction*>("scorefxns", scorefxn_name);
	} else {
		stringstream error_msg;
		error_msg
			<< "The " << type_name() << " reporter requires a 'scorefxn' tag:" << endl
			<< endl
			<< "    <feature name=" << type_name() <<" scorefxn=(name_of_score_function) />" << endl;
		throw utility::excn::EXCN_RosettaScriptsOption(error_msg.str());
	}

	string const definition_type(
		tag->getOption<string>("definition_type", "energy"));
	if(definition_type == "energy"){
		definition_type_ = hbdef_ENERGY;
	} else if(definition_type == "AHdist"){
		definition_type_ = hbdef_AHDIST;
	} else {
		stringstream error_msg;
		error_msg
			<< "The hbond definition type '" << definition_type << "' is not recognized." << endl
			<< "Available hbond definition types are:" << endl
			<< "	'energy' => A polar-polar contact is an hbond when energy is below the definition_threshold." << endl
			<< "  'AHdist' => A polar-polar contact is an hbond when the Acceptor-Hydrogen distance is less than the definition_threshold." << endl;
		throw utility::excn::EXCN_RosettaScriptsOption(error_msg.str());
	}

	definition_threshold_ =
		tag->getOption<Real>("definition_threshold", 0);

}

Size
HBondFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	HBondSet hbond_set;

	// assert pose.update_residue_neighbors() has been called:
	runtime_assert(
		 !pose.conformation().structure_moved() &&
		 pose.energies().residue_neighbors_updated());


	if(definition_type_ == hbdef_ENERGY){
		hbond_set.setup_for_residue_pair_energies( pose, false, false );
	} else if(definition_type_ == hbdef_AHDIST){
		fill_hbond_set_by_AHdist_threshold(pose, definition_threshold_, hbond_set);
	} else {
		utility_exit_with_message("Unrecognized hbond definition type.");
	}

	TR << "Number of hydrogen bonds found: " << hbond_set.nhbonds() << endl;

	Real const probe_radius_s(1.0), probe_radius_m(1.4), probe_radius_l(2.0);
	AtomID_Map< Real > atom_sasa_s, atom_sasa_m, atom_sasa_l;
	vector1< Real > residue_sasa_s, residue_sasa_m, residue_sasa_l;

	if(pose.is_fullatom()){
		calc_per_atom_sasa(pose, atom_sasa_s, residue_sasa_s, probe_radius_s);
		calc_per_atom_sasa(pose, atom_sasa_m, residue_sasa_m, probe_radius_m);
		calc_per_atom_sasa(pose, atom_sasa_l, residue_sasa_l, probe_radius_l);
	}

	AtomID_Map< vector1<HBondCOP> > site_partners;
	initialize_atomid_map(site_partners, pose);

	AtomID_Map<Real>site_hbond_energies;
	initialize_atomid_map(site_hbond_energies, pose);

	for (Size i = 1; i<= hbond_set.nhbonds(); i++) {
		HBondCOP hbond(hbond_set.hbond(i));
		if(!relevant_residues[hbond->don_res()] ||
			!relevant_residues[hbond->acc_res()]) continue;

		site_partners(hbond->don_res(),hbond->don_hatm()).push_back(hbond);
		site_hbond_energies(hbond->don_res(),hbond->don_hatm()) += hbond->energy()/2;

		site_partners(hbond->acc_res(),hbond->acc_atm()).push_back(hbond);
		site_hbond_energies(hbond->acc_res(),hbond->acc_atm()) += hbond->energy()/2;

	}

	Size site_id(0);
	AtomID_Map< Size > site_ids;
	core::pose::initialize_atomid_map(site_ids, pose);
	for( Size resNum =1; resNum <= pose.n_residue(); ++resNum ){
		if(!relevant_residues[resNum]) continue;

		Residue const & res(pose.residue(resNum));
		// donor sites
		for ( AtomIndices::const_iterator
			atmNum  = res.Hpos_polar().begin(),
			atmNume = res.Hpos_polar().end(); atmNum != atmNume; ++atmNum ) {

			site_id++;
			insert_site_row(pose, struct_id, site_id, resNum, *atmNum, true /*is donor*/, db_session);
			site_ids(resNum,*atmNum) = site_id;
			insert_site_pdb_row(pose, resNum, *atmNum, res.atom_base(*atmNum), struct_id, site_id, db_session);
			insert_site_environment_row(pose, resNum, *atmNum, struct_id, site_id, atom_sasa_s, atom_sasa_m, atom_sasa_l, site_partners, site_hbond_energies, db_session);
			insert_site_atoms_row(pose, resNum, *atmNum, struct_id, site_id, db_session);

			vector1< HBondCOP >::iterator partners_begin( site_partners(resNum, *atmNum).begin());
			vector1< HBondCOP >::iterator partners_end( site_partners(resNum, *atmNum).end() );

			sort(partners_begin, partners_end, HBond::hbond_energy_comparer);

		}
		// acceptor sites
		for( AtomIndices::const_iterator
			atmNum  = res.accpt_pos().begin(),
			atmNume = res.accpt_pos().end(); atmNum != atmNume; ++atmNum ) {
			site_id++;
			insert_site_row(pose, struct_id, site_id, resNum, *atmNum, false /*is not donor*/, db_session);
			site_ids(resNum,*atmNum) = site_id;
			insert_site_pdb_row(pose, resNum, *atmNum, *atmNum, struct_id, site_id, db_session);
			insert_site_environment_row(pose, resNum, *atmNum, struct_id, site_id, atom_sasa_s, atom_sasa_m, atom_sasa_l, site_partners, site_hbond_energies, db_session);
			insert_site_atoms_row(pose, resNum, *atmNum, struct_id, site_id, db_session);

			sort(site_partners(resNum,*atmNum).begin(),
				site_partners(resNum,*atmNum).end(),
				HBond::hbond_energy_comparer);
		}
	}

	for (Size hbond_id = 1; hbond_id <= hbond_set.nhbonds(); hbond_id++) {
		HBond const & hbond = hbond_set.hbond( hbond_id );
		if(!relevant_residues[hbond.don_res()] ||
		  !relevant_residues[hbond.acc_res()]) continue;

		insert_hbond_row(hbond, struct_id, hbond_id, site_ids, site_partners, db_session);
		insert_hbond_geom_coords(pose, hbond_set.hbond_options(), hbond, struct_id, hbond_id, db_session);
		insert_hbond_lennard_jones_row(pose, hbond, struct_id, hbond_id, db_session);
		insert_hbond_dehydron_row(pose, hbond, struct_id, hbond_id, db_session);
	}
	return 0;
}

void
HBondFeatures::insert_site_row(
	Pose const & pose,
	StructureID struct_id,
	Size site_id,
	Size resNum,
	Size atmNum,
	bool is_donor,
	sessionOP db_session
){


	Size chain( pose.chain(resNum) );
	string const & resType( pose.residue_type(resNum).name() );
	string atmType;


	string HBChemType;
	if (is_donor){
		Size batmNum = pose.residue(resNum).atom_base( atmNum );
		HBDonChemType hb_don_chem_type =
			get_hb_don_chem_type( batmNum, pose.residue(resNum) );
		HBChemType = HBondTypeManager::name_from_don_chem_type(hb_don_chem_type);
		atmType = pose.residue(resNum).atom_type(batmNum).name();
	} else {
		HBAccChemType hb_acc_chem_type =
			get_hb_acc_chem_type( atmNum, pose.residue(resNum) );
		HBChemType = HBondTypeManager::name_from_acc_chem_type(hb_acc_chem_type);
		atmType = pose.residue(resNum).atom_type(atmNum).name();
	}

	std::string statement_string = "INSERT INTO hbond_sites (struct_id, site_id, resNum, HBChemType, atmNum, is_donor, chain, resType, atmType) VALUES (?,?,?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,site_id);
	stmt.bind(3,resNum);
	stmt.bind(4,HBChemType);
	stmt.bind(5,atmNum);
	stmt.bind(6,is_donor);
	stmt.bind(7,chain);
	stmt.bind(8,resType);
	stmt.bind(9,atmType);
	basic::database::safely_write_to_database(stmt);

}

void
HBondFeatures::insert_site_pdb_row(
	Pose const & pose,
	Size resNum,
	Size,
  Size heavy_atmNum,
	StructureID struct_id,
	Size site_id,
	sessionOP db_session
){
	if(!pose.pdb_info()) return; //eg if this is a silent file structure

	string const pdb_chain(1,pose.pdb_info()->chain(resNum));
	int const pdb_resNum( pose.pdb_info()->number(resNum) );
	string const pdb_iCode(1,pose.pdb_info()->icode(resNum));
	Real const pdb_heavy_atom_temperature(
		pose.pdb_info()->temperature(resNum,heavy_atmNum) );
	Real const pdb_heavy_atom_occupancy(
		pose.pdb_info()->occupancy(resNum, heavy_atmNum) );

	std::string statement_string = "INSERT INTO hbond_sites_pdb (struct_id, site_id, chain, resNum, iCode, heavy_atom_temperature, heavy_atom_occupancy) VALUES (?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,site_id);
	stmt.bind(3,pdb_chain);
	stmt.bind(4,pdb_resNum);
	stmt.bind(5,pdb_iCode);
	stmt.bind(6,pdb_heavy_atom_temperature);
	stmt.bind(7,pdb_heavy_atom_occupancy);
	basic::database::safely_write_to_database(stmt);
}


void
HBondFeatures::insert_site_environment_row(
	Pose const & pose,
	Size resNum,
	Size atmNum,
	StructureID struct_id,
	Size site_id,
	AtomID_Map< Real > const & atom_sasa_s,
	AtomID_Map< Real > const & atom_sasa_m,
	AtomID_Map< Real > const & atom_sasa_l,
	AtomID_Map< vector1<HBondCOP> > const & site_partners,
	AtomID_Map< Real > const & site_hbond_energies,
	sessionOP db_session
){

	Real const hbond_energy (site_hbond_energies(resNum, atmNum) );
	Size const num_hbonds(site_partners(resNum,atmNum).size() );

	string stmt_string("INSERT INTO hbond_site_environment (struct_id, site_id, sasa_r100, sasa_r140, sasa_r200, hbond_energy, num_hbonds) VALUES (?,?,?,?,?,?,?);");
	statement stmt(basic::database::safely_prepare_statement(stmt_string, db_session));

	stmt.bind(1, struct_id);
	stmt.bind(2, site_id);

	if(pose.is_fullatom()){
		stmt.bind(3, atom_sasa_s[AtomID(atmNum, resNum)]);
		stmt.bind(4, atom_sasa_m[AtomID(atmNum, resNum)]);
		stmt.bind(5, atom_sasa_l[AtomID(atmNum, resNum)]);
	} else {
		stmt.bind_null(3);
		stmt.bind_null(4);
		stmt.bind_null(5);
	}

	stmt.bind(6, hbond_energy);
	stmt.bind(7, num_hbonds);
	basic::database::safely_write_to_database(stmt);

}

void
HBondFeatures::insert_site_atoms_row(
	Pose const & pose,
	Size resNum,
	Size atmNum,
	StructureID struct_id,
	Size site_id,
	sessionOP db_session
){
	Residue const & rsd( pose.residue(resNum) );


	Real const atm_x( rsd.atom(atmNum).xyz().x() );
	Real const atm_y( rsd.atom(atmNum).xyz().y() );
	Real const atm_z( rsd.atom(atmNum).xyz().z() );
	Real const base_x( rsd.atom(rsd.atom_base(atmNum)).xyz().x() );
	Real const base_y( rsd.atom(rsd.atom_base(atmNum)).xyz().y() );
	Real const base_z( rsd.atom(rsd.atom_base(atmNum)).xyz().z() );
	Real const bbase_x( rsd.atom(rsd.atom_base(rsd.atom_base(atmNum))).xyz().x() );
	Real const bbase_y( rsd.atom(rsd.atom_base(rsd.atom_base(atmNum))).xyz().y() );
	Real const bbase_z( rsd.atom(rsd.atom_base(rsd.atom_base(atmNum))).xyz().z() );

	bool has_base2( rsd.abase2(atmNum) );
	Real base2_x=0, base2_y=0, base2_z =0;
	if( has_base2 ){
		base2_x = rsd.atom(rsd.abase2(atmNum)).xyz().x();
		base2_y = rsd.atom(rsd.abase2(atmNum)).xyz().y();
		base2_z = rsd.atom(rsd.abase2(atmNum)).xyz().z();
	}


	statement stmt = (*db_session)
		<< "INSERT INTO hbond_site_atoms (struct_id, site_id, atm_x, atm_y, atm_z, base_x, base_y, base_z, bbase_x, bbase_y, bbase_z, base2_x, base2_y, base2_z) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);"
		<< struct_id
		<< site_id
		<< atm_x << atm_y << atm_z
		<< base_x << base_y << base_z
		<< bbase_x << bbase_y << bbase_z;
	if( has_base2 ){
		stmt << base2_x << base2_y << base2_z;
	} else {
		stmt.bind_null();
		stmt.bind_null();
		stmt.bind_null();
	}
	basic::database::safely_write_to_database(stmt);

}

void
HBondFeatures::insert_hbond_row(
	HBond const & hbond,
	StructureID struct_id,
	Size hbond_id,
	AtomID_Map< Size > const & site_ids,        // This is for ranking hbonds
	AtomID_Map< vector1<HBondCOP> > const & site_partners,
	sessionOP db_session
){
	ASSERT_ONLY( bool found_don_partner( false ); )

	//Zero if unique
	//i in 1 through n if ith lowest energy hbond made with this site
	Size donRank=0, accRank=0;
	vector1< HBondCOP > don_partners(
		site_partners(hbond.don_res(), hbond.don_hatm()));
	if (don_partners.size() > 1){
		donRank++;
	}
	foreach(HBondCOP candidate_hbond, don_partners){
		if(hbond == *candidate_hbond){
			ASSERT_ONLY( found_don_partner = true; )
			break;
		} else {
			++donRank;
		}
	}
	assert(found_don_partner);

	ASSERT_ONLY( bool found_acc_partner(false); )
	vector1< HBondCOP > acc_partners(
		site_partners(hbond.acc_res(), hbond.acc_atm()));

	if (acc_partners.size() > 1){
		accRank++;
	}

	foreach(HBondCOP candidate_hbond, acc_partners){
		if(hbond == *candidate_hbond){
			ASSERT_ONLY( found_acc_partner = true; )
			break;
		} else {
			++accRank;
		}
	}
	assert(found_acc_partner);

	Size const don_id( site_ids(hbond.don_res(), hbond.don_hatm()) );
	Size const acc_id( site_ids(hbond.acc_res(), hbond.acc_atm()) );
	Size const HBEvalType( hbond.eval_type() );
	Real const energy( hbond.energy() );
	Real const envWeight( hbond.weight() );
	Real const score_weight(
		hb_eval_type_weight(hbond.eval_type(), scfxn_->weights(), false /*intra_res*/ ));

	statement stmt = (*db_session)
		<< "INSERT INTO hbonds (struct_id, hbond_id, don_id, acc_id, HBEvalType, energy, envWeight, score_weight, donRank, accRank) VALUES (?,?,?,?,?,?,?,?,?,?);"
		<< struct_id
		<< hbond_id
		<< don_id
		<< acc_id
		<< HBEvalType
		<< energy
		<< envWeight
		<< score_weight
		<< donRank
		<< accRank;
	basic::database::safely_write_to_database(stmt);
}

void
HBondFeatures::insert_hbond_geom_coords(
	Pose const & pose,
	HBondOptions const & hbond_options,
	HBond const & hbond,
	StructureID struct_id,
	Size hbond_id,
	sessionOP db_session
){
	const Residue acc_res( pose.residue(hbond.acc_res()));
	const Vector Axyz(acc_res.atom(hbond.acc_atm()).xyz());
	const Vector Bxyz(acc_res.atom(acc_res.atom_base(hbond.acc_atm())).xyz());
	const Vector B2xyz(acc_res.atom(acc_res.abase2(hbond.acc_atm())).xyz());
	const Residue don_res( pose.residue(hbond.don_res()));
	const Vector Hxyz(don_res.atom(hbond.don_hatm()).xyz());
	const Vector Dxyz(don_res.atom(don_res.atom_base(hbond.don_hatm())).xyz());


	// Paraphrase hb_energy_deriv (when called with atomic coordinates):
	Vector HDunit( Dxyz - Hxyz );
	Real const HDdist2( HDunit.length_squared() );
	Real const inv_HDdist = 1.0f / std::sqrt( HDdist2 );
	HDunit *= inv_HDdist;
	Vector BAunit;
	Vector PBxyz; // pseudo base atom
	Hybridization acc_hybrid( get_hbe_acc_hybrid( hbond.eval_type() ) );
	make_hbBasetoAcc_unitvector(hbond_options, acc_hybrid, Axyz, Bxyz, B2xyz, PBxyz, BAunit);

	// Paraphrase hb_energy_deriv (when called with coords/vectors)
	Vector AH = Hxyz - Axyz;
	Real const AHdist2( AH.length_squared() );
	Real const AHdist = std::sqrt(AHdist2);
	Real const inv_AHdist = 1.0f / AHdist;
	Vector AHunit = AH * inv_AHdist;

	Real const cosAHD = dot( AHunit, HDunit );
	Real const cosBAH = dot( BAunit, AHunit );

	float const chi(dihedral_radians(B2xyz, Bxyz, Axyz, Hxyz));

	statement stmt = (*db_session)
		<< "INSERT INTO hbond_geom_coords VALUES (?,?,?,?,?,?);"
		<< struct_id
		<< hbond_id
		<< AHdist
		<< cosBAH
		<< cosAHD
		<< chi;
	basic::database::safely_write_to_database(stmt);

}

// right now this just reports the lennard jones values for hbond atoms
void
HBondFeatures::insert_hbond_lennard_jones_row(
	Pose const & pose,
	HBond const & hbond,
	StructureID struct_id,
	Size hbond_id,
	sessionOP db_session
){

	Residue const & don_res(pose.residue(hbond.don_res()));
	Residue const & acc_res(pose.residue(hbond.acc_res()));
	Size const don_hatmNum(hbond.don_hatm());
	Size const acc_atmNum(hbond.acc_atm());
	Size const don_datmNum(don_res.atom_base(don_hatmNum));
	Size const acc_batmNum(acc_res.atom_base(acc_atmNum));


	Real bb_dummy, dsq_dummy;
	Real don_acc_atrE, don_acc_repE, don_acc_solv;
	Real don_acc_base_atrE, don_acc_base_repE, don_acc_base_solv;
	Real h_acc_atrE, h_acc_repE, h_acc_solv;
	Real h_acc_base_atrE, h_acc_base_repE, h_acc_base_solv;

	if ( !(scfxn_->energy_method_options().analytic_etable_evaluation()) ) {
		TableLookupEtableEnergy const etable_energy(
			*ScoringManager::get_instance()->etable(
				scfxn_->energy_method_options().etable_type() ),
			scfxn_->energy_method_options() );

		etable_energy.atom_pair_energy(
			don_res.atom(don_datmNum), acc_res.atom(acc_atmNum),
			/*weight*/ 1,
			don_acc_atrE, don_acc_repE, don_acc_solv,
			bb_dummy, dsq_dummy );

		etable_energy.atom_pair_energy(
			don_res.atom(don_datmNum), acc_res.atom(acc_batmNum),
			/*weight*/ 1,
			don_acc_base_atrE, don_acc_base_repE, don_acc_base_solv,
			bb_dummy, dsq_dummy );

		etable_energy.atom_pair_energy(
			don_res.atom(don_hatmNum), acc_res.atom(acc_atmNum),
			/*weight*/ 1,
			h_acc_atrE, h_acc_repE, h_acc_solv,
			bb_dummy, dsq_dummy );

		etable_energy.atom_pair_energy(
			don_res.atom(don_hatmNum), acc_res.atom(acc_batmNum),
			/*weight*/ 1,
			h_acc_base_atrE, h_acc_base_repE, h_acc_base_solv,
			bb_dummy, dsq_dummy );
	} else {
		AnalyticEtableEnergy const etable_energy(
			*ScoringManager::get_instance()->etable(
				scfxn_->energy_method_options().etable_type() ),
			scfxn_->energy_method_options() );

		etable_energy.atom_pair_energy(
			don_res.atom(don_datmNum), acc_res.atom(acc_atmNum),
			/*weight*/ 1,
			don_acc_atrE, don_acc_repE, don_acc_solv,
			bb_dummy, dsq_dummy );

		etable_energy.atom_pair_energy(
			don_res.atom(don_datmNum), acc_res.atom(acc_batmNum),
			/*weight*/ 1,
			don_acc_base_atrE, don_acc_base_repE, don_acc_base_solv,
			bb_dummy, dsq_dummy );

		etable_energy.atom_pair_energy(
			don_res.atom(don_hatmNum), acc_res.atom(acc_atmNum),
			/*weight*/ 1,
			h_acc_atrE, h_acc_repE, h_acc_solv,
			bb_dummy, dsq_dummy );

		etable_energy.atom_pair_energy(
			don_res.atom(don_hatmNum), acc_res.atom(acc_batmNum),
			/*weight*/ 1,
			h_acc_base_atrE, h_acc_base_repE, h_acc_base_solv,
			bb_dummy, dsq_dummy );
	}

	don_acc_atrE *= (*scfxn_)[ fa_atr ];
	don_acc_repE *= (*scfxn_)[ fa_rep ];
	don_acc_solv *= (*scfxn_)[ fa_sol ];

	don_acc_base_atrE *= (*scfxn_)[ fa_atr ];
	don_acc_base_repE *= (*scfxn_)[ fa_rep ];
	don_acc_base_solv *= (*scfxn_)[ fa_sol ];

	h_acc_atrE *= (*scfxn_)[ fa_atr ];
	h_acc_repE *= (*scfxn_)[ fa_rep ];
	h_acc_solv *= (*scfxn_)[ fa_sol ];

	h_acc_base_atrE *= (*scfxn_)[ fa_atr ];
	h_acc_base_repE *= (*scfxn_)[ fa_rep ];
	h_acc_base_solv *= (*scfxn_)[ fa_sol ];


	statement stmt = (*db_session)
		<< "INSERT INTO hbond_lennard_jones (struct_id, hbond_id, don_acc_atrE, don_acc_repE, don_acc_solv, don_acc_base_atrE, don_acc_base_repE, don_acc_base_solv, h_acc_atrE, h_acc_repE, h_acc_solv, h_acc_base_atrE, h_acc_base_repE, h_acc_base_solv) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);"
		<< struct_id
		<< hbond_id
		<< don_acc_atrE
		<< don_acc_repE
		<< don_acc_solv
		<< don_acc_base_atrE
		<< don_acc_base_repE
		<< don_acc_base_solv
		<< h_acc_atrE
		<< h_acc_repE
		<< h_acc_solv
		<< h_acc_base_atrE
		<< h_acc_base_repE
		<< h_acc_base_solv;
	basic::database::safely_write_to_database(stmt);

}

// This follows the definition of the dehydron described in:
//﻿1. Fernández A, Scheraga H a. Insufficiently dehydrated hydrogen bonds as determinants of protein interactions. Proceedings of the National Academy of Sciences of the United States of America. 2003;100(1):113-8. Available at: http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=140898&tool=pmcentrez&rendertype=abstract.

// The idea is that if the wrapping number is low, then the hydrogen
// bond is "underwrapped". Underwrapped hydrogen bonds are 'sticky'
// because, when they are wrapped, the hydrogen bond becomes harder to
// break. I'm not sure how this makes them sticky, but it is at least
// something that can be easily measured so why not?
void
HBondFeatures::insert_hbond_dehydron_row(
	Pose const & pose,
	HBond const & hbond,
    StructureID struct_id,
	Size hbond_id,
	sessionOP db_session
){


	Real const wrapping_radius(6.5);
	Size wrapping_count(0);

	Residue const & don_res(pose.residue(hbond.don_res()));
	Residue const & acc_res(pose.residue(hbond.acc_res()));

	Vector const & don_res_ca_xyz(don_res.xyz("CA"));
	Vector const & acc_res_ca_xyz(acc_res.xyz("CA"));


	TenANeighborGraph const & tenA(pose.energies().tenA_neighbor_graph());

	// For each neighboring residue to the donor
	for(EdgeListConstIterator
		ni = tenA.get_node(hbond.don_res())->const_edge_list_begin(),
		ni_end = tenA.get_node(hbond.don_res())->const_edge_list_end();
		ni != ni_end; ++ni){
		Residue const & nbr_res(pose.residue((*ni)->get_other_ind(hbond.don_res())));

		// sum all CH_n groups with in the wrapping radius of the c-alpha
		// atom of the donor residue.
		for(Size atm_i = 1; atm_i <= nbr_res.nheavyatoms(); ++atm_i){
			if(nbr_res.type().atom_type(atm_i).element() == "C" &&
			 don_res_ca_xyz.distance(nbr_res.xyz(atm_i)) <= wrapping_radius){
				wrapping_count++;
			}
		}
	}

	// for each neighboring residue to the acceptor
	for(EdgeListConstIterator
		ni = tenA.get_node(hbond.don_res())->const_edge_list_begin(),
		ni_end = tenA.get_node(hbond.don_res())->const_edge_list_end();
		ni != ni_end; ++ni){
		Residue const & nbr_res(pose.residue((*ni)->get_other_ind(hbond.don_res())));

		// sum all CH_n groups with the wrapping radius of the c-alpha
		// atom of the acceptor but not within wrapping radius of the
		// c-alpha atom of the donor. This prevents double counting
		// wrapping non-polar groups that are in both wrapping spheres.
		for(Size atm_i = 1; atm_i <= nbr_res.nheavyatoms(); ++atm_i){
			if(nbr_res.type().atom_type(atm_i).element() == "C" &&
				don_res_ca_xyz.distance(nbr_res.xyz(atm_i)) > wrapping_radius &&
				acc_res_ca_xyz.distance(nbr_res.xyz(atm_i)) <= wrapping_radius){
				wrapping_count++;
			}
		}
	}

	statement stmt = (*db_session)
		<< "INSERT INTO hbond_dehydrons (struct_id, hbond_id, wrapping_count) VALUES (?,?,?);"
		<< struct_id
		<< hbond_id
		<< wrapping_count;
	basic::database::safely_write_to_database(stmt);

}




} // namesapce
} // namespace
