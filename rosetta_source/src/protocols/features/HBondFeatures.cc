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
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <cmath>
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
using core::scoring::etable::EtableEnergy;
using core::scoring::hbonds::HBDonChemType;
using core::scoring::hbonds::HBAccChemType;
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

string
HBondFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS hbond_chem_types(\n"
		"	chem_type TEXT,\n"
		"	label TEXT,\n"
		"	PRIMARY KEY(chem_type));\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_NONE', 'aNONE');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_PBA', 'aPBA: bb');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_CXA', 'aCXA: n,q');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_CXL', 'aCXL: d,e');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_IMD', 'aIMD: h');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_IME', 'aIME: h');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_AHX', 'aAHX: t');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_HXL', 'aHXL: s,t');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_PCA_DNA', 'aPCA_DNA: O{1,2}P');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_PES_DNA', 'aPES_DNA: O{3,5}*');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_RRI_DNA', 'aRRI_DNA: O4*');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_PCA_RNA', 'aPCA_RNA: O{1,2}P');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_PES_RNA', 'aPES_RNA: O{3,5}*');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_RRI_RNA', 'aRRI_RNA: O4*');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_H2O', 'aH2O');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_GENERIC_SP2BB', 'aGEN: sp2 bb');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_GENERIC_SP2SC', 'aGEN: sp2 sc');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_GENERIC_SP3BB', 'aGEN: sp3 bb');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_GENERIC_SP3SC', 'aGEN: sp3 sc');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_GENERIC_RINGBB', 'aGEN: ring bb');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbacc_GENERIC_RINGSC', 'aGEN: ring sc');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_NONE', 'dNONE');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_PBA', 'dPBA: bb');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_CXA', 'dCXA: n,q');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_IMD', 'dIMD: h');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_IME', 'dIME: h');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_IND', 'dIND: w');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_AMO', 'dAMO: k');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_GDE', 'dGDE: r');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_GDH', 'dGDH: r');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_AHX', 'dAHX: s,t');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_HXL', 'dHXL: y');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_H2O', 'dH2O');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_GENERIC_BB', 'dGEN: bb');\n"
		"	INSERT OR IGNORE INTO hbond_chem_types VALUES('hbdon_GENERIC_SC', 'dGEN: sc');\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS hbond_sites (\n"
		"	struct_id INTEGER,\n"
		"	site_id INTEGER,\n"
		"	resNum INTEGER,\n"
		"	atmNum INTEGER,\n"
		"	is_donor BOOLEAN,\n"
		"	chain INTEGER,\n"
		"	resType TEXT,\n"
		"	atmType TEXT,\n"
		"	HBChemType TEXT,\n"
		"	FOREIGN KEY(struct_id, resNum)\n"
		"		REFERENCES residues(struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY(HBChemType)\n"
		"		REFERENCES hbond_chem_types(chem_type)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, site_id));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS hbond_sites_pdb (\n"
		"	struct_id INTEGER,\n"
		"	site_id INTEGER,\n"
		"	chain TEXT,\n"
		"	resNum INTEGER,\n"
		"	iCode TEXT,\n"
		"	heavy_atom_temperature REAL,\n"
		"	heavy_atom_occupancy REAL,\n"
		"	FOREIGN KEY(struct_id, site_id)\n"
		"		REFERENCES hbond_sites(struct_id, site_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, site_id));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS hbond_site_environment (\n"
		"	struct_id INTEGER,\n"
		"	site_id INTEGER,\n"
		"	sasa_r100 REAL,\n"
		"	sasa_r140 REAL,\n"
		"	sasa_r200 REAL,\n"
		"	hbond_energy REAL,\n"
		"	num_hbonds INTEGER,\n"
		"	FOREIGN KEY(struct_id, site_id)\n"
		"		REFERENCES hbond_sites(struct_id, site_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, site_id));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS hbond_site_atoms (\n"
		"	struct_id INTEGER,\n"
		"	site_id INTEGER,\n"
		"	atm_x REAL,\n"
		"	atm_y REAL,\n"
		"	atm_z REAL,\n"
		"	base_x REAL,\n"
		"	base_y REAL,\n"
		"	base_z REAL,\n"
		"	bbase_x REAL,\n"
		"	bbase_y REAL,\n"
		"	bbase_z REAL,\n"
		"	base2_x REAL,\n"
		"	base2_y REAL,\n"
		"	base2_z REAL,\n"
		"	FOREIGN KEY(struct_id, site_id)\n"
		"		REFERENCES hbond_sites(struct_id, site_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, site_id));\n"

		"\n"
		"CREATE TABLE IF NOT EXISTS hbonds (\n"
		"	struct_id INTEGER,\n"
		"	hbond_id INTEGER,\n"
		"	don_id INTEGER,\n"
		"	acc_id INTEGER,\n"
		"	HBEvalType INTEGER,\n"
		"	energy REAL,\n"
		"	envWeight REAL,\n"
		"	score_weight REAL,\n"
		"	donRank INTEGER,\n"
		"	accRank INTEGER,\n"
		"	FOREIGN KEY(struct_id, don_id)\n"
		"		REFERENCES hbond_sites(struct_id, site_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY(struct_id, acc_id)\n"
		"		REFERENCES hbond_sites(struct_id, site_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, hbond_id));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS hbond_lennard_jones (\n"
		"	struct_id INTEGER,\n"
		"	hbond_id INTEGER,\n"
		"	don_acc_atrE REAL,\n"
		"	don_acc_repE REAL,\n"
		"	don_acc_solv REAL,\n"
		"	don_acc_base_atrE REAL,\n"
		"	don_acc_base_repE REAL,\n"
		"	don_acc_base_solv REAL,\n"
		"	h_acc_atrE REAL,\n"
		"	h_acc_repE REAL,\n"
		"	h_acc_solv REAL,\n"
		"	h_acc_base_atrE REAL,\n"
		"	h_acc_base_repE REAL,\n"
		"	h_acc_base_solv REAL,\n"
		"	FOREIGN KEY(struct_id, hbond_id)\n"
		"		REFERENCES hbonds(struct_id, hbond_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, hbond_id));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS hbond_geom_coords (\n"
		"	struct_id INTEGER,\n"
		"	hbond_id INTEGER,\n"
		"	AHdist REAL,\n"
		"	cosBAH REAL,\n"
		"	cosAHD REAL,\n"
		"	chi REAL,\n"
		"	FOREIGN KEY (struct_id, hbond_id)\n"
		"		REFERENCES hbonds (struct_id, hbond_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, hbond_id));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS hbond_dehydrons (\n"
		"	struct_id INTEGER,\n"
		"	hbond_id INTEGER,\n"
		"	wrapping_count INTEEGER,\n"
		"	FOREIGN KEY(struct_id, hbond_id)\n"
		"		REFERENCES hbonds (struct_id, hbond_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, hbond_id));\n";

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
		utility_exit_with_message(error_msg.str());
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
		utility_exit_with_message(error_msg.str());
	}

	definition_threshold_ =
		tag->getOption<Real>("definition_threshold", 0);

}

Size
HBondFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
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

	Real const probe_radius_s(1.0);
	AtomID_Map< Real > atom_sasa_s;
	vector1< Real > residue_sasa_s;
	calc_per_atom_sasa(pose, atom_sasa_s, residue_sasa_s, probe_radius_s);

	Real const probe_radius_m(1.4);
	AtomID_Map< Real > atom_sasa_m;
	vector1< Real > residue_sasa_m;
	calc_per_atom_sasa(pose, atom_sasa_m, residue_sasa_m, probe_radius_m);

	Real const probe_radius_l(2.0);
	AtomID_Map< Real > atom_sasa_l;
	vector1< Real > residue_sasa_l;
	calc_per_atom_sasa(pose, atom_sasa_l, residue_sasa_l, probe_radius_l);

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
		insert_hbond_geom_coords(pose, hbond, struct_id, hbond_id, db_session);
		insert_hbond_lennard_jones_row(pose, hbond, struct_id, hbond_id, db_session);
		insert_hbond_dehydron_row(pose, hbond, struct_id, hbond_id, db_session);
	}
	return 0;
}

void
HBondFeatures::insert_site_row(
	Pose const & pose,
	Size struct_id,
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

	std::string statement_string = "INSERT INTO hbond_sites VALUES (?,?,?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,site_id);
	stmt.bind(3,resNum);
	stmt.bind(4,atmNum);
	stmt.bind(5,is_donor);
	stmt.bind(6,chain);
	stmt.bind(7,resType);
	stmt.bind(8,atmType);
	stmt.bind(9,HBChemType);
	basic::database::safely_write_to_database(stmt);

}

void
HBondFeatures::insert_site_pdb_row(
	Pose const & pose,
	Size resNum,
	Size atmNum,
  Size heavy_atmNum,
	Size struct_id,
	Size site_id,
	sessionOP db_session
){
	if(!pose.pdb_info()) return; //eg if this is a silent file structure

	string const pdb_chain(1,pose.pdb_info()->chain(resNum));
	Size const pdb_resNum( pose.pdb_info()->number(resNum) );
	string const pdb_iCode(1,pose.pdb_info()->icode(resNum));
	Real const pdb_heavy_atom_temperature(
		pose.pdb_info()->temperature(resNum,heavy_atmNum) );
	Real const pdb_heavy_atom_occupancy(
		pose.pdb_info()->occupancy(resNum, heavy_atmNum) );

	std::string statement_string = "INSERT INTO hbond_sites_pdb VALUES (?,?,?,?,?,?,?);";
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
	Size struct_id,
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

	statement stmt = (*db_session)
		<< "INSERT INTO hbond_site_environment VALUES (?,?,?,?,?,?,?);"
		<< struct_id
		<< site_id
		<< atom_sasa_s[AtomID(atmNum, resNum)]
		<< atom_sasa_m[AtomID(atmNum, resNum)]
		<< atom_sasa_l[AtomID(atmNum, resNum)]
		<< hbond_energy
		<< num_hbonds;
	basic::database::safely_write_to_database(stmt);

}

void
HBondFeatures::insert_site_atoms_row(
	Pose const & pose,
	Size resNum,
	Size atmNum,
	Size struct_id,
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
		<< "INSERT INTO hbond_site_atoms VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);"
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

Size
HBondFeatures::insert_hbond_row(
	HBond const & hbond,
	Size struct_id,
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
		hb_eval_type_weight(hbond.eval_type(), scfxn_->weights()));

	statement stmt = (*db_session)
		<< "INSERT INTO hbonds VALUES (?,?,?,?,?,?,?,?,?,?);"
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
	return stmt.last_insert_id();
}

void
HBondFeatures::insert_hbond_geom_coords(
	Pose const & pose,
	HBond const & hbond,
	Size struct_id,
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
	make_hbBasetoAcc_unitvector(acc_hybrid, Axyz, Bxyz, B2xyz, PBxyz, BAunit);

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
	Size struct_id,
	Size hbond_id,
	sessionOP db_session
){

	Residue const & don_res(pose.residue(hbond.don_res()));
	Residue const & acc_res(pose.residue(hbond.acc_res()));
	Size const don_hatmNum(hbond.don_hatm());
	Size const acc_atmNum(hbond.acc_atm());
	Size const don_datmNum(don_res.atom_base(don_hatmNum));
	Size const acc_batmNum(acc_res.atom_base(acc_atmNum));

	EtableEnergy const etable_energy(
		*ScoringManager::get_instance()->etable(
			scfxn_->energy_method_options().etable_type() ),
		scfxn_->energy_method_options() );

	Real bb_dummy, dsq_dummy;

	Real don_acc_atrE, don_acc_repE, don_acc_solv;
	etable_energy.atom_pair_energy(
		don_res.atom(don_datmNum), acc_res.atom(acc_atmNum),
		/*weight*/ 1,
		don_acc_atrE, don_acc_repE, don_acc_solv,
		bb_dummy, dsq_dummy );
	don_acc_atrE *= (*scfxn_)[ fa_atr ];
	don_acc_repE *= (*scfxn_)[ fa_rep ];
	don_acc_solv *= (*scfxn_)[ fa_sol ];

	Real don_acc_base_atrE, don_acc_base_repE, don_acc_base_solv;
	etable_energy.atom_pair_energy(
		don_res.atom(don_datmNum), acc_res.atom(acc_batmNum),
		/*weight*/ 1,
		don_acc_base_atrE, don_acc_base_repE, don_acc_base_solv,
		bb_dummy, dsq_dummy );
	don_acc_base_atrE *= (*scfxn_)[ fa_atr ];
	don_acc_base_repE *= (*scfxn_)[ fa_rep ];
	don_acc_base_solv *= (*scfxn_)[ fa_sol ];

	Real h_acc_atrE, h_acc_repE, h_acc_solv;
	etable_energy.atom_pair_energy(
		don_res.atom(don_hatmNum), acc_res.atom(acc_atmNum),
		/*weight*/ 1,
		h_acc_atrE, h_acc_repE, h_acc_solv,
		bb_dummy, dsq_dummy );
	h_acc_atrE *= (*scfxn_)[ fa_atr ];
	h_acc_repE *= (*scfxn_)[ fa_rep ];
	h_acc_solv *= (*scfxn_)[ fa_sol ];

	Real h_acc_base_atrE, h_acc_base_repE, h_acc_base_solv;
	etable_energy.atom_pair_energy(
		don_res.atom(don_hatmNum), acc_res.atom(acc_batmNum),
		/*weight*/ 1,
		h_acc_base_atrE, h_acc_base_repE, h_acc_base_solv,
		bb_dummy, dsq_dummy );
	h_acc_base_atrE *= (*scfxn_)[ fa_atr ];
	h_acc_base_repE *= (*scfxn_)[ fa_rep ];
	h_acc_base_solv *= (*scfxn_)[ fa_sol ];


	statement stmt = (*db_session)
		<< "INSERT INTO hbond_lennard_jones VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);"
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
	Size struct_id,
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
		<< "INSERT INTO hbond_dehydrons VALUES (?,?,?);"
		<< struct_id
		<< hbond_id
		<< wrapping_count;
	basic::database::safely_write_to_database(stmt);

}




} // namesapce
} // namespace
