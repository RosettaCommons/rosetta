// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/InterfaceFeatures.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/features/InterfaceFeatures.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

//Core Headers
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>

//Protocols Headers
#include <protocols/features/util.hh>

//Basic Headers
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/sql_utils.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

#include <cppdb/frontend.h>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/InterfaceFeaturesCreator.hh>

static basic::Tracer TR( "protocols.InterfaceFeatures" );

namespace protocols {
namespace features {

using utility::vector1;
using cppdb::statement;
using cppdb::result;
using std::string;
using namespace core::scoring;
using namespace protocols::analysis; // for total, side1, etc.

InterfaceFeatures::InterfaceFeatures() :
	FeaturesReporter()
{
	interface_analyzer_ = nullptr;
	scorefxn_ = core::scoring::get_score_function();
	set_defaults();
}

InterfaceFeatures::InterfaceFeatures(core::scoring::ScoreFunctionCOP scorefxn) :
	FeaturesReporter()
{
	interface_analyzer_ = nullptr;
	scorefxn_ = scorefxn;
	set_defaults();
}

void
InterfaceFeatures::set_defaults(){
	pack_separated_ = true;
	pack_together_ = false;
	dSASA_cutoff_ = 100;
	compute_packstat_ = true;
}

void
InterfaceFeatures::set_pack_separated(const bool pack_separated){
	pack_separated_ = pack_separated;
}

void
InterfaceFeatures::set_pack_together(const bool pack_together) {
	pack_together_ = pack_together;
}

void
InterfaceFeatures::set_dSASA_cutoff(core::Real dSASA_cutoff){
	dSASA_cutoff_ = dSASA_cutoff;
}

void
InterfaceFeatures::set_scorefxn(ScoreFunctionOP scorefxn){
	scorefxn_ = scorefxn;
}


void
InterfaceFeatures::set_compute_packstat(const bool compute_packstat) {
	compute_packstat_ = compute_packstat;
}

void
InterfaceFeatures::set_interface_chains(vector1<std::string> const & interfaces){
	interfaces_ = interfaces;
}

InterfaceFeatures::~InterfaceFeatures()= default;

// XRW TEMP std::string
// XRW TEMP InterfaceFeatures::type_name() const {
// XRW TEMP  return "InterfaceFeatures";
// XRW TEMP }

void
InterfaceFeatures::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & /*data*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/)
{
	pack_separated_ = tag->getOption<bool>("pack_separated", true);
	pack_together_ = tag->getOption<bool>("pack_together", false);
	dSASA_cutoff_ = tag->getOption<core::Real>("dSASA_cutoff", 100);
	compute_packstat_ = tag->getOption<bool>("compute_packstat", true);

	if ( tag->hasOption("scorefxn") ) {
		string const scorefxn_name(tag->getOption<std::string>("scorefxn"));
		scorefxn_ = data.get_ptr<core::scoring::ScoreFunction>("scorefxns", scorefxn_name);
	}

	if ( tag->hasOption("interfaces") && tag->hasOption("interface") ) {
		utility_exit_with_message("Cannot specify both interface and interfaces option for InterfaceFeatures Reporter");
	}

	if ( tag->hasOption("interface") ) {
		std::string interface = tag->getOption<std::string>("interface");
		if ( interface.find(",") != std::string::npos ) { utility_exit_with_message("Only one interface should be specified using the interface option");}
		interfaces_.push_back(interface);
	}

	if ( tag->hasOption("interfaces") ) {
		std::string interfaces = tag->getOption<std::string>("interfaces");
		interfaces_ = utility::string_split_multi_delim(interfaces, ":,'`~+*&|;.");//Why not?
	}


}

void
InterfaceFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	write_interface_schema_to_db(db_session);
	write_interface_residues_schema_to_db(db_session);
	write_interface_side_schema_to_db(db_session);
}


utility::vector1< std::string >
InterfaceFeatures::features_reporter_dependencies() const{
	utility::vector1< std::string > dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("StructureFeatures");
	return dependencies;
}

core::Size
InterfaceFeatures::report_features(
	core::pose::Pose const & pose,
	const utility::vector1<bool>& relevant_residues,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session) {

	if ( interfaces_.empty() ) {
		make_interface_combos(pose, interfaces_);
	}

	for ( core::Size i = 1; i <= interfaces_.size(); ++i ) {

		std::string interface = interfaces_[i];
		report_all_interface_features(pose, relevant_residues, struct_id, db_session, interface, interface);
	}
	return 0;
}

void
InterfaceFeatures::report_all_interface_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & relevant_residues,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	std::string const & interface,
	std::string const & db_interface)
{
	TR << "reporting features for: "<< interface << std::endl;
	//Check to make sure interface/chain definition is solid.
	if ( ! interface.find('_') || !db_interface.find('_') ) {
		utility_exit_with_message("Unrecognized interface: "+interface+" must have side1 and side2, ex: LH_A or L_H to calculate interface data");
	}

	if ( !chains_exist_in_pose(pose, interface) ) {
		TR <<"All chains do not exist in the given pose.  Skipping interface: " << interface << std::endl;
		return;
	}

	vector1<std::string> interface_sides = utility::string_split(db_interface, '_');
	std::string interface_side1 = interface_sides[1];
	std::string interface_side2 = interface_sides[2];
	interface_analyzer_ = protocols::analysis::InterfaceAnalyzerMoverOP( new protocols::analysis::InterfaceAnalyzerMover(interface, true, scorefxn_, compute_packstat_, pack_together_, pack_separated_) );
	interface_analyzer_->set_use_centroid_dG(false); //Uses score3 by default without rg - used in zinc homodimer design .  Used for clash detection
	interface_analyzer_->apply_const(pose);

	if ( interface_analyzer_->get_interface_delta_sasa() < dSASA_cutoff_ ) {
		TR << "Interface dSASA lower than set cutoff value of "<< dSASA_cutoff_ << " : not including data in database..." << std::endl;
		return;
	}
	report_interface_features(pose, struct_id, db_session, interface_side1, interface_side2);
	report_interface_residue_features(pose, relevant_residues, struct_id, db_session, interface_side1, interface_side2);

	report_interface_side_features(pose, struct_id, db_session, interface_side1, interface_side2, total, "total");
	report_interface_side_features(pose, struct_id, db_session, interface_side1, interface_side2, side1, "side1");
	report_interface_side_features(pose, struct_id, db_session, interface_side1, interface_side2, side2,  "side2");
}


void
InterfaceFeatures::make_interface_combos(const core::pose::Pose& pose, vector1< std::string > & interfaces) {


	TR << "making interface combos"<<std::endl;
	vector1< std::string > chain_combos;
	vector1< std::string > size_combos;

	std::string chains = get_all_pose_chains(pose);
	get_length_combos(chains, size_combos);

	for ( core::Size i = 1; i <= size_combos.size(); ++i ) {
		get_all_string_combos(size_combos[i], "", chain_combos);
	}

	//Remove duplicates
	std::sort(chain_combos.begin(), chain_combos.end());
	chain_combos.erase(std::unique(chain_combos.begin(), chain_combos.end()), chain_combos.end());

	TR << std::endl;
	for ( core::Size i = 1; i<= chain_combos.size(); ++i ) {
		for ( core::Size s = 1; s <= chain_combos[i].length()-1; ++s ) {
			std::string dock_chains = chain_combos[i].substr(0, s)+"_"+chain_combos[i].substr(s);
			if ( ! interface_exists(interfaces, dock_chains) ) {
				interfaces.push_back(dock_chains);
			}
		}
	}
	std::sort(interfaces.begin(), interfaces.end());
	for ( core::Size i = 1; i<= interfaces.size(); ++i ) {
		TR << "added interface:  " + interfaces[i] << std::endl;
	}

}

void
InterfaceFeatures::get_length_combos(std::string current, vector1<std::string> & sizes) const {
	if ( current.length() >=2 ) {
		sizes.push_back(current);
		if ( current.length() == 2 ) {
			return;
		}
	}

	for ( std::size_t k = 1; k <= current.length(); ++k ) {
		std::string new_string(current);
		//char chain = current.at(k-1);
		new_string.erase(k-1, 1); //Remove the chain from the interface
		get_length_combos(new_string, sizes);
	}
}


void
InterfaceFeatures::get_all_string_combos(std::string& interface, std::string current, vector1<std::string>& chains) const {

	if ( interface.length() == current.length() ) {
		chains.push_back(current);
		//TR << "Chain added: "<<current << std::endl;
		return;
	}
	for ( core::Size i = 1; i<=interface.length(); ++i ) {

		std::string letter = interface.substr(i-1, 1);
		if ( current.find(letter) != std::string::npos ) {
			continue;
		}

		std::string new_current = current + letter;
		get_all_string_combos(interface, new_current, chains);
	}
}

bool
InterfaceFeatures::chains_exist_in_pose(core::pose::Pose const & pose, std::string const & interface) const {
	utility::vector1<std::string> interfaceSP = utility::string_split(interface, '_');
	std::string chains = interfaceSP[1]+interfaceSP[2];

	bool has_all_chains = true;
	for ( core::Size i = 1; i<=chains.length(); ++i ) {
		char chain = chains.at(i-1);
		if ( ! core::pose::has_chain(chain, pose) ) {
			has_all_chains = false;
			break;
		}
	}
	return has_all_chains;
}

std::string
InterfaceFeatures::get_all_pose_chains(core::pose::Pose const & pose){

	std::string chains = "";
	for ( core::Size i = 1; i <= pose.conformation().num_chains(); ++i ) {
		char chain_c = core::pose::get_chain_from_chain_id(i, pose);
		std::string chain = utility::to_string(chain_c);
		chains = chains+chain;
	}
	TR <<"Pose chains: "<<chains << std::endl;
	return chains;
}


bool
InterfaceFeatures::interface_exists(vector1<std::string> & interfaces, std::string & dock_chains) const {
	utility::vector1<std::string> newSP = utility::string_split(dock_chains, '_');

	//Sort so order of chains doesn't matter (LH_AB = HL_BA)
	std::sort(newSP[1].begin(), newSP[1].end());
	std::sort(newSP[2].begin(), newSP[2].end());

	for ( core::Size i = 1; i <= interfaces.size(); ++i ) {
		utility::vector1<std::string> oldSP = utility::string_split(interfaces[i], '_');
		std::sort(oldSP[1].begin(), oldSP[1].end());
		std::sort(oldSP[2].begin(), oldSP[2].end());


		if ( oldSP[1] == newSP[2] && oldSP[2] == newSP[1] ) {
			return true; //LH_AB == AB_LH
		} else if ( oldSP[1] == newSP[1] && oldSP[2] == newSP[2] ) {
			return true; // LH_AB == LH_AB
		} else {
			continue;
		}

	}
	return false;
}

void
InterfaceFeatures::report_interface_residue_features(
	const core::pose::Pose& pose,
	const utility::vector1<bool>& relevant_residues,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	std::string const & chains_side1,
	std::string const & chains_side2) const
{
	std::map<protocols::analysis::InterfaceRegion, std::string> regions;
	regions[side1] = "side1";
	regions[side2] = "side2";
	regions[total] = "total";

	protocols::analysis::PerResidueInterfaceData interface_data = interface_analyzer_->get_all_per_residue_data();
	protocols::analysis::InterfaceData all_data = interface_analyzer_->get_all_data();
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( check_relevant_residues(relevant_residues, i) && interface_data.interface_residues[i] ) {
			protocols::analysis::InterfaceRegion side;
			//Check the side.  Can only be either side1 or side2
			if ( all_data.interface_residues[side1][i] ) {
				side = side1;
			} else {
				side = side2;
			}
			write_interface_residue_data_row_to_db(struct_id, db_session, chains_side1, chains_side2, regions[side],i, interface_data);
		}
	}
}

void
InterfaceFeatures::write_interface_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column interface("interface", DbDataTypeOP( new DbText(255) ));
	Column chains_side1("chains_side1", DbDataTypeOP( new DbText(255) ));
	Column chains_side2("chains_side2", DbDataTypeOP( new DbText(255) ));

	Columns primary_keys;
	primary_keys.push_back(struct_id);
	primary_keys.push_back(interface);

	Column nchains_side1("nchains_side1", DbDataTypeOP( new DbInteger() ));
	Column nchains_side2("nchains_side2", DbDataTypeOP( new DbInteger() ));
	Column dG_cross("dG_cross", DbDataTypeOP( new DbReal() ));
	Column dG_cross_dev_dSASAx100("dG_cross_dev_dSASAx100", DbDataTypeOP( new DbReal() ));
	Column dG_separated("dG", DbDataTypeOP( new DbReal() ));
	Column dG_separated_dev_dSASAx100("dG_dev_dSASAx100", DbDataTypeOP( new DbReal() ));
	Column dSASA_hphobic("dSASA_hphobic", DbDataTypeOP( new DbReal() ));
	Column dSASA_int("dSASA", DbDataTypeOP( new DbReal() ));
	Column dSASA_polar("dSASA_polar", DbDataTypeOP( new DbReal() ));
	Column delta_unsatHbonds("delta_unsatHbonds", DbDataTypeOP( new DbReal() ));
	Column hbond_E_fraction("hbond_E_fraction", DbDataTypeOP( new DbReal() ));
	Column nres_all("nres_all", DbDataTypeOP( new DbInteger() ));
	Column nres_int("nres_int", DbDataTypeOP( new DbInteger() ));
	Column packstat("packstat", DbDataTypeOP( new DbReal() ));
	Column sc_value("sc_value", DbDataTypeOP( new DbReal() ));
	//Column side1_normalized("side1_normalized", new DbReal());
	//Column side1_score("side1_complexed_interface_score", new DbReal());
	//Column side1_nres("side1_nres", new DbReal());
	//Column side2_normalized("side2_normalized", new DbReal());
	//Column side2_nres("side2_nres", new DbReal());
	Column complex_normalized("complex_normalized", DbDataTypeOP( new DbReal() ));
	//Column side2_score("side2_complexed_interface_score", new DbReal());

	ForeignKey foreign_key(struct_id, "structures", "struct_id", true);
	PrimaryKey primary_key(primary_keys);

	Schema table("interfaces", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(chains_side1);
	table.add_column(chains_side2);
	table.add_column(nchains_side1);
	table.add_column(nchains_side2);
	table.add_column(dSASA_int);
	table.add_column(dSASA_hphobic);
	table.add_column(dSASA_polar);
	table.add_column(dG_separated);
	table.add_column(dG_cross);
	table.add_column(dG_separated_dev_dSASAx100);
	table.add_column(dG_cross_dev_dSASAx100);
	table.add_column(delta_unsatHbonds);
	table.add_column(hbond_E_fraction);
	table.add_column(sc_value);
	table.add_column(packstat);
	table.add_column(nres_int);
	table.add_column(nres_all);
	table.add_column(complex_normalized);
	table.write(db_session);
}

void
InterfaceFeatures::report_interface_features(
	const core::pose::Pose& pose,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	std::string const & chains_side1,
	std::string const & chains_side2) const
{
	using namespace protocols::analysis;

	std::string stmt_string = "INSERT INTO interfaces ("
		"struct_id,"
		"interface,"
		"chains_side1,"
		"chains_side2,"
		"nchains_side1,"
		"nchains_side2,"
		"dSASA,"
		"dSASA_hphobic,"
		"dSASA_polar,"
		"dG,"
		"dG_cross,"
		"dG_dev_dSASAx100,"
		"dG_cross_dev_dSASAx100,"
		"delta_unsatHbonds,"
		"hbond_E_fraction,"
		"sc_value,"
		"packstat,"
		"nres_int,"
		"nres_all,"
		"complex_normalized)"
		" VALUES "+get_question_mark_string(20);

	statement stmt(basic::database::safely_prepare_statement(stmt_string, db_session));

	protocols::analysis::InterfaceData data = interface_analyzer_->get_all_data();
	protocols::analysis::PerResidueInterfaceData res_data = interface_analyzer_->get_all_per_residue_data();

	// strings are passed to bind by reference
	// we need the lifetime of the combined string to exist at least until safely_write_to_database()
	std::string partners( chains_side1+"_"+chains_side2 );

	core::Size i = 0;
	stmt.bind(i+=1, struct_id);
	stmt.bind(i+=1, partners);
	stmt.bind(i+=1, chains_side1);
	stmt.bind(i+=1, chains_side2);
	stmt.bind(i+=1, chains_side1.length());
	stmt.bind(i+=1, chains_side2.length());
	stmt.bind(i+=1, data.dSASA[total]);
	stmt.bind(i+=1, data.dhSASA[total]);
	stmt.bind(i+=1, data.dSASA[total] - data.dhSASA[total]);
	stmt.bind(i+=1, data.dG[total]);
	stmt.bind(i+=1, data.crossterm_interface_energy);
	stmt.bind(i+=1, data.dG_dSASA_ratio*100);
	stmt.bind(i+=1,data.crossterm_interface_energy_dSASA_ratio*100);
	stmt.bind(i+=1,data.delta_unsat_hbonds);
	stmt.bind(i+=1,data.hbond_E_fraction);
	stmt.bind(i+=1,data.sc_value);
	stmt.bind(i+=1,data.packstat);
	stmt.bind(i+=1,data.interface_nres[total]);
	stmt.bind(i+=1,pose.size());
	stmt.bind(i+=1,data.complex_total_energy[total]/pose.size());
	basic::database::safely_write_to_database(stmt);
}

void
InterfaceFeatures::write_interface_side_schema_to_db(utility::sql_database::sessionOP db_session) const {
	using namespace basic::database::schema_generator;
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column interface("interface", DbDataTypeOP( new DbText(255) ));
	Column chains_side1("chains_side1", DbDataTypeOP( new DbText(255) ));
	Column chains_side2("chains_side2", DbDataTypeOP( new DbText(255) ));
	Column side("side", DbDataTypeOP( new DbText(255) ));
	Columns primary_keys;
	primary_keys.push_back(struct_id);
	primary_keys.push_back(interface);
	primary_keys.push_back(side);

	PrimaryKey primary_key(primary_keys);

	ForeignKey foreign_key(struct_id, "structures", "struct_id", true);

	Column complexed_interface_score("energy_int", DbDataTypeOP( new DbReal() ));
	Column separated_interface_score("energy_sep", DbDataTypeOP( new DbReal() ));
	Column dSASA("dSASA", DbDataTypeOP( new DbReal() ));
	Column dSASA_bb("dSASA_bb", DbDataTypeOP( new DbReal() ));
	Column dSASA_sc("dSASA_sc", DbDataTypeOP( new DbReal() ));

	Column dhSASA("dhSASA", DbDataTypeOP( new DbReal() ));
	Column dhSASA_bb("dhSASA_bb", DbDataTypeOP( new DbReal() ));
	Column dhSASA_sc("dhSASA_sc", DbDataTypeOP( new DbReal() ));
	Column dhSASA_rel("dhSASA_rel_by_charge", DbDataTypeOP( new DbReal() ));

	Column dG("dG", DbDataTypeOP( new DbReal() ));
	Column interface_nres("interface_nres", DbDataTypeOP( new DbInteger() ));

	Column aromatic_fraction("aromatic_fraction", DbDataTypeOP( new DbReal() ));
	Column aromatic_dSASA_fraction("aromatic_dSASA_fraction", DbDataTypeOP( new DbReal() ));
	Column aromatic_energy_fraction("aromatic_dG_fraction", DbDataTypeOP( new DbReal() ));

	Column interface_to_surface_fraction("interface_to_surface_fraction", DbDataTypeOP( new DbReal() ));

	Column ss_sheet_fraction("ss_sheet_fraction", DbDataTypeOP( new DbReal() ));
	Column ss_helix_fraction("ss_helix_fraction", DbDataTypeOP( new DbReal() ));
	Column ss_loop_fraction("ss_loop_fraction", DbDataTypeOP( new DbReal() ));

	Column avg_per_res_energy_int("avg_per_residue_energy_int", DbDataTypeOP( new DbReal() ));
	Column avg_per_res_energy_sep("avg_per_residue_energy_sep", DbDataTypeOP( new DbReal() ));
	Column avg_per_res_energy_dG("avg_per_residue_energy_dG", DbDataTypeOP( new DbReal() ));
	Column avg_per_res_SASA_int("avg_per_residue_SASA_int", DbDataTypeOP( new DbReal() ));
	Column avg_per_res_SASA_sep("avg_per_residue_SASA_sep", DbDataTypeOP( new DbReal() ));
	Column avg_per_res_dSASA("avg_per_residue_dSASA", DbDataTypeOP( new DbReal() ));

	Schema table("interface_sides", primary_keys);
	table.add_foreign_key(foreign_key);
	table.add_column(chains_side1);
	table.add_column(chains_side2);
	table.add_column(interface_nres);
	table.add_column(dSASA);
	table.add_column(dSASA_sc);
	table.add_column(dhSASA);
	table.add_column(dhSASA_sc);
	table.add_column(dhSASA_rel);
	table.add_column(dG);
	table.add_column(complexed_interface_score);
	table.add_column(separated_interface_score);
	table.add_column(avg_per_res_energy_dG);
	table.add_column(avg_per_res_energy_int);
	table.add_column(avg_per_res_energy_sep);
	table.add_column(avg_per_res_dSASA);
	table.add_column(avg_per_res_SASA_int);
	table.add_column(avg_per_res_SASA_sep);
	table.add_column(aromatic_fraction);
	table.add_column(aromatic_dSASA_fraction);
	table.add_column(aromatic_energy_fraction);
	table.add_column(interface_to_surface_fraction);
	table.add_column(ss_sheet_fraction);
	table.add_column(ss_helix_fraction);
	table.add_column(ss_loop_fraction);

	table.write(db_session);

}

void
InterfaceFeatures::report_interface_side_features(
	core::pose::Pose const &,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	std::string const & chains_side1,
	std::string const & chains_side2,
	protocols::analysis::InterfaceRegion region,
	std::string const & region_string) const
{

	std::string stmt_string = "INSERT INTO interface_sides ("
		"struct_id,"
		"interface,"
		"side,"
		"chains_side1,"
		"chains_side2,"
		"interface_nres,"
		"dSASA,"
		"dSASA_sc,"
		"dhSASA,"
		"dhSASA_sc,"
		"dhSASA_rel_by_charge,"
		"dG,"
		"energy_int,"
		"energy_sep,"
		"avg_per_residue_energy_dG,"
		"avg_per_residue_energy_int,"
		"avg_per_residue_energy_sep,"
		"avg_per_residue_dSASA,"
		"avg_per_residue_SASA_int,"
		"avg_per_residue_SASA_sep,"
		"aromatic_fraction,"
		"aromatic_dSASA_fraction,"
		"aromatic_dG_fraction,"
		"interface_to_surface_fraction,"
		"ss_sheet_fraction,"
		"ss_helix_fraction,"
		"ss_loop_fraction) VALUES "+get_question_mark_string(27);

	statement stmt(basic::database::safely_prepare_statement(stmt_string, db_session));

	protocols::analysis::InterfaceData data = interface_analyzer_->get_all_data();
	protocols::analysis::PerResidueInterfaceData res_data = interface_analyzer_->get_all_per_residue_data();

	// strings are passed to bind by reference
	// we need the lifetime of the combined string to exist at least until safely_write_to_database()
	std::string partners( chains_side1+"_"+chains_side2 );

	core::Size i = 0;
	stmt.bind(i+=1, struct_id);
	stmt.bind(i+=1, partners);
	stmt.bind(i+=1, region_string);
	stmt.bind(i+=1, chains_side1);
	stmt.bind(i+=1, chains_side2);
	stmt.bind(i+=1, data.interface_nres[region]);
	stmt.bind(i+=1, data.dSASA[region]);
	stmt.bind(i+=1, data.dSASA_sc[region]);
	stmt.bind(i+=1, data.dhSASA[region]);
	stmt.bind(i+=1, data.dhSASA_sc[region]);
	stmt.bind(i+=1, data.dhSASA_rel_by_charge[region]);
	stmt.bind(i+=1, data.dG[region]);
	stmt.bind(i+=1, data.complexed_interface_score[region]);
	stmt.bind(i+=1, data.separated_interface_score[region]);
	stmt.bind(i+=1, res_data.regional_avg_per_residue_dG[region]);
	stmt.bind(i+=1, res_data.regional_avg_per_residue_energy_int[region]);
	stmt.bind(i+=1, res_data.regional_avg_per_residue_energy_sep[region]);
	stmt.bind(i+=1, res_data.regional_avg_per_residue_dSASA[region]);
	stmt.bind(i+=1, res_data.regional_avg_per_residue_SASA_int[region]);
	stmt.bind(i+=1, res_data.regional_avg_per_residue_SASA_sep[region]);
	if ( data.interface_nres[region] > 0 ) {
		stmt.bind(i+=1, data.aromatic_nres[region]/(core::Real)data.interface_nres[region]);
	} else {
		stmt.bind(i+=1, 0.0);
	}

	stmt.bind(i+=1, data.aromatic_dSASA_fraction[region]);
	stmt.bind(i+=1, data.aromatic_dG_fraction[region]);
	stmt.bind(i+=1, data.interface_to_surface_fraction[region]);
	stmt.bind(i+=1, data.ss_sheet_nres[region]/(core::Real)data.interface_nres[region]);
	stmt.bind(i+=1, data.ss_helix_nres[region]/(core::Real)data.interface_nres[region]);
	stmt.bind(i+=1, data.ss_loop_nres[region]/(core::Real)data.interface_nres[region]);

	basic::database::safely_write_to_database(stmt);
}

void
InterfaceFeatures::write_interface_residues_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column interface("interface", DbDataTypeOP( new DbText(255) ));
	Column chains_side1("chains_side1", DbDataTypeOP( new DbText(255) ));
	Column chains_side2("chains_side2", DbDataTypeOP( new DbText(255) ));
	Column resnum("resNum", DbDataTypeOP( new DbInteger() ));

	Columns primary_keys;
	primary_keys.push_back(struct_id);
	primary_keys.push_back(interface);
	primary_keys.push_back(resnum);
	PrimaryKey primary_key(primary_keys);

	ForeignKey struct_foreign_key(struct_id, "structures", "struct_id", true);

	Columns foreign_keys;
	foreign_keys.push_back(struct_id);
	foreign_keys.push_back(resnum);

	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey res_foreign_key(foreign_keys, "residues", reference_columns, true);

	Column side("side", DbDataTypeOP( new DbText(255) ));
	Column SASA_sep("SASA_sep", DbDataTypeOP( new DbReal() ));
	Column SASA_int("SASA_int", DbDataTypeOP( new DbReal() ));
	Column dSASA("dSASA", DbDataTypeOP( new DbReal() ));
	Column dSASA_bb("dSASA_bb", DbDataTypeOP( new DbReal() ));
	Column dSASA_sc("dSASA_sc", DbDataTypeOP( new DbReal() ));
	Column dhSASA("dhSASA", DbDataTypeOP( new DbReal() ));
	Column dhSASA_bb("dhSASA_bb", DbDataTypeOP( new DbReal() ));
	Column dhSASA_sc("dhSASA_sc", DbDataTypeOP( new DbReal() ));
	Column dhSASA_rel("dhSASA_rel_by_charge", DbDataTypeOP( new DbReal() ));

	Column rel_dSASA_fraction("relative_dSASA_fraction", DbDataTypeOP( new DbReal() ));
	Column energy_sep("energy_sep", DbDataTypeOP( new DbReal() ));
	Column energy_int("energy_int", DbDataTypeOP( new DbReal() ));
	Column dG("dG", DbDataTypeOP( new DbReal() ));

	Schema table("interface_residues", primary_keys);
	table.add_foreign_key(res_foreign_key);
	table.add_foreign_key(struct_foreign_key);
	table.add_column(chains_side1);
	table.add_column(chains_side2);
	table.add_column(side);
	table.add_column(dSASA);
	table.add_column(dSASA_sc);
	table.add_column(dhSASA);
	table.add_column(dhSASA_sc);
	table.add_column(dhSASA_rel);
	table.add_column(SASA_int);
	table.add_column(SASA_sep);
	table.add_column(rel_dSASA_fraction);
	table.add_column(dG);
	table.add_column(energy_int);
	table.add_column(energy_sep);
	table.write(db_session);
}

void
InterfaceFeatures::write_interface_residue_data_row_to_db(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	std::string const & chains_side1,
	std::string const & chains_side2,
	std::string const & side,
	core::Size const resnum,
	protocols::analysis::PerResidueInterfaceData const & interface_data) const
{


	std::string stmnt_string = "INSERT INTO interface_residues ("
		"struct_id,"
		"interface,"
		"resNum,"
		"chains_side1,"
		"chains_side2,"
		"side,"
		"dSASA,"
		"dSASA_sc,"
		"dhSASA,"
		"dhSASA_sc,"
		"dhSASA_rel_by_charge,"
		"SASA_int,"
		"SASA_sep,"
		"relative_dSASA_fraction,"
		"dG,"
		"energy_int,"
		"energy_sep) VALUES "+get_question_mark_string(17);

	statement stmnt(basic::database::safely_prepare_statement(stmnt_string, db_session));

	// strings are passed to bind by reference
	// we need the lifetime of the combined string to exist at least until safely_write_to_database()
	std::string partners( chains_side1+"_"+chains_side2 );

	core::Size i = 0;
	stmnt.bind(i+=1, struct_id);
	stmnt.bind(i+=1, partners);
	stmnt.bind(i+=1, resnum);
	stmnt.bind(i+=1, chains_side1);
	stmnt.bind(i+=1, chains_side2);
	stmnt.bind(i+=1, side);
	stmnt.bind(i+=1, interface_data.dSASA[resnum]);
	stmnt.bind(i+=1, interface_data.dSASA_sc[resnum]);
	stmnt.bind(i+=1, interface_data.dhSASA[resnum]);
	stmnt.bind(i+=1, interface_data.dhSASA_sc[resnum]);
	stmnt.bind(i+=1, interface_data.dhSASA_rel_by_charge[resnum]);
	stmnt.bind(i+=1, interface_data.complexed_sasa[resnum]);
	stmnt.bind(i+=1, interface_data.separated_sasa[resnum]);
	if ( interface_data.separated_sasa[resnum]>0 ) {
		stmnt.bind(i+=1, interface_data.dSASA[resnum]/interface_data.separated_sasa[resnum]);
	} else {
		stmnt.bind(i+=1, 0);
	}
	stmnt.bind(i+=1, interface_data.dG[resnum]);
	stmnt.bind(i+=1,interface_data.complexed_energy[resnum]);
	stmnt.bind(i+=1, interface_data.separated_energy[resnum]);
	basic::database::safely_write_to_database(stmnt);

}

std::string InterfaceFeatures::type_name() const {
	return class_name();
}

std::string InterfaceFeatures::class_name() {
	return "InterfaceFeatures";
}

void InterfaceFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;

	std::string interface_pattern = "[" + chr_chains_nonrepeated() + "]+_[" + chr_chains_nonrepeated() + "]+";
	//interface_features_string
	//interface_features_list
	XMLSchemaRestriction str_res;
	str_res.name( "interface_features_string" );
	str_res.base_type( xs_string );
	str_res.add_restriction( xsr_pattern, interface_pattern );
	xsd.add_top_level_element( str_res );
	XMLSchemaRestriction list_res;
	list_res.name( "interface_features_list" );
	list_res.base_type( xs_string );
	//list_res.add_restriction( xsr_pattern, interface_pattern + "([:,'`~+*&|;.]" + interface_pattern + ")*" );
	// AMW: temp can't have &
	list_res.add_restriction( xsr_pattern, interface_pattern + "([:,'`~+*|;.]" + interface_pattern + ")*" );
	xsd.add_top_level_element( list_res );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "pack_separated", xsct_rosetta_bool, "Pack the structures separately", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "pack_together", xsct_rosetta_bool, "Pack the structures together", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "dSASA_cutoff", xsct_real, "Cutoff for buried solvent accessible surface area to ignore reporting most values", "100" )
		+ XMLSchemaAttribute::attribute_w_default( "compute_packstat", xsct_rosetta_bool, "Compute packstat score for interface?", "true" )
		+ XMLSchemaAttribute( "scorefxn", xs_string, "Score function to use when evaluating interface")
		+ XMLSchemaAttribute( "interface", "interface_features_string", "Specify the interface as ChainsSide1_ChainsSide2 (e.g. AB_C)" )
		+ XMLSchemaAttribute( "interfaces", "interface_features_list", "Provide a comma-separated list of interfaces." );

	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "FeaturesReporter wrapper for InterfaceAnalyzer", attlist );

}

std::string InterfaceFeaturesCreator::type_name() const {
	return InterfaceFeatures::class_name();
}

protocols::features::FeaturesReporterOP
InterfaceFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new InterfaceFeatures );
}

void InterfaceFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterfaceFeatures::provide_xml_schema( xsd );
}


} //features
} //protocols
