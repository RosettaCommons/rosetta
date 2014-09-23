// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <core/types.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/scoring/ScoringManager.fwd.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

#include <cppdb/frontend.h>

#include <boost/lexical_cast.hpp>

using std::string;
using utility::vector1;

namespace core {
namespace scoring {
namespace methods {

EnergyMethodOptions::EnergyMethodOptions():
	// hard-wired default, but you can set this with etable_type( string )
	etable_type_(FA_STANDARD_DEFAULT),
	analytic_etable_evaluation_( false ),
	atom_vdw_atom_type_set_name_(chemical::CENTROID), // can be set, see below
	unfolded_energies_type_( UNFOLDED_SCORE12 ),
	exclude_protein_protein_fa_elec_(false), // rosetta++ defaulted to true!
	exclude_monomer_fa_elec_(false),
	elec_max_dis_(5.5),
	elec_min_dis_(1.5),
	elec_die_(10.0),
	elec_no_dis_dep_die_(false),
	smooth_fa_elec_( false ),
	exclude_DNA_DNA_(true), // rosetta++ default
	intrares_elec_correction_scale_( 0.0 ),
	hbond_options_(hbonds::HBondOptionsOP( new hbonds::HBondOptions() )),
	etable_options_(core::scoring::etable::EtableOptionsOP( new core::scoring::etable::EtableOptions() )),
	cst_max_seq_sep_( std::numeric_limits<core::Size>::max() ),
	cartbonded_len_(-1.0),
	cartbonded_ang_(-1.0),
	cartbonded_tors_(-1.0),
	cartbonded_proton_(-1.0),
	cartbonded_improper_(-1.0),
	cartbonded_linear_(false),
	pb_bound_tag_("bound"),
	pb_unbound_tag_("unbound"),
	bond_angle_residue_type_param_set_(/* NULL */)
{
	initialize_from_options();
}

void EnergyMethodOptions::initialize_from_options() {
	analytic_etable_evaluation_ = basic::options::option[ basic::options::OptionKeys::score::analytic_etable_evaluation ];
	elec_max_dis_ = basic::options::option[basic::options::OptionKeys::score::elec_max_dis ]();
	elec_min_dis_ = basic::options::option[basic::options::OptionKeys::score::elec_min_dis ]();
	elec_die_ = basic::options::option[ basic::options::OptionKeys::score::elec_die ]();
	elec_no_dis_dep_die_ = basic::options::option[ basic::options::OptionKeys::score::elec_r_option ]();
	smooth_fa_elec_ = basic::options::option[ basic::options::OptionKeys::score::smooth_fa_elec ]();
	exclude_DNA_DNA_ = basic::options::option[basic::options::OptionKeys::dna::specificity::exclude_dna_dna]; // adding because this parameter should absolutely be false for any structure with DNA in it and it doesn't seem to be read in via the weights file method, so now it's an option - sthyme
	intrares_elec_correction_scale_ = basic::options::option[ basic::options::OptionKeys::score::intrares_elec_correction_scale ]();
}

/// copy constructor
EnergyMethodOptions::EnergyMethodOptions(EnergyMethodOptions const & src)
	: ReferenceCount( src )
{
	*this = src;
}

EnergyMethodOptions::~EnergyMethodOptions() {}

/// copy operator
EnergyMethodOptions const &
EnergyMethodOptions::operator=(EnergyMethodOptions const & src) {
	if ( this != &src ) {
		etable_type_ = src.etable_type_;
		analytic_etable_evaluation_ = src.analytic_etable_evaluation_;
		atom_vdw_atom_type_set_name_ = src.atom_vdw_atom_type_set_name_;
		unfolded_energies_type_ = src.unfolded_energies_type_;
		method_weights_ = src.method_weights_;
		ss_weights_ = src.ss_weights_;
		exclude_protein_protein_fa_elec_ = src.exclude_protein_protein_fa_elec_;
		exclude_monomer_fa_elec_ = src.exclude_monomer_fa_elec_;
		elec_max_dis_ = src.elec_max_dis_;
		elec_min_dis_ = src.elec_min_dis_;
		elec_die_ = src.elec_die_;
		elec_no_dis_dep_die_ = src.elec_no_dis_dep_die_;
		smooth_fa_elec_ = src.smooth_fa_elec_;
		exclude_DNA_DNA_ = src.exclude_DNA_DNA_;
		intrares_elec_correction_scale_ = src.intrares_elec_correction_scale_;
		hbond_options_ = hbonds::HBondOptionsOP( new hbonds::HBondOptions( *(src.hbond_options_) ) );
		etable_options_ = core::scoring::etable::EtableOptionsOP( new etable::EtableOptions( *(src.etable_options_) ) );
		cst_max_seq_sep_ = src.cst_max_seq_sep_;
		bond_angle_central_atoms_to_score_ = src.bond_angle_central_atoms_to_score_;
		bond_angle_residue_type_param_set_ = src.bond_angle_residue_type_param_set_;
		cartbonded_len_ = src.cartbonded_len_;
		cartbonded_ang_ = src.cartbonded_ang_;
		cartbonded_tors_ = src.cartbonded_tors_;
		cartbonded_proton_ = src.cartbonded_proton_;
		cartbonded_linear_ = src.cartbonded_linear_;
		pb_bound_tag_ = src.pb_bound_tag_;
		pb_unbound_tag_ = src.pb_unbound_tag_;
	}
	return *this;
}

string const &
EnergyMethodOptions::etable_type() const {
	return etable_type_;
}

void
EnergyMethodOptions::etable_type(string const & type ) {
	etable_type_ = type;
}

string const &
EnergyMethodOptions::unfolded_energies_type() const {
	return unfolded_energies_type_;
}

void
EnergyMethodOptions::unfolded_energies_type(string const & type ) {
	unfolded_energies_type_ = type;
}

bool
EnergyMethodOptions::exclude_protein_protein_fa_elec() const {
	return exclude_protein_protein_fa_elec_;
}

void
EnergyMethodOptions::exclude_protein_protein_fa_elec( bool const setting ) {
	exclude_protein_protein_fa_elec_ = setting;
}

bool
EnergyMethodOptions::exclude_monomer_fa_elec() const {
	return exclude_monomer_fa_elec_;
}

void
EnergyMethodOptions::exclude_monomer_fa_elec( bool const setting ) {
	exclude_monomer_fa_elec_ = setting;
}

core::Real
EnergyMethodOptions::elec_max_dis() const {
	return elec_max_dis_;
}

void
EnergyMethodOptions::elec_max_dis( core::Real const setting ) {
	elec_max_dis_ = setting;
}

core::Real
EnergyMethodOptions::elec_min_dis() const {
	return elec_min_dis_;
}

void
EnergyMethodOptions::elec_min_dis( core::Real const setting ) {
	elec_min_dis_ = setting;
}

core::Real
EnergyMethodOptions::elec_die() const {
	return elec_die_;
}

void
EnergyMethodOptions::elec_die( core::Real const setting ) {
	elec_die_ = setting;
}

bool
EnergyMethodOptions::elec_no_dis_dep_die() const {
	return elec_no_dis_dep_die_;
}

void
EnergyMethodOptions::elec_no_dis_dep_die( bool const setting ) {
	elec_no_dis_dep_die_ = setting;
}

bool
EnergyMethodOptions::smooth_fa_elec() const {
	return smooth_fa_elec_;
}

void
EnergyMethodOptions::smooth_fa_elec( bool setting )
{
	smooth_fa_elec_ = setting;
}


bool
EnergyMethodOptions::exclude_DNA_DNA() const {
	runtime_assert( hbond_options_->exclude_DNA_DNA() == exclude_DNA_DNA_ );
	return exclude_DNA_DNA_;
}

void
EnergyMethodOptions::exclude_DNA_DNA( bool const setting ) {
	exclude_DNA_DNA_ = setting;
	hbond_options_->exclude_DNA_DNA( setting );
}

core::Real
EnergyMethodOptions::intrares_elec_correction_scale() const {
	return intrares_elec_correction_scale_;
}

void
EnergyMethodOptions::intrares_elec_correction_scale( core::Real const setting ) {
  intrares_elec_correction_scale_ = setting;
}

hbonds::HBondOptions const &
EnergyMethodOptions::hbond_options() const {
	return *hbond_options_;
}

hbonds::HBondOptions &
EnergyMethodOptions::hbond_options() {
	return *hbond_options_;
}

void
EnergyMethodOptions::hbond_options( hbonds::HBondOptions const & opts ) {
 	hbond_options_ = hbonds::HBondOptionsOP( new hbonds::HBondOptions( opts ) );
}

etable::EtableOptions const &
EnergyMethodOptions::etable_options() const {
	return *etable_options_;
}

etable::EtableOptions &
EnergyMethodOptions::etable_options() {
	return *etable_options_;
}

void
EnergyMethodOptions::etable_options( etable::EtableOptions const & opts ) {
	etable_options_ = core::scoring::etable::EtableOptionsOP( new etable::EtableOptions( opts ) );
}

std::string const &
EnergyMethodOptions::pb_bound_tag() const {
	return pb_bound_tag_;
}
std::string &
EnergyMethodOptions::pb_bound_tag() {
	return pb_bound_tag_;
}
void
EnergyMethodOptions::pb_bound_tag( std::string const & tag ) {
	 pb_bound_tag_ = tag;
}
std::string const &
EnergyMethodOptions::pb_unbound_tag() const {
	return pb_unbound_tag_;
}
std::string &
EnergyMethodOptions::pb_unbound_tag() {
	return pb_unbound_tag_;
}
void
EnergyMethodOptions::pb_unbound_tag( std::string const & tag ) {
	pb_unbound_tag_ = tag;
}

/// @brief  This is used in the construction of the VDW_Energy's AtomVDW object
string const &
EnergyMethodOptions::atom_vdw_atom_type_set_name() const {
	return atom_vdw_atom_type_set_name_;
}

void
EnergyMethodOptions::atom_vdw_atom_type_set_name(string const & setting ) {
	atom_vdw_atom_type_set_name_ = setting;
}

void
EnergyMethodOptions::set_strand_strand_weights(
	int ss_lowstrand,
	int ss_cutoff
) {
	ss_weights_.set_ss_cutoff(ss_cutoff);
	ss_weights_.set_ss_lowstrand(ss_lowstrand);
}

SecondaryStructureWeights const &
EnergyMethodOptions::secondary_structure_weights() const {
	return ss_weights_;
}

SecondaryStructureWeights &
EnergyMethodOptions::secondary_structure_weights() {
	return ss_weights_;
}

bool
EnergyMethodOptions::has_method_weights( ScoreType const & type ) const {
	return ( method_weights_.find( type ) != method_weights_.end() );
}

vector1< Real > const &
EnergyMethodOptions::method_weights( ScoreType const & type ) const {
	MethodWeights::const_iterator it( method_weights_.find( type ) );
	if ( it == method_weights_.end() ) {
		utility_exit_with_message( "EnergyMethodOptions::method_weights do not exist: " +
			name_from_score_type(type));
	}
	return it->second;
}

void
EnergyMethodOptions::set_method_weights(
	ScoreType const & type,
	vector1< Real > const & wts)
{
	method_weights_[ type ] = wts;
}

core::Size
EnergyMethodOptions::cst_max_seq_sep() const {
	return cst_max_seq_sep_;
}

void
EnergyMethodOptions::cst_max_seq_sep( Size const setting ) {
	cst_max_seq_sep_ = setting;
}

/// deprecated
vector1<string> const &
EnergyMethodOptions::bond_angle_central_atoms_to_score() const {
	if (bond_angle_residue_type_param_set_) {
		return bond_angle_residue_type_param_set_->central_atoms_to_score();
	}
	return bond_angle_central_atoms_to_score_;
}

/// deprecated
void
EnergyMethodOptions::bond_angle_central_atoms_to_score(vector1<string> const & atom_names) {
	bond_angle_central_atoms_to_score_ = atom_names;
	if (bond_angle_residue_type_param_set_) {
		bond_angle_residue_type_param_set_->central_atoms_to_score(atom_names);
	}
}

core::scoring::mm::MMBondAngleResidueTypeParamSetOP
EnergyMethodOptions::bond_angle_residue_type_param_set() {
	return bond_angle_residue_type_param_set_;
}

core::scoring::mm::MMBondAngleResidueTypeParamSetCOP
EnergyMethodOptions::bond_angle_residue_type_param_set() const {
	return bond_angle_residue_type_param_set_;
}

void
EnergyMethodOptions::bond_angle_residue_type_param_set(core::scoring::mm::MMBondAngleResidueTypeParamSetOP param_set) {
	bond_angle_residue_type_param_set_ = param_set;
}

/// used inside ScoreFunctionInfo::operator==
bool
operator==( EnergyMethodOptions const & a, EnergyMethodOptions const & b ) {

	return ( ( a.etable_type_ == b.etable_type_ ) &&
		( a.analytic_etable_evaluation_ == b.analytic_etable_evaluation_ ) &&
		( a.atom_vdw_atom_type_set_name_ == b.atom_vdw_atom_type_set_name_ ) &&
		( a.unfolded_energies_type_ == b.unfolded_energies_type_ ) &&
		( a.method_weights_ == b.method_weights_ ) &&
		( a.ss_weights_ == b.ss_weights_ ) &&
		( a.exclude_protein_protein_fa_elec_ == b.exclude_protein_protein_fa_elec_ ) &&
		( a.exclude_monomer_fa_elec_ == b.exclude_monomer_fa_elec_ ) &&
		( a.elec_max_dis_ == b.elec_max_dis_ ) &&
		( a.elec_min_dis_ == b.elec_min_dis_ ) &&
		( a.elec_die_ == b.elec_die_ ) &&
		( a.elec_no_dis_dep_die_ == b.elec_no_dis_dep_die_ ) &&
		( a.smooth_fa_elec_ == b.smooth_fa_elec_ ) &&
		( a.exclude_DNA_DNA_ == b.exclude_DNA_DNA_ ) &&
		( * (a.hbond_options_) == * (b.hbond_options_) ) &&
		( * (a.etable_options_) == * (b.etable_options_) ) &&
		( a.intrares_elec_correction_scale_ == b.intrares_elec_correction_scale_ ) &&
		( a.cst_max_seq_sep_ == b.cst_max_seq_sep_ ) &&
		( a.cartbonded_len_ == b.cartbonded_len_ ) &&
		( a.cartbonded_ang_ == b.cartbonded_ang_ ) &&
		( a.cartbonded_tors_ == b.cartbonded_tors_ ) &&
		( a.cartbonded_proton_ == b.cartbonded_proton_ ) &&
		( a.cartbonded_linear_ == b.cartbonded_linear_ ) &&
		( a.bond_angle_central_atoms_to_score_ == b.bond_angle_central_atoms_to_score_ ) &&
		( a.bond_angle_residue_type_param_set_ == b.bond_angle_residue_type_param_set_ ) &&
		( a.pb_bound_tag_ == b.pb_bound_tag_ ) &&
					 ( a.pb_unbound_tag_ == b.pb_unbound_tag_ ) )
		;
}

/// used inside ScoreFunctionInfo::operator==
bool
operator!=( EnergyMethodOptions const & a, EnergyMethodOptions const & b ) {
	return !( a == b );
}

void
EnergyMethodOptions::show( std::ostream & out ) const {
	if ( etable_type_.size() ) out << "EnergyMethodOptions::show: etable_type: " << etable_type_ <<'\n';
	out << "analytic_etable_evaluation: " << analytic_etable_evaluation_ << '\n';
	for ( MethodWeights::const_iterator it=method_weights_.begin(), ite = method_weights_.end(); it != ite; ++it ) {
		out << "EnergyMethodOptions::show: method_weights: " << it->first;
		for ( Size i=1; i<= it->second.size(); ++i ) {
			out << ' ' << it->second[i];
		}
		out << '\n';
	}
	out << "EnergyMethodOptions::show: unfolded_energies_type: " << unfolded_energies_type_ << std::endl;
	out << "EnergyMethodOptions::show: atom_vdw_atom_type_set_name: " << atom_vdw_atom_type_set_name_ << std::endl;
	out << "EnergyMethodOptions::show: exclude_protein_protein_fa_elec: "
			<< (exclude_protein_protein_fa_elec_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: exclude_monomer_fa_elec: "
			<< (exclude_monomer_fa_elec_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: elec_max_dis: " << elec_max_dis_ << std::endl;
	out << "EnergyMethodOptions::show: elec_min_dis: " << elec_min_dis_ << std::endl;
	out << "EnergyMethodOptions::show: elec_die: " << elec_die_ << std::endl;
	out << "EnergyMethodOptions::show: elec_no_dis_dep_die: "
			<< (elec_no_dis_dep_die_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: smooth_fa_elec: " << ( smooth_fa_elec_ ? "true" : "false" ) << std::endl;
	out << "EnergyMethodOptions::show: exclude_DNA_DNA: "
			<< (exclude_DNA_DNA_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: cst_max_seq_sep: " << cst_max_seq_sep_ << std::endl;
	out << "EnergyMethodOptions::show: pb_bound_tag: " << pb_bound_tag_ << std::endl;
	out << "EnergyMethodOptions::show: pb_unbound_tag: " << pb_unbound_tag_ << std::endl;
	out << "EnergyMethodOptions::show: bond_angle_central_atoms_to_score:";
	if (bond_angle_residue_type_param_set_) {
		out << "setting ignored";
	} else {
		for ( Size i=1; i <= bond_angle_central_atoms_to_score_.size(); ++i ) {
			out << " \"" << bond_angle_central_atoms_to_score_[i] << "\"";
		}
	}
	out << std::endl;
	out << "EnergyMethodOptions::show: bond_angle_residue_type_param_set: "
			<< (bond_angle_residue_type_param_set_ ? "in use" : "none") << std::endl;
	if (bond_angle_residue_type_param_set_) {
		out << "  central_atoms_to_score:";
		if (!bond_angle_residue_type_param_set_->central_atoms_to_score().size()) out << "all";
		for ( Size i=1; i <= bond_angle_residue_type_param_set_->central_atoms_to_score().size(); ++i ) {
			out << " \"" << bond_angle_residue_type_param_set_->central_atoms_to_score()[i] << "\"";
		}
		out << std::endl;
		out << "  use_residue_type_theta0: "
				<< (bond_angle_residue_type_param_set_->use_residue_type_theta0() ? "true" : "false") << std::endl;
	}
	out << *hbond_options_;
}

std::ostream& operator<<(std::ostream & out, EnergyMethodOptions const & options) {
	options.show( out );
	return out;
}

void
EnergyMethodOptions::write_score_function_method_options_table_schema(
	utility::sql_database::sessionOP db_session
) {
	using namespace basic::database::schema_generator;

	Column batch_id("batch_id", DbDataTypeOP( new DbInteger() ), true);
	Column score_function_name("score_function_name", DbDataTypeOP( new DbText(255) ), true);
	Column option_key("option_key", DbDataTypeOP( new DbText(255) ), true);
	Column option_value("option_value", DbDataTypeOP( new DbText(255) ), true);

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(batch_id);
	pkey_cols.push_back(score_function_name);
	pkey_cols.push_back(option_key);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(batch_id);
	vector1< string > reference_columns;
	reference_columns.push_back("batch_id");
	ForeignKey foreign_key(foreign_key_columns, "batches", reference_columns, true);


	Schema table("score_function_method_options", PrimaryKey(pkey_cols));
	table.add_foreign_key(foreign_key);
	table.add_column(option_value);

	table.write(db_session);

}


void
EnergyMethodOptions::insert_score_function_method_options_rows(
	Size batch_id,
	std::string const & score_function_name,
	utility::sql_database::sessionOP db_session
) const {

	vector1< std::string > option_keys;
	vector1< std::string > option_values;
	if(etable_type_.size()){
		option_keys.push_back("etable_type");
		option_values.push_back(etable_type_);
	}
	option_keys.push_back("analytic_etable_evaluation");
	option_values.push_back(analytic_etable_evaluation_ ? "1" : "0");

	option_keys.push_back("atom_vdw_atom_type_set_name");
	option_values.push_back(atom_vdw_atom_type_set_name_);

	option_keys.push_back("unfolded_energies_type");
	option_values.push_back(unfolded_energies_type_);

	option_keys.push_back("exclude_protein_protein_fa_elec");
	option_values.push_back(exclude_protein_protein_fa_elec_ ? "1" : "0");

	option_keys.push_back("exclude_monomer_fa_elec");
	option_values.push_back(exclude_monomer_fa_elec_ ? "1" : "0");

	option_keys.push_back("elec_max_dis");
	option_values.push_back(boost::lexical_cast<std::string>(elec_max_dis_));

	option_keys.push_back("elec_min_dis");
	option_values.push_back(boost::lexical_cast<std::string>(elec_min_dis_));

	option_keys.push_back("elec_die");
	option_values.push_back(boost::lexical_cast<std::string>(elec_die_));

	option_keys.push_back("elec_no_dis_dep_die");
	option_values.push_back(elec_no_dis_dep_die_ ? "1" : "0");

	option_keys.push_back("exclude_DNA_DNA");
	option_values.push_back(exclude_DNA_DNA_ ? "1" : "0");

	option_keys.push_back("intrares_elec_correction_scale");
	option_values.push_back(boost::lexical_cast<std::string>(intrares_elec_correction_scale_));

	option_keys.push_back("cst_max_seq_sep");
	option_values.push_back(boost::lexical_cast<std::string>(cst_max_seq_sep_));

	option_keys.push_back("cartbonded_len");
	option_values.push_back(boost::lexical_cast<std::string>(cartbonded_len_));

	option_keys.push_back("cartbonded_ang");
	option_values.push_back(boost::lexical_cast<std::string>(cartbonded_ang_));

	option_keys.push_back("cartbonded_tors");
	option_values.push_back(boost::lexical_cast<std::string>(cartbonded_tors_));

	option_keys.push_back("cartbonded_proton");
	option_values.push_back(boost::lexical_cast<std::string>(cartbonded_proton_));

	option_keys.push_back("cartbonded_linear");
	option_values.push_back(cartbonded_linear_ ? "1" : "0");

	string statement_string;
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3:
		statement_string = "INSERT OR IGNORE INTO score_function_method_options (batch_id, score_function_name, option_key, option_value) VALUES (?,?,?,?);";
		break;
	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres:
		statement_string = "INSERT IGNORE INTO score_function_method_options (batch_id, score_function_name, option_key, option_value) VALUES (?,?,?,?);";
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}

	cppdb::statement stmt(
		basic::database::safely_prepare_statement(statement_string, db_session));

	for(Size i=1; i <= option_keys.size(); ++i){
		stmt.bind(1, batch_id);
		stmt.bind(2, score_function_name);
		stmt.bind(3, option_keys[i]);
		stmt.bind(4, option_values[i]);
		basic::database::safely_write_to_database(stmt);
	}
}



}
}
}
