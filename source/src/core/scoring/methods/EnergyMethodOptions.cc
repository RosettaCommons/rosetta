// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/scoring/methods/FreeDOF_Options.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/scoring/ScoringManager.fwd.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/unfolded_state.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

#include <cppdb/frontend.h>

#include <boost/lexical_cast.hpp>

using std::string;
using utility::vector1;

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace methods {

/// @details Delegating constructor using the global option collection
EnergyMethodOptions::EnergyMethodOptions():
	EnergyMethodOptions( basic::options::option )
{
}

///
EnergyMethodOptions::EnergyMethodOptions( utility::options::OptionCollection const & options ) :
	aa_composition_setup_files_(),
	netcharge_setup_files_(),
	aspartimide_penalty_value_(25.0),
	// hard-wired default, but you can set this with etable_type( string )
	atom_vdw_atom_type_set_name_(chemical::CENTROID), // can be set, see below
	unfolded_energies_type_( UNFOLDED_SCORE12 ),
	split_unfolded_label_type_(SPLIT_UNFOLDED_MM),
	split_unfolded_value_type_(SPLIT_UNFOLDED_BOLTZ),
	exclude_protein_protein_fa_elec_(false), // rosetta++ defaulted to true!
	exclude_RNA_RNA_fa_elec_(false),
	exclude_monomer_fa_elec_(false),
	elec_max_dis_(5.5),
	elec_min_dis_(/*1.5*/1.6),
	elec_die_(10.0),
	elec_no_dis_dep_die_(false),
	smooth_fa_elec_( /*false*/true ),
	elec_sigmoidal_die_(/*false*/true),
	elec_sigmoidal_D_(/*78.0*/80.0),
	elec_sigmoidal_D0_(/*2.0*/6.0),
	elec_sigmoidal_S_(/*0.36*/0.4),
	grpelec_fade_type_( "false" ),
	grpelec_fade_param1_( 1.0 ),
	grpelec_fade_param2_( 1.0 ),
	grpelec_fade_hbond_( true ),
	grp_cpfxn_( true ),
	elec_group_file_( "" ),
	grpelec_context_dependent_( false ),
	use_polarization_(true),
	use_gen_kirkwood_(true),
	protein_dielectric_( 1.0 ),
	water_dielectric_( 78.3 ),
	exclude_DNA_DNA_(true),
	exclude_intra_res_protein_(true), // rosetta++ default
	put_intra_into_total_(false ),
	geom_sol_interres_path_distance_cutoff_( 0 ), // rosetta++ default -- should be 4.
	geom_sol_intrares_path_distance_cutoff_( 7 ), // originally implemented for RNA base/phosphate.
	eval_intrares_elec_ST_only_( false ),
	envsmooth_zero_negatives_( false ),
	hbond_options_(hbonds::HBondOptionsOP( new hbonds::HBondOptions( options ) )),
	etable_options_(core::scoring::etable::EtableOptionsOP( new core::scoring::etable::EtableOptions( options ) )),
	rna_options_( new rna::RNA_EnergyMethodOptions() ),
	free_dof_options_( new methods::FreeDOF_Options() ),
	cst_max_seq_sep_( std::numeric_limits<core::Size>::max() ),
	cartbonded_len_(-1.0),
	cartbonded_ang_(-1.0),
	cartbonded_tors_(-1.0),
	cartbonded_proton_(-1.0),
	cartbonded_improper_(-1.0),
	cartbonded_linear_(false),
	pb_bound_tag_("bound"),
	pb_unbound_tag_("unbound"),
	ordered_wat_penalty_(1.5),
	ordered_pt_wat_penalty_(2.15),
	symmetric_gly_tables_(false),
	loop_close_use_6D_potential_(false),
	fa_stack_base_all_(false),

	voids_penalty_energy_containing_cones_cutoff_(6),
	voids_penalty_energy_cone_dotproduct_cutoff_(0.1),
	voids_penalty_energy_cone_distance_cutoff_(8.0),
	voids_penalty_energy_voxel_size_(0.5),
	voids_penalty_energy_voxel_grid_padding_(1.0),
	voids_penalty_energy_disabled_except_during_packing_(true),

	bond_angle_residue_type_param_set_(/* NULL */)
{
	initialize_from_options( options );
	method_weights_[ free_res ] = utility::vector1< Real >();
}

/// copy constructor
EnergyMethodOptions::EnergyMethodOptions(EnergyMethodOptions const & src)
: ReferenceCount( src )
{
	*this = src;
}

EnergyMethodOptions::~EnergyMethodOptions() = default;

/// copy operator
EnergyMethodOptions &
EnergyMethodOptions::operator = (EnergyMethodOptions const & src) {
	if ( this != &src ) {
		aa_composition_setup_files_ = src.aa_composition_setup_files_;
		netcharge_setup_files_ = src.netcharge_setup_files_;
		aspartimide_penalty_value_ = src.aspartimide_penalty_value_;
		atom_vdw_atom_type_set_name_ = src.atom_vdw_atom_type_set_name_;
		unfolded_energies_type_ = src.unfolded_energies_type_;
		split_unfolded_label_type_=src.split_unfolded_label_type_;
		split_unfolded_value_type_=src.split_unfolded_value_type_;
		method_weights_ = src.method_weights_;
		ss_weights_ = src.ss_weights_;
		exclude_protein_protein_fa_elec_ = src.exclude_protein_protein_fa_elec_;
		exclude_RNA_RNA_fa_elec_ = src.exclude_RNA_RNA_fa_elec_;
		exclude_monomer_fa_elec_ = src.exclude_monomer_fa_elec_;
		elec_max_dis_ = src.elec_max_dis_;
		elec_min_dis_ = src.elec_min_dis_;
		elec_die_ = src.elec_die_;
		elec_no_dis_dep_die_ = src.elec_no_dis_dep_die_;
		elec_sigmoidal_die_ = src.elec_sigmoidal_die_;
		elec_sigmoidal_D_ = src.elec_sigmoidal_D_;
		elec_sigmoidal_D0_ = src.elec_sigmoidal_D0_;
		elec_sigmoidal_S_ = src.elec_sigmoidal_S_;
		smooth_fa_elec_ = src.smooth_fa_elec_;
		grpelec_fade_type_ = src.grpelec_fade_type_;
		grpelec_fade_param1_ = src.grpelec_fade_param1_;
		grpelec_fade_param2_ = src.grpelec_fade_param2_;
		grpelec_fade_hbond_  = src.grpelec_fade_hbond_;
		grp_cpfxn_  = src.grp_cpfxn_;
		elec_group_file_ = src.elec_group_file_;
		grpelec_context_dependent_ = src.grpelec_context_dependent_;
		use_polarization_ = src.use_polarization_;
		use_gen_kirkwood_ = src.use_gen_kirkwood_;
		protein_dielectric_ = src.protein_dielectric_;
		water_dielectric_ = src.water_dielectric_;
		exclude_DNA_DNA_ = src.exclude_DNA_DNA_;
		exclude_intra_res_protein_ = src.exclude_intra_res_protein_;
		put_intra_into_total_ = src.put_intra_into_total_;
		geom_sol_interres_path_distance_cutoff_ = src.geom_sol_interres_path_distance_cutoff_;
		geom_sol_intrares_path_distance_cutoff_ = src.geom_sol_intrares_path_distance_cutoff_;
		eval_intrares_elec_ST_only_ = src.eval_intrares_elec_ST_only_;
		envsmooth_zero_negatives_ = src.envsmooth_zero_negatives_;
		hbond_options_ = hbonds::HBondOptionsOP( new hbonds::HBondOptions( *(src.hbond_options_) ) );
		etable_options_ = core::scoring::etable::EtableOptionsOP( new etable::EtableOptions( *(src.etable_options_) ) );
		rna_options_ = rna::RNA_EnergyMethodOptionsOP( new rna::RNA_EnergyMethodOptions( *(src.rna_options_) ) );
		free_dof_options_ = methods::FreeDOF_OptionsOP( new methods::FreeDOF_Options( *(src.free_dof_options_) ) );
		cst_max_seq_sep_ = src.cst_max_seq_sep_;
		bond_angle_central_atoms_to_score_ = src.bond_angle_central_atoms_to_score_;
		bond_angle_residue_type_param_set_ = src.bond_angle_residue_type_param_set_;
		cartbonded_len_ = src.cartbonded_len_;
		cartbonded_ang_ = src.cartbonded_ang_;
		cartbonded_tors_ = src.cartbonded_tors_;
		cartbonded_proton_ = src.cartbonded_proton_;
		cartbonded_improper_ = src.cartbonded_improper_;
		cartbonded_linear_ = src.cartbonded_linear_;
		pb_bound_tag_ = src.pb_bound_tag_;
		pb_unbound_tag_ = src.pb_unbound_tag_;
		fastdens_perres_weights_ = src.fastdens_perres_weights_;
		ordered_wat_penalty_ = src.ordered_wat_penalty_;
		ordered_pt_wat_penalty_ = src.ordered_pt_wat_penalty_;
		symmetric_gly_tables_ = src.symmetric_gly_tables_;
		loop_close_use_6D_potential_ = src.loop_close_use_6D_potential_;
		fa_stack_base_all_ = src.fa_stack_base_all_;
		voids_penalty_energy_containing_cones_cutoff_ = src.voids_penalty_energy_containing_cones_cutoff_;
		voids_penalty_energy_cone_dotproduct_cutoff_ = src.voids_penalty_energy_cone_dotproduct_cutoff_;
		voids_penalty_energy_cone_distance_cutoff_ = src.voids_penalty_energy_cone_distance_cutoff_;
		voids_penalty_energy_voxel_size_ = src.voids_penalty_energy_voxel_size_;
		voids_penalty_energy_voxel_grid_padding_ = src.voids_penalty_energy_voxel_grid_padding_;
		voids_penalty_energy_disabled_except_during_packing_ = src.voids_penalty_energy_disabled_except_during_packing_;
	}
	return *this;
}

void EnergyMethodOptions::initialize_from_options() {
	initialize_from_options( basic::options::option );
}

/// @details If you add an option to the EnergyMethodOptions that is initialized from the input option collection,
/// make sure you add that option key to the list_options_read function below.
void EnergyMethodOptions::initialize_from_options( utility::options::OptionCollection const & options ) {
	utility::vector1 < std::string > emptyvector;

	hbond_options_->initialize_from_options( options );
	etable_options_->initialize_from_options( options );

	aa_composition_setup_files_ = (options[ basic::options::OptionKeys::score::aa_composition_setup_file ].user() ? options[ basic::options::OptionKeys::score::aa_composition_setup_file ]() : emptyvector);
	netcharge_setup_files_ = ( options[ basic::options::OptionKeys::score::netcharge_setup_file ].user() ? options[ basic::options::OptionKeys::score::netcharge_setup_file ]() : emptyvector);
	aspartimide_penalty_value_ = options[ basic::options::OptionKeys::score::aspartimide_penalty_value ]();
	elec_max_dis_ = options[basic::options::OptionKeys::score::elec_max_dis ]();
	elec_min_dis_ = options[basic::options::OptionKeys::score::elec_min_dis ]();
	elec_die_ = options[ basic::options::OptionKeys::score::elec_die ]();
	elec_no_dis_dep_die_ = options[ basic::options::OptionKeys::score::elec_r_option ]();
	elec_sigmoidal_die_ = options[ basic::options::OptionKeys::score::elec_sigmoidal_die ]();
	elec_sigmoidal_D_ = options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D ]();
	elec_sigmoidal_D0_ = options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D0 ]();
	elec_sigmoidal_S_ = options[ basic::options::OptionKeys::score::elec_sigmoidal_die_S ]();
	smooth_fa_elec_ = options[ basic::options::OptionKeys::score::smooth_fa_elec ]();
	grpelec_fade_type_ = options[ basic::options::OptionKeys::score::grpelec_fade_type ]();
	grpelec_fade_param1_ = options[ basic::options::OptionKeys::score::grpelec_fade_param1 ]();
	grpelec_fade_param2_ = options[ basic::options::OptionKeys::score::grpelec_fade_param2 ]();
	grpelec_fade_hbond_ = options[ basic::options::OptionKeys::score::grpelec_fade_hbond ]();
	grp_cpfxn_ = options[ basic::options::OptionKeys::score::grp_cpfxn ]();
	elec_group_file_ = options[ basic::options::OptionKeys::score::elec_group_file ]();
	grpelec_context_dependent_ = options[ basic::options::OptionKeys::score::grpelec_context_dependent ]();
	use_polarization_= options[ basic::options::OptionKeys::score::use_polarization]();
	use_gen_kirkwood_= options[ basic::options::OptionKeys::score::use_gen_kirkwood]();
	protein_dielectric_= options[ basic::options::OptionKeys::score::protein_dielectric]();
	water_dielectric_= options[ basic::options::OptionKeys::score::water_dielectric]();
	exclude_DNA_DNA_ = options[ basic::options::OptionKeys::dna::specificity::exclude_dna_dna](); // adding because this parameter should absolutely be false for any structure with DNA in it and it doesn't seem to be read in via the weights file method, so now it's an option - sthyme
	exclude_intra_res_protein_ = !options[ basic::options::OptionKeys::score::include_intra_res_protein]();
	put_intra_into_total( options[ basic::options::OptionKeys::score::put_intra_into_total]() );
	geom_sol_interres_path_distance_cutoff_ = options[ basic::options::OptionKeys::score::geom_sol_interres_path_distance_cutoff]();
	geom_sol_intrares_path_distance_cutoff_ = options[ basic::options::OptionKeys::score::geom_sol_intrares_path_distance_cutoff]();
	eval_intrares_elec_ST_only_ = options[ basic::options::OptionKeys::score::eval_intrares_elec_ST_only]();
	envsmooth_zero_negatives_ = options[ basic::options::OptionKeys::score::envsmooth_zero_negatives ]();
	symmetric_gly_tables_ = options[ basic::options::OptionKeys::score::symmetric_gly_tables ]();
	loop_close_use_6D_potential_ = options[ basic::options::OptionKeys::score::loop_close::use_6D_potential ]();
	fa_stack_base_all_ = !options[ basic::options::OptionKeys::score::fa_stack_base_base_only ]();

	runtime_assert_string_msg( options[basic::options::OptionKeys::score::voids_penalty_energy_containing_cones_cutoff]() >= 0, "Error in EnergyMethodOptions::initialize_from_options: The -score:voids_penalty_energy_containing_cones_cutoff flag must be set to a non-negative integer." );
	voids_penalty_energy_containing_cones_cutoff_ = options[basic::options::OptionKeys::score::voids_penalty_energy_containing_cones_cutoff]();
	voids_penalty_energy_cone_dotproduct_cutoff_ = options[basic::options::OptionKeys::score::voids_penalty_energy_cone_dotproduct_cutoff]();
	voids_penalty_energy_cone_distance_cutoff_ = options[basic::options::OptionKeys::score::voids_penalty_energy_cone_distance_cutoff]();
	voids_penalty_energy_voxel_size_ = options[basic::options::OptionKeys::score::voids_penalty_energy_voxel_size]();
	voids_penalty_energy_voxel_grid_padding_ = options[basic::options::OptionKeys::score::voids_penalty_energy_voxel_grid_padding]();
	voids_penalty_energy_disabled_except_during_packing_ = options[basic::options::OptionKeys::score::voids_penalty_energy_disabled_except_during_packing]();

	// check to see if the unfolded state command line options are set by the user
	if ( options[ basic::options::OptionKeys::unfolded_state::unfolded_energies_file].user() ) {
		unfolded_energies_type_ = UNFOLDED_SPLIT_USER_DEFINED;
	}

	if ( options[ basic::options::OptionKeys::unfolded_state::split_unfolded_energies_file].user() ) {
		split_unfolded_label_type_ = SPLIT_UNFOLDED_USER_DEFINED;
		split_unfolded_value_type_ = SPLIT_UNFOLDED_USER_DEFINED;
	}

	if ( options[ basic::options::OptionKeys::edensity::sc_scaling ].user() ) {
		fastdens_perres_weights_.resize( core::chemical::num_canonical_aas, options[ basic::options::OptionKeys::edensity::sc_scaling ]() );
	}

	ordered_wat_penalty_ = options[ basic::options::OptionKeys::corrections::water::ordered_wat_penalty ]();
	ordered_pt_wat_penalty_ = options[ basic::options::OptionKeys::corrections::water::ordered_pt_wat_penalty ]();

}

/// @details If you add a read to an option in initialize_from_options above, you must
/// update this function.
void
EnergyMethodOptions::list_options_read( utility::options::OptionKeyList & read_options )
{
	hbonds::HBondOptions::list_options_read( read_options );
	etable::EtableOptions::list_options_read( read_options );
	read_options
		+ basic::options::OptionKeys::dna::specificity::exclude_dna_dna
		+ basic::options::OptionKeys::edensity::sc_scaling
		+ basic::options::OptionKeys::score::aa_composition_setup_file
		+ basic::options::OptionKeys::score::aspartimide_penalty_value
		+ basic::options::OptionKeys::score::elec_die
		+ basic::options::OptionKeys::score::elec_group_file
		+ basic::options::OptionKeys::score::elec_max_dis
		+ basic::options::OptionKeys::score::elec_min_dis
		+ basic::options::OptionKeys::score::elec_r_option
		+ basic::options::OptionKeys::score::elec_sigmoidal_die
		+ basic::options::OptionKeys::score::elec_sigmoidal_die_D
		+ basic::options::OptionKeys::score::elec_sigmoidal_die_D0
		+ basic::options::OptionKeys::score::elec_sigmoidal_die_S
		+ basic::options::OptionKeys::score::envsmooth_zero_negatives
		+ basic::options::OptionKeys::score::eval_intrares_elec_ST_only
		+ basic::options::OptionKeys::score::geom_sol_interres_path_distance_cutoff
		+ basic::options::OptionKeys::score::geom_sol_intrares_path_distance_cutoff
		+ basic::options::OptionKeys::score::grp_cpfxn
		+ basic::options::OptionKeys::score::grpelec_context_dependent
		+ basic::options::OptionKeys::score::grpelec_fade_hbond
		+ basic::options::OptionKeys::score::grpelec_fade_param1
		+ basic::options::OptionKeys::score::grpelec_fade_param2
		+ basic::options::OptionKeys::score::grpelec_fade_type
		+ basic::options::OptionKeys::score::include_intra_res_protein
		+ basic::options::OptionKeys::score::netcharge_setup_file
		+ basic::options::OptionKeys::score::protein_dielectric
		+ basic::options::OptionKeys::score::put_intra_into_total
		+ basic::options::OptionKeys::score::smooth_fa_elec
		+ basic::options::OptionKeys::score::symmetric_gly_tables
		+ basic::options::OptionKeys::score::loop_close::use_6D_potential
		+ basic::options::OptionKeys::score::fa_stack_base_base_only
		+ basic::options::OptionKeys::score::use_gen_kirkwood
		+ basic::options::OptionKeys::score::use_polarization
		+ basic::options::OptionKeys::score::voids_penalty_energy_containing_cones_cutoff
		+ basic::options::OptionKeys::score::voids_penalty_energy_cone_distance_cutoff
		+ basic::options::OptionKeys::score::voids_penalty_energy_cone_dotproduct_cutoff
		+ basic::options::OptionKeys::score::voids_penalty_energy_voxel_grid_padding
		+ basic::options::OptionKeys::score::voids_penalty_energy_voxel_size
		+ basic::options::OptionKeys::score::voids_penalty_energy_disabled_except_during_packing
		+ basic::options::OptionKeys::score::water_dielectric
		+ basic::options::OptionKeys::unfolded_state::split_unfolded_energies_file
		+ basic::options::OptionKeys::unfolded_state::unfolded_energies_file
		+ basic::options::OptionKeys::corrections::water::ordered_wat_penalty
		+ basic::options::OptionKeys::corrections::water::ordered_pt_wat_penalty;
}


string const &
EnergyMethodOptions::etable_type() const {
	return etable_options_->etable_type;
}

void
EnergyMethodOptions::etable_type(string const & type ) {
	etable_options_->etable_type = type;
}

bool
EnergyMethodOptions::analytic_etable_evaluation() const { return etable_options_->analytic_etable_evaluation; }

void
EnergyMethodOptions::analytic_etable_evaluation( bool setting ) { etable_options_->analytic_etable_evaluation = setting; }

string const &
EnergyMethodOptions::unfolded_energies_type() const {
	return unfolded_energies_type_;
}

void
EnergyMethodOptions::unfolded_energies_type(string const & type ) {
	unfolded_energies_type_ = type;
}

string const &
EnergyMethodOptions::split_unfolded_label_type() const {
	return split_unfolded_label_type_;
}

void
EnergyMethodOptions::split_unfolded_label_type(string const & label_type) {
	split_unfolded_label_type_=label_type;
}

string const &
EnergyMethodOptions::split_unfolded_value_type() const {
	return split_unfolded_value_type_;
}

void
EnergyMethodOptions::split_unfolded_value_type(string const & value_type) {
	split_unfolded_value_type_=value_type;
}

bool
EnergyMethodOptions::exclude_protein_protein_fa_elec() const {
	return exclude_protein_protein_fa_elec_;
}

bool
EnergyMethodOptions::exclude_RNA_RNA_fa_elec() const {
	return exclude_RNA_RNA_fa_elec_;
}

void
EnergyMethodOptions::exclude_protein_protein_fa_elec( bool const setting ) {
	exclude_protein_protein_fa_elec_ = setting;
}

void
EnergyMethodOptions::exclude_RNA_RNA_fa_elec( bool const setting ) {
	exclude_RNA_RNA_fa_elec_ = setting;
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
EnergyMethodOptions::elec_sigmoidal_die() const {
	return elec_sigmoidal_die_;
}

void
EnergyMethodOptions::elec_sigmoidal_die( bool const setting ) {
	elec_sigmoidal_die_ = setting;
}

void
EnergyMethodOptions::elec_sigmoidal_die_params(Real &D, Real &D0, Real &S) const {
	D = elec_sigmoidal_D_;
	D0 = elec_sigmoidal_D0_;
	S = elec_sigmoidal_S_;
}

void
EnergyMethodOptions::set_elec_sigmoidal_die_params( Real D, Real D0, Real S ) {
	elec_sigmoidal_D_ = D;
	elec_sigmoidal_D0_ = D0;
	elec_sigmoidal_S_ = S;
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

std::string
EnergyMethodOptions::grpelec_fade_type() const {
	return grpelec_fade_type_;
}

void
EnergyMethodOptions::grpelec_fade_type( std::string setting )
{
	grpelec_fade_type_ = setting;
}

core::Real
EnergyMethodOptions::grpelec_fade_param1() const {
	return grpelec_fade_param1_;
}

void
EnergyMethodOptions::grpelec_fade_param1( core::Real setting )
{
	grpelec_fade_param1_ = setting;
}

core::Real
EnergyMethodOptions::grpelec_fade_param2() const {
	return grpelec_fade_param2_;
}

void
EnergyMethodOptions::grpelec_fade_param2( core::Real setting )
{
	grpelec_fade_param2_ = setting;
}

bool
EnergyMethodOptions::grpelec_fade_hbond() const {
	return grpelec_fade_hbond_;
}

void
EnergyMethodOptions::grpelec_fade_hbond( bool setting )
{
	grpelec_fade_hbond_ = setting;
}

bool
EnergyMethodOptions::grp_cpfxn() const {
	return grp_cpfxn_;
}

void
EnergyMethodOptions::grp_cpfxn( bool setting )
{
	grp_cpfxn_ = setting;
}

std::string
EnergyMethodOptions::elec_group_file() const {
	return elec_group_file_;
}

void
EnergyMethodOptions::elec_group_file( std::string setting )
{
	elec_group_file_ = setting;
}

bool
EnergyMethodOptions::grpelec_context_dependent() const {
	return grpelec_context_dependent_;
}

void
EnergyMethodOptions::grpelec_context_dependent( bool setting )
{
	grpelec_context_dependent_ = setting;
}

bool
EnergyMethodOptions::use_polarization() const {
	return use_polarization_;
}

void
EnergyMethodOptions::use_polarization( bool const setting ) {
	use_polarization_ = setting;
}

bool
EnergyMethodOptions::use_gen_kirkwood() const {
	return use_gen_kirkwood_;
}

void
EnergyMethodOptions::use_gen_kirkwood( bool const setting ) {
	use_gen_kirkwood_ = setting;
}

Real
EnergyMethodOptions::protein_dielectric() const {
	return protein_dielectric_;
}

void
EnergyMethodOptions::protein_dielectric( Real const setting ) {
	protein_dielectric_ = setting;
}

Real
EnergyMethodOptions::water_dielectric() const {
	return water_dielectric_;
}

void
EnergyMethodOptions::water_dielectric( Real const setting ) {
	water_dielectric_ = setting;
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

bool
EnergyMethodOptions::exclude_intra_res_protein() const {
	runtime_assert( hbond_options_->exclude_intra_res_protein() == exclude_intra_res_protein_ );
	return exclude_intra_res_protein_;
}

void
EnergyMethodOptions::exclude_intra_res_protein( bool const setting ) {
	exclude_intra_res_protein_ = setting;
	hbond_options_->exclude_intra_res_protein( setting );
}

bool
EnergyMethodOptions::put_intra_into_total() const {
	runtime_assert( hbond_options_->put_intra_into_total() == put_intra_into_total_ );
	return put_intra_into_total_;
}

void
EnergyMethodOptions::put_intra_into_total( bool const setting ) {
	put_intra_into_total_ = setting;
	hbond_options_->put_intra_into_total( setting );
}

core::Size
EnergyMethodOptions::geom_sol_interres_path_distance_cutoff() const {
	return geom_sol_interres_path_distance_cutoff_;
}

void
EnergyMethodOptions::geom_sol_interres_path_distance_cutoff( core::Size const setting ) {
	geom_sol_interres_path_distance_cutoff_ = setting;
}

core::Size
EnergyMethodOptions::geom_sol_intrares_path_distance_cutoff() const {
	return geom_sol_intrares_path_distance_cutoff_;
}

void
EnergyMethodOptions::geom_sol_intrares_path_distance_cutoff( core::Size const setting ) {
	geom_sol_intrares_path_distance_cutoff_ = setting;
}

bool
EnergyMethodOptions::eval_intrares_elec_ST_only() const {
	return eval_intrares_elec_ST_only_;
}

void
EnergyMethodOptions::eval_intrsres_elec_ST_only( bool setting ) {
	eval_intrares_elec_ST_only_ = setting;
}

bool
EnergyMethodOptions::envsmooth_zero_negatives() const
{
	return envsmooth_zero_negatives_;
}

void
EnergyMethodOptions::envsmooth_zero_negatives( bool const setting )
{
	envsmooth_zero_negatives_ = setting;
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

rna::RNA_EnergyMethodOptions const &
EnergyMethodOptions::rna_options() const {
	return *rna_options_;
}

rna::RNA_EnergyMethodOptions &
EnergyMethodOptions::rna_options() {
	return *rna_options_;
}

void
EnergyMethodOptions::rna_options( rna::RNA_EnergyMethodOptions const & opts ) {
	rna_options_ = rna::RNA_EnergyMethodOptionsOP( new rna::RNA_EnergyMethodOptions( opts ) );
}

methods::FreeDOF_Options const &
EnergyMethodOptions::free_dof_options() const {
	return *free_dof_options_;
}

methods::FreeDOF_Options &
EnergyMethodOptions::free_dof_options() {
	return *free_dof_options_;
}

void
EnergyMethodOptions::free_dof_options( methods::FreeDOF_Options const & opts ) {
	free_dof_options_ = methods::FreeDOF_OptionsOP( new methods::FreeDOF_Options( opts ) );
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

/// @brief Should glyceine's Ramachandran and P_AA_PP tables be symmetrized (e.g. for scoring in a mixed D/L context)?
/// @details Default false.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
EnergyMethodOptions::symmetric_gly_tables() const {
	return symmetric_gly_tables_;
}

/// @brief Set whether glyceine's Ramachandran and P_AA_PP tables should be symmetrized (e.g. for scoring in a mixed D/L context).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
EnergyMethodOptions::symmetric_gly_tables( bool const setting ) {
	symmetric_gly_tables_ = setting;
}

/// @brief Get the number of cones in which a voxel must lie in order for that voxel to be considered
/// to be "buried".
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::Size EnergyMethodOptions::voids_penalty_energy_containing_cones_cutoff() const { return voids_penalty_energy_containing_cones_cutoff_; }

/// @brief Get the cone distance cutoff for the voids penalty energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::Real EnergyMethodOptions::voids_penalty_energy_cone_distance_cutoff() const { return voids_penalty_energy_cone_distance_cutoff_; }

/// @brief Get the cone dot product cutoff for the voids penalty energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::Real EnergyMethodOptions::voids_penalty_energy_cone_dotproduct_cutoff() const { return voids_penalty_energy_cone_dotproduct_cutoff_; }

/// @brief Get the voxel grid padding for the voids penalty energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::Real EnergyMethodOptions::voids_penalty_energy_voxel_grid_padding() const { return voids_penalty_energy_voxel_grid_padding_; }

/// @brief Get the voxel size for the voids penalty energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::Real EnergyMethodOptions::voids_penalty_energy_voxel_size() const { return voids_penalty_energy_voxel_size_; }

/// @brief Get whether we're prohibiting evaluation of the voids_penalty score term outside
/// of the context of the packer.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool EnergyMethodOptions::voids_penalty_energy_disabled_except_during_packing() const { return voids_penalty_energy_disabled_except_during_packing_; }

/// @brief Set the number of cones in which a voxel must lie in order for that voxel to be considered
/// to be "buried".
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void EnergyMethodOptions::voids_penalty_energy_containing_cones_cutoff( core::Size const setting ) {
	voids_penalty_energy_containing_cones_cutoff_ = setting;
}

/// @brief Set the cone distance cutoff for the voids penalty energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void EnergyMethodOptions::voids_penalty_energy_cone_distance_cutoff( core::Real const &setting ) {
	runtime_assert_string_msg( setting > 0, "Error in EnergyMethodOptions::voids_penalty_energy_cone_distance_cutoff(): The setting must be positive!" );
	voids_penalty_energy_cone_distance_cutoff_ = setting;
}

/// @brief Set the cone dot product cutoff for the voids penalty energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void EnergyMethodOptions::voids_penalty_energy_cone_dotproduct_cutoff( core::Real const &setting ) {
	runtime_assert_string_msg( -1.0 <= setting && setting <= 1.0, "Error in EnergyMethodOptions::voids_penalty_energy_cone_dotproduct_cutoff(): The setting must be between -1.0 and 1.0." );
	voids_penalty_energy_cone_dotproduct_cutoff_ = setting;
}

/// @brief Set the voxel grid padding for the voids penalty energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void EnergyMethodOptions::voids_penalty_energy_voxel_grid_padding( core::Real const &setting ) {
	runtime_assert_string_msg( setting > 0, "Error in EnergyMethodOptions::voids_penalty_energy_voxel_grid_padding(): The setting must be positive!" );
	voids_penalty_energy_voxel_grid_padding_ = setting;
}

/// @brief Set the voxel size for the voids penalty energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void EnergyMethodOptions::voids_penalty_energy_voxel_size( core::Real const &setting ) {
	runtime_assert_string_msg( setting > 0, "Error in EnergyMethodOptions::voids_penalty_energy_voxel_size(): The setting must be positive!" );
	voids_penalty_energy_voxel_size_ = setting;
}

/// @brief Set whether we're prohibiting evaluation of the voids_penalty score term outside
/// of the context of the packer.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void EnergyMethodOptions::voids_penalty_energy_disabled_except_during_packing( bool const setting ) { voids_penalty_energy_disabled_except_during_packing_ = setting; }

bool
EnergyMethodOptions::loop_close_use_6D_potential() const {
	return loop_close_use_6D_potential_;
}

void
EnergyMethodOptions::loop_close_use_6D_potential( bool const setting ) {
	loop_close_use_6D_potential_ = setting;
}

bool
EnergyMethodOptions::fa_stack_base_all() const {
	return fa_stack_base_all_;
}

void
EnergyMethodOptions::fa_stack_base_all( bool const setting ) {
	fa_stack_base_all_ = setting;
}


utility::vector1< core::Real > const &
EnergyMethodOptions::get_density_sc_scale_byres() const {
	return fastdens_perres_weights_;
}

void
EnergyMethodOptions::set_density_sc_scale_byres(core::chemical::AA aa, core::Real newscscale) {
	runtime_assert ( aa <= core::chemical::num_canonical_aas );
	if ( fastdens_perres_weights_.size() == 0 ) {
		fastdens_perres_weights_.resize( core::chemical::num_canonical_aas, 1.0 );
	}
	fastdens_perres_weights_[(int)aa] = newscscale;
}

void
EnergyMethodOptions::set_density_sc_scale_byres(core::Real newscscale) {
	fastdens_perres_weights_ = utility::vector1< core::Real >( core::chemical::num_canonical_aas , newscscale );
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
	auto it( method_weights_.find( type ) );
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
	if ( bond_angle_residue_type_param_set_ ) {
		return bond_angle_residue_type_param_set_->central_atoms_to_score();
	}
	return bond_angle_central_atoms_to_score_;
}

/// deprecated
void
EnergyMethodOptions::bond_angle_central_atoms_to_score(vector1<string> const & atom_names) {
	bond_angle_central_atoms_to_score_ = atom_names;
	if ( bond_angle_residue_type_param_set_ ) {
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

core::Real
EnergyMethodOptions::ordered_wat_penalty() const {
	return ordered_wat_penalty_;
}

core::Real
EnergyMethodOptions::ordered_pt_wat_penalty() const {
	return ordered_pt_wat_penalty_;
}

/// used inside ScoreFunctionInfo::operator==
bool
operator==( EnergyMethodOptions const & a, EnergyMethodOptions const & b ) {

	return (
		( a.aa_composition_setup_files_ == b.aa_composition_setup_files_ ) &&
		( a.netcharge_setup_files_ == b.netcharge_setup_files_ ) &&
		( a.aspartimide_penalty_value_ == b.aspartimide_penalty_value_ ) &&
		( a.atom_vdw_atom_type_set_name_ == b.atom_vdw_atom_type_set_name_ ) &&
		( a.unfolded_energies_type_ == b.unfolded_energies_type_ ) &&
		( a.split_unfolded_label_type_ == b.split_unfolded_label_type_ ) &&
		( a.split_unfolded_value_type_ == b.split_unfolded_value_type_ ) &&
		( a.method_weights_ == b.method_weights_ ) &&
		( a.ss_weights_ == b.ss_weights_ ) &&
		( a.exclude_protein_protein_fa_elec_ == b.exclude_protein_protein_fa_elec_ ) &&
		( a.exclude_RNA_RNA_fa_elec_ == b.exclude_RNA_RNA_fa_elec_ ) &&
		( a.exclude_monomer_fa_elec_ == b.exclude_monomer_fa_elec_ ) &&
		( a.elec_max_dis_ == b.elec_max_dis_ ) &&
		( a.elec_min_dis_ == b.elec_min_dis_ ) &&
		( a.elec_die_ == b.elec_die_ ) &&
		( a.elec_no_dis_dep_die_ == b.elec_no_dis_dep_die_ ) &&
		( a.elec_sigmoidal_die_ == b.elec_sigmoidal_die_ ) &&
		( a.elec_sigmoidal_D_ == b.elec_sigmoidal_D_ ) &&
		( a.elec_sigmoidal_D0_ == b.elec_sigmoidal_D0_ ) &&
		( a.elec_sigmoidal_S_ == b.elec_sigmoidal_S_ ) &&
		( a.smooth_fa_elec_ == b.smooth_fa_elec_ ) &&
		( a.grpelec_fade_type_ == b.grpelec_fade_type_ ) &&
		( a.grpelec_fade_param1_ == b.grpelec_fade_param1_ ) &&
		( a.grpelec_fade_param2_ == b.grpelec_fade_param2_ ) &&
		( a.grpelec_fade_hbond_ == b.grpelec_fade_hbond_ ) &&
		( a.grp_cpfxn_ == b.grp_cpfxn_ ) &&
		( a.elec_group_file_ == b.elec_group_file_ ) &&
		( a.grpelec_context_dependent_ == b.grpelec_context_dependent_ ) &&
		( a.use_polarization_ == b.use_polarization_ ) &&
		( a.use_gen_kirkwood_ == b.use_gen_kirkwood_ ) &&
		( a.protein_dielectric_ == b.protein_dielectric_ ) &&
		( a.water_dielectric_ == b.water_dielectric_ ) &&
		( a.exclude_DNA_DNA_ == b.exclude_DNA_DNA_ ) &&
		( a.exclude_intra_res_protein_ == b.exclude_intra_res_protein_ ) &&
		( a.put_intra_into_total_ == b.put_intra_into_total_ ) &&
		( a.geom_sol_interres_path_distance_cutoff_ == b.geom_sol_interres_path_distance_cutoff_ ) &&
		( a.geom_sol_intrares_path_distance_cutoff_ == b.geom_sol_intrares_path_distance_cutoff_ ) &&
		( a.eval_intrares_elec_ST_only_ == b.eval_intrares_elec_ST_only_ ) &&
		( a.envsmooth_zero_negatives_ == b.envsmooth_zero_negatives_ ) &&
		( * (a.hbond_options_) == * (b.hbond_options_) ) &&
		( * (a.etable_options_) == * (b.etable_options_) ) &&
		( * (a.rna_options_) == * (b.rna_options_) ) &&
		( * (a.free_dof_options_) == * (b.free_dof_options_) ) &&
		( a.cst_max_seq_sep_ == b.cst_max_seq_sep_ ) &&
		( a.cartbonded_len_ == b.cartbonded_len_ ) &&
		( a.cartbonded_ang_ == b.cartbonded_ang_ ) &&
		( a.cartbonded_tors_ == b.cartbonded_tors_ ) &&
		( a.cartbonded_proton_ == b.cartbonded_proton_ ) &&
		( a.cartbonded_improper_ == b.cartbonded_improper_ ) &&
		( a.cartbonded_linear_ == b.cartbonded_linear_ ) &&
		( a.bond_angle_central_atoms_to_score_ == b.bond_angle_central_atoms_to_score_ ) &&
		( a.bond_angle_residue_type_param_set_ == b.bond_angle_residue_type_param_set_ ) &&
		( a.pb_bound_tag_ == b.pb_bound_tag_ ) &&
		( a.pb_unbound_tag_ == b.pb_unbound_tag_ ) &&
		( a.ordered_wat_penalty_ == b.ordered_wat_penalty_ ) &&
		( a.ordered_pt_wat_penalty_ == b.ordered_pt_wat_penalty_ ) &&
		( a.voids_penalty_energy_containing_cones_cutoff_ == b.voids_penalty_energy_containing_cones_cutoff_ ) &&
		( a.voids_penalty_energy_cone_distance_cutoff_ == b.voids_penalty_energy_cone_distance_cutoff_ ) &&
		( a.voids_penalty_energy_cone_dotproduct_cutoff_ == b.voids_penalty_energy_cone_dotproduct_cutoff_ ) &&
		( a.voids_penalty_energy_voxel_grid_padding_ == b.voids_penalty_energy_voxel_grid_padding_ ) &&
		( a.voids_penalty_energy_voxel_size_ == b.voids_penalty_energy_voxel_size_ ) &&
		( a.voids_penalty_energy_disabled_except_during_packing_ == b.voids_penalty_energy_disabled_except_during_packing_ ) &&
		( a.fastdens_perres_weights_ == b.fastdens_perres_weights_ ) );
}

/// used inside ScoreFunctionInfo::operator==
bool
operator!=( EnergyMethodOptions const & a, EnergyMethodOptions const & b ) {
	return !( a == b );
}

void
EnergyMethodOptions::show( std::ostream & out ) const {
	out << "EnergyMethodOptions::show: aa_composition_setup_files: ";
	for ( core::Size i=1, imax=aa_composition_setup_file_count(); i<=imax; ++i ) {
		out << aa_composition_setup_file(i);
		if ( i<imax ) {
			out << ", ";
		} else {
			out << ".";
		}
	}
	out << std::endl;

	out << "EnergyMethodOptions::show: netcharge_setup_files: ";
	for ( core::Size i(1), imax(netcharge_setup_file_count()); i<=imax; ++i ) {
		out << netcharge_setup_file(i);
		if ( i < imax ) out << ", ";
		else out << ".";
	}
	out << std::endl;

	out << "EnergyMethodOptions::show: aspartimide_penalty_value: " << aspartimide_penalty_value() << std::endl;

	if ( etable_options_->etable_type.size() ) out << "EnergyMethodOptions::show: etable_type: " << etable_options_->etable_type <<'\n';
	out << "analytic_etable_evaluation: " << etable_options_->analytic_etable_evaluation << '\n';
	for ( auto const & method_weight : method_weights_ ) {
		out << "EnergyMethodOptions::show: method_weights: " << method_weight.first;
		for ( Size i=1; i<= method_weight.second.size(); ++i ) {
			out << ' ' << method_weight.second[i];
		}
		out << '\n';
	}
	out << "EnergyMethodOptions::show: unfolded_energies_type: " << unfolded_energies_type_ << std::endl;
	out << "EnergyMethodOptions::show: split_unfolded_label_type: " << split_unfolded_label_type_ << std::endl;
	out << "EnergyMethodOptions::show: split_unfolded_value_type: " << split_unfolded_value_type_ << std::endl;
	out << "EnergyMethodOptions::show: atom_vdw_atom_type_set_name: " << atom_vdw_atom_type_set_name_ << std::endl;
	out << "EnergyMethodOptions::show: exclude_protein_protein_fa_elec: "
		<< (exclude_protein_protein_fa_elec_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: exclude_RNA_RNA_fa_elec: "
		<< (exclude_RNA_RNA_fa_elec_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: exclude_monomer_fa_elec: "
		<< (exclude_monomer_fa_elec_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: elec_max_dis: " << elec_max_dis_ << std::endl;
	out << "EnergyMethodOptions::show: elec_min_dis: " << elec_min_dis_ << std::endl;
	out << "EnergyMethodOptions::show: elec_die: " << elec_die_ << std::endl;
	out << "EnergyMethodOptions::show: elec_no_dis_dep_die: "
		<< (elec_no_dis_dep_die_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: elec_sigmoidal_die: "
		<< (elec_sigmoidal_die_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: elec_sigmoidal_D: " << elec_sigmoidal_D_ << std::endl;
	out << "EnergyMethodOptions::show: elec_sigmoidal_D0: " << elec_sigmoidal_D0_ << std::endl;
	out << "EnergyMethodOptions::show: elec_sigmoidal_S: " << elec_sigmoidal_S_ << std::endl;
	out << "EnergyMethodOptions::show: smooth_fa_elec: " << ( smooth_fa_elec_ ? "true" : "false" ) << std::endl;
	out << "EnergyMethodOptions::show: grpelec_fade_type: " << grpelec_fade_type_ << std::endl;
	out << "EnergyMethodOptions::show: grpelec_fade_param1: " << grpelec_fade_param1_ << std::endl;
	out << "EnergyMethodOptions::show: grpelec_fade_param2: " << grpelec_fade_param2_ << std::endl;
	out << "EnergyMethodOptions::show: grpelec_fade_hbond: " << grpelec_fade_hbond_ << std::endl;
	out << "EnergyMethodOptions::show: grp_cpfxn: " << grp_cpfxn_ << std::endl;
	out << "EnergyMethodOptions::show: elec_group_file: " << elec_group_file_ << std::endl;
	out << "EnergyMethodOptions::show: grpelec_context_dependent: " << grpelec_context_dependent_ << std::endl;
	out << "EnergyMethodOptions::show: use_polarization: "
		<< (use_polarization_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: use_gen_kirkwood: "
		<< (use_gen_kirkwood_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: protein_dielectric: " << protein_dielectric_ << std::endl;
	out << "EnergyMethodOptions::show: water_dielectric: " << water_dielectric_ << std::endl;
	out << "EnergyMethodOptions::show: exclude_DNA_DNA: "
		<< (exclude_DNA_DNA_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: exclude_intra_res_protein: "
		<< (exclude_intra_res_protein_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: put_intra_into_total: "
		<< (put_intra_into_total_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: geom_sol_interres_path_distance_cutoff: "
		<< (geom_sol_interres_path_distance_cutoff_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: geom_sol_intrares_path_distance_cutoff: "
		<< (geom_sol_intrares_path_distance_cutoff_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: eval_intrares_elec_ST_only: "
		<< ( eval_intrares_elec_ST_only_ ? "true" : "false" ) << std::endl;
	out << "EnergyMethodOptions::show: envsmooth_zero_negatives: "
		<< (envsmooth_zero_negatives_ ? "true" : "false") << std::endl;
	out << "EnergyMethodOptions::show: cst_max_seq_sep: " << cst_max_seq_sep_ << std::endl;
	out << "EnergyMethodOptions::show: pb_bound_tag: " << pb_bound_tag_ << std::endl;
	out << "EnergyMethodOptions::show: pb_unbound_tag: " << pb_unbound_tag_ << std::endl;
	out << "EnergyMethodOptions::show: ordered_wat_penalty: " << ordered_wat_penalty_ << std::endl;
	out << "EnergyMethodOptions::show: ordered_pt_wat_penalty: " << ordered_pt_wat_penalty_ << std::endl;
	out << "EnergyMethodOptions::show: voids_penalty_energy_containing_cones_cutoff_:" << voids_penalty_energy_containing_cones_cutoff_ << std::endl;
	out << "EnergyMethodOptions::show: voids_penalty_energy_cone_distance_cutoff_: " << voids_penalty_energy_cone_distance_cutoff_ << std::endl;
	out << "EnergyMethodOptions::show: voids_penalty_energy_cone_dotproduct_cutoff_: " << voids_penalty_energy_cone_dotproduct_cutoff_ << std::endl;
	out << "EnergyMethodOptions::show: voids_penalty_energy_voxel_grid_padding_: " << voids_penalty_energy_voxel_grid_padding_ << std::endl;
	out << "EnergyMethodOptions::show: voids_penalty_energy_voxel_size_: " << voids_penalty_energy_voxel_size_<< std::endl;
	out << "EnergyMethodOptions::show: voids_penalty_energy_disabled_except_during_packing_: " << (voids_penalty_energy_disabled_except_during_packing_ ? "TRUE" : "FALSE")  << std::endl;
	out << "EnergyMethodOptions::show: bond_angle_central_atoms_to_score:";
	if ( bond_angle_residue_type_param_set_ ) {
		out << "setting ignored";
	} else {
		for ( Size i=1; i <= bond_angle_central_atoms_to_score_.size(); ++i ) {
			out << " \"" << bond_angle_central_atoms_to_score_[i] << "\"";
		}
	}
	out << std::endl;
	out << "EnergyMethodOptions::show: bond_angle_residue_type_param_set: "
		<< (bond_angle_residue_type_param_set_ ? "in use" : "none") << std::endl;
	if ( bond_angle_residue_type_param_set_ ) {
		out << "  central_atoms_to_score:";
		if ( !bond_angle_residue_type_param_set_->central_atoms_to_score().size() ) out << "all";
		for ( Size i=1; i <= bond_angle_residue_type_param_set_->central_atoms_to_score().size(); ++i ) {
			out << " \"" << bond_angle_residue_type_param_set_->central_atoms_to_score()[i] << "\"";
		}
		out << std::endl;
		out << "  use_residue_type_theta0: "
			<< (bond_angle_residue_type_param_set_->use_residue_type_theta0() ? "true" : "false") << std::endl;
	}
	out << *hbond_options_;
	out << *rna_options_;
	out << *free_dof_options_;
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
	if ( etable_options_->etable_type.size() ) {
		option_keys.push_back("etable_type");
		option_values.push_back(etable_options_->etable_type);
	}
	option_keys.push_back("analytic_etable_evaluation");
	option_values.push_back( etable_options_->analytic_etable_evaluation ? "1" : "0" );

	option_keys.push_back("atom_vdw_atom_type_set_name");
	option_values.push_back(atom_vdw_atom_type_set_name_);

	option_keys.push_back("unfolded_energies_type");
	option_values.push_back(unfolded_energies_type_);

	option_keys.push_back("split_unfolded_label_type");
	option_values.push_back(split_unfolded_label_type_);

	option_keys.push_back("split_unfolded_value_type");
	option_values.push_back(split_unfolded_value_type_);

	option_keys.push_back("exclude_protein_protein_fa_elec");
	option_values.push_back(exclude_protein_protein_fa_elec_ ? "1" : "0");

	option_keys.push_back("exclude_RNA_RNA_fa_elec");
	option_values.push_back(exclude_RNA_RNA_fa_elec_ ? "1" : "0");

	option_keys.push_back("exclude_monomer_fa_elec");
	option_values.push_back(exclude_monomer_fa_elec_ ? "1" : "0");

	option_keys.push_back("elec_max_dis");
	option_values.push_back(boost::lexical_cast<std::string>(elec_max_dis_));

	option_keys.push_back("elec_min_dis");
	option_values.push_back(boost::lexical_cast<std::string>(elec_min_dis_));

	option_keys.push_back("elec_die");
	option_values.push_back(boost::lexical_cast<std::string>(elec_die_));

	option_keys.push_back("elec_sigmoidal_die");
	option_values.push_back(elec_sigmoidal_die_ ? "1" : "0");

	option_keys.push_back("elec_sigmoidal_D");
	option_values.push_back(boost::lexical_cast<std::string>(elec_sigmoidal_D_));

	option_keys.push_back("elec_sigmoidal_D0");
	option_values.push_back(boost::lexical_cast<std::string>(elec_sigmoidal_D0_));

	option_keys.push_back("elec_sigmoidal_S");
	option_values.push_back(boost::lexical_cast<std::string>(elec_sigmoidal_S_));

	option_keys.push_back("elec_no_dis_dep_die");
	option_values.push_back(elec_no_dis_dep_die_ ? "1" : "0");

	option_keys.push_back("use_polarization");
	option_values.push_back(use_polarization_ ? "1" : "0");

	option_keys.push_back("use_gen_kirkwood");
	option_values.push_back(use_gen_kirkwood_ ? "1" : "0");

	option_keys.push_back("protein_dielectric");
	option_values.push_back(boost::lexical_cast<std::string>(protein_dielectric_));

	option_keys.push_back("water_dielectric");
	option_values.push_back(boost::lexical_cast<std::string>(water_dielectric_));

	option_keys.push_back("exclude_DNA_DNA");
	option_values.push_back(exclude_DNA_DNA_ ? "1" : "0");

	option_keys.push_back("exclude_intra_res_protein");
	option_values.push_back(exclude_intra_res_protein_ ? "1" : "0");

	option_keys.push_back("put_intra_into_total");
	option_values.push_back(put_intra_into_total_ ? "1" : "0");

	option_keys.push_back("geom_sol_interres_path_distance_cutoff");
	option_values.push_back(geom_sol_interres_path_distance_cutoff_ ? "1" : "0");

	option_keys.push_back("geom_sol_intrares_path_distance_cutoff");
	option_values.push_back(geom_sol_intrares_path_distance_cutoff_ ? "1" : "0");

	option_keys.push_back("eval_intrares_elec_ST_only");
	option_values.push_back(eval_intrares_elec_ST_only_ ? "1" : "0");

	option_keys.push_back("envsmooth_zero_negatives");
	option_values.push_back(envsmooth_zero_negatives_ ? "1" : "0");

	option_keys.push_back("cst_max_seq_sep");
	option_values.push_back(boost::lexical_cast<std::string>(cst_max_seq_sep_));

	option_keys.push_back("cartbonded_len");
	option_values.push_back(boost::lexical_cast<std::string>(cartbonded_len_));

	option_keys.push_back("cartbonded_ang");
	option_values.push_back(boost::lexical_cast<std::string>(cartbonded_ang_));

	option_keys.push_back("cartbonded_tors");
	option_values.push_back(boost::lexical_cast<std::string>(cartbonded_tors_));

	option_keys.push_back("cartbonded_improper");
	option_values.push_back(boost::lexical_cast<std::string>(cartbonded_improper_));

	option_keys.push_back("cartbonded_proton");
	option_values.push_back(boost::lexical_cast<std::string>(cartbonded_proton_));

	option_keys.push_back("cartbonded_linear");
	option_values.push_back(cartbonded_linear_ ? "1" : "0");

	option_keys.push_back("ordered_wat_penalty");
	option_values.push_back(boost::lexical_cast<std::string>(ordered_wat_penalty_));

	option_keys.push_back("ordered_pt_wat_penalty");
	option_values.push_back(boost::lexical_cast<std::string>(ordered_wat_penalty_));

	string statement_string;
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		statement_string = "INSERT OR IGNORE INTO score_function_method_options (batch_id, score_function_name, option_key, option_value) VALUES (?,?,?,?);";
		break;
	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres :
		statement_string = "INSERT IGNORE INTO score_function_method_options (batch_id, score_function_name, option_key, option_value) VALUES (?,?,?,?);";
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}

	cppdb::statement stmt(
		basic::database::safely_prepare_statement(statement_string, db_session));

	for ( Size i=1; i <= option_keys.size(); ++i ) {
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

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::methods::EnergyMethodOptions::save( Archive & arc ) const {
	arc( CEREAL_NVP( aa_composition_setup_files_ ) );
	arc( CEREAL_NVP( netcharge_setup_files_ ) );
	arc( CEREAL_NVP( aspartimide_penalty_value_ ) ); //core::Real
	arc( CEREAL_NVP( atom_vdw_atom_type_set_name_ ) ); // std::string
	arc( CEREAL_NVP( unfolded_energies_type_ ) ); // std::string
	arc( CEREAL_NVP( split_unfolded_label_type_ ) ); // std::string
	arc( CEREAL_NVP( split_unfolded_value_type_ ) ); // std::string
	arc( CEREAL_NVP( method_weights_ ) ); // MethodWeights
	arc( CEREAL_NVP( ss_weights_ ) ); // class core::scoring::SecondaryStructureWeights
	arc( CEREAL_NVP( exclude_protein_protein_fa_elec_ ) ); // _Bool
	arc( CEREAL_NVP( exclude_RNA_RNA_fa_elec_ ) ); // _Bool
	arc( CEREAL_NVP( exclude_monomer_fa_elec_ ) ); // _Bool
	arc( CEREAL_NVP( elec_max_dis_ ) ); // core::Real
	arc( CEREAL_NVP( elec_min_dis_ ) ); // core::Real
	arc( CEREAL_NVP( elec_die_ ) ); // core::Real
	arc( CEREAL_NVP( elec_no_dis_dep_die_ ) ); // _Bool
	arc( CEREAL_NVP( smooth_fa_elec_ ) ); // _Bool
	arc( CEREAL_NVP( elec_sigmoidal_die_ ) ); // _Bool
	arc( CEREAL_NVP( elec_sigmoidal_D_ ) ); // Real
	arc( CEREAL_NVP( elec_sigmoidal_D0_ ) ); // Real
	arc( CEREAL_NVP( elec_sigmoidal_S_ ) ); // Real
	arc( CEREAL_NVP( grpelec_fade_type_ ) ); // std::string
	arc( CEREAL_NVP( grpelec_fade_param1_ ) ); // core::Real
	arc( CEREAL_NVP( grpelec_fade_param2_ ) ); // core::Real
	arc( CEREAL_NVP( grpelec_fade_hbond_ ) ); // _Bool
	arc( CEREAL_NVP( grp_cpfxn_ ) ); // _Bool
	arc( CEREAL_NVP( elec_group_file_ ) ); // std::string
	arc( CEREAL_NVP( grpelec_context_dependent_ ) ); // _Bool
	arc( CEREAL_NVP( use_polarization_ ) ); // _Bool
	arc( CEREAL_NVP( use_gen_kirkwood_ ) ); // _Bool
	arc( CEREAL_NVP( protein_dielectric_ ) );
	arc( CEREAL_NVP( water_dielectric_ ) );
	arc( CEREAL_NVP( exclude_DNA_DNA_ ) ); // _Bool
	arc( CEREAL_NVP( exclude_intra_res_protein_ ) ); // _Bool
	arc( CEREAL_NVP( put_intra_into_total_ ) ); // _Bool
	arc( CEREAL_NVP( geom_sol_interres_path_distance_cutoff_ ) ); // core::Size
	arc( CEREAL_NVP( geom_sol_intrares_path_distance_cutoff_ ) ); // core::Size
	arc( CEREAL_NVP( eval_intrares_elec_ST_only_ ) ); // _Bool
	arc( CEREAL_NVP( envsmooth_zero_negatives_ ) ); // _Bool
	arc( CEREAL_NVP( hbond_options_ ) ); // hbonds::HBondOptionsOP
	arc( CEREAL_NVP( etable_options_ ) ); // core::scoring::etable::EtableOptionsOP
	arc( CEREAL_NVP( rna_options_ ) ); // rna::RNA_EnergyMethodOptionsOP
	arc( CEREAL_NVP( free_dof_options_ ) ); // methods::FreeDOF_OptionsOP
	arc( CEREAL_NVP( cst_max_seq_sep_ ) ); // core::Size
	arc( CEREAL_NVP( cartbonded_len_ ) ); // core::Real
	arc( CEREAL_NVP( cartbonded_ang_ ) ); // core::Real
	arc( CEREAL_NVP( cartbonded_tors_ ) ); // core::Real
	arc( CEREAL_NVP( cartbonded_proton_ ) ); // core::Real
	arc( CEREAL_NVP( cartbonded_improper_ ) ); // core::Real
	arc( CEREAL_NVP( cartbonded_linear_ ) ); // _Bool
	arc( CEREAL_NVP( pb_bound_tag_ ) ); // std::string
	arc( CEREAL_NVP( pb_unbound_tag_ ) ); // std::string
	arc( CEREAL_NVP( fastdens_perres_weights_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( symmetric_gly_tables_ ) ); // _Bool
	arc( CEREAL_NVP( loop_close_use_6D_potential_ ) ); // _Bool
	arc( CEREAL_NVP( fa_stack_base_all_ ) ); // _Bool
	arc( CEREAL_NVP( voids_penalty_energy_containing_cones_cutoff_ ) ); // core::Size
	arc( CEREAL_NVP( voids_penalty_energy_cone_distance_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( voids_penalty_energy_cone_dotproduct_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( voids_penalty_energy_voxel_grid_padding_ ) ); // core::Real
	arc( CEREAL_NVP( voids_penalty_energy_voxel_size_ ) ); // core::Real
	arc( CEREAL_NVP( voids_penalty_energy_disabled_except_during_packing_ ) ); //bool
	arc( CEREAL_NVP( bond_angle_central_atoms_to_score_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( bond_angle_residue_type_param_set_ ) ); // core::scoring::mm::MMBondAngleResidueTypeParamSetOP
	arc( CEREAL_NVP( ordered_pt_wat_penalty_ ) );
	arc( CEREAL_NVP( ordered_wat_penalty_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::methods::EnergyMethodOptions::load( Archive & arc ) {
	arc( aa_composition_setup_files_ );
	arc( netcharge_setup_files_ );
	arc( aspartimide_penalty_value_ ); // core::Real
	arc( atom_vdw_atom_type_set_name_ ); // std::string
	arc( unfolded_energies_type_ ); // std::string
	arc( split_unfolded_label_type_ ); // std::string
	arc( split_unfolded_value_type_ ); // std::string
	arc( method_weights_ ); // MethodWeights
	arc( ss_weights_ ); // class core::scoring::SecondaryStructureWeights
	arc( exclude_protein_protein_fa_elec_ ); // _Bool
	arc( exclude_RNA_RNA_fa_elec_ ); // _Bool
	arc( exclude_monomer_fa_elec_ ); // _Bool
	arc( elec_max_dis_ ); // core::Real
	arc( elec_min_dis_ ); // core::Real
	arc( elec_die_ ); // core::Real
	arc( elec_no_dis_dep_die_ ); // _Bool
	arc( smooth_fa_elec_ ); // _Bool
	arc( elec_sigmoidal_die_ ); // _Bool
	arc( elec_sigmoidal_D_ ); // Real
	arc( elec_sigmoidal_D0_ ); // Real
	arc( elec_sigmoidal_S_ ); // Real
	arc( grpelec_fade_type_ ); // std::string
	arc( grpelec_fade_param1_ ); // core::Real
	arc( grpelec_fade_param2_ ); // core::Real
	arc( grpelec_fade_hbond_ ); // _Bool
	arc( grp_cpfxn_ ); // _Bool
	arc( elec_group_file_ ); // std::string
	arc( grpelec_context_dependent_ ); // _Bool
	arc( use_polarization_ ); // _Bool
	arc( use_gen_kirkwood_ ); // _Bool
	arc( protein_dielectric_ );
	arc( water_dielectric_ );
	arc( exclude_DNA_DNA_ ); // _Bool
	arc( exclude_intra_res_protein_ ); // _Bool
	arc( put_intra_into_total_ ); // _Bool
	arc( geom_sol_interres_path_distance_cutoff_ ); // core::Size
	arc( geom_sol_intrares_path_distance_cutoff_ ); // core::Size
	arc( eval_intrares_elec_ST_only_ ); // _Bool
	arc( envsmooth_zero_negatives_ ); // _Bool
	arc( hbond_options_ ); // hbonds::HBondOptionsOP
	arc( etable_options_ ); // core::scoring::etable::EtableOptionsOP
	arc( rna_options_ ); // rna::RNA_EnergyMethodOptionsOP
	arc( free_dof_options_ ); // methods::FreeDOF_OptionsOP
	arc( cst_max_seq_sep_ ); // core::Size
	arc( cartbonded_len_ ); // core::Real
	arc( cartbonded_ang_ ); // core::Real
	arc( cartbonded_tors_ ); // core::Real
	arc( cartbonded_proton_ ); // core::Real
	arc( cartbonded_improper_ ); // core::Real
	arc( cartbonded_linear_ ); // _Bool
	arc( pb_bound_tag_ ); // std::string
	arc( pb_unbound_tag_ ); // std::string
	arc( fastdens_perres_weights_ ); // utility::vector1<core::Real>
	arc( symmetric_gly_tables_ ); // _Bool
	arc( loop_close_use_6D_potential_ ); // _Bool
	arc( fa_stack_base_all_ ); // _Bool
	arc( voids_penalty_energy_containing_cones_cutoff_ ); // core::Size
	arc( voids_penalty_energy_cone_distance_cutoff_ ); // core::Real
	arc( voids_penalty_energy_cone_dotproduct_cutoff_ ); // core::Real
	arc( voids_penalty_energy_voxel_grid_padding_ ); // core::Real
	arc( voids_penalty_energy_voxel_size_ ); // core::Real
	arc( voids_penalty_energy_disabled_except_during_packing_ ); //bool
	arc( bond_angle_central_atoms_to_score_ ); // utility::vector1<std::string>
	arc( bond_angle_residue_type_param_set_ ); // core::scoring::mm::MMBondAngleResidueTypeParamSetOP
	arc( ordered_pt_wat_penalty_ );
	arc( ordered_wat_penalty_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::methods::EnergyMethodOptions );
CEREAL_REGISTER_TYPE( core::scoring::methods::EnergyMethodOptions )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_methods_EnergyMethodOptions )
#endif // SERIALIZATION
