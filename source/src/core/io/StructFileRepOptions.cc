// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/StructFileRepOptions.cc
/// @brief  Definitions for StructFileRepOptions.
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Unit headers
#include <core/io/StructFileRepOptions.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>

// Tracer instance for this file
static basic::Tracer TR( "core.io.StructFileRepOptions" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace io {

StructFileRepOptions::StructFileRepOptions() {
	init_from_options( basic::options::option );
}

StructFileRepOptions::StructFileRepOptions( utility::options::OptionCollection const & options )
{
	init_from_options( options );
}

StructFileRepOptions::~StructFileRepOptions() = default;

/// @brief Copy this object and return an owning pointer to the copy.
StructFileRepOptionsOP
StructFileRepOptions::clone() const {
	return StructFileRepOptionsOP( new StructFileRepOptions(*this) );
}


void StructFileRepOptions::parse_my_tag( utility::tag::TagCOP tag )
{
	set_check_if_residues_are_Ntermini( tag->getOption< std::string >( "Ntermini", "ALL" ) );
	set_check_if_residues_are_Ctermini( tag->getOption< std::string >( "Ctermini", "ALL" ) );
	set_skip_connect_info( tag->getOption< bool >( "skip_connect_info", false ) );
	set_connect_info_cutoff( tag->getOption< Real >( "connect_info_cutoff", 0.0 ) );
	set_do_not_autoassign_SS( tag->getOption< bool >( "do_not_autoassign_SS", false ) ); //NOTE to later person writing provide_xml_schema: this option does nothing if output_secondary_structure is not set true
	set_exit_if_missing_heavy_atoms( tag->getOption< bool >( "exit_if_missing_heavy_atoms", false ) );
	set_fold_tree_io( tag->getOption< bool >( "fold_tree_io", false ) );
	set_ignore_unrecognized_res( tag->getOption< bool >( "ignore_unrecognized_res", false ) );
	set_ignore_sugars( ! tag->getOption< bool >( "include_sugars", false ) );
	set_ignore_waters( tag->getOption< bool >( "ignore_waters", false ) );
	set_ignore_zero_occupancy( tag->getOption< bool >( "ignore_zero_occupancy", true ) );
	set_guarantee_no_DNA( tag->getOption< bool >( "guarantee_no_DNA", true ) );
	set_keep_input_protonation_state( tag->getOption< bool >( "keep_input_protonation_state", false ) );
	set_preserve_header( tag->getOption< bool >( "preserve_header", false ) );
	set_preserve_crystinfo( tag->getOption< bool >( "preserve_crystinfo", false ) );
	set_missing_dens_as_jump( tag->getOption< bool >( "missing_dens_as_jump", false ) );
	set_no_chainend_ter( tag->getOption< bool >( "no_chainend_ter", false ) );
	set_no_detect_pseudobonds( tag->getOption< bool >("no_detect_pseudobonds", false) );
	set_no_output_cen( tag->getOption< bool >( "no_output_cen", false ) );
	set_normalize_to_thk( tag->getOption< bool >( "normalize_to_thk", false ) );
	set_output_secondary_structure( tag->getOption< bool >( "output_secondary_structure", false ) );
	set_output_torsions( tag->getOption< bool >( "output_torsions", false ) );
	set_output_virtual( tag->getOption< bool >( "output_virtual", false ) );
	set_output_virtual_zero_occ( tag->getOption< bool >( "output_virtual_zero_occ", false ) );
	set_output_only_asymmetric_unit( tag->getOption< bool >( "output_only_asymmetric_unit", false) );
	set_pdb_comments( tag->getOption< bool >( "pdb_comments", false ) );
	set_pdb_parents( tag->getOption< bool >( "pdb_parents", false ) );
	set_per_chain_renumbering( tag->getOption< bool >( "per_chain_renumbering", false ) );
	set_randomize_missing_coords( tag->getOption< bool >( "randomize_missing_coords", false ) );
	set_remember_unrecognized_res( tag->getOption< bool >( "remember_unrecognized_res", false ) );
	set_remember_unrecognized_water( tag->getOption< bool >( "remember_unrecognized_water", false ) );
	set_renumber_pdb( tag->getOption< bool >( "renumber_pdb", false ) );
	set_suppress_zero_occ_pdb_output( tag->getOption< bool >( "suppress_zero_occ_pdb_output", false ) );
	set_auto_detect_glycan_connections( tag->getOption< bool >( "auto_detect_glycan_connections", false) );
	set_write_glycan_pdb_codes( tag->getOption< bool >( "write_glycan_pdb_codes", false) );
	set_output_alternate_atomids( tag->getOption< bool >( "output_alternate_atomids", false) );
	set_output_ligands_as_separate_chains( tag->getOption< bool >( "output_ligands_as_separate_chains", false ) );
	set_maintain_links( tag->getOption< bool >( "maintain_links", false) );
	set_max_bond_length( tag->getOption< core::Real >( "max_bond_length", 1.6) );
	set_min_bond_length( tag->getOption< core::Real >( "min_bond_length", 1.3) );
	set_use_pdb_format_HETNAM_records( tag->getOption< bool >( "use_pdb_format_HETNAM_records", false ) );
	set_write_pdb_link_records( tag->getOption <bool >( "write_pdb_link_records", false ) );
	set_write_pdb_parametric_info( tag->getOption< bool >("write_pdb_parametric_info", true) );
	set_write_all_connect_info( tag->getOption< bool >("write_all_connect_info", false) );
	set_write_seqres_records( tag->getOption< bool >("write_seqres_records", 0) );
	set_chains_whose_residues_are_separate_chemical_entities(
		tag->getOption< std::string >( "treat_residues_in_these_chains_as_separate_chemical_entities", " " ) );

	set_residues_for_atom_name_remapping(
		utility::string_split( tag->getOption< std::string >( "remap_pdb_atom_names_for", "" ), ',') );

	set_show_all_fixes( tag->getOption< bool >( "show_all_fixes", false ) );
	set_constraints_from_link_records( tag->getOption< bool >( "constraints_from_link_records", false ) );
}

std::string StructFileRepOptions::type() const { return "file_data_options"; }


// accessors
//std::string StructFileRepOptions::check_if_residues_are_termini() const { return check_if_residues_are_termini_; }
std::string const & StructFileRepOptions::check_if_residues_are_Ntermini() const { return check_if_residues_are_Ntermini_; }
std::string const & StructFileRepOptions::check_if_residues_are_Ctermini() const { return check_if_residues_are_Ctermini_; }
bool StructFileRepOptions::skip_connect_info() const { return skip_connect_info_; }
core::Real StructFileRepOptions::connect_info_cutoff() const { return connect_info_cutoff_; }
bool StructFileRepOptions::do_not_autoassign_SS() const { return do_not_autoassign_SS_;}
bool StructFileRepOptions::exit_if_missing_heavy_atoms() const { return exit_if_missing_heavy_atoms_; }
bool StructFileRepOptions::fold_tree_io() const { return fold_tree_io_; }
bool StructFileRepOptions::ignore_unrecognized_res() const { return ignore_unrecognized_res_; }
bool StructFileRepOptions::ignore_sugars() const { return ignore_sugars_; }
bool StructFileRepOptions::ignore_waters() const { return ignore_waters_; }
bool StructFileRepOptions::ignore_zero_occupancy() const { return ignore_zero_occupancy_; }
bool StructFileRepOptions::guarantee_no_DNA() const { return guarantee_no_DNA_; }
bool StructFileRepOptions::keep_input_protonation_state() const { return keep_input_protonation_state_; }
bool StructFileRepOptions::preserve_header() const { return preserve_header_; }
bool StructFileRepOptions::preserve_crystinfo() const { return preserve_crystinfo_; }
bool StructFileRepOptions::missing_dens_as_jump() const { return missing_dens_as_jump_; }
bool StructFileRepOptions::no_chainend_ter() const { return no_chainend_ter_; }
bool StructFileRepOptions::no_detect_pseudobonds() const { return no_detect_pseudobonds_; }
bool StructFileRepOptions::no_output_cen() const { return no_output_cen_; }
bool StructFileRepOptions::normalize_to_thk() const { return normalize_to_thk_; }
bool StructFileRepOptions::output_secondary_structure() const { return output_secondary_structure_; }
bool StructFileRepOptions::output_torsions() const { return output_torsions_; }
bool StructFileRepOptions::output_virtual() const { return output_virtual_; }
bool StructFileRepOptions::output_virtual_zero_occ() const { return output_virtual_zero_occ_; }
bool StructFileRepOptions::pdb_parents() const { return pdb_parents_; }
bool StructFileRepOptions::per_chain_renumbering() const { return per_chain_renumbering_; }
bool StructFileRepOptions::randomize_missing_coords() const { return randomize_missing_coords_; }
bool StructFileRepOptions::remember_unrecognized_res() const { return remember_unrecognized_res_; }
bool StructFileRepOptions::remember_unrecognized_water() const { return remember_unrecognized_water_; }
bool StructFileRepOptions::renumber_pdb() const { return renumber_pdb_; }
bool StructFileRepOptions::read_only_ATOM_entries() const { return read_only_ATOM_entries_; }
bool StructFileRepOptions::suppress_zero_occ_pdb_output() const { return suppress_zero_occ_pdb_output_; }
bool StructFileRepOptions::auto_detect_glycan_connections() const { return auto_detect_glycan_connections_; }
bool StructFileRepOptions::write_glycan_pdb_codes() const { return write_glycan_pdb_codes_; }
bool StructFileRepOptions::output_alternate_atomids() const { return output_alternate_atomids_; }
bool StructFileRepOptions::output_ligands_as_separate_chains() const { return output_ligands_as_separate_chains_; }
bool StructFileRepOptions::maintain_links() const { return maintain_links_; }
core::Real StructFileRepOptions::max_bond_length() const { return max_bond_length_; }
core::Real StructFileRepOptions::min_bond_length() const { return min_bond_length_; }
bool StructFileRepOptions::use_pdb_format_HETNAM_records() const { return use_pdb_format_HETNAM_records_; }
bool StructFileRepOptions::write_pdb_link_records() const { return write_pdb_link_records_; }
bool StructFileRepOptions::write_pdb_parametric_info() const { return write_pdb_parametric_info_; }
bool StructFileRepOptions::write_all_connect_info() const { return write_all_connect_info_; }
bool StructFileRepOptions::write_seqres_records() const { return write_seqres_records_; }
std::string const & StructFileRepOptions::chains_whose_residues_are_separate_chemical_entities() const { return chains_whose_residues_are_separate_chemical_entities_; }
utility::vector1<std::string> const & StructFileRepOptions::residues_for_atom_name_remapping() const { return residues_for_atom_name_remapping_; }
bool StructFileRepOptions::pdb_comments() const { return pdb_comments_; }
bool StructFileRepOptions::show_all_fixes() const { return show_all_fixes_; }
bool StructFileRepOptions::constraints_from_link_records() const { return constraints_from_link_records_; }
bool StructFileRepOptions::output_pose_cache() const { return output_pose_cache_data_; }
bool StructFileRepOptions::output_pose_energies_table() const { return output_pose_energies_table_; }
bool StructFileRepOptions::output_only_asymmetric_unit() const { return output_only_asymmetric_unit_; }
// mutators

void StructFileRepOptions::set_check_if_residues_are_Ntermini( std::string const & check_if_residues_are_Ntermini )
{ check_if_residues_are_Ntermini_ = check_if_residues_are_Ntermini; }

void StructFileRepOptions::set_check_if_residues_are_Ctermini( std::string const & check_if_residues_are_Ctermini )
{ check_if_residues_are_Ctermini_ = check_if_residues_are_Ctermini; }

void StructFileRepOptions::set_skip_connect_info( bool const skip_connect_info )
{ skip_connect_info_ = skip_connect_info; }

void StructFileRepOptions::set_connect_info_cutoff( core::Real const & connect_info_cutoff )
{ connect_info_cutoff_ = connect_info_cutoff; }

void StructFileRepOptions::set_do_not_autoassign_SS( bool const do_not_autoassign_SS)
{ do_not_autoassign_SS_ = do_not_autoassign_SS; }

void StructFileRepOptions::set_exit_if_missing_heavy_atoms( bool const exit_if_missing_heavy_atoms )
{ exit_if_missing_heavy_atoms_ = exit_if_missing_heavy_atoms; }

void StructFileRepOptions::set_fold_tree_io( bool const fold_tree_io )
{ fold_tree_io_ = fold_tree_io; }


void StructFileRepOptions::set_ignore_unrecognized_res( bool const ignore_unrecognized_res )
{
	ignore_unrecognized_res_ = ignore_unrecognized_res;
}

void StructFileRepOptions::set_ignore_sugars( bool const setting )
{ ignore_sugars_ = setting; }

void StructFileRepOptions::set_ignore_waters( bool const ignore_waters )
{ ignore_waters_ = ignore_waters; }

void StructFileRepOptions::set_ignore_zero_occupancy( bool const ignore_zero_occupancy )
{ ignore_zero_occupancy_ = ignore_zero_occupancy; }

void StructFileRepOptions::set_guarantee_no_DNA( bool const guarantee_no_DNA )
{ guarantee_no_DNA_ = guarantee_no_DNA; }

void StructFileRepOptions::set_keep_input_protonation_state( bool const keep_input_protonation_state )
{ keep_input_protonation_state_ = keep_input_protonation_state; }

void StructFileRepOptions::set_preserve_header( bool const preserve_header )
{ preserve_header_ = preserve_header; }

void StructFileRepOptions::set_preserve_crystinfo( bool const preserve_crystinfo )
{ preserve_crystinfo_ = preserve_crystinfo; }

void StructFileRepOptions::set_missing_dens_as_jump( bool const missing_dens_as_jump )
{ missing_dens_as_jump_ = missing_dens_as_jump; }

void StructFileRepOptions::set_no_chainend_ter( bool const no_chainend_ter )
{ no_chainend_ter_ = no_chainend_ter; }

void StructFileRepOptions::set_no_detect_pseudobonds( bool const setting )
{ no_detect_pseudobonds_ = setting; }

void StructFileRepOptions::set_no_output_cen( bool const no_output_cen )
{ no_output_cen_ = no_output_cen; }

void StructFileRepOptions::set_normalize_to_thk( bool const normalize_to_thk )
{ normalize_to_thk_ = normalize_to_thk; }

void StructFileRepOptions::set_output_secondary_structure( bool const output_secondary_structure )
{ output_secondary_structure_ = output_secondary_structure; }

void StructFileRepOptions::set_output_torsions( bool const output_torsions )
{ output_torsions_ = output_torsions; }

void StructFileRepOptions::set_output_virtual( bool const output_virtual )
{ output_virtual_ = output_virtual; }

void StructFileRepOptions::set_output_virtual_zero_occ( bool const output_virtual_zero_occ )
{ output_virtual_zero_occ_ = output_virtual_zero_occ; }

void StructFileRepOptions::set_output_only_asymmetric_unit( bool const output_only_asymmetric_unit )
{ output_only_asymmetric_unit_ = output_only_asymmetric_unit; }

void StructFileRepOptions::set_pdb_comments( bool const pdb_comments )
{ pdb_comments_ = pdb_comments; }

void StructFileRepOptions::set_pdb_parents( bool const pdb_parents )
{ pdb_parents_ = pdb_parents; }

void StructFileRepOptions::set_per_chain_renumbering( bool const per_chain_renumbering )
{ per_chain_renumbering_ = per_chain_renumbering; }

void StructFileRepOptions::set_randomize_missing_coords( bool const randomize_missing_coords )
{ randomize_missing_coords_ = randomize_missing_coords; }

void StructFileRepOptions::set_remember_unrecognized_res( bool const remember_unrecognized_res )
{ remember_unrecognized_res_ = remember_unrecognized_res; }

void StructFileRepOptions::set_remember_unrecognized_water( bool const remember_unrecognized_water )
{ remember_unrecognized_water_ = remember_unrecognized_water; }

void StructFileRepOptions::set_renumber_pdb( bool const setting )
{ renumber_pdb_ = setting; }

void StructFileRepOptions::set_read_only_ATOM_entries( bool const setting )
{ read_only_ATOM_entries_ = setting; }

void StructFileRepOptions::set_suppress_zero_occ_pdb_output( bool const setting )
{ suppress_zero_occ_pdb_output_ = setting; }

void StructFileRepOptions::set_auto_detect_glycan_connections( bool const auto_detect_glycan_connections )
{ auto_detect_glycan_connections_ = auto_detect_glycan_connections; }

void StructFileRepOptions::set_write_glycan_pdb_codes( bool const write_glycan_pdb_codes )
{ write_glycan_pdb_codes_ = write_glycan_pdb_codes; }

void StructFileRepOptions::set_output_alternate_atomids( bool const output_alternate_atomids )
{ output_alternate_atomids_ = output_alternate_atomids; }

void StructFileRepOptions::set_output_ligands_as_separate_chains( bool const output_ligands_as_separate_chains )
{ output_ligands_as_separate_chains_ = output_ligands_as_separate_chains; }

void StructFileRepOptions::set_maintain_links( bool const maintain_links )
{ maintain_links_ = maintain_links; }

void StructFileRepOptions::set_max_bond_length( core::Real const max_bond_length )
{ max_bond_length_ = max_bond_length; }

void StructFileRepOptions::set_min_bond_length( core::Real const min_bond_length )
{ min_bond_length_ = min_bond_length; }

void StructFileRepOptions::set_use_pdb_format_HETNAM_records( bool const setting )
{ use_pdb_format_HETNAM_records_ = setting; }

void StructFileRepOptions::set_write_pdb_link_records( bool const setting )
{ write_pdb_link_records_ = setting; }

void StructFileRepOptions::set_write_pdb_parametric_info( bool const setting )
{ write_pdb_parametric_info_ = setting; }

void StructFileRepOptions::set_write_all_connect_info( bool const setting )
{ write_all_connect_info_ = setting; }

void StructFileRepOptions::set_write_seqres_records(bool const setting )
{ write_seqres_records_ = setting; }

void StructFileRepOptions::set_chains_whose_residues_are_separate_chemical_entities( std::string const & chains_whose_residues_are_separate_chemical_entities )
{ chains_whose_residues_are_separate_chemical_entities_ = chains_whose_residues_are_separate_chemical_entities; }


void StructFileRepOptions::set_residues_for_atom_name_remapping(utility::vector1<std::string> const & setting) {
	residues_for_atom_name_remapping_ = setting;
}

void StructFileRepOptions::set_show_all_fixes( bool const setting ) { show_all_fixes_ = setting; }
void StructFileRepOptions::set_constraints_from_link_records( bool const setting ) { constraints_from_link_records_ = setting; }

void StructFileRepOptions::set_output_pose_energies_table( bool const setting ) { output_pose_energies_table_ = setting; }
void StructFileRepOptions::set_output_pose_cache_data( bool const setting ) { output_pose_cache_data_ = setting; }

/// @details List all of the options (by option key) that are read in the init_from_options function.
void
StructFileRepOptions::list_options_read( utility::options::OptionKeyList & read_options )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	read_options
		+ in::Ntermini
		+ in::Ctermini
		+ inout::skip_connect_info
		+ inout::connect_info_cutoff
		+ out::file::do_not_autoassign_SS
		+ run::exit_if_missing_heavy_atoms
		+ inout::fold_tree_io
		+ in::include_sugars
		+ in::ignore_unrecognized_res
		+ in::ignore_waters
		+ run::ignore_zero_occupancy
		+ in::guarantee_no_DNA
		+ pH::keep_input_protonation_state
		+ run::preserve_header
		+ in::preserve_crystinfo
		+ in::missing_density_to_jump
		+ out::file::no_chainend_ter
		+ in::file::no_detect_pseudobonds
		+ out::file::no_output_cen
		+ mp::output::normalize_to_thk
		+ out::file::output_secondary_structure
		+ out::file::output_torsions
		+ out::file::output_virtual
		+ out::file::output_virtual_zero_occ
		+ out::file::output_only_asymmetric_unit
		+ out::file::pdb_comments
		+ out::file::pdb_parents
		+ out::file::per_chain_renumbering
		+ run::randomize_missing_coords
		+ in::remember_unrecognized_res
		+ in::remember_unrecognized_water
		+ out::file::renumber_pdb
		+ out::file::suppress_zero_occ_pdb_output
		+ out::file::use_pdb_format_HETNAM_records
		+ out::file::write_pdb_link_records
		+ in::file::treat_residues_in_these_chains_as_separate_chemical_entities
		+ in::file::remap_pdb_atom_names_for
		+ out::file::write_pdb_parametric_info
		+ inout::write_all_connect_info
		+ inout::output_alternate_atomids
		+ out::file::output_ligands_as_separate_chains
		+ in::show_all_fixes
		+ in::constraints_from_link_records
		+ out::file::output_pose_energies_table
		+ out::file::output_pose_cache_data
		+ out::file::output_pose_fold_tree
		+ out::file::write_glycan_pdb_codes
		+ inout::output_alternate_atomids
		+ in::maintain_links
		+ in::auto_detect_glycan_connections
		+ in::max_bond_length
		+ in::min_bond_length
		+ in::read_only_ATOM_entries
		+ out::file::write_seqres_records;

}

/// @details The StructFileRepOptions object can be used to initialize a Pose within
/// the ResourceManager and elsewhere. The attributes (aka options) listed here are
/// those read in the parse_my_tag method.
void
StructFileRepOptions::append_schema_attributes( utility::tag::AttributeList & attributes )
{
	using namespace utility::tag;
	typedef XMLSchemaAttribute Attr;

	attributes
		+ Attr::attribute_w_default( "Ntermini", xs_string, "TO DO", "ALL" )
		+ Attr::attribute_w_default( "Ctermini", xs_string, "TO DO", "ALL" )
		+ Attr::attribute_w_default( "skip_connect_info", xsct_rosetta_bool, "TO DO", "0" )
		+ Attr::attribute_w_default( "connect_info_cutoff", xsct_real, "TO DO", "0.0" )
		+ Attr::attribute_w_default( "do_not_autoassign_SS", xsct_rosetta_bool, "Requires output_secondary_structure to be set to true",  "false" )
		+ Attr::attribute_w_default( "exit_if_missing_heavy_atoms", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "fold_tree_io", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "ignore_unrecognized_res", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "include_sugars", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "ignore_waters", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "ignore_zero_occupancy", xsct_rosetta_bool, "TO DO",  "1" )
		+ Attr::attribute_w_default( "guarantee_no_DNA", xsct_rosetta_bool, "TO DO",  "1" )
		+ Attr::attribute_w_default( "keep_input_protonation_state", xsct_rosetta_bool, "TO DO",   "0" )
		+ Attr::attribute_w_default( "preserve_header", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "preserve_crystinfo", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "missing_dens_as_jump", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "no_chainend_ter", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "no_output_cen", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "normalize_to_thk", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "output_ligands_as_separate_chains", xsct_rosetta_bool, "Forces Rosetta to reassign .pdb file chainIDs to ligands sharing chain IDs with peptide chains.", "false" )
		+ Attr::attribute_w_default( "output_secondary_structure", xsct_rosetta_bool, "TO DO",  "false" )
		+ Attr::attribute_w_default( "output_torsions", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "output_virtual", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "output_virtual_zero_occ", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "output_only_asymmetric_unit", xsct_rosetta_bool, "If symmetrical pose, only output PDB coord of master subunit.", "false")
		+ Attr::attribute_w_default( "pdb_comments", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "pdb_parents", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "per_chain_renumbering", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "randomize_missing_coords", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "remember_unrecognized_res", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "remember_unrecognized_water", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "renumber_pdb", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "suppress_zero_occ_pdb_output", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "auto_detect_glycan_connections", xsct_rosetta_bool, "TO DO", "0" )
		+ Attr::attribute_w_default( "write_glycan_pdb_codes", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "maintain_links", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "max_bond_length", xsct_real, "TO DO",  "1.6" )
		+ Attr::attribute_w_default( "min_bond_length", xsct_real, "TO DO", "1.3" )
		+ Attr::attribute_w_default( "use_pdb_format_HETNAM_records", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "write_pdb_link_records", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "write_pdb_parametric_info", xsct_rosetta_bool, "TO DO", "1")
		+ Attr::attribute_w_default( "write_all_connect_info", xsct_rosetta_bool, "TO DO", "0" )
		+ Attr::attribute_w_default( "treat_residues_in_these_chains_as_separate_chemical_entities", xs_string, "TO DO", " " )
		+ Attr::attribute_w_default( "remap_pdb_atom_names_for", xs_string, "Comma separated list of atom names to remap to other names", "" )
		+ Attr::attribute_w_default( "show_all_fixes", xsct_rosetta_bool, "TO DO",  "0" )
		+ Attr::attribute_w_default( "constraints_from_link_records", xsct_rosetta_bool, "TO DO", "0" );

}


/// @details Note that if you should add a new option to this function, that you should
/// 1. add a data member to the class (of course!)
/// 2. add an accessor for that data member (of course!)
/// 3. add a mutator for that data member such that it can be set programmatically
/// 4. add the option key to the OptionKeyList in list_options_read
/// 5. add a read to the (XML) Tag object in the parse_my_tag function
/// 6. list the XML attribute (aka "Tag option") that you read from in parse_my_tag to the provide_xml_schema function.
/// 7. add to operator < and operator ==
/// 8. add to init_from_options (here)
///
/// The usability of Rosetta is improved by continued vigilance of its developers to
/// provide a well documented and accessible code, and the process of reading in structures
/// is vital to Rosetta. Do your part and follow in the careful efforts of those who
/// came before you.
void StructFileRepOptions::init_from_options( utility::options::OptionCollection const & options )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_check_if_residues_are_Ntermini( options[ in::Ntermini ].value());
	set_check_if_residues_are_Ctermini( options[ in::Ctermini ].value());
	set_skip_connect_info( options[ inout::skip_connect_info ].value());
	set_connect_info_cutoff( options[ inout::connect_info_cutoff ].value());
	set_do_not_autoassign_SS( options [ OptionKeys::out::file::do_not_autoassign_SS ].value());
	set_exit_if_missing_heavy_atoms( options[ run::exit_if_missing_heavy_atoms ].value());
	set_fold_tree_io( options[ inout::fold_tree_io ]() );

	set_ignore_sugars( ! options[ in::include_sugars ]() );
	set_ignore_unrecognized_res( options[ in::ignore_unrecognized_res ]());
	set_ignore_waters( options[ in::ignore_waters ]() );

	set_ignore_zero_occupancy( options[ run::ignore_zero_occupancy ]());
	set_guarantee_no_DNA( options[ in::guarantee_no_DNA].value());
	set_keep_input_protonation_state( options[ pH::keep_input_protonation_state ]());
	set_preserve_header( options[ run::preserve_header ].value());
	set_preserve_crystinfo( options[ in::preserve_crystinfo ]() );
	set_missing_dens_as_jump( options[ in::missing_density_to_jump ]() );
	set_no_chainend_ter( options[ OptionKeys::out::file::no_chainend_ter ]() );
	set_no_detect_pseudobonds( options[OptionKeys::in::file::no_detect_pseudobonds]() );
	set_no_output_cen( options[ OptionKeys::out::file::no_output_cen ]() );
	set_normalize_to_thk( options[ OptionKeys::mp::output::normalize_to_thk ]() );
	set_output_secondary_structure( options[ OptionKeys::out::file::output_secondary_structure ]() );
	set_output_torsions( options[ OptionKeys::out::file::output_torsions ]() );
	set_output_virtual( options[ OptionKeys::out::file::output_virtual ]() );
	set_output_virtual_zero_occ( options[ OptionKeys::out::file::output_virtual_zero_occ ]() );
	set_output_only_asymmetric_unit( options[ OptionKeys::out::file::output_only_asymmetric_unit]() );
	set_pdb_comments( options[ OptionKeys::out::file::pdb_comments ].value() );
	set_pdb_parents( options[ OptionKeys::out::file::pdb_parents ].value() );
	set_per_chain_renumbering( options[ OptionKeys::out::file::per_chain_renumbering ].value() );
	set_randomize_missing_coords( options[ run::randomize_missing_coords ]());
	set_remember_unrecognized_res( options[ in::remember_unrecognized_res ]());
	set_remember_unrecognized_water( options[ in::remember_unrecognized_water ]() );
	set_renumber_pdb( options[ OptionKeys::out::file::renumber_pdb ].value() );
	set_suppress_zero_occ_pdb_output( options[ OptionKeys::out::file::suppress_zero_occ_pdb_output ]() );
	set_auto_detect_glycan_connections( options[ in::auto_detect_glycan_connections ]() );
	set_write_glycan_pdb_codes( options[ out::file::write_glycan_pdb_codes ]() );
	set_output_alternate_atomids( options[ inout::output_alternate_atomids ]() );
	set_output_ligands_as_separate_chains( options[ out::file::output_ligands_as_separate_chains ]() );
	set_maintain_links( options[ in::maintain_links ]() );
	set_max_bond_length( options[ in::max_bond_length ]() );
	set_min_bond_length( options[ in::min_bond_length ]() );
	set_read_only_ATOM_entries( options[ in::read_only_ATOM_entries ]() );
	set_use_pdb_format_HETNAM_records( options[ OptionKeys::out::file::use_pdb_format_HETNAM_records ]() );
	set_write_pdb_link_records( options[ out::file::write_pdb_link_records ]() );
	set_write_seqres_records( options[out::file::write_seqres_records ]() );
	set_chains_whose_residues_are_separate_chemical_entities( options[ in::file::treat_residues_in_these_chains_as_separate_chemical_entities].user_or(""));
	if ( options[ in::file::remap_pdb_atom_names_for ].active() ) {
		set_residues_for_atom_name_remapping( options[ in::file::remap_pdb_atom_names_for ] );
	}
	set_write_pdb_parametric_info( options[out::file::write_pdb_parametric_info]() );
	set_write_all_connect_info( options[inout::write_all_connect_info]() );
	set_show_all_fixes( options[ in::show_all_fixes ]() );
	set_constraints_from_link_records( options[ in::constraints_from_link_records ]() );

	set_output_pose_energies_table(options[ OptionKeys::out::file::output_pose_energies_table ]()); //JD2 only?
	set_output_pose_cache_data(options[ OptionKeys::out::file::output_pose_cache_data ]()); //JD2 only?
	set_fold_tree_io( options[ OptionKeys::out::file::output_pose_fold_tree]());
}

bool
StructFileRepOptions::operator == ( StructFileRepOptions const & other ) const
{
	if ( check_if_residues_are_Ntermini_                        != other.check_if_residues_are_Ntermini_                       ) return false;
	if ( check_if_residues_are_Ctermini_                        != other.check_if_residues_are_Ctermini_                      ) return false;
	if ( skip_connect_info_                                     != other.skip_connect_info_                                    ) return false;
	if ( connect_info_cutoff_                                   != other.connect_info_cutoff_                                  ) return false;
	if ( do_not_autoassign_SS_                                  != other.do_not_autoassign_SS_                                 ) return false;
	if ( exit_if_missing_heavy_atoms_                           != other.exit_if_missing_heavy_atoms_                          ) return false;
	if ( fold_tree_io_                                          != other.fold_tree_io_                                         ) return false;
	if ( ignore_unrecognized_res_                               != other.ignore_unrecognized_res_                              ) return false;
	if ( ignore_sugars_                                         != other.ignore_sugars_                                        ) return false;
	if ( ignore_waters_                                         != other.ignore_waters_                                        ) return false;
	if ( ignore_zero_occupancy_                                 != other.ignore_zero_occupancy_                                ) return false;
	if ( keep_input_protonation_state_                          != other.keep_input_protonation_state_                         ) return false;
	if ( preserve_header_                                       != other.preserve_header_                                      ) return false;
	if ( preserve_crystinfo_                                    != other.preserve_crystinfo_                                   ) return false;
	if ( missing_dens_as_jump_                                  != other.missing_dens_as_jump_                                 ) return false;
	if ( no_chainend_ter_                                       != other.no_chainend_ter_                                      ) return false;
	if ( no_detect_pseudobonds_        != other.no_detect_pseudobonds_       ) return false;
	if ( no_output_cen_                                         != other.no_output_cen_                                        ) return false;
	if ( normalize_to_thk_                                      != other.normalize_to_thk_                                     ) return false;
	if ( output_secondary_structure_                            != other.output_secondary_structure_                           ) return false;
	if ( output_torsions_                                       != other.output_torsions_                                      ) return false;
	if ( output_virtual_                                        != other.output_virtual_                                       ) return false;
	if ( output_virtual_zero_occ_                               != other.output_virtual_zero_occ_                              ) return false;
	if ( output_only_asymmetric_unit_                            != other.output_only_asymmetric_unit_                              ) return false;
	if ( pdb_parents_                                           != other.pdb_parents_                                          ) return false;
	if ( per_chain_renumbering_                                 != other.per_chain_renumbering_                                ) return false;
	if ( randomize_missing_coords_                              != other.randomize_missing_coords_                             ) return false;
	if ( remember_unrecognized_res_                             != other.remember_unrecognized_res_                            ) return false;
	if ( remember_unrecognized_water_                           != other.remember_unrecognized_water_                          ) return false;
	if ( renumber_pdb_                                          != other.renumber_pdb_                                         ) return false;
	if ( suppress_zero_occ_pdb_output_                          != other.suppress_zero_occ_pdb_output_                         ) return false;
	if ( use_pdb_format_HETNAM_records_                         != other.use_pdb_format_HETNAM_records_                        ) return false;
	if ( write_pdb_link_records_                                != other.write_pdb_link_records_                               ) return false;
	if ( write_seqres_records_                                  != other.write_seqres_records_                                 ) return false;
	if ( write_pdb_parametric_info_                             != other.write_pdb_parametric_info_                            ) return false;
	if ( write_all_connect_info_                                != other.write_all_connect_info_                               ) return false;
	if ( chains_whose_residues_are_separate_chemical_entities_  != other.chains_whose_residues_are_separate_chemical_entities_ ) return false;
	if ( residues_for_atom_name_remapping_                      != other.residues_for_atom_name_remapping_                     ) return false;
	if ( pdb_comments_                                          != other.pdb_comments_                                         ) return false;
	if ( show_all_fixes_                                        != other.show_all_fixes_                                       ) return false;
	if ( constraints_from_link_records_                         != other.constraints_from_link_records_                        ) return false;
	if ( auto_detect_glycan_connections_                        != other.auto_detect_glycan_connections_                       ) return false;
	if ( write_glycan_pdb_codes_                                != other.write_glycan_pdb_codes_                               ) return false;
	if ( output_alternate_atomids_                              != other.output_alternate_atomids_                             ) return false;
	if ( output_ligands_as_separate_chains_                     != other.output_ligands_as_separate_chains_                    ) return false;
	if ( maintain_links_                                        != other.maintain_links_                                       ) return false;
	if ( max_bond_length_                                       != other.max_bond_length_                                      ) return false;
	if ( min_bond_length_                                       != other.min_bond_length_                                      ) return false;

	return true;
}

/// @note This function makes no sense.  If ANY option is less than that in other, this returns true. --VKM 24 Jan 2018.
bool
StructFileRepOptions::operator < ( StructFileRepOptions const & other ) const
{
	if ( check_if_residues_are_Ntermini_                        <  other.check_if_residues_are_Ntermini_                       ) return true;
	if ( check_if_residues_are_Ntermini_                        != other.check_if_residues_are_Ntermini_                       ) return false;
	if ( check_if_residues_are_Ctermini_                        <  other.check_if_residues_are_Ctermini_                       ) return true;
	if ( check_if_residues_are_Ctermini_                        != other.check_if_residues_are_Ctermini_                       ) return false;
	if ( skip_connect_info_                                     <  other.skip_connect_info_                                    ) return true;
	if ( skip_connect_info_                                     != other.skip_connect_info_                                    ) return false;
	if ( connect_info_cutoff_                                   <  other.connect_info_cutoff_                                  ) return true;
	if ( connect_info_cutoff_                                   != other.connect_info_cutoff_                                  ) return false;
	if ( do_not_autoassign_SS_                                  < other.do_not_autoassign_SS_                                 ) return true;
	if ( do_not_autoassign_SS_                                  != other.do_not_autoassign_SS_                                 ) return false;
	if ( exit_if_missing_heavy_atoms_                           <  other.exit_if_missing_heavy_atoms_                          ) return true;
	if ( exit_if_missing_heavy_atoms_                           != other.exit_if_missing_heavy_atoms_                          ) return false;
	if ( fold_tree_io_                                          <  other.fold_tree_io_                                         ) return true;
	if ( fold_tree_io_                                          != other.fold_tree_io_                                         ) return false;
	if ( ignore_unrecognized_res_                               <  other.ignore_unrecognized_res_                              ) return true;
	if ( ignore_unrecognized_res_                               != other.ignore_unrecognized_res_                              ) return false;
	if ( ignore_sugars_                                         <  other.ignore_sugars_                                        ) return true;
	if ( ignore_sugars_                                         != other.ignore_sugars_                                        ) return false;
	if ( ignore_waters_                                         <  other.ignore_waters_                                        ) return true;
	if ( ignore_waters_                                         != other.ignore_waters_                                        ) return false;
	if ( ignore_zero_occupancy_                                 <  other.ignore_zero_occupancy_                                ) return true;
	if ( ignore_zero_occupancy_                                 != other.ignore_zero_occupancy_                                ) return false;
	if ( keep_input_protonation_state_                          <  other.keep_input_protonation_state_                         ) return true;
	if ( keep_input_protonation_state_                          != other.keep_input_protonation_state_                         ) return false;
	if ( preserve_header_                                       <  other.preserve_header_                                      ) return true;
	if ( preserve_header_                                       != other.preserve_header_                                      ) return false;
	if ( preserve_crystinfo_                                    <  other.preserve_crystinfo_                                   ) return true;
	if ( preserve_crystinfo_                                    != other.preserve_crystinfo_                                   ) return false;
	if ( missing_dens_as_jump_                                  <  other.missing_dens_as_jump_                                 ) return true;
	if ( missing_dens_as_jump_                                  != other.missing_dens_as_jump_                                 ) return false;
	if ( no_chainend_ter_                                       <  other.no_chainend_ter_                                      ) return true;
	if ( no_chainend_ter_                                       != other.no_chainend_ter_                                      ) return false;
	if ( no_detect_pseudobonds_        <  other.no_detect_pseudobonds_       ) return true;
	if ( no_detect_pseudobonds_        !=  other.no_detect_pseudobonds_       ) return false;
	if ( no_output_cen_                                         <  other.no_output_cen_                                        ) return true;
	if ( no_output_cen_                                         != other.no_output_cen_                                        ) return false;
	if ( normalize_to_thk_                                      <  other.normalize_to_thk_                                     ) return true;
	if ( normalize_to_thk_                                      != other.normalize_to_thk_                                     ) return false;
	if ( output_secondary_structure_                            <  other.output_secondary_structure_                           ) return true;
	if ( output_secondary_structure_                            != other.output_secondary_structure_                           ) return false;
	if ( output_torsions_                                       <  other.output_torsions_                                      ) return true;
	if ( output_torsions_                                       != other.output_torsions_                                      ) return false;
	if ( output_virtual_                                        <  other.output_virtual_                                       ) return true;
	if ( output_virtual_                                        != other.output_virtual_                                       ) return false;
	if ( output_virtual_zero_occ_                               <  other.output_virtual_zero_occ_                              ) return true;
	if ( output_virtual_zero_occ_                               != other.output_virtual_zero_occ_                              ) return false;
	if ( output_only_asymmetric_unit_                               <  other.output_only_asymmetric_unit_                              ) return true;
	if ( output_only_asymmetric_unit_                               != other.output_only_asymmetric_unit_                              ) return false;
	if ( pdb_parents_                                           <  other.pdb_parents_                                          ) return true;
	if ( pdb_parents_                                           != other.pdb_parents_                                          ) return false;
	if ( per_chain_renumbering_                                 <  other.per_chain_renumbering_                                ) return true;
	if ( per_chain_renumbering_                                 != other.per_chain_renumbering_                                ) return false;
	if ( randomize_missing_coords_                              <  other.randomize_missing_coords_                             ) return true;
	if ( randomize_missing_coords_                              != other.randomize_missing_coords_                             ) return false;
	if ( remember_unrecognized_res_                             <  other.remember_unrecognized_res_                            ) return true;
	if ( remember_unrecognized_res_                             != other.remember_unrecognized_res_                            ) return false;
	if ( remember_unrecognized_water_                           <  other.remember_unrecognized_water_                          ) return true;
	if ( remember_unrecognized_water_                           != other.remember_unrecognized_water_                          ) return false;
	if ( renumber_pdb_                                          <  other.renumber_pdb_                                         ) return true;
	if ( renumber_pdb_                                          != other.renumber_pdb_                                         ) return false;
	if ( suppress_zero_occ_pdb_output_                          <  other.suppress_zero_occ_pdb_output_                         ) return true;
	if ( suppress_zero_occ_pdb_output_                          != other.suppress_zero_occ_pdb_output_                         ) return false;
	if ( use_pdb_format_HETNAM_records_                         <  other.use_pdb_format_HETNAM_records_                        ) return true;
	if ( use_pdb_format_HETNAM_records_                         != other.use_pdb_format_HETNAM_records_                        ) return false;
	if ( write_pdb_link_records_                                <  other.write_pdb_link_records_                               ) return true;
	if ( write_pdb_link_records_                                != other.write_pdb_link_records_                               ) return false;
	if ( write_seqres_records_                                  <  other.write_seqres_records_                                 ) return true;
	if ( write_seqres_records_                                  != other.write_seqres_records_                                 ) return false;
	if ( write_pdb_parametric_info_                             <  other.write_pdb_parametric_info_                            ) return true;
	if ( write_pdb_parametric_info_                             != other.write_pdb_parametric_info_                            ) return false;
	if ( write_all_connect_info_                                <  other.write_all_connect_info_                               ) return true;
	if ( write_all_connect_info_                                != other.write_all_connect_info_                               ) return false;
	if ( chains_whose_residues_are_separate_chemical_entities_  <  other.chains_whose_residues_are_separate_chemical_entities_ ) return true;
	if ( chains_whose_residues_are_separate_chemical_entities_  != other.chains_whose_residues_are_separate_chemical_entities_ ) return false;
	if ( residues_for_atom_name_remapping_                      <  other.residues_for_atom_name_remapping_                     ) return true;
	if ( residues_for_atom_name_remapping_                      != other.residues_for_atom_name_remapping_                     ) return false;
	if ( pdb_comments_                                          <  other.pdb_comments_                                         ) return true;
	if ( pdb_comments_                                          != other.pdb_comments_                                         ) return false;
	if ( show_all_fixes_                                        <  other.show_all_fixes_                                       ) return true;
	if ( show_all_fixes_                                        != other.show_all_fixes_                                       ) return false;
	if ( constraints_from_link_records_                         <  other.constraints_from_link_records_                        ) return true;
	//if ( constraints_from_link_records_                         != other.constraints_from_link_records_                        ) return false;
	if ( auto_detect_glycan_connections_                        <  other.auto_detect_glycan_connections_                       ) return true;
	if ( auto_detect_glycan_connections_                        != other.auto_detect_glycan_connections_                       ) return false;
	if ( write_glycan_pdb_codes_                                <  other.write_glycan_pdb_codes_                               ) return true;
	if ( write_glycan_pdb_codes_                                != other.write_glycan_pdb_codes_                               ) return false;
	if ( output_alternate_atomids_                              <  other.output_alternate_atomids_                             ) return true;
	if ( output_alternate_atomids_                              != other.output_alternate_atomids_                             ) return false;
	if ( output_ligands_as_separate_chains_                     <  other.output_ligands_as_separate_chains_                    ) return true;
	if ( output_ligands_as_separate_chains_                     != other.output_ligands_as_separate_chains_                    ) return false;
	if ( maintain_links_                                        <  other.maintain_links_                                       ) return true;
	if ( maintain_links_                                        != other.maintain_links_                                       ) return false;
	if ( max_bond_length_                                       <  other.max_bond_length_                                      ) return true;
	if ( max_bond_length_                                       != other.max_bond_length_                                      ) return false;
	if ( min_bond_length_                                       <  other.min_bond_length_                                      ) return true;
	if ( min_bond_length_                                       != other.min_bond_length_                                      ) return false;
	return false;
}


} // namespace io
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::io::StructFileRepOptions::save( Archive & arc ) const {
	arc( CEREAL_NVP( check_if_residues_are_Ntermini_ ) ); // std::string
	arc( CEREAL_NVP( check_if_residues_are_Ctermini_ ) ); // std::string
	arc( CEREAL_NVP( skip_connect_info_ ) ); // _Bool
	arc( CEREAL_NVP( connect_info_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( do_not_autoassign_SS_ ) ); // _Bool
	arc( CEREAL_NVP( exit_if_missing_heavy_atoms_ ) ); // _Bool
	arc( CEREAL_NVP( fold_tree_io_ ) ); // _Bool
	arc( CEREAL_NVP( ignore_unrecognized_res_ ) ); // _Bool
	arc( CEREAL_NVP( ignore_sugars_ ) ); // _Bool
	arc( CEREAL_NVP( ignore_waters_ ) ); // _Bool
	arc( CEREAL_NVP( ignore_zero_occupancy_ ) ); // _Bool
	arc( CEREAL_NVP( guarantee_no_DNA_ ) ); // _Bool
	arc( CEREAL_NVP( keep_input_protonation_state_ ) ); // _Bool
	arc( CEREAL_NVP( preserve_header_ ) ); // _Bool
	arc( CEREAL_NVP( preserve_crystinfo_ ) ); // _Bool
	arc( CEREAL_NVP( missing_dens_as_jump_ ) ); // _Bool
	arc( CEREAL_NVP( no_chainend_ter_ ) ); // _Bool
	arc( CEREAL_NVP( no_detect_pseudobonds_ ) ); // _Bool
	arc( CEREAL_NVP( no_output_cen_ ) ); // _Bool
	arc( CEREAL_NVP( normalize_to_thk_ ) ); // _Bool
	arc( CEREAL_NVP( output_secondary_structure_ ) ); // _Bool
	arc( CEREAL_NVP( output_torsions_ ) ); // _Bool
	arc( CEREAL_NVP( output_virtual_ ) ); // _Bool
	arc( CEREAL_NVP( output_virtual_zero_occ_ ) ); // _Bool
	arc( CEREAL_NVP( output_only_asymmetric_unit_ ) ); // _Bool
	arc( CEREAL_NVP( pdb_parents_ ) ); // _Bool
	arc( CEREAL_NVP( per_chain_renumbering_ ) ); // _Bool
	arc( CEREAL_NVP( randomize_missing_coords_ ) ); // _Bool
	arc( CEREAL_NVP( remember_unrecognized_res_ ) ); // _Bool
	arc( CEREAL_NVP( remember_unrecognized_water_ ) ); // _Bool
	arc( CEREAL_NVP( renumber_pdb_ ) ); // _Bool
	arc( CEREAL_NVP( read_only_ATOM_entries_ ) ); // _Bool
	arc( CEREAL_NVP( suppress_zero_occ_pdb_output_ ) ); // _Bool
	arc( CEREAL_NVP( auto_detect_glycan_connections_ ) ); // _Bool
	arc( CEREAL_NVP( write_glycan_pdb_codes_ ) ); // _Bool
	arc( CEREAL_NVP( output_alternate_atomids_ ) ); // _Bool
	arc( CEREAL_NVP( output_ligands_as_separate_chains_ ) );  // _Bool
	arc( CEREAL_NVP( maintain_links_ ) ); // _Bool
	arc( CEREAL_NVP( max_bond_length_ ) ); // core::Real
	arc( CEREAL_NVP( min_bond_length_ ) ); // core::Real
	arc( CEREAL_NVP( use_pdb_format_HETNAM_records_ ) ); // _Bool
	arc( CEREAL_NVP( write_pdb_link_records_ ) ); // _Bool
	arc( CEREAL_NVP( write_pdb_parametric_info_ ) ); // _Bool
	arc( CEREAL_NVP( write_all_connect_info_ ) ); // _Bool
	arc( CEREAL_NVP( write_seqres_records_ ) ); // _Bool
	arc( CEREAL_NVP( chains_whose_residues_are_separate_chemical_entities_ ) ); // std::string
	arc( CEREAL_NVP( residues_for_atom_name_remapping_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( pdb_comments_ ) ); // _Bool
	arc( CEREAL_NVP( show_all_fixes_ ) ); // _Bool
	arc( CEREAL_NVP( constraints_from_link_records_ ) ); // _Bool
	arc( CEREAL_NVP( output_pose_energies_table_ ) ); // _Bool
	arc( CEREAL_NVP( output_pose_cache_data_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::io::StructFileRepOptions::load( Archive & arc ) {
	arc( check_if_residues_are_Ntermini_ ); // std::string
	arc( check_if_residues_are_Ctermini_ ); // std::string
	arc( skip_connect_info_ ); // _Bool
	arc( connect_info_cutoff_ ); // core::Real
	arc( do_not_autoassign_SS_ ); // _Bool
	arc( exit_if_missing_heavy_atoms_ ); // _Bool
	arc( fold_tree_io_ ); // _Bool
	arc( ignore_unrecognized_res_ ); // _Bool
	arc( ignore_sugars_ ); // _Bool
	arc( ignore_waters_ ); // _Bool
	arc( ignore_zero_occupancy_ ); // _Bool
	arc( guarantee_no_DNA_ ); // _Bool
	arc( keep_input_protonation_state_ ); // _Bool
	arc( preserve_header_ ); // _Bool
	arc( preserve_crystinfo_ ); // _Bool
	arc( missing_dens_as_jump_ ); // _Bool
	arc( no_chainend_ter_ ); // _Bool
	arc( no_detect_pseudobonds_ ); // _Bool
	arc( no_output_cen_ ); // _Bool
	arc( normalize_to_thk_ ); // _Bool
	arc( output_secondary_structure_ ); // _Bool
	arc( output_torsions_ ); // _Bool
	arc( output_virtual_ ); // _Bool
	arc( output_virtual_zero_occ_ ); // _Bool
	arc( output_only_asymmetric_unit_ ); // _Bool
	arc( pdb_parents_ ); // _Bool
	arc( per_chain_renumbering_ ); // _Bool
	arc( randomize_missing_coords_ ); // _Bool
	arc( remember_unrecognized_res_ ); // _Bool
	arc( remember_unrecognized_water_ ); // _Bool
	arc( renumber_pdb_ ); // _Bool
	arc( read_only_ATOM_entries_ ); // _Bool
	arc( suppress_zero_occ_pdb_output_ ); // _Bool
	arc( auto_detect_glycan_connections_ ); // _Bool
	arc( write_glycan_pdb_codes_ ); // _Bool
	arc( output_alternate_atomids_ ); // _Bool
	arc( output_ligands_as_separate_chains_ );  // _Bool
	arc( maintain_links_ ); // _Bool
	arc( max_bond_length_ ); // core::Real
	arc( min_bond_length_ ); // core::Real
	arc( use_pdb_format_HETNAM_records_ ); // _Bool
	arc( write_pdb_link_records_ ); // _Bool
	arc( write_pdb_parametric_info_ ); // _Bool
	arc( write_all_connect_info_ ); // _Bool
	arc( write_seqres_records_ ); // _Bool
	arc( chains_whose_residues_are_separate_chemical_entities_ ); // std::string
	arc( residues_for_atom_name_remapping_ ); // utility::vector1<std::string>
	arc( pdb_comments_ ); // _Bool
	arc( show_all_fixes_ ); // _Bool
	arc( constraints_from_link_records_ ); // _Bool
	arc( output_pose_energies_table_ ); // _Bool
	arc( output_pose_cache_data_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::io::StructFileRepOptions );
CEREAL_REGISTER_TYPE( core::io::StructFileRepOptions )

CEREAL_REGISTER_DYNAMIC_INIT( core_io_StructFileRepOptions )
#endif // SERIALIZATION
