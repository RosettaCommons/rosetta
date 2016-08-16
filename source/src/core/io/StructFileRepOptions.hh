// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/pdb/file_data_options.hh
/// @brief  Declarations for StructFileRepOptions.
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_io_StructFileRepOptions_HH
#define INCLUDED_core_io_StructFileRepOptions_HH


// Unit headers
#include <core/io/StructFileRepOptions.fwd.hh>

// Core headers
#include <core/types.hh>

// Basic headers
#include <basic/resource_manager/ResourceOptions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>


namespace core {
namespace io {

class StructFileRepOptions : public basic::resource_manager::ResourceOptions
{
public:
	/// @brief Constructor that takes default values from the global OptionCollection object, basic::options::option.
	StructFileRepOptions();
	/// @brief Constructor that takes default values from a provided OptionCollection object
	StructFileRepOptions( utility::options::OptionCollection const & options );

	virtual ~StructFileRepOptions();

	/// @brief Copy this object and return an owning pointer to the copy.
	virtual StructFileRepOptionsOP clone() const;

	virtual void parse_my_tag( utility::tag::TagCOP tag );

	virtual std::string type() const;

	// accessors
	std::string const & check_if_residues_are_Ntermini() const;
	std::string const & check_if_residues_are_Ctermini() const;
	bool skip_connect_info() const;
	core::Real connect_info_cutoff() const;
	bool exit_if_missing_heavy_atoms() const;
	bool fold_tree_io() const;
	bool ignore_unrecognized_res() const;
	bool ignore_sugars() const;
	bool ignore_waters() const;
	bool ignore_zero_occupancy() const;
	bool guarantee_no_DNA() const;
	bool keep_input_protonation_state() const;
	bool preserve_header() const;
	bool preserve_crystinfo() const;
	bool missing_dens_as_jump() const;
	bool no_chainend_ter() const;
	bool no_output_cen() const;
	bool normalize_to_thk() const;
	bool output_torsions() const;
	bool output_virtual() const;
	bool output_virtual_zero_occ() const;
	bool pdb_parents() const;
	bool per_chain_renumbering() const;
	bool randomize_missing_coords() const;
	bool remember_unrecognized_res() const;
	bool remember_unrecognized_water() const;
	bool renumber_pdb() const;
	bool suppress_zero_occ_pdb_output() const;
	bool use_pdb_format_HETNAM_records() const;
	bool write_pdb_link_records() const;
	bool write_pdb_parametric_info() const;
	bool write_all_connect_info() const;
	std::string const & chains_whose_residues_are_separate_chemical_entities() const;
	utility::vector1<std::string> const & residues_for_atom_name_remapping() const;
	bool pdb_comments() const;
	bool show_all_fixes() const;
	bool constraints_from_link_records() const;

	// mutators
	void set_check_if_residues_are_Ntermini( std::string const & check_if_residues_are_Ntermini );
	void set_check_if_residues_are_Ctermini( std::string const & check_if_residues_are_Ctermini );
	void set_skip_connect_info( bool const skip_connect_info );
	void set_connect_info_cutoff( core::Real const & connect_info_cutoff );
	void set_exit_if_missing_heavy_atoms( bool const exit_if_missing_heavy_atoms );
	void set_fold_tree_io( bool const fold_tree_io );
	void set_ignore_unrecognized_res( bool const ignore_unrecognized_res );
	void set_ignore_sugars( bool const setting );
	void set_ignore_waters( bool const ignore_waters );
	void set_ignore_zero_occupancy( bool const ignore_zero_occupancy );
	void set_guarantee_no_DNA( bool const guarantee_no_DNA );
	void set_keep_input_protonation_state( bool const keep_input_protonation_state );
	void set_preserve_header( bool const preserve_header );
	void set_preserve_crystinfo( bool const preserve_crystinfo );
	void set_missing_dens_as_jump( bool const missing_dens_as_jump );
	void set_no_chainend_ter( bool const no_chainend_ter );
	void set_no_output_cen( bool const no_output_cen );
	void set_normalize_to_thk( bool const normalize_to_thk );
	void set_output_torsions( bool const output_torsions );
	void set_output_virtual( bool const output_virtual );
	void set_output_virtual_zero_occ( bool const output_virtual_zero_occ );
	void set_pdb_parents( bool const pdb_parents );
	void set_per_chain_renumbering( bool const per_chain_renumbering );
	void set_randomize_missing_coords( bool const randomize_missing_coords );
	void set_remember_unrecognized_res( bool const remember_unrecognized_res );
	void set_remember_unrecognized_water( bool const remember_unrecognized_water );
	void set_renumber_pdb( bool const setting );
	void set_suppress_zero_occ_pdb_output( bool const setting );
	void set_use_pdb_format_HETNAM_records( bool const setting );
	void set_write_pdb_link_records( bool const setting );
	void set_write_pdb_parametric_info( bool const setting );
	void set_write_all_connect_info( bool const setting );
	void set_chains_whose_residues_are_separate_chemical_entities( std::string const & setting );
	void set_residues_for_atom_name_remapping( utility::vector1<std::string> const & setting );
	void set_pdb_comments( bool const pdb_comments );
	void set_show_all_fixes( bool setting );
	void set_constraints_from_link_records( bool setting );

	/// @brief Declare the list of options that are read in the process of reading a PDB (or SDF) and converting
	/// it into a Pose.
	static
	void
	list_options_read( utility::options::OptionKeyList & read_options );

	/// @brief Describe the XML Schema for this ResourceOption object
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options( utility::options::OptionCollection const & options );

private:
	std::string check_if_residues_are_Ntermini_; // DEFAULT "ALL" chains. String of 1-letter chains to apply termini to
	std::string check_if_residues_are_Ctermini_; // DEFAULT "ALL" chains. String of 1-letter chains to apply termini to
	bool skip_connect_info_;
	core::Real connect_info_cutoff_;
	bool exit_if_missing_heavy_atoms_;
	bool fold_tree_io_;
	bool ignore_unrecognized_res_;
	bool ignore_sugars_;
	bool ignore_waters_;
	bool ignore_zero_occupancy_;
	bool guarantee_no_DNA_;
	bool keep_input_protonation_state_;
	bool preserve_header_;
	bool preserve_crystinfo_;
	bool missing_dens_as_jump_;
	bool no_chainend_ter_;
	bool no_output_cen_;
	bool normalize_to_thk_;
	bool output_torsions_;
	bool output_virtual_;
	bool output_virtual_zero_occ_;
	bool pdb_parents_;
	bool per_chain_renumbering_;
	bool randomize_missing_coords_;
	bool remember_unrecognized_res_;
	bool remember_unrecognized_water_;
	bool renumber_pdb_;
	bool suppress_zero_occ_pdb_output_;
	bool use_pdb_format_HETNAM_records_;
	bool write_pdb_link_records_;
	bool write_pdb_parametric_info_;
	bool write_all_connect_info_;

	std::string chains_whose_residues_are_separate_chemical_entities_; //treat_residues_in_these_chains_as_separate_chemical_entities
	/// @brief Three letter codes of residues for which to allow atom renaming.
	utility::vector1< std::string > residues_for_atom_name_remapping_;

	bool pdb_comments_;
	bool show_all_fixes_;
	bool constraints_from_link_records_;
};

} // namespace io
} // namespace core

#endif // INCLUDED_core_io_StructFileRepOptions_HH
