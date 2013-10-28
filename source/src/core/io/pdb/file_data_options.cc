// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data_options.cc
/// @brief  Definitions for FileDataOptions.
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Unit headers
#include <core/io/pdb/file_data_options.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>

// Utility headers
#include <utility/tag/Tag.hh>

namespace core {
namespace io {
namespace pdb {

FileDataOptions::FileDataOptions() { init_from_options(); }

FileDataOptions::~FileDataOptions() {}

void FileDataOptions::parse_my_tag( utility::tag::TagCOP tag )
{
	//set_check_if_residues_are_termini( tag->getOption< std::string >( "termini", "ALL" ));
	set_check_if_residues_are_Ntermini( tag->getOption< std::string >( "Ntermini", "ALL" ));
	set_check_if_residues_are_Ctermini( tag->getOption< std::string >( "Ctermini", "ALL" ));
	set_exit_if_missing_heavy_atoms( tag->getOption< bool >( "exit_if_missing_heavy_atoms", 0 ));
	set_ignore_unrecognized_res( tag->getOption< bool >( "ignore_unrecognized_res", 0 ));
	set_ignore_waters( tag->getOption< bool >( "ignore_waters", 0 ));
	set_ignore_zero_occupancy( tag->getOption< bool >( "ignore_zero_occupancy", 1 ));
	set_keep_input_protonation_state( tag->getOption< bool >( "keep_input_protonation_state", 0 ));
	set_preserve_header( tag->getOption< bool >( "preserve_header", 0 ));
	set_preserve_crystinfo( tag->getOption< bool >( "preserve_crystinfo", 0 ));
	set_missing_dens_as_jump( tag->getOption< bool >( "missing_dens_as_jump", 0 ));
	set_randomize_missing_coords( tag->getOption< bool >( "randomize_missing_coords", 0 ));
	set_remember_unrecognized_res( tag->getOption< bool >( "remember_unrecognized_res", 0 ));
	set_remember_unrecognized_water( tag->getOption< bool >( "remember_unrecognized_water", 0 ));
	set_write_pdb_link_records(tag->getOption<bool>("write_pdb_link_records", 0));

	set_chains_whose_residues_are_separate_chemical_entities(
		tag->getOption< std::string >( "treat_residues_in_these_chains_as_separate_chemical_entities", " " ));
}

std::string FileDataOptions::type() const { return "file_data_options"; }


// accessors
//std::string FileDataOptions::check_if_residues_are_termini() const { return check_if_residues_are_termini_; }
std::string FileDataOptions::check_if_residues_are_Ntermini() const { return check_if_residues_are_Ntermini_; }
std::string FileDataOptions::check_if_residues_are_Ctermini() const { return check_if_residues_are_Ctermini_; }
bool FileDataOptions::exit_if_missing_heavy_atoms() const { return exit_if_missing_heavy_atoms_; }
bool FileDataOptions::ignore_unrecognized_res() const { return ignore_unrecognized_res_; }
bool FileDataOptions::ignore_waters() const { return ignore_waters_; }
bool FileDataOptions::ignore_zero_occupancy() const { return ignore_zero_occupancy_; }
bool FileDataOptions::keep_input_protonation_state() const { return keep_input_protonation_state_; }
bool FileDataOptions::preserve_header() const { return preserve_header_; }
bool FileDataOptions::preserve_crystinfo() const { return preserve_crystinfo_; }
bool FileDataOptions::missing_dens_as_jump() const { return missing_dens_as_jump_; }
bool FileDataOptions::randomize_missing_coords() const { return randomize_missing_coords_; }
bool FileDataOptions::remember_unrecognized_res() const { return remember_unrecognized_res_; }
bool FileDataOptions::remember_unrecognized_water() const { return remember_unrecognized_water_; }
bool FileDataOptions::write_pdb_link_records() const {return write_pdb_link_records_;}
std::string const & FileDataOptions::chains_whose_residues_are_separate_chemical_entities() const { return chains_whose_residues_are_separate_chemical_entities_; }

// mutators

void FileDataOptions::set_check_if_residues_are_Ntermini( std::string check_if_residues_are_Ntermini )
{ check_if_residues_are_Ntermini_ = check_if_residues_are_Ntermini; }

void FileDataOptions::set_check_if_residues_are_Ctermini( std::string check_if_residues_are_Ctermini )
{ check_if_residues_are_Ctermini_ = check_if_residues_are_Ctermini; }

void FileDataOptions::set_exit_if_missing_heavy_atoms( bool exit_if_missing_heavy_atoms )
{ exit_if_missing_heavy_atoms_ = exit_if_missing_heavy_atoms; }

void FileDataOptions::set_ignore_unrecognized_res( bool ignore_unrecognized_res )
{ ignore_unrecognized_res_ = ignore_unrecognized_res; }

void FileDataOptions::set_ignore_waters( bool ignore_waters )
{ ignore_waters_ = ignore_waters; }

void FileDataOptions::set_ignore_zero_occupancy( bool ignore_zero_occupancy )
{ ignore_zero_occupancy_ = ignore_zero_occupancy; }

void FileDataOptions::set_keep_input_protonation_state( bool keep_input_protonation_state )
{ keep_input_protonation_state_ = keep_input_protonation_state; }

void FileDataOptions::set_preserve_header( bool preserve_header )
{ preserve_header_ = preserve_header; }

void FileDataOptions::set_preserve_crystinfo( bool preserve_crystinfo )
{ preserve_crystinfo_ = preserve_crystinfo; }

void FileDataOptions::set_missing_dens_as_jump( bool missing_dens_as_jump )
{ missing_dens_as_jump_ = missing_dens_as_jump; }

void FileDataOptions::set_randomize_missing_coords( bool randomize_missing_coords )
{ randomize_missing_coords_ = randomize_missing_coords; }

void FileDataOptions::set_remember_unrecognized_res( bool remember_unrecognized_res )
{ remember_unrecognized_res_ = remember_unrecognized_res; }

void FileDataOptions::set_remember_unrecognized_water( bool remember_unrecognized_water )
{ remember_unrecognized_water_ = remember_unrecognized_water; }

void FileDataOptions::set_write_pdb_link_records(bool setting)
{
	write_pdb_link_records_ = setting;
}

void FileDataOptions::set_chains_whose_residues_are_separate_chemical_entities( std::string const & chains_whose_residues_are_separate_chemical_entities )
{ chains_whose_residues_are_separate_chemical_entities_ = chains_whose_residues_are_separate_chemical_entities; }


void FileDataOptions::init_from_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_check_if_residues_are_Ntermini( option[ in::Ntermini ].value());
	set_check_if_residues_are_Ctermini( option[ in::Ctermini ].value());
	set_exit_if_missing_heavy_atoms( option[ run::exit_if_missing_heavy_atoms ].value());
	set_ignore_unrecognized_res( option[ in::ignore_unrecognized_res ]());
	set_ignore_waters( option[ in::ignore_waters ]());
	set_ignore_zero_occupancy( option[ run::ignore_zero_occupancy ]());
	set_keep_input_protonation_state( option[ pH::keep_input_protonation_state ]());
	set_preserve_header( option[ run::preserve_header ].value());
	set_preserve_crystinfo( option[ in::preserve_crystinfo ]() );
	set_missing_dens_as_jump( option[ in::missing_density_to_jump ]() );
	set_randomize_missing_coords( option[ run::randomize_missing_coords ]());
	set_remember_unrecognized_res( option[ in::remember_unrecognized_res ]());
	set_remember_unrecognized_water( option[ in::remember_unrecognized_water ]());
	set_write_pdb_link_records(option[out::file::write_pdb_link_records]());
	set_chains_whose_residues_are_separate_chemical_entities( option[ in::file::treat_residues_in_these_chains_as_separate_chemical_entities].user_or(""));
}

} // namespace pdb
} // namespace io
} // namespace core
