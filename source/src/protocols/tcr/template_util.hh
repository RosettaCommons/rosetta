// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/template_util.hh
/// @brief template identification utility functions for protocols/tcr/TCRmodel
/// @author Ragul (ragul@umd.edu)

#ifndef INCLUDED_protocols_tcr_template_util_hh
#define INCLUDED_protocols_tcr_template_util_hh

#include <protocols/tcr/TCRmodel.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/TCRloopRefine.hh>
#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/sequence/ScoringScheme.hh>

namespace protocols {
namespace tcr {

/// @brief initialize template files for all TCR segments
/// @details template files for CDR sequences were choosen by CDR sequence length
void initialize_template_db_files(TCRseqInfo::tcrsegs const &atcr, TCRseqInfo::tcrsegs const &btcr, TCRmodel::tcrtmplts &atmplt, TCRmodel::tcrtmplts &btmplt, std::string const &tcr_seq_db);

/// @brief check and update if template files based on sequence length exists on the template
utility::vector1<std::string> update_template_files(utility::vector1<std::string> const &file_names, std::string const &tcr_seq_db);

/// @brief initialize list of template structures to ignore
/// @details ignore list can be created by user provide list with flag -ignore_list
void initialize_template_ignore_list(TCRseqInfo::tcrsegs const &atcr, TCRseqInfo::tcrsegs const &btcr, std::list< std::set< std::string > > &ignore_lists, core::Real const &blast_cutoff, std::string const &tcr_seq_db);

/// @brief create template blacklist by runing blast with cut-off value
std::set< std::string > create_ignore_list(std::string const query, std::string const blast_db, core::Real const &blast_cutoff);

/// @brief find best scoring template from a fasta sequence file
/// @details template choosen by best alignment score b/w input query and seq from template file
std::pair<core::Real, std::string> find_template(std::string const &query, std::string const &db, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss);

/// @brief find templates from one or more template fasta sequence files
/// @details calls find_template for each file in the list
std::pair<core::Real, std::string> find_template_from_multiple_input_db(std::string const &query, utility::vector1< std::string > const &multidb, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss);

/// @brief Check if a particular structure should be ignored for consideration as template
/// @details checks if a template id is present in the ignore list
bool ignore_template(std::string const &tmpltid, std::list< std::set< std::string > > const &ignore_lists);

/// @brief get the CDR template structure from user provided PDB ID
void setup_cdr_template_from_pdbid(TCRmodel::tmpltinfo &currtmplt, core::Size const cdrstart, core::Size const cdrend, std::string const &tpdb_path);

/// @brief read the template structure from user provided PDB file
void setup_template_from_pdbfile(TCRmodel::tmpltinfo &currtmplt);

/// @brief Search templates for a CDR loop sequence in the template database
void setup_cdr_template_from_db(std::string const &seq, TCRmodel::tmpltinfo &currtmplt, core::Size const cdrstart, core::Size const cdrend, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss);

/// @brief CDR template identification and grafting of template structure
/// @brief parent function to get the CDR template from PDB file or PDB ID or template databse
void setup_cdr_template(std::string const &seq, TCRmodel::tmpltinfo &currtmplt, core::Size const cdrstart, core::Size const cdrend, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss);

/// @brief helper method to read the template structure from user provided PDB id
/// @details for framework/germline template structures
void setup_fw_template_from_pdbid(TCRmodel::tmpltinfo &currtmplt, TCRseqInfo::tcrposi const &aho_pos, std::string const &tpdb_path);

/// @brief Search template structures in the template database
/// @details for Framework/Germline templates
void setup_fw_template_from_db(std::string const &seq, TCRseqInfo::tcrposi const &aho_pos, TCRmodel::tmpltinfo &currtmplt, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss);

/// @brief Framework/germline template identification and grafting of template structure
/// @details parent function to get the FE/GM template from PDB file or PDB ID or template databse
void setup_fw_template(std::string const &seq, TCRseqInfo::tcrposi const &aho_pos, TCRmodel::tmpltinfo &currtmplt, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss);

/// @brief orientation template identification and grafting from template db
std::pair<std::string, std::string> get_orientation_template_from_db(std::string const &fra, std::string const &frb, TCRmodel::tmpltinfo const &oritmplt, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &ss);

/// @brief parent function for orientation template identification
/// @details checks for user provided template pdb id or structure before searching for templates
void setup_orientation_template(TCRseqInfo::tcrsegs const &atcr, TCRseqInfo::tcrsegs const &btcr, TCRmodel::tmpltinfo &oria_tmplt, TCRmodel::tmpltinfo &orib_tmplt, TCRseqInfo::tcrposi const &aaho_posi, TCRseqInfo::tcrposi const &baho_posi, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &ss);

/// @brief dump all templates found for various tcr segments
void dump_templates(TCRmodel::tcrtmplts const &atmplt, TCRmodel::tcrtmplts const &btmplt, bool const use_gma, bool const use_gmb);

} //protocols
} //tcr


#endif //protocols/tcr_template_util_hh

