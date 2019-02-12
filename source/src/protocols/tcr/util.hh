// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/util.hh
/// @brief utility functions for protocols/tcr
/// @author Ragul (ragul@umd.edu)

#ifndef INCLUDED_protocols_tcr_util_hh
#define INCLUDED_protocols_tcr_util_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/tcr/TCRmodel.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/Pose.hh>
#include <utility/file/FileName.hh>
#include <basic/database/open.hh>
#include <core/sequence/MatrixScoringScheme.hh>

namespace protocols {
namespace tcr {

/// @brief initialize aho numbers
/// @details Aho number are standard numbering scheme for TCR structures
/// @details template structures are renumbered by Aho numbers, grafting of TCr segments are based on these numbering scheme
void initialize_aho_numbers(TCRseqInfo::tcrposi &a_aho_pos, TCRseqInfo::tcrposi &b_aho_pos);

/// @brief calculate percent identity between 2 sequences
core::Real calculate_identity_score( std::string const &query, std::string const &content);

/// @brief calculte alignment score b/w 2 sequences of same length
/// @details PAM30 pssm matrix is used for scoring the AA pairs
core::Real score_alignment( std::string const &query, std::string const &content, core::sequence::ScoringSchemeOP const &tcr_ss );

/// @brief repack the tcr model pose
void repack_tcr_model(core::pose::Pose &model_pose);

/// @brief minimize the tcr model pose
void minimize_tcr_model(core::pose::Pose &model_pose);

/// @brief relax the tcr model pose
void relax_tcr_model(core::pose::Pose &model_pose);

/// @brief function to initialize and regular expressions for parsing achain TCR sequence
/// @details different functions for achain and bchain due to different regex used
/// @details due to compiler issues (reqd adv compilers to compile regex) not passing around regex
void assign_achain_CDRs_using_REGEX( std::string const &seq, TCRseqInfo::tcrsegs &parsedseqs, TCRseqInfo::tcrposi &pose_posi );

/// @brief regex intialization for parsing bchain TCR sequence
/// @details bchian TCR segments were assigned using regular expression
void assign_bchain_CDRs_using_REGEX( std::string const &seq, TCRseqInfo::tcrsegs &parsedseqs, TCRseqInfo::tcrposi &pose_posi );

/// @brief TCR segments were assigned using start and end position residue numbers
void assign_CDRs_using_numbers( std::string const &seq, antibody::grafting::CDR_Bounds const &cdr1seqpos, antibody::grafting::CDR_Bounds const &cdr2seqpos, antibody::grafting::CDR_Bounds const &cdr3seqpos, antibody::grafting::CDR_Bounds const &cap, TCRseqInfo::tcrsegs &parsedseqs, TCRseqInfo::tcrposi &pose_posi );

/// @brief adjust pose position for partner chain
/// @details add partner chain length to all values
void adjust_position_for_chain(TCRseqInfo::tcrposi &pose_posi, core::Size const &chainlen );

/// @brief TCR segments were assigned using anarci program and Aho numbering
void assign_CDRs_using_anarci( std::string const &seq, std::string const &anarci_path, TCRseqInfo::tcrposi &a_aho_pos, TCRseqInfo::tcrposi &pose_posi, TCRseqInfo::tcrsegs &parsedseqs );

/// @brief check if query/target sequence are same from a input template db fasta file
/// @details used for deciding whether germline template can be used instead of framework template
/// @details for germline templates only CDR3 loop in grafted into germline template structure
bool check_seq_match( std::string const &query, std::string const &db, std::list< std::set< std::string > > const &ignore_lists );

/// @brief check if query/target sequence are same from multiple input template db fasta file
/// @details calls check_seq_match for each input db file in the list
bool check_seq_match_from_multiple_input_db( std::string const &query, utility::vector1< std::string > const &multidb, std::list< std::set< std::string > > const &ignore_lists );

/// @brief find cdr positions (pose numbers) in input pose with AHO numbering scheme
utility::vector1< std::pair<core::Size, core::Size> > get_cdr_pdb_pos_from_input_pose( core::pose::Pose const &inpose, TCRseqInfo::tcrposi const &selepos);

/// @brief run blast to find templates
/// @details used for blacklisting pdb's based on cutoff
void run_blast( std::string const &query, std::string const &blast_db);

/// @brief convert string to numbers used in CDR_bounds
antibody::grafting::CDR_Bounds string_to_CDRbounds(std::string const &position );


} //protocols
} //tcr

#endif // __ANTIBODY_GRAFTING__

#endif //protocols/tcr_util_hh

