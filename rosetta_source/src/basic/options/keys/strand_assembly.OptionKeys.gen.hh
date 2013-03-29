// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/strand_assembly.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_strand_assembly_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_strand_assembly_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace strand_assembly { extern BooleanOptionKey const strand_assembly; }
namespace strand_assembly { extern IntegerOptionKey const min_num_strands_to_deal; }
namespace strand_assembly { extern IntegerOptionKey const max_num_strands_to_deal; }
namespace strand_assembly { extern BooleanOptionKey const extract_native_only; }
namespace strand_assembly { extern BooleanOptionKey const extract_sandwich; }
namespace strand_assembly { extern IntegerOptionKey const min_res_in_strand; }
namespace strand_assembly { extern IntegerOptionKey const max_res_in_strand; }
namespace strand_assembly { extern RealOptionKey const min_CA_CA_dis; }
namespace strand_assembly { extern RealOptionKey const max_CA_CA_dis; }
namespace strand_assembly { extern RealOptionKey const min_O_N_dis; }
namespace strand_assembly { extern RealOptionKey const max_O_N_dis; }
namespace strand_assembly { extern RealOptionKey const min_C_O_N_angle; }
namespace strand_assembly { extern RealOptionKey const min_sheet_dis; }
namespace strand_assembly { extern RealOptionKey const max_sheet_dis; }
namespace strand_assembly { extern RealOptionKey const min_sheet_torsion; }
namespace strand_assembly { extern RealOptionKey const max_sheet_torsion; }
namespace strand_assembly { extern RealOptionKey const min_sheet_torsion_cen_res; }
namespace strand_assembly { extern RealOptionKey const max_sheet_torsion_cen_res; }
namespace strand_assembly { extern RealOptionKey const min_sheet_angle; }
namespace strand_assembly { extern RealOptionKey const max_sheet_angle; }
namespace strand_assembly { extern IntegerOptionKey const min_num_strands_in_sheet; }
namespace strand_assembly { extern RealOptionKey const min_inter_sheet_dis_CA_CA; }
namespace strand_assembly { extern RealOptionKey const max_inter_sheet_dis_CA_CA; }
namespace strand_assembly { extern RealOptionKey const min_shortest_dis_sidechain_inter_sheet; }
namespace strand_assembly { extern BooleanOptionKey const write_chain_B_resnum; }
namespace strand_assembly { extern BooleanOptionKey const write_phi_psi; }
namespace strand_assembly { extern IntegerOptionKey const max_starting_loop_size; }
namespace strand_assembly { extern IntegerOptionKey const max_ending_loop_size; }
namespace strand_assembly { extern BooleanOptionKey const no_helix_in_pdb; }
namespace strand_assembly { extern IntegerOptionKey const max_helix_in_extracted_sw_loop; }
namespace strand_assembly { extern IntegerOptionKey const max_E_in_extracted_sw_loop; }
namespace strand_assembly { extern BooleanOptionKey const exclude_sandwich_that_is_linked_w_same_direction_strand; }
namespace strand_assembly { extern RealOptionKey const max_inter_strand_angle_to_not_be_same_direction_strands; }
namespace strand_assembly { extern RealOptionKey const max_abs_inter_strand_dihedral_to_not_be_same_direction_strands; }
namespace strand_assembly { extern IntegerOptionKey const max_num_sw_per_pdb; }
namespace strand_assembly { extern StringOptionKey const check_N_to_C_direction_by; }
namespace strand_assembly { extern BooleanOptionKey const do_not_connect_sheets_by_loops; }
namespace strand_assembly { extern RealOptionKey const check_canonicalness_cutoff; }
namespace strand_assembly { extern BooleanOptionKey const count_AA_with_direction; }
namespace strand_assembly { extern RealOptionKey const inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets; }
namespace strand_assembly { extern BooleanOptionKey const exclude_desinated_pdbs; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
