// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/modeler/rna/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_legacy_modeler_rna_util_HH
#define INCLUDED_protocols_stepwise_legacy_modeler_rna_util_HH

#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

bool
check_can_prepend( utility::vector1< Size > const & seq_num_list );

bool
check_can_append( utility::vector1< Size > const & seq_num_list );


void
output_pair_size( std::pair < Size, Size > const & pair_size, std::ostream & outstream = std::cout );

void
output_pair_size( utility::vector1 < std::pair < Size, Size > > const & pair_size_vector, std::string const & output_string, std::ostream & outstream = std::cout, Size const spacing = 40 );

void output_is_prepend_map( std::string const & tag, std::map< Size, bool > const & my_map, Size const max_seq_num, std::ostream & outstream = std::cout, Size const tag_spacing = 40 );

void
output_bool_list( std::string const & tag, utility::vector1< bool > const & bool_list, std::ostream & outstream = std::cout, Size const spacing = 40 );

void
output_bool_list( std::string const & tag, utility::vector1< Size > const & size_list, std::ostream & outstream = std::cout, Size const spacing = 40 );

void
output_size_list( std::string const & tag, utility::vector1< Size > const & size_list, std::ostream & outstream = std::cout, Size const spacing = 40 );

void
sort_pair_list( utility::vector1< std::pair < Size, Size > > pair_list );

void
output_fold_tree_info( kinematics::FoldTree const & fold_tree, std::string const & pose_name, std::ostream & outstream = std::cout );

void
output_fold_tree_info( pose::Pose const & pose, std::string pose_name, std::ostream & outstream = std::cout );


Real
full_length_rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, utility::vector1 < Size > const & residue_list, std::string const & full_sequence, bool const verbose, bool const ignore_virtual_atom );


void
print_backbone_torsions( pose::Pose const & pose, Size five_prime_chainbreak );

core::Size
setup_chain_break_jump_point( core::pose::Pose & pose,
															core::Size const moving_res,
															core::Size const reference_res );

void
remove_chain_break_jump_point( core::pose::Pose & pose,
															 core::Size const moving_res,
															 core::Size const reference_res );

Size
setup_bulge_jump_point( pose::Pose & pose, Size const & moving_base, Size const & reference_base, bool verbose = false );

void
apply_rotamer( pose::Pose & pose, utility::vector1< protocols::stepwise::modeler::rna::Torsion_Info >  const & rotamer_list );

bool
is_same_sugar_pucker( pose::Pose const & current_pose, pose::Pose const & cluster_center_pose, Size const seq_num );

void
setup_simple_fold_tree( pose::Pose & pose );

void
import_pose_from_silent_file(
    pose::Pose & import_pose,
    std::string const & silent_file,
    std::string const & input_tag );

std::string
get_tag_from_pdb_filename( std::string const & pdb_filename );

void
move_jump_atom_to_base(
    kinematics::FoldTree & fold_tree,
    std::string const & working_sequence );

void
print_WorkingParameters_info( protocols::stepwise::modeler::working_parameters::StepWiseWorkingParametersOP const & WP, std::string const & WP_name, std::ostream & outstream = std::cout, bool const is_simple_full_length_WP = false );

void
print_WorkingParameters_info( protocols::stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP const & const_WP, std::string const & WP_name, std::ostream & outstream = std::cout, bool const is_simple_full_length_WP = false );

void
set_nucleotide_to_A_form( pose::Pose & pose, Size const seq_num );

std::string
path_basename( std::string const & full_path );

} //rna
} //modeler
} //legacy
} //stepwise
} //protocols

#endif
