// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_ResidueInfo.hh
/// @brief
/// @detailed
///
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_ResidueInfo_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_ResidueInfo_HH


#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <string>
#include <map>


namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {

	struct Residue_info{
 	  std::string name;
  	core::Size seq_num; //Full_pose_sequence_num
	};

	void
	print_torsion_info( core::pose::Pose const & pose, core::Size const seq_num, core::Size const rna_torsion_number, std::string const type );

	utility::vector1 < Residue_info >
	Convert_rebuild_residue_string_to_list( std::string const & rebuild_residue_string );

	void
	output_residue_struct( Residue_info const & residue );

	std::string
	get_one_letter_name( std::string const & three_letter_name );

	std::string
	get_three_letter_name( std::string const & one_letter_name );

	core::Size
	get_max_seq_num_from_res_map( std::map< core::Size, core::Size > const & my_map );


	void
	output_res_map( std::map< core::Size, core::Size > const & my_map, core::Size const max_seq_num );

	// Undefined, commenting out to fix PyRosetta build  void output_is_prepend_map(std::map< core::Size, bool > const & my_map, core::Size const max_seq_num);

	void
	output_residue_list( utility::vector1< Residue_info > residue_list );

	utility::vector1< Residue_info >
	get_residue_list_from_fasta( std::string const full_fasta_sequence );

	Residue_info
	get_residue_from_seq_num( core::Size const & seq_num, utility::vector1 < Residue_info > const & residue_list );

	bool
	Contain_residue_at_seq_num( core::Size seq_num, utility::vector1 < Residue_info > const & residue_list );

	utility::vector1 < utility::vector1 < Residue_info > >
	Create_strand_list( utility::vector1 < Residue_info > const & residue_list );

	utility::vector1 < Residue_info >
	Set_Difference( utility::vector1 < Residue_info > const & residue_list_1, utility::vector1 < Residue_info > const & residue_list_2 );

	utility::vector1 < Residue_info >
	Set_Union( utility::vector1 < Residue_info > const & residue_list_1, utility::vector1 < Residue_info > const & residue_list_2 );

	bool
	residue_list_sort_criterion( Residue_info residue_info_1, Residue_info residue_info_2 );

	void
	sort_residue_list( utility::vector1< Residue_info > & residue_list );


} //rna
} //enumerate
} //stepwise
} //protocols

#endif
