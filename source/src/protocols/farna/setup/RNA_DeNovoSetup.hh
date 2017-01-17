// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/setup/RNA_DeNovoSetup.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_setup_RNA_DeNovoSetup_HH
#define INCLUDED_protocols_farna_setup_RNA_DeNovoSetup_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/farna/setup/RNA_DeNovoSetup.fwd.hh>
#include <protocols/farna/setup/RNA_DeNovoParameters.fwd.hh>
#include <core/pose/rna/RNA_SecStruct.fwd.hh>
#include <protocols/farna/options/RNA_DeNovoProtocolOptions.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/types.hh>
#include <core/pose/rna/BasePair.hh>

namespace protocols {
namespace farna {
namespace setup {

class RNA_DeNovoSetup: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_DeNovoSetup();

	//destructor
	~RNA_DeNovoSetup();

public:

	void
	initialize_from_command_line();

	protocols::farna::options::RNA_DeNovoProtocolOptionsCOP options() const { return  options_; }
	protocols::farna::setup::RNA_DeNovoParametersCOP rna_params() const { return rna_params_; }
	utility::vector1< core::pose::PoseOP > const & refine_pose_list() const { return refine_pose_list_; }
	core::pose::PoseOP pose() { return pose_; }
	core::pose::PoseCOP native_pose() { return native_pose_; }

private:

	void
	de_novo_setup_from_command_line();

	void
	de_novo_setup_from_command_line_legacy();

	void
	setup_refine_pose_list();

	utility::vector1< core::pose::PoseOP >
	get_refine_pose_list( std::string const & input_silent_file,
		std::pair< utility::vector1< int >, utility::vector1< char > > const & output_res_and_chain,
		core::chemical::ResidueTypeSetCOP rsd_set ) const;

	utility::vector1< core::Size >
	working_res_map( utility::vector1< core::Size > const & vec,
		utility::vector1< core::Size > const & working_res,
		bool const leave_out_last_working_residue = false ) const;

	std::string
	working_res_map( std::string const & seq,
		utility::vector1< core::Size > const & working_res,
		bool const annotations_in_brackets = true ) const;

	core::pose::rna::RNA_SecStruct
	working_res_map( core::pose::rna::RNA_SecStruct const & seq,
		utility::vector1< core::Size > const & working_res ) const;

	void
	get_seq_and_resnum( std::string const & pdb,
		std::string & seq,
		utility::vector1< int > & resnum,
		utility::vector1< char > & chain ) const;

	std::string
	get_silent_seq( std::string const & silent_file ) const;

	std::pair< utility::vector1< int >, utility::vector1< char > >
	get_silent_resnum( std::string const & silent_file ) const;

	bool
	already_listed_in_obligate_pair( utility::vector1< core::Size > const & new_pair,
		utility::vector1< core::Size > const & obligate_pair ) const;

	bool
	already_listed_in_obligate_pair( utility::vector1< Size > const & new_pair,
		utility::vector1< Size > const & obligate_pair,
		utility::vector1< Size > const & obligate_pair_explicit ) const;

	void
	update_working_obligate_pairs_with_stems(
		utility::vector1< core::pose::rna::BasePair > & working_obligate_pairs,
		utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > const & working_stems,
		utility::vector1< core::Size > const & working_input_res ) const;

private:

	protocols::farna::options::RNA_DeNovoProtocolOptionsOP options_;
	protocols::farna::setup::RNA_DeNovoParametersOP rna_params_;
	core::chemical::ResidueTypeSetCOP rsd_set_;
	core::pose::PoseOP pose_, native_pose_;

	utility::vector1< core::pose::PoseOP > refine_pose_list_; // does this belong here?

};

} //setup
} //farna
} //protocols

#endif
