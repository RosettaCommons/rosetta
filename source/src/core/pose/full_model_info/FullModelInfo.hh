// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/FullModelInfo.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

#ifndef INCLUDED_core_pose_full_model_info_FullModelInfo_hh
#define INCLUDED_core_pose_full_model_info_FullModelInfo_hh

#include <core/types.hh>

// Project headers
#ifdef WIN32
#include <core/pose/Pose.hh>
#endif
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pose/full_model_info/FullModelInfo.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/pose/full_model_info/SubMotifInfo.hh>


// C++
#include <string>
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace full_model_info {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of all information related to how a subpose 'fits in' to global modeling scheme.
///
/// @details  See FullModelParameters for more information on the global modeling scheme, including
///           full_sequence, etc.
///
class FullModelInfo: public basic::datacache::CacheableData  {

public:


	FullModelInfo( std::string const & full_sequence,
		utility::vector1< Size > const & cutpoint_open_in_full_model,
		utility::vector1< Size > const & res_numbers_in_pose );

	FullModelInfo( FullModelParametersCOP full_model_parameters );

	FullModelInfo( pose::Pose & pose );

	FullModelInfo( FullModelInfo const & src );

	~FullModelInfo();

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new FullModelInfo( *this ) );
	}

	FullModelInfoOP
	clone_info() const
	{
		return FullModelInfoOP( new FullModelInfo( *this ) );
	}

	// properties of full model.
	FullModelParametersCOP full_model_parameters() const;
	void set_full_model_parameters( FullModelParametersCOP setting );

	std::string const & full_sequence() const;
	utility::vector1< int > const & conventional_numbering() const;
	utility::vector1< char > const & conventional_chains() const;
	utility::vector1< Size > const & cutpoint_open_in_full_model() const;
	utility::vector1< Size > const & dock_domain_map() const;
	utility::vector1< Size > const & fixed_domain_map() const;
	utility::vector1< Size > const & input_domain_map() const;
	utility::vector1< Size > const & extra_minimize_res() const;
	utility::vector1< Size > const & preferred_root_res() const;
	utility::vector1< Size > const & sample_res() const;
	utility::vector1< Size > const & working_res() const;
	utility::vector1< Size > const & calc_rms_res() const;
	utility::vector1< Size > const & rna_syn_chi_res() const;
	utility::vector1< Size > const & rna_anti_chi_res() const;
	utility::vector1< Size > const & rna_terminal_res() const;
	utility::vector1< Size > const & rna_block_stack_above_res() const;
	utility::vector1< Size > const & rna_block_stack_below_res() const;
	utility::vector1< Size > const & rna_bulge_res() const;
	utility::vector1< Size > const & alignment_anchor_res() const;
	utility::vector1< Size > const & rna_sample_sugar_res() const;
	utility::vector1< std::pair< Size, Size > > jump_pairs() const;
	utility::vector1< Size > const & jump_res_map() const;

	void clear_res_list();

	void clear_other_pose_list();

	utility::vector1< core::pose::PoseOP > const & other_pose_list() const { return other_pose_list_; }

	utility::vector1< Size > const & res_list() const { return res_list_; }

	utility::vector1< Size > full_to_sub( utility::vector1< Size > const & res_in_full_model_numbering ) const;

	Size full_to_sub( Size const res_in_full_model_numbering ) const;

	utility::vector1< Size > sub_to_full( utility::vector1< Size > const & res ) const;

	std::map< Size, Size > full_to_sub() const;

	Size sub_to_full( Size const & res ) const;

	void add_other_pose( core::pose::PoseOP pose );


	void set_other_pose_list( utility::vector1< pose::PoseOP > const & setting );

	void set_res_list( utility::vector1< Size > const & res_list ){ res_list_ = res_list; }

	void remove_other_pose_at_idx( Size const idx );

	Size find_index_in_other_pose_list( pose::Pose const & pose ) const;

	Size get_idx_for_other_pose_with_residues( utility::vector1< Size > const & input_res ) const;

	Size get_idx_for_other_pose_with_residue( Size const input_res ) const;

	Size get_idx_for_other_pose( pose::Pose const & pose ) const;

	utility::vector1< Size > chains_in_full_model() const;

	Size size() const;

	// submotif accessors/modifiers
	utility::vector1< SubMotifInfoOP > & submotif_info_list() const { return submotif_info_list_; }

	bool is_a_submotif() const;

	bool is_a_submotif_seed() const;

	bool is_a_submotif( utility::vector1< Size > const & res_list, bool const & check_other_poses = false ) const;

	bool is_a_submotif_seed( utility::vector1< Size > const & res_list, bool const & check_other_poses = false ) const;

	SubMotifInfoOP submotif_info( utility::vector1< Size > const & res_list, bool const & check_other_poses = false ) const;

	bool in_a_submotif( utility::vector1< Size > const & res_list, bool const & check_other_poses = false ) const;

	bool in_a_submotif_seed( utility::vector1< Size > const & res_list, bool const & check_other_poses = false ) const;

	SubMotifInfoOP submotif_info_containing_residues( utility::vector1< Size > const & res_list, bool const & check_other_poses = false ) const;

	void update_submotif_info_list();

	void add_submotif_info(
		utility::vector1< Size > const & res_list,
		std::string const & tag,
		bool const & seed = false );

	void add_submotif_info( utility::vector1< SubMotifInfoOP > submotif_info_list );

	void add_submotif_info( SubMotifInfoOP submotif_info_op );

	void delete_submotif_info(
		utility::vector1< Size > const & res_list,
		std::string const & tag,
		bool const & seed = false );

	void delete_submotif_info( utility::vector1< SubMotifInfoOP > submotif_info_list );

	void delete_submotif_info( SubMotifInfoOP submotif_info_op );

	void show_submotif_info_list( bool const & for_all_poses = false ) const;

private:

	utility::vector1< Size >
	get_cutpoint_open_from_pdb_info( pose::Pose const & pose ) const;

	void
	get_sequence_with_gaps_filled_with_n( pose::Pose const & pose,
		std::string & sequence,
		utility::vector1< Size > & full_numbering ) const;

private:

	// residues that go with pose. In principle, this is redundant with PDBInfo,
	// but not in eventual case where user has favorite numbering/chain scheme.
	utility::vector1< Size > res_list_;

	// what's known about this pose and any neighbors in a "PoseTree"
	utility::vector1< core::pose::PoseOP > other_pose_list_;

	// FullModelParameters: properties of full model.
	// Note that this holds full_sequence, cutpoint_open_in_full_model,
	//  conventional numbering scheme, etc.
	// And it should not change as pieces are added/deleted from the pose!
	FullModelParametersCOP full_model_parameters_;

	// SubMotifInfo: information about submotifs in poses;
	mutable utility::vector1< SubMotifInfoOP > submotif_info_list_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	FullModelInfo();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};

FullModelInfo const &
const_full_model_info( pose::Pose const & pose );

FullModelInfo &
nonconst_full_model_info( core::pose::Pose & pose );

// basically an alias to nonconst_full_model_info.
FullModelInfo const &
make_sure_full_model_info_is_setup( pose::Pose & pose );

bool
full_model_info_defined( pose::Pose const & pose );

void
set_full_model_info( pose::Pose & pose, FullModelInfoOP & full_model_info );

void
update_full_model_info_from_pose( pose::Pose & pose );

} //full_model_info
} //pose
} //core
#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_full_model_info_FullModelInfo )
#endif // SERIALIZATION


#endif
