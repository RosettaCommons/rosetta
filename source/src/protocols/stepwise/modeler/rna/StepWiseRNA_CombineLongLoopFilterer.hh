// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseRNA_CombineLongLoopFilterer.hh
/// @brief
/// @details
///
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)
/// @author Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_CombineLongLoopFilterer_hh
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_CombineLongLoopFilterer_hh

#include <protocols/stepwise/modeler/rna/util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>

#include <string>
#include <map>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class Combine_Tags_Info{

public:

	Combine_Tags_Info():
		side_one_tag( "" ),
		side_two_tag( "" ),
		//combine_score( 999999999999.99 ) //Feb 12, 2012 This might lead to server-test error at R47200
		combine_score( 999999.9 ) //Feb 12, 2012
	{
	}

	~Combine_Tags_Info(){};

public:

	std::string side_one_tag;
	std::string side_two_tag;
	core::Real combine_score;
};


class Filterer_Count{

public:

	Filterer_Count():
		total_count( 0 ),
		score_cut_count( 0 ),
		chain_closable_geometry_screen( 0 ),
		filter_for_previous_contact( 0 ),
		filter_for_previous_clash( 0 ),
		filter_for_moving_res_contact( 0 )
	{
	}

	~Filterer_Count(){};

public:

	core::Size total_count;
	core::Size score_cut_count;
	core::Size chain_closable_geometry_screen;
	core::Size filter_for_previous_contact;
	core::Size filter_for_previous_clash;
	core::Size filter_for_moving_res_contact;
};


class StepWiseRNA_CombineLongLoopFilterer: public utility::pointer::ReferenceCount {
public:

	//constructor!
	StepWiseRNA_CombineLongLoopFilterer( working_parameters::StepWiseWorkingParametersCOP const & working_parameters, bool const combine_helical_silent_file );

	//destructor -- necessary?
	virtual ~StepWiseRNA_CombineLongLoopFilterer();

	/////////////////////////////////////////////////////////////////////////


	void
	set_silent_files_in( utility::vector1< std::string > const & setting ){ silent_files_in_ = setting; } //Only called if COPY_DOF is true

	void
	set_parin_favorite_output( bool const setting ){ parin_favorite_output_ = setting; }

	void
	filter();

	void
	set_output_filename( std::string const & setting ){ output_filename_ = setting; }

	//void
	//set_score_diff_cut(core::Real const setting){score_diff_cut_=setting;}

	void
	set_filter_for_previous_contact( core::Real const setting ){ filter_for_previous_contact_ = setting; }

	void
	set_filter_for_previous_clash( core::Real const setting ){ filter_for_previous_clash_ = setting; }

	void
	set_undercount_sugar_rotamers( bool const setting ){ undercount_sugar_rotamers_ = setting; }

	void set_max_decoys( core::Size const & setting ){ max_decoys_ = setting; }

private:

	void
	figure_out_appended_and_prepended_res_list();

	void
	figure_out_last_appended_and_last_prepended_res();


	utility::vector1< core::pose::PoseOP >
	convert_silent_file_to_pose_data_list( core::import_pose::pose_stream::SilentFilePoseInputStreamOP & silent_file_stream, core::Size const pose_list_id );


	bool
	previously_builded_res_VDW_filter( core::pose::PoseOP const & side_ONE_pose_data,
		core::pose::PoseOP const & side_TWO_pose_data,
		core::Real const overlap_dist_cutoff,
		core::Size const num_atom_contacts_cutoff );

	bool
	previously_builded_res_contact_filter( core::pose::PoseOP const & side_ONE_pose_data, core::pose::PoseOP const & side_TWO_pose_data );

	bool
	previously_builded_res_clash_filter( core::pose::PoseOP const & side_ONE_pose_data, core::pose::PoseOP const & side_TWO_pose_data );

	bool
	moving_res_contact_filter( core::pose::PoseOP const & side_ONE_pose_data, core::pose::PoseOP const & side_TWO_pose_data );

	void
	align_all_pose( utility::vector1< core::pose::PoseOP > const & side_ONE_pose_data_list,
		utility::vector1< core::pose::PoseOP > const & side_TWO_pose_data_list );


	void
	do_some_filtering();

	bool
	pass_all_filters( core::pose::PoseOP const & side_ONE_pose_data, core::pose::PoseOP const & side_TWO_pose_data );

	void
	setup_silent_file_stream();

	void
	figure_out_NUM_pose_list();

	void
	setup_tag_to_source_map();

	//  bool
	//  score_sort_criterion(Combine_Tags_Info tag_info_1, Combine_Tags_Info tag_info_2);

	void
	sort_Combine_Tags_Info( utility::vector1< Combine_Tags_Info > & combine_tags_info_list );


	std::string
	get_parent_tag( utility::vector1< std::string > const & tag_token ) const;

	bool
	is_virt_sample_sugar_tag( std::string const & tag, utility::vector1< std::string > const & tag_token ) const;

	bool
	is_sibling_sugar_rotamer_pose( std::string const & curr_tag, std::string const & prev_tag, std::map< std::string, std::string > const & tag_to_source_map ) const;


private:

	//  utility::vector1< utility::vector1< Size > > input_res_vectors_;

	Filterer_Count filterer_count_;

	core::chemical::ResidueTypeSetCAP rsd_set_;
	utility::vector1< std::string > silent_files_in_;
	core::import_pose::pose_stream::SilentFilePoseInputStreamOP silent_file_stream_ONE_;
	core::import_pose::pose_stream::SilentFilePoseInputStreamOP silent_file_stream_TWO_;
	working_parameters::StepWiseWorkingParametersCOP const working_parameters_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool verbose_;
	bool parin_favorite_output_;
	bool filter_for_previous_contact_;
	bool filter_for_previous_clash_;
	bool undercount_sugar_rotamers_;
	bool const filter_for_chain_closable_geometry_;
	bool const filter_for_moving_res_contact_;
	bool const moving_res_to_base_contact_only_;

	core::Size total_input_struct_pair_;
	core::Size pass_screen_struct_pair_;
	core::Size input_pose_ONE_last_appended_res_;
	core::Size input_pose_TWO_last_prepended_res_;
	utility::vector1< core::Size > input_pose_ONE_appended_res_list_;
	utility::vector1< core::Size > input_pose_TWO_prepended_res_list_;

	std::map< core::Size, core::Size > full_to_input_res_map_ONE_;
	std::map< core::Size, core::Size > full_to_input_res_map_TWO_;
	std::string output_filename_;
	core::Real best_combine_score_;
	core::Real worst_combine_score_;
	//core::Real score_diff_cut_;
	core::Real const contact_dist_cutoff_;
	core::Real const clash_dist_cutoff_;
	core::Size const num_contact_cutoff_;
	core::Size const num_clash_cutoff_;

	core::Size const max_pose_data_list_size_;
	core::Size side_ONE_NUM_pose_list_;
	core::Size side_TWO_NUM_pose_list_;
	core::Size side_ONE_pose_list_id_;
	core::Size side_TWO_pose_list_id_;
	core::Real moving_res_contact_dist_cutoff_;
	utility::vector1< Combine_Tags_Info > filterered_combine_tag_info_list_;

	Size max_decoys_;
	bool combine_helical_silent_file_;

	std::map< std::string, std::string > tag_to_source_map_ONE_;
	std::map< std::string, std::string > tag_to_source_map_TWO_;


};

} //rna
} //modeler
} //stepwise
} //protocols

#endif
