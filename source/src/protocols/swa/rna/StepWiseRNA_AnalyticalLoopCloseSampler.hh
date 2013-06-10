// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software res and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_AnalyticalLoopCloseSampler.hh
/// @brief Alternative SWA sampling using Analytical Loop Closure
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_AnalyticalLoopCloseSampler_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_AnalyticalLoopCloseSampler_HH

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator_Wrapper.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBase_Sampler_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_Bin_Screener.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>


namespace protocols {
namespace swa {
namespace rna {

//	typedef std::map< std::string, core::pose::PoseOP > PoseList;

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWiseRNA_AnalyticalLoopCloseSampler: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseRNA_AnalyticalLoopCloseSampler ( StepWiseRNA_JobParametersCOP & job_parameters_ );

	//destructor -- necessary?
	~StepWiseRNA_AnalyticalLoopCloseSampler();

	/// @brief Apply the minimizer to one pose
	virtual void apply ( core::pose::Pose & pose_to_visualize );

	/// @brief Each derived class must specify its name.  The class name.
	virtual std::string get_name() const {
		return "StepWiseRNA_AnalyticalLoopCloseSampler";
	}

	void
	set_centroid_screen ( bool const setting ) {
		centroid_screen_ = setting;
	}

	void
	set_base_pairing_only_screen ( bool const setting ) {
		base_pairing_only_screen_ = setting;
	}

	void
	set_VDW_atr_rep_screen ( bool const setting ) {
		VDW_atr_rep_screen_ = setting;
	}


	void
	set_silent_file ( std::string const & setting );

	void
	set_output_filename ( std::string const & output_filename );

	void
	set_scorefxn ( core::scoring::ScoreFunctionOP const & scorefxn );

	void
	set_fast ( bool const & setting );

	void
	set_medium_fast ( bool const & setting );

	void
	set_native_rmsd_screen ( bool const & setting );

	void
	set_native_screen_rmsd_cutoff ( core::Real const & setting ) {
		native_screen_rmsd_cutoff_ = setting ;
	}

	void
	set_verbose ( bool const & setting );

	void
	set_o2star_screen ( bool const & setting );

	void
	set_cluster_rmsd ( core::Real const & setting );

	// Undefined, commenting out to fix PyRosetta build  void set_allow_bulge_at_chainbreak ( bool const & setting );

	utility::vector1< pose_data_struct2 > &
	get_pose_data_list();

	core::io::silent::SilentFileDataOP & silent_file_data();

	void output_pose_data_list ( std::string const final_sampler_output_silent_file ) const;

	void set_num_pose_kept ( core::Size const & num_pose_kept );

	void
	set_base_centroid_screener ( StepWiseRNA_BaseCentroidScreenerOP & screener );

	void
	set_distinguish_pucker ( bool const & setting ) {		distinguish_pucker_ = setting ;	}

	void
	set_finer_sampling_at_chain_closure ( bool const & setting ) {		finer_sampling_at_chain_closure_ = setting;	}

	void
	set_PBP_clustering_at_chain_closure ( bool const & setting ) {		PBP_clustering_at_chain_closure_ = setting;	}

	void
	set_user_input_VDW_bin_screener ( StepWiseRNA_VDW_Bin_ScreenerOP const & user_input_VDW_bin_screener );

	void
	set_allow_syn_pyrimidine( bool const & setting ) {		allow_syn_pyrimidine_ =setting;	}

	void
	set_extra_syn_chi_rotamer ( bool const & setting ) {		extra_syn_chi_rotamer_ = setting;	}

	void
	set_extra_anti_chi_rotamer ( bool const & setting ) {		extra_anti_chi_rotamer_ = setting;	}

	void
	set_use_phenix_geo ( bool const & setting ) {		use_phenix_geo_ = setting;	}

	void
	set_choose_random ( bool const & setting ) {		choose_random_ = setting;	}

	void
	set_force_centroid_interaction ( bool const & setting ) {		force_centroid_interaction_ = setting;	}

private:

	void
	initialize_scorefunctions();

	void
	standard_sampling ( core::pose::Pose & pose, utility::vector1< pose_data_struct2 > & pose_data_list );

	void
	copy_suite_torsion ( core::pose::Pose & pose, core::pose::Pose const & ref_pose, core::Size const & suite_num );

	void
	get_base_atr_rep_score ( core::pose::Pose const & pose, core::Real & base_atr_score, core::Real & base_rep_score );

	bool
	Full_atom_van_der_Waals_screening_REPLICATE (
	  core::pose::Pose & current_pose_screen,
	  core::Real const & base_rep_score,
	  core::Real const &  base_atr_score,
	  core::Real & delta_rep_score,
	  core::Real & delta_atr_score,
	  core::Size const & gap_size,
	  bool const & Is_internal );


	bool
	Full_atom_van_der_Waals_screening (
	  core::pose::Pose & current_pose_screen,
	  core::Real const & base_rep_score,
	  core::Real const &  base_atr_score,
	  core::Real & delta_rep_score,
	  core::Real & delta_atr_score,
	  core::Size const & gap_size,
	  bool const & Is_internal );

	void
	initialize_o2star_packer_task ( core::pose::Pose const & pose );

	void
	sample_o2star_hydrogen ( core::pose::Pose & pose , core::pose::Pose & pose_with_original_HO2star_torsion );

	core::Real
	Pose_selection_by_full_score ( utility::vector1< pose_data_struct2 >& pose_data_list, core::pose::Pose & current_pose, std::string const & tag );


	void
	Update_pose_data_list ( std::string const & tag, utility::vector1< pose_data_struct2 > & pose_data_list, core::pose::Pose const & current_pose, core::Real const & current_score ) const;

	void
	cluster_pose_data_list ( utility::vector1< pose_data_struct2 > & pose_data_list );

	std::string
	create_tag ( std::string const prestring, StepWiseRNA_RotamerGenerator_WrapperOP const & rotamer_generator ) const;

	std::string //silly function to convert to real to string
	create_torsion_value_string ( core::Real const & torsion_value ) const;

	std::string
	create_rotamer_string ( core::pose::Pose const & pose ) const;

	/////////////////////////////////////function related to sampling/setup virtual sugar /////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


private:

	StepWiseRNA_JobParametersCOP job_parameters_; //need to use the full_to_sub map...should convert to const style.. Parin Feb 28, 2010

	core::io::silent::SilentFileDataOP sfd_;
	utility::vector1< pose_data_struct2 > pose_data_list_;

	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP atr_rep_screening_scorefxn_;
	core::scoring::ScoreFunctionOP sampling_scorefxn_;
	core::scoring::ScoreFunctionOP o2star_pack_scorefxn_;
	core::scoring::ScoreFunctionOP chainbreak_scorefxn_;

	utility::vector1 < core::Size > working_rmsd_res_;

	SillyCountStruct count_data_;
	std::string silent_file_;
	std::string output_filename_;
	core::Real const bin_size_; /*ALWAYS 20!!!*/
	core::Real const rep_cutoff_;
	core::Size num_pose_kept_;
	core::Size const multiplier_;
	core::Real cluster_rmsd_;
	bool verbose_;
	bool native_rmsd_screen_;
	core::Real native_screen_rmsd_cutoff_;

	bool o2star_screen_;
	core::pack::task::PackerTaskOP o2star_pack_task_;
	bool fast_;
	bool medium_fast_;

	StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;

	bool centroid_screen_;
	bool base_pairing_only_screen_;
	bool VDW_atr_rep_screen_;
	bool distinguish_pucker_;
	core::Real current_score_cutoff_;
	bool finer_sampling_at_chain_closure_;
	bool PBP_clustering_at_chain_closure_;
	bool extra_anti_chi_rotamer_;
	bool extra_syn_chi_rotamer_;
	bool allow_syn_pyrimidine_;
	bool use_phenix_geo_;
	bool choose_random_;
	bool force_centroid_interaction_;

	StepWiseRNA_VDW_Bin_ScreenerOP user_input_VDW_bin_screener_;


};

}
} //swa
} // protocols

#endif
