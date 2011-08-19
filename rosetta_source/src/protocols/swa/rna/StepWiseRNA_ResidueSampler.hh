// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software res and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Parin Sripakdeevong
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_ResidueSampler_hh
#define INCLUDED_protocols_swa_rna_StepWiseRNA_ResidueSampler_hh

//#include <numeric/xyzMatrix.hh>
//#include <numeric/xyzVector.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/moves/GreenPacker.fwd.hh>
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
class StepWiseRNA_ResidueSampler: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseRNA_ResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters_ );

	//destructor -- necessary?
	~StepWiseRNA_ResidueSampler();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void
	set_silent_file( std::string const & setting );

	void
	set_output_filename( std::string const & output_filename);

	void
	set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	void
	set_fast( bool const & setting );

	void
	set_verbose( bool const & setting );

	void
	set_o2star_screen( bool const & setting );

	void
	set_allow_bulge_at_chainbreak( bool const & setting );

	utility::vector1< pose_data_struct2 > &
	get_pose_data_list();

	core::io::silent::SilentFileDataOP & silent_file_data();

	void output_pose_data_list( std::string const & silent_file ) const;

	void
	set_base_centroid_screener( StepWiseRNA_BaseCentroidScreenerOP & screener );

private:

	void
	initialize_scorefunctions();

	PuckerState
	Get_residue_pucker_state(core::pose::Pose const & pose, core::Size const seq_num) const;

	void
	Copy_CCD_torsions(core::pose::Pose & pose, core::pose::Pose const & template_pose) const;

	bool
	Chain_break_screening(
				core::pose::Pose & chain_break_screening_pose,
				core::scoring::ScoreFunctionOP const & constraint_scorefxn );

	/////////////////////////////////////////////////////////////////////////////Will create a RNA_rotamer class for these function/////////////////
	void
	get_full_rotamers_exclude_syn_chi( utility::vector1< utility::vector1 <utility::vector1 <core::Real > > >& backbone_rotamers_groups, PuckerState const & pucker1, PuckerState const & pucker2) const;

	void
	get_rotamers(
		utility::vector1< utility::vector1 <utility::vector1 <core::Real > > >&  res_rotamers_groups,
		PuckerState const & pucker1,
		PuckerState const & pucker2 ) const;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	Build_dinucleotide( core::pose::Pose & pose);

	void
	Build_single_nucleotide( core::pose::Pose & pose);

	bool
	Full_atom_van_der_Waals_screening(
		core::pose::Pose & current_pose_screen,
		core::Real const & base_rep_score,
		core::Real const &  base_atr_score,
		core::Real & delta_atr_score,
		core::Real & delta_rep_score,
		bool const apply_fa_atr_cut = true
	);

	void
	initialize_o2star_packer_task( core::pose::Pose const & pose );

	void
	initialize_o2star_green_packer( core::pose::Pose & pose );

	void
	sample_o2star_hydrogen( core::pose::Pose & pose );

	core::Real
	Pose_selection_by_full_score(
		utility::vector1< pose_data_struct2 >& pose_data_list,
		core::pose::Pose & current_pose,
		/*utility::vector1<core::Real> const & current_rotamer, */
		std::string const & tag,
		core::scoring::ScoreFunctionOP const & scorefxn
	) const;

	void
	apply_bulge_variant( core::pose::Pose & pose, bool & bulge_added, core::Real const & delta_atr_score );

	void
	Update_pose_data_list(std::string const & tag, utility::vector1< pose_data_struct2 > & pose_data_list, core::pose::Pose const & current_pose, core::Real const & current_score) const;

	void
	cluster_pose_data_list(utility::vector1< pose_data_struct2 > & pose_data_list) const;


	//Best way to make this robust is to tokenize the tag.
	std::string
	create_tag(std::string const & prestring, Size const & group_rotamer, Size const & subgroup_rotamer, std::string const & old_tag);

	std::string
	create_tag(std::string const & prestring, Size const bulge_rotamer_ID, Size const & group_rotamer, Size const & subgroup_rotamer, std::string const & old_tag);


private:

	StepWiseRNA_JobParametersCOP job_parameters_;

	core::io::silent::SilentFileDataOP sfd_;
	utility::vector1< pose_data_struct2 > pose_data_list_;

	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP atr_rep_screening_scorefxn_;
	core::scoring::ScoreFunctionOP chainbreak_scorefxn_;
	core::scoring::ScoreFunctionOP sampling_scorefxn_;
	core::scoring::ScoreFunctionOP o2star_pack_scorefxn_;

	SillyCountStruct count_data_;
	std::string silent_file_;
	std::string output_filename_;
	core::Real const bin_size_; /*ALWAYS 20!!!*/
	core::Real const rep_cutoff_;
	core::Size num_pose_kept_;
	core::Size const multiplier_;
	core::Real const cluster_rmsd_;
	bool verbose_;
	bool const native_rmsd_screen_;
	core::Real const native_screen_rmsd_cutoff_;

	bool o2star_screen_;
	core::pack::task::PackerTaskOP o2star_pack_task_;
	protocols::moves::GreenPackerOP o2star_green_packer_;
	bool const use_green_packer_;
	bool allow_bulge_at_chainbreak_;
	bool fast_;

	StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;

};

}
} //swa
} // protocols

#endif
