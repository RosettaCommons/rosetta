// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_PoseSetup.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_StepWiseProteinPoseSetup_HH
#define INCLUDED_protocols_stepwise_StepWiseProteinPoseSetup_HH

#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinPoseSetup.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

#include <string>
#include <map>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::protein;


namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace protein {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWiseProteinPoseSetup: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseProteinPoseSetup( utility::vector1< core::Size > const & moving_res_list,
		std::string const & desired_sequence,
		utility::vector1< InputStreamWithResidueInfoOP > & input_streams_with_residue_info,
		utility::vector1< core::Size > const & cutpoint_open,
		utility::vector1< core::Size > const & cutpoint_closed );

	//destructor -- necessary?
	~StepWiseProteinPoseSetup();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;


	stepwise::modeler::working_parameters::StepWiseWorkingParametersOP & working_parameters();

	void
	set_fixed_res(utility::vector1 < core::Size > const & fixed_res );

	void
	set_extra_minimize_res( utility::vector1 < core::Size > const & setting ){  extra_minimize_res_ = setting; }

	void
	set_jump_res(utility::vector1 < core::Size > const & jump_res );

	void
	set_virtual_res(utility::vector1 < core::Size > const & set_virtual_res_list);

	void
	set_terminal_res(utility::vector1 < core::Size > const & terminal_res );

	void
	set_superimpose_res(utility::vector1 < core::Size > const & superimpose_res );

	void
	set_calc_rms_res(utility::vector1 < core::Size > const & calc_rms_res );

	void
	set_bulge_res(utility::vector1 < core::Size > const & bulge_res );

	void
	set_bridge_res(utility::vector1 < core::Size > const & bridge_res );

	void
	set_parin_favorite_output( bool const & setting){ parin_favorite_output_=setting ; }

	void
	set_rsd_set( core::chemical::ResidueTypeSetCAP & rsd_set );

	void
	set_cst_file( std::string const & cst_file );

	void
	set_disulfide_file( std::string const & disulfide_file );

	void
	set_align_file( std::string const & align_file );

	void
	set_add_peptide_plane_variants( bool const & setting );

	void
	set_remove_nterminus_variant( bool const & setting ){ remove_nterminus_variant_ = setting; }

	void
	set_remove_cterminus_variant( bool const & setting ){ remove_cterminus_variant_ = setting; }

	void
	align_pose( core::pose::Pose & pose ) const;

	bool
	ready_to_align() const;

	void
	set_dump( bool const dump );

	void
	set_add_virt_res( bool const setting );

	void
	set_secstruct( std::string const & secstruct );


private:

	void
	Import_pose( core::pose::Pose & import_pose, InputStreamWithResidueInfoOP & stream );


	void
	check_moving_res_in_chain( Size const & start_chain, Size const & end_chain,
		Size const & num_chains, Size & which_chain_has_moving_res  );

	void
	figure_out_working_sequence_and_mapping();

	void
	setup_secstruct( core::pose::Pose & pose ) const;

	void
	figure_out_jump_partners();

	Size
	which_chain( Size const & i );

	bool
	already_connected( std::pair< Size, Size > const & potential_chain_partner,
		utility::vector1< std::pair< Size, Size > > const & chain_partners ) const;

	bool
	already_connected( Size const start_chain,
		Size const stop_chain,
		utility::vector1< std::pair< Size, Size > > const & chain_partners,
		utility::vector1< bool > already_checked ) const;


	void
	figure_out_cuts();

	void
	make_pose( core::pose::Pose & pose );

	void
	make_full_pose( core::pose::Pose & pose );

	void
	setup_constraints( core::pose::Pose & pose );

	void
	setup_disulfides( core::pose::Pose & pose );

	void
	initialize_pose_from_streams( core::pose::Pose & pose );

	void
	figure_out_Prepend_Internal( core::pose::Pose const & pose );

	void
	figure_out_partition_definition( core::pose::Pose const & pose );

	void
	figure_out_gap_size_and_first_chain_break_res();

	void
	reroot_fold_tree( core::pose::Pose & pose );

	void
	apply_cutpoint_variants( core::pose::Pose & pose ) const;

	void
	check_close_chain_break( core::pose::Pose const & pose ) const;

	void
	apply_bulge_variants( core::pose::Pose & pose ) const;

	void
	apply_terminus_variants_at_protein_rna_boundaries( core::pose::Pose & pose ) const;

	//  void
	//  slice_native();

	//Replacement for slice_native(); Parin Jan 29, 2010
	void
	setup_working_native_pose();

	void
	get_working_pose( core::pose::Pose const & pose, core::pose::Pose & working_pose );

	void
	align_poses( core::pose::Pose & pose );

	utility::vector1< Size > const
	apply_full_to_sub_mapping( utility::vector1< Size > & res_vector) const;

	void
	apply_virtual_phosphate_variants( core::pose::Pose & pose ) const;

	void
	apply_peptide_plane_variants_OLD( core::pose::Pose & pose ) const;

	void
	apply_peptide_plane_variants( core::pose::Pose & pose ) const;

	void
	add_terminal_res_repulsion( core::pose::Pose & pose ) const;

	void
	apply_virtual_res_variant(core::pose::Pose & pose ) const;

	void
	initialize_phi_psi_offsets( core::pose::Pose const & pose );

	void
	save_phi_psi_offsets(
		core::pose::Pose const & start_pose,
		utility::vector1< core::Size > const & input_res,
		utility::vector1< core::Size > const & slice_res );

	void
	fix_phi_psi_offsets( core::pose::Pose & pose ) const;

	void
	copy_rna_chi( core::pose::Pose & pose,
		core::pose::Pose const & import_pose,
		utility::vector1< core::Size > const & input_res,
		utility::vector1< core::Size > const & slice_res );


	void
	check_superimpose_res( core::pose::Pose const & pose );

	bool
	is_working_cutpoint_closed( Size const res, std::map< Size, Size > & full_to_sub ) const;

	std::string
	get_stepwise_jump_atom( core::conformation::Residue const & rsd );

	void
	add_aa_virt_rsd_as_root( core::pose::Pose & pose);

	void
	setup_full_model_info( core::pose::Pose & pose ) const;

private:

	utility::vector1< core::Size > const & moving_res_list_;
	std::string const desired_sequence_;
	core::chemical::ResidueTypeSetCAP rsd_set_;

	utility::vector1< InputStreamWithResidueInfoOP > input_streams_with_residue_info_;

	utility::vector1< Size > const cutpoint_open_;
	utility::vector1< Size > const cutpoint_closed_;
	ObjexxFCL::FArray1D< bool > is_cutpoint_;
	std::string secstruct_;
	utility::vector1< Size > fixed_res_;
	utility::vector1< Size > virtual_res_list_;
	utility::vector1< Size > terminal_res_;
	utility::vector1< Size > superimpose_res_;
	utility::vector1< Size > calc_rms_res_;
	utility::vector1< Size > bulge_res_;
	utility::vector1< Size > jump_res_;
	utility::vector1< Size > bridge_res_;

	stepwise::modeler::working_parameters::StepWiseWorkingParametersOP working_parameters_;

	ObjexxFCL::FArray1D< core::Real > phi_offsets_;
	ObjexxFCL::FArray1D< core::Real > psi_offsets_;

	utility::vector1< std::pair< core::Size, core::Size > > jump_partners_;
	utility::vector1< core::Size > cuts_;
	utility::vector1< core::Size > extra_minimize_res_;


	bool virtualize_5prime_phosphates_;
	bool add_peptide_plane_variants_;
	bool remove_nterminus_variant_;
	bool remove_cterminus_variant_;
	bool parin_favorite_output_;
	bool add_virt_res_;

	std::string cst_file_;
	std::string disulfide_file_;
	core::scoring::constraints::ConstraintSetOP cst_set_;
	core::pose::PoseOP working_native_pose;

	std::string align_file_;
	bool ready_to_align_;
	core::pose::PoseOP working_align_pose_;
	core::id::AtomID_Map< core::id::AtomID > alignment_atom_id_map_;
	bool dump_;

};


} //protein
} //modeler
} //legacy
} //stepwise
} //protocols

#endif
