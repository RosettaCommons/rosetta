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
/// @detailed
///
/// @author Parin Sripakdeevong
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_SWA_PoseSetup_hh
#define INCLUDED_protocols_swa_SWA_PoseSetup_hh

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <string>
#include <map>


namespace protocols {
namespace swa {
namespace rna {

	//	typedef std::map< std::string, core::pose::PoseOP > PoseList;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseRNA_PoseSetup: public protocols::moves::Mover {
  public:

    //constructor!
		StepWiseRNA_PoseSetup(StepWiseRNA_JobParametersOP & job_parameters);

    //destructor -- necessary?
    ~StepWiseRNA_PoseSetup();

		/////////////////////////////////////////////////////////////////////////

	  virtual void apply( core::pose::Pose & pose);    

		virtual std::string get_name() const;

		void
		set_input_tags( utility::vector1< std::string > const & setting){ input_tags_=setting; } //Only called if COPY_DOF is true

		void
		set_silent_files_in( utility::vector1< std::string > const & setting){ silent_files_in_= setting; } //Only called if COPY_DOF is true

		void
		set_parin_favorite_output( bool const setting){ parin_favorite_output_=setting; }

		void
		set_bulge_res( utility::vector1 < core::Size > const & bulge_res ){ bulge_res_ = bulge_res; }

		void
		set_virtual_res( utility::vector1 < core::Size > const & virtual_res_list){ virtual_res_list_ = virtual_res_list; } 	

		void
		set_native_virtual_res( utility::vector1 < core::Size > const & native_virtual_res_list ){ native_virtual_res_list_ = native_virtual_res_list; } //Parin Mar 22, 2010

		void
		set_copy_DOF( bool const setting){ copy_DOF_=setting;} //Parin Mar 29, 2010

		void
		set_verbose( bool const setting){ verbose_=setting;}
		
		void
		setup_native_pose( core::pose::Pose const & pose );

		//void
		//set_sampler_native_rmsd_screen( bool const setting) { sampler_native_rmsd_screen_= setting; }

		void
		set_FARFAR_start_pdb( std::string const setting){ FARFAR_start_pdb_= setting; }

		void
		set_rebuild_bulge_mode( bool const setting){ rebuild_bulge_mode_=setting;}

		void
		set_output_pdb( bool const setting){ output_pdb_=setting;}

  	private:

		void
		Import_pose( Size const & i, core::pose::Pose & import_pose) const; //Only called if COPY_DOF is true

		void
		make_pose( core::pose::Pose & pose ); //Only called if COPY_DOF is true

		void
		read_input_pose_and_copy_dofs( core::pose::Pose & pose ); //Only called if COPY_DOF is true

		void
		apply_cutpoint_variants( core::pose::Pose & pose , core::pose::Pose & pose_without_cutpoints);

		void
		check_close_chain_break( core::pose::Pose const & pose ) const;

		void
		apply_bulge_variants( core::pose::Pose & pose ) const;

		void
		apply_virtual_phosphate_variants( core::pose::Pose & pose ) const;

		void
		add_terminal_res_repulsion( core::pose::Pose & pose ) const;

 		void
 		apply_virtual_res_variant(core::pose::Pose & pose ) const;

		void
		correctly_copy_HO2star_positions( core::pose::Pose & full_pose , utility::vector1<core::pose::Pose> const & start_pose_list);

		core::Real
		get_nearest_dist_to_O2star( core::Size const O2star_seq_num, core::pose::Pose const & input_pose, utility::vector1< core::Size > const input_res_list , utility::vector1< core::Size > const & common_res_list);

		//void 
		//ensure_idealize_bond_length_bond_angle_at_cutpoint( core::pose::Pose & working_pose);

		void
		add_protonated_H1_adenosine_variants( core::pose::Pose & pose ) const;

		void
		verify_protonated_H1_adenosine_variants( core::pose::Pose & pose ) const;

	private:

//		utility::vector1< utility::vector1< Size > > input_res_vectors_;

		core::chemical::ResidueTypeSetCAP rsd_set_;
		utility::vector1< std::string > input_tags_;
		utility::vector1< std::string > silent_files_in_;
		StepWiseRNA_JobParametersOP job_parameters_;
		bool copy_DOF_;
		bool verbose_;

		bool parin_favorite_output_;


		utility::vector1< Size > bulge_res_;
		utility::vector1< Size > virtual_res_list_;
		utility::vector1< Size > native_virtual_res_list_;

		//bool sampler_native_rmsd_screen_;
		std::string FARFAR_start_pdb_;
		bool rebuild_bulge_mode_; 
		bool output_pdb_;

	};

}
} //swa
} // protocols

#endif
