// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_Minimizer_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_Minimizer_HH

//#include <numeric/xyzMatrix.hh>
//#include <numeric/xyzVector.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_BinScreener.hh>     //Feb 02, 2012: Need this to pass rosetta_tools/python_cc_reader/test_all_headers_compile.py
#include <protocols/swa/rna/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <core/conformation/Residue.hh>
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
  class StepWiseRNA_Minimizer: public protocols::moves::Mover {
  public:


    //constructor!
    StepWiseRNA_Minimizer(
            utility::vector1 < pose_data_struct2 > const & pose_data_list,
						StepWiseRNA_JobParametersCOP & job_parameters );

    //destructor -- necessary?
    ~StepWiseRNA_Minimizer();

    /// @brief Apply the minimizer to one pose
    virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

    void
		set_silent_file( std::string const & silent_file );

		void
		set_move_map_list( utility::vector1 < core::kinematics::MoveMap > const & move_map_list );

		void
		set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

		core::io::silent::SilentFileDataOP & silent_file_data();

		void
		set_base_centroid_screener( StepWiseRNA_BaseCentroidScreenerOP & screener );

		void
		set_native_edensity_score_cutoff( core::Real const & setting ); //Fang's electron density code

		void
		set_verbose( bool const setting ){ verbose_ = setting; }

		void
		set_centroid_screen( bool const setting ){ centroid_screen_ = setting; } //For testing purposes.

		void
		set_perform_o2star_pack( bool const setting ){ perform_o2star_pack_ = setting; } //For testing purposes.

		void
		set_output_before_o2star_pack( bool const setting ){ output_before_o2star_pack_ = setting; } //For testing purposes.

		void
		set_rm_virt_phosphate( bool const setting ){ rm_virt_phosphate_ = setting; }

		void
		set_perform_minimize( bool const setting ){ perform_minimize_ = setting; }

		void
		set_num_pose_minimize( Size const setting ){ num_pose_minimize_ = setting; }

		void
		set_minimize_and_score_sugar( bool const setting ){ minimize_and_score_sugar_ = setting; }

		void
		set_native_rmsd_screen( bool const setting ){ native_screen_ = setting; }

		void
		set_native_screen_rmsd_cutoff( core::Real const setting ){ native_screen_rmsd_cutoff_ = setting; }

		void
		set_user_input_VDW_bin_screener( StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener ){ user_input_VDW_bin_screener_ = user_input_VDW_bin_screener; }

		void
		set_rename_tag( bool const setting ){ rename_tag_ = setting; }

		void
		set_output_minimized_pose_data_list( bool const setting ){ output_minimized_pose_data_list_ = setting; }

  private:

		utility::vector1 < core::kinematics::MoveMap >
		Get_default_movemap( core::pose::Pose const & pose ) const;

		void
		Figure_out_moving_residues( core::kinematics::MoveMap & mm, core::pose::Pose const & pose ) const;

		bool
		pass_all_pose_screens( core::pose::Pose & pose, std::string const in_tag, core::io::silent::SilentFileData & silent_file_data ) const;

		void
		Freeze_sugar_torsions( core::kinematics::MoveMap & mm, Size const nres ) const;

		void
		output_empty_minimizer_silent_file() const; //FEB 09, 2012

		void
		output_pose_data_wrapper( std::string & tag, char tag_first_char, core::pose::Pose & pose, core::io::silent::SilentFileData & silent_file_data, std::string const out_silent_file ) const;

		core::Size
		figure_out_actual_five_prime_chain_break_res( core::pose::Pose const & pose ) const;

		bool
		native_edensity_score_screener( core::pose::Pose & pose, core::pose::Pose & native_pose ) const; //Fang's electron density code

		void
		output_minimized_pose_data_list();

	private:

		utility::vector1 < pose_data_struct2 > const pose_data_list_;
		utility::vector1 < pose_data_struct2 > minimized_pose_data_list_;
		StepWiseRNA_JobParametersCOP job_parameters_;

		core::io::silent::SilentFileDataOP sfd_;
		utility::vector1 < core::kinematics::MoveMap > move_map_list_;
		core::scoring::ScoreFunctionOP scorefxn_;
		std::string silent_file_;
		bool verbose_;
		bool screen_verbose_;
		bool native_screen_;
		core::Real native_screen_rmsd_cutoff_;

		bool perform_electron_density_screen_, rm_virt_phosphate_; //Fang's electron density code
		core::Real native_edensity_score_cutoff_; //Fang's electron density code

		bool centroid_screen_; //for testing purposes

		bool perform_o2star_pack_; //Jan 19, 2012
		bool output_before_o2star_pack_; //for testing purposes

		bool perform_minimize_; //Parin Mar 12, 2012

		core::Size num_pose_minimize_;
		bool minimize_and_score_sugar_;
		bool rename_tag_;
		bool output_minimized_pose_data_list_;


		std::map< core::id::AtomID, core::id::AtomID > pose_to_native_map_;

		utility::vector1< core::Size > fixed_res_;

		StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;

		StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener_;


  };

}
} //swa
} // protocols

#endif
