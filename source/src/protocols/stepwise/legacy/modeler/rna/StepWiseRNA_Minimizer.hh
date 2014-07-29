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
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_Minimizer_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_Minimizer_HH

#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_Minimizer.fwd.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>     //Feb 02, 2012: Need this to pass rosetta_tools/python_cc_reader/test_all_headers_compile.py
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>
#include <protocols/toolbox/AllowInsert.hh>

using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::rna;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

	//	typedef std::map< std::string, core::pose::PoseOP > PoseList;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseRNA_Minimizer: public protocols::moves::Mover {
  public:


    //constructor!
    StepWiseRNA_Minimizer(
													utility::vector1 < core::pose::PoseOP > const & pose_list,
													stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP & working_parameters );

    //destructor -- necessary?
    ~StepWiseRNA_Minimizer();

    /// @brief Apply the minimizer to one pose
    virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

    void
		set_silent_file( std::string const & silent_file );

		void
		set_move_map( core::kinematics::MoveMapOP move_map_ );

		void
		set_scorefxn( core::scoring::ScoreFunctionCOP const & scorefxn );

		core::io::silent::SilentFileDataOP & silent_file_data();

		void
		set_base_centroid_checker( checker::RNA_BaseCentroidCheckerOP & checker );

		void
		set_user_input_VDW_bin_checker( checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker ){ user_input_VDW_bin_checker_ = user_input_VDW_bin_checker; }


		void
		output_pose_wrapper( std::string const & tag,
															core::pose::Pose & pose,
															std::string const & out_silent_file ) const;
		void
		set_working_extra_minimize_res( utility::vector1< core::Size > setting ){ working_extra_minimize_res_ = setting; }

		void
		set_allow_variable_bond_geometry( bool const setting ){ allow_variable_bond_geometry_ = setting; }

		void
		set_vary_bond_geometry_frequency( core::Real const setting ) { vary_bond_geometry_frequency_ = setting; }

		void
		set_allow_insert( toolbox::AllowInsertOP setting ){ allow_insert_ = setting; }

		void
		set_options( protocols::stepwise::modeler::options::StepWiseModelerOptionsOP options );

		utility::vector1< core::pose::PoseOP >
		const & minimized_pose_list() const { return minimized_pose_list_; }

  private:

		utility::vector1 < core::Size >
		get_working_moving_res( Size const & nres ) const;

		core::kinematics::MoveMapOP
		get_default_movemap( core::pose::Pose const & pose );

		bool
		pass_all_pose_screens( core::pose::Pose & pose, std::string const in_tag, core::io::silent::SilentFileData & silent_file_data ) const;

		void
		freeze_sugar_torsions( core::kinematics::MoveMap & mm, Size const nres ) const;

		void
		output_empty_minimizer_silent_file() const; //FEB 09, 2012

		void
		output_pose_wrapper( std::string & tag, char tag_first_char, core::pose::Pose & pose, core::io::silent::SilentFileData & silent_file_data, std::string const out_silent_file ) const;

		bool
		native_edensity_score_checker( core::pose::Pose & pose, core::pose::Pose & native_pose ) const; //Fang's electron density code

		void
		output_minimized_pose_list();

		void
		output_parameters();

	private:

		utility::vector1 < core::pose::PoseOP > const pose_list_;
		utility::vector1 < core::pose::PoseOP > minimized_pose_list_;
		stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters_;

		core::io::silent::SilentFileDataOP sfd_;
		kinematics::MoveMapOP move_map_;
		core::scoring::ScoreFunctionCOP scorefxn_;
		std::string silent_file_;
		bool perform_electron_density_screen_;

		protocols::stepwise::modeler::options::StepWiseModelerOptionsOP options_;

		std::map< core::id::AtomID, core::id::AtomID > pose_to_native_map_;

		utility::vector1< core::Size > fixed_res_;

		checker::RNA_BaseCentroidCheckerOP base_centroid_checker_;
		checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker_;

		core::Real vary_bond_geometry_frequency_;

		bool allow_variable_bond_geometry_;

		utility::vector1< core::Size > working_extra_minimize_res_;

		toolbox::AllowInsertOP allow_insert_;

		core::Real original_geometry_weight_;


  };

} //rna
} //modeler
} //legacy
} //stepwise
} //protocols

#endif
