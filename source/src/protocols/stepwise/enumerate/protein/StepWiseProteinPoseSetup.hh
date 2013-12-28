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
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinPoseSetup_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinPoseSetup_HH

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray1D.hh>

#include <string>
#include <map>

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

	typedef std::map< std::string, core::pose::PoseOP > PoseList;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseProteinPoseSetup: public protocols::moves::Mover {
  public:

    //constructor!
    StepWiseProteinPoseSetup( std::string const & desired_sequence,
											 utility::vector1< std::string > const & start_tags,
											 utility::vector1< std::string > const & silent_files_in );


    //destructor -- necessary?
    ~StepWiseProteinPoseSetup();

    /// @brief Apply the minimizer to one pose
    virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;


    utility::vector1< core::Size > const & moving_residues() const;

		void
		set_n_terminus( bool const setting ){ n_terminus_ = setting; }

		void
		set_c_terminus( bool const setting ){ c_terminus_ = setting; }

		void
		set_sample_junction( bool const setting ){ sample_junction_ = setting; }

		void
		set_add_peptide_plane( bool const & setting );

		//		void
		//		set_rsd_set( core::chemical::ResidueTypeSetCAP rsd_set );

  private:


		void
		match_specific_variants( core::pose::Pose const & pose, core::Size const & pose_res,
														 core::pose::Pose & scratch_pose, core::Size const & scratch_pose_res,
														 utility::vector1< core::chemical::VariantType > const & variant_types ) const;


		void
		match_end_variants( core::pose::Pose const & pose, core::pose::Pose & scratch_pose, core::Size const & start_res, core::Size const & end_res ) const;

		void
		add_end_variants( core::pose::Pose & pose );

		Size
		figure_out_nested_positions(
																std::string const & inside_sequence,
																std::string const & desired_sequence ) const;

		void
		initialize_from_scratch( core::pose::Pose & pose );

		void
		prepend_residues( core::pose::Pose & pose );

		void
		append_residues( core::pose::Pose & pose, core::pose::Pose const & start_pose );

	private:

    std::string const desired_sequence_;
		utility::vector1< std::string > const start_tags_;
		utility::vector1< std::string > const silent_files_in_;
		core::chemical::ResidueTypeSetCAP rsd_set_;

		Size const nres_;

		ObjexxFCL::FArray1D< bool > input_residue_array_;
		ObjexxFCL::FArray1D< bool > junction_residue_array_;
		ObjexxFCL::FArray1D< bool > moving_residue_array_;

		bool sample_junction_;
    utility::vector1< core::Size > moving_residues_;

		bool n_terminus_;
		bool c_terminus_;
		//bool fullatom_;
		bool add_peptide_plane_;
		bool verbose_;

  };

} //protein
} //enumerate
} //stepwise
} //protocols

#endif
