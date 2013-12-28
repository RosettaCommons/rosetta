// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_Clusterer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_RigidBodySampler_HH
#define INCLUDED_protocols_stepwise_RigidBodySampler_HH

#include <protocols/stepwise/enumerate/legacy/RigidBodySampler.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <core/types.hh>

#include <map>

using namespace core;
typedef numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace legacy {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class RigidBodySampler: public utility::pointer::ReferenceCount {
  public:

    //constructor!
		RigidBodySampler( utility::vector1< Size > const & fixed_res,
											utility::vector1< Size > const & moving_res );

    //destructor -- necessary?
    virtual ~RigidBodySampler();

		void
		do_the_sampling( pose::Pose & pose );

		void
		save_silent_struct( pose::Pose & pose, std::string const tag );

		void
		output_results( utility::io::ozstream & out );

		void
		output_histogram( utility::io::ozstream & out );

		void
		output_silent_file( std::string const silent_file, bool const write_score_only = false );


		void set_native_pose( pose::PoseOP native_pose ){ native_pose_ = native_pose; }

		void set_n_sample_alpha_full_range( Size const setting );
		void set_n_sample_cosbeta_full_range( Size const setting );
		void set_n_sample_gamma_full_range( Size const setting );

		void set_alpha_min( Real const setting ){ alpha_min_ = setting;}
		void set_alpha_max( Real const setting ){ alpha_max_ = setting;}
		void set_alpha_increment( Real const setting ){ alpha_increment_ = setting;}

		void set_cosbeta_min( Real const setting ){ cosbeta_min_ = setting;}
		void set_cosbeta_max( Real const setting ){ cosbeta_max_ = setting;}
		void set_cosbeta_increment( Real const setting ){ cosbeta_increment_ = setting;}

		void set_gamma_min( Real const setting ){ gamma_min_ = setting;}
		void set_gamma_max( Real const setting ){ gamma_max_ = setting;}
		void set_gamma_increment( Real const setting ){ gamma_increment_ = setting;}

		void set_translation_sample( Real const box_size, Real const xyz_increment );

		void set_x_min( Real const setting ){ x_min_ = setting;}
		void set_x_max( Real const setting ){ x_max_ = setting;}
		void set_x_increment( Real const setting ){ x_increment_ = setting;}

		void set_y_min( Real const setting ){ y_min_ = setting;}
		void set_y_max( Real const setting ){ y_max_ = setting;}
		void set_y_increment( Real const setting ){ y_increment_ = setting;}

		void set_z_min( Real const setting ){ z_min_ = setting;}
		void set_z_max( Real const setting ){ z_max_ = setting;}
		void set_z_increment( Real const setting ){ z_increment_ = setting;}

		void set_rmsd_cutoff( Real const setting ){ rmsd_cutoff_ = setting;}

		void		force_coplanar();
		void		force_antiparallel();
		void		force_parallel();

		void set_score_cutoff( Real const setting ){ score_cutoff_ = setting;}
		void set_score_function( core::scoring::ScoreFunctionOP setting ){ scorefxn_ = setting;}
		core::scoring::ScoreFunctionOP  score_function();

		void	set_silent_file_data(  core::io::silent::SilentFileDataOP sfd );
		core::io::silent::SilentFileDataOP silent_file_data(){ return sfd_; }

		void		set_contact_cutoff( Real const setting );
		void		set_min_num_contacts( Size const setting );
		void		set_steric_dist_cutoff( Real const setting );
		void		set_min_hbonds( Size const value ){ min_hbonds_ = value;}
		void		set_fa_rep_cutoff( Real const value ){ fa_rep_cutoff_ = value;}
		void	  set_o2prime_trials( bool const setting );
		void	  set_ignore_o2prime_hbonds_in_filter( bool const setting ){ ignore_o2prime_hbonds_in_filter_ = setting;}
		void	  set_assign_WC_edges( bool const setting ){ assign_WC_edges_ = setting;}


		void
		assign_WC_edges_to_base_pair12( pose::Pose & pose, io::silent::SilentStruct & s );

		void
		apply_input_samples( pose::Pose & pose,
												 std::string const rigid_body_sample_file );

		void
		apply_rigid_body_settings( pose::Pose & pose, pose::Pose const & pose_start,
															 Real const alpha,
															 Real const beta,
															 Real const gamma,
															 Real const x,
															 Real const y,
															 Real const z );

		void
		initialize_reference_axes_and_centroid( conformation::Residue const & rsd );


	private:
		void
		initialize_counters();

		void
		search_rotations_and_translations( pose::Pose & pose );

		void
		search_translations( pose::Pose & pose,
												 pose::Pose const & pose_to_translate );


		void
		setup_heavy_atoms( pose::Pose const & pose,
											 utility::vector1< Vector > & pose_atoms,
											 utility::vector1< Size > const & subset_res );

		bool
		check_contact( Vector const & translation,
									 utility::vector1< Vector > const & moving_atoms,
									 utility::vector1< Vector > const & partner_atoms
									 );

		bool
		check_steric_overlap( Vector const & translation,
													utility::vector1< Vector > const & moving_atoms,
													utility::vector1< Vector > const & partner_atoms
													);

		bool
		check_o2prime_needs_optimization( pose::Pose const & pose );

		bool
		check_num_hbonds( pose::Pose & pose );

		bool
		check_fa_rep( pose::Pose & pose );

		void
		save_rigid_body_settings( Real const energy );

		void
		figure_out_reference_energy( pose::Pose & pose );

  private:

		utility::vector1< Size > const fixed_res_;
		utility::vector1< Size > const moving_res_;
		numeric::xyzMatrix< Real > reference_axes_;
		numeric::xyzVector< Real > reference_centroid_;

		core::io::silent::SilentFileDataOP sfd_;

		core::scoring::ScoreFunctionOP scorefxn_, o2prime_pack_scorefxn_;
		bool o2prime_trials_, ignore_o2prime_hbonds_in_filter_, assign_WC_edges_;

		Real alpha_, beta_, gamma_, delx_, dely_, delz_;

		Real alpha_min_, alpha_max_, alpha_increment_;
		Real cosbeta_min_, cosbeta_max_, cosbeta_increment_;
		Real gamma_min_, gamma_max_, gamma_increment_;
		Real x_min_, x_max_, x_increment_;
		Real y_min_, y_max_, y_increment_;
		Real z_min_, z_max_, z_increment_;
		Real score_cutoff_, best_energy_, reference_energy_;

		Size count_total_, count_good_, count_no_contact_,	count_clash_;
		Size min_hbonds_;
		Real fa_rep_cutoff_;

		utility::vector1< utility::vector1< Real > > all_rigid_body_settings_save_;

		Real CONTACT_CUTOFF_squared_, STERIC_DIST_CUTOFF_squared_;
		Size MIN_NUM_CONTACTS_;

		Real rmsd_cutoff_;

		pose::PoseOP native_pose_;

  };

} //legacy
} //enumerate
} //stepwise
} //protocols

#endif
