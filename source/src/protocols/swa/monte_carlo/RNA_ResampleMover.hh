// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/RNA_ResampleMover.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_monte_carlo_RNA_ResampleMover_HH
#define INCLUDED_protocols_swa_monte_carlo_RNA_ResampleMover_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/swa/monte_carlo/RNA_ResampleMover.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace swa {
namespace monte_carlo {

	class RNA_ResampleMover: public utility::pointer::ReferenceCount {

	public:

		//constructor
		RNA_ResampleMover(	core::scoring::ScoreFunctionOP scorefxn );

		//destructor
		~RNA_ResampleMover();

	public:

		bool
		apply( core::pose::Pose & pose, std::string & move_type );


		void set_just_min_after_mutation_frequency( core::Real const & setting ){ just_min_after_mutation_frequency_ = setting; }
		core::Real just_min_after_mutation_frequency() const{ return just_min_after_mutation_frequency_; }

		void set_allow_internal_moves( bool const & setting ){ allow_internal_moves_ = setting; }
		bool allow_internal_moves() const{ return allow_internal_moves_; }

		void set_num_random_samples( core::Size const & setting ){ num_random_samples_ = setting; }
		core::Size num_random_samples() const{ return num_random_samples_; }

		void set_use_phenix_geo( bool const & setting ){ use_phenix_geo_ = setting; }
		bool use_phenix_geo() const{ return use_phenix_geo_; }

		void set_erraser( bool const & setting ){ erraser_ = setting; }
		bool erraser() const{ return erraser_; }

		void set_minimize_single_res( bool const & setting ){ minimize_single_res_ = setting; }
		bool minimize_single_res() const{ return minimize_single_res_; }

	private:

		core::scoring::ScoreFunctionOP scorefxn_;

		core::Real just_min_after_mutation_frequency_;
		bool allow_internal_moves_;
		core::Size num_random_samples_;
		bool use_phenix_geo_;
		bool erraser_;
		bool minimize_single_res_;

	};

} //monte_carlo
} //swa
} //protocols

#endif
