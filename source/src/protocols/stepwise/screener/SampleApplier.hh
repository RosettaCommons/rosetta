// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/SampleApplier.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_SampleApplier_HH
#define INCLUDED_protocols_stepwise_screener_SampleApplier_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/SampleApplier.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	class SampleApplier: public StepWiseScreener {

	public:

		//constructor
		SampleApplier();

		//constructor
		SampleApplier( pose::Pose & pose,
									 bool const apply_residue_alternative_sampler = true );

		//destructor
		~SampleApplier();

	public:

		virtual
		bool
		check_screen(){ return true; }

		virtual
		void
		get_update( rotamer_sampler::RotamerBaseOP sampler );

		virtual
		std::string
		name() const { return "SampleApplier"; }

		virtual
		StepWiseScreenerType
		type() const { return SAMPLE_APPLIER; }

		void
		apply_mover( moves::CompositionMoverOP mover, Size const i, Size const j );

		pose::Pose & pose(){ return pose_; }

		void set_apply_residue_alternative_sampler_( bool const setting ){ apply_residue_alternative_sampler_ = setting; }

	private:

		pose::Pose & pose_;

		core::conformation::ResidueOP moving_rsd_at_origin; // only in use for rigid-body sampling.
		utility::vector1< core::conformation::ResidueOP > moving_rsd_at_origin_list; // only in use for rigid-body sampling.
		bool apply_residue_alternative_sampler_;
	};

} //screener
} //stepwise
} //protocols

#endif
