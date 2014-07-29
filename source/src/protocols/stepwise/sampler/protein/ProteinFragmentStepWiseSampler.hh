// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/protein/ProteinFragmentStepWiseSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_protein_ProteinFragmentStepWiseSampler_HH
#define INCLUDED_protocols_sampler_protein_ProteinFragmentStepWiseSampler_HH

#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <protocols/stepwise/sampler/protein/ProteinFragmentStepWiseSampler.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace protein {

	class ProteinFragmentStepWiseSampler: public protocols::stepwise::sampler::StepWiseSamplerSized {

	public:

		//constructor
		ProteinFragmentStepWiseSampler( std::string const frag_file,
														utility::vector1< core::Size > const & slice_res,
														utility::vector1< core::Size > const & moving_residues	 );


		//destructor
		~ProteinFragmentStepWiseSampler();

	public:

		/// @brief Get the total number of rotamers in sampler
		virtual core::Size size() const;

		/// @brief Apply the i-th rotamer to pose
		virtual void apply( core::pose::Pose &, core::Size const );

		/// @brief Name of the class
		virtual std::string get_name() const { return "ProteinFragmentStepWiseSampler"; }

		/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
		virtual StepWiseSamplerType type() const { return PROTEIN_FRAGMENT; }

	private:

		void
		initialize( std::string const frag_file,
								utility::vector1< core::Size > const & slice_res,
								utility::vector1< core::Size > const & moving_residues	 );

	private:

		core::Size insert_pos_;
		core::fragment::FrameOP frame_;

	};

} //protein
} //sampler
} //stepwise
} //protocols

#endif
