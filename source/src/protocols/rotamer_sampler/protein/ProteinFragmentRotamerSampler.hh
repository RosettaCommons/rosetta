// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/protein/ProteinFragmentRotamerSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_protein_ProteinFragmentRotamerSampler_HH
#define INCLUDED_protocols_rotamer_sampler_protein_ProteinFragmentRotamerSampler_HH

#include <protocols/rotamer_sampler/RotamerSamplerSized.hh>
#include <protocols/rotamer_sampler/protein/ProteinFragmentRotamerSampler.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace rotamer_sampler {
namespace protein {

	class ProteinFragmentRotamerSampler: public protocols::rotamer_sampler::RotamerSamplerSized {

	public:

		//constructor
		ProteinFragmentRotamerSampler( std::string const frag_file,
														utility::vector1< core::Size > const & slice_res,
														utility::vector1< core::Size > const & moving_residues	 );


		//destructor
		~ProteinFragmentRotamerSampler();

	public:

		/// @brief Get the total number of rotamers in sampler
		virtual core::Size size() const;

		/// @brief Apply the i-th rotamer to pose
		virtual void apply( core::pose::Pose &, core::Size const );

		/// @brief Name of the class
		virtual std::string get_name() const { return "ProteinFragmentRotamerSampler"; }

		/// @brief Type of class (see enum in RotamerSamplerTypes.hh)
		virtual RotamerSamplerType type() const { return PROTEIN_FRAGMENT; }

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
} //rotamer_sampler
} //protocols

#endif
