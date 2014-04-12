// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/NoOpRotamer.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_NoOpRotamer_HH
#define INCLUDED_protocols_rotamer_sampler_NoOpRotamer_HH

#include <protocols/rotamer_sampler/RotamerSized.hh>
#include <protocols/rotamer_sampler/NoOpRotamer.fwd.hh>

namespace protocols {
namespace rotamer_sampler {

	class NoOpRotamer: public rotamer_sampler::RotamerSized {

	public:

		//constructor
		NoOpRotamer();

		//destructor
		~NoOpRotamer();

	public:

		/// @brief Get the total number of rotamers in sampler
		virtual core::Size size() const{ return 1;}

		/// @brief Apply the i-th rotamer to pose
		virtual void apply( core::pose::Pose&, core::Size const ){}

		/// @brief Name of the class
		virtual std::string get_name() const { return "NoOpRotamer"; }

		/// @brief Type of class (see enum in RotamerTypes.hh)
		virtual RotamerType type() const { return NO_OP; }

	};

} //rotamer_sampler
} //protocols

#endif
