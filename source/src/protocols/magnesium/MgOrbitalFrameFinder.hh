// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/magnesium/MgOrbitalFrameFinder.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_MgOrbitalFrameFinder_HH
#define INCLUDED_protocols_magnesium_MgOrbitalFrameFinder_HH

#include <protocols/moves/Mover.hh>
#include <protocols/magnesium/MgOrbitalFrameFinder.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <numeric/UniformRotationSampler.fwd.hh>

namespace protocols {
namespace magnesium {

class MgOrbitalFrameFinder: public moves::Mover {

public:
	//constructor
	MgOrbitalFrameFinder();

	//destructor
	~MgOrbitalFrameFinder() override;

public:

	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override{ return "MgOrbitalFrameFinder"; }

private:

	void
	determine_mg_orbital_frame( core::pose::Pose & pose,
		core::Size const i /* mg2+ res number*/ );

	void
	point_orbitals_to_closest_ligands( core::pose::Pose & pose,
		core::Size const i /* mg2+ res number*/,
		utility::vector1< core::id::AtomID > const & ligands );

	core::Real
	get_orbital_frame_score_upon_rotation( utility::vector1< core::Vector > const & r_lig,
		utility::vector1< core::Vector > const & v_lig,
		numeric::xyzMatrix< core::Real > const & R );
	void
	sample_orbital_frame( core::pose::Pose & pose,
		core::Size const i /* mg2+ res number*/,
		utility::vector1< core::id::AtomID > const & ligands );

private:

	//bool legacy_mode_;
	numeric::UniformRotationSamplerCOP urs_;
};

} //magnesium
} //protocols

#endif
