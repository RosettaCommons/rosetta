// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/AddCavitiesMover.hh
///
/// @brief
/// @author


#ifndef INCLUDED_protocols_simple_moves_AddCavitiesMover_hh
#define INCLUDED_protocols_simple_moves_AddCavitiesMover_hh

#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/CavityBall.hh>
#include <core/scoring/packstat/compute_sasa.hh>

#include <protocols/moves/Mover.hh>


#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {


class AddCavitiesMover : public protocols::moves::Mover {
public:

	AddCavitiesMover(
		core::Size max_cav  = 100,
		core::Real min_size = 1.2,
		core::Size min_nb   = 150,
		core::Real min_sep  = 3.0
	);

	void
	clear_suckers(
		core::pose::Pose & pose
	);

	void
	add_suckers(
		core::pose::Pose & pose
	);

	virtual
	void
	apply(
		core::pose::Pose & pose
	);
	virtual std::string get_name() const;

	core::scoring::packstat::CavBalls
	get_cavities(
		core::pose::Pose const & pose,
		core::Real nbdis,
		int nbcount,
		core::Real minsep
	);

protected:

	core::id::AtomID
	get_closest_heavy_atom(
		core::pose::Pose & pose,
		numeric::xyzVector<core::Real> xyz
	);

	core::conformation::ResidueOP
	get_suck_res();

	core::Size max_cav_;
	core::Real min_size_;
	core::Size min_nb_;
	core::Real min_sep_;
	core::scoring::packstat::SasaOptions opts;

};


} // end namespace simple_moves
} // end namespace protocols

#endif
