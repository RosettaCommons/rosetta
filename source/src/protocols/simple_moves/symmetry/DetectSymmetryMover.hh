// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/symmetry/DetectSymmetryMover.hh
/// @brief Automatical detection and setup of the symmetry machinery from an asymetric pose made of symmetric chains. Only works with cyclic simmetries.
/// @author Javier Castellanos ( javiercv@uw.edu )

#ifndef _INCLUDED_protocols_simple_moves_symmetry_DetectSymmetryMover_hh_
#define _INCLUDED_protocols_simple_moves_symmetry_DetectSymmetryMover_hh_

#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyz.functions.hh>

namespace protocols {
namespace simple_moves {
namespace symmetry {

class DetectSymmetry : public protocols::moves::Mover {
public:
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef numeric::xyzMatrix< core::Real > xyzMatrix;
	typedef numeric::xyzVector< core::Real > xyzVector;

public:

	DetectSymmetry();
	DetectSymmetry(core::Real subunit_tolerance, core::Real plane_tolerance);

	virtual void apply(Pose & pose);

	virtual std::string get_name() const {return "DetectSymmetry";}

	moves::MoverOP clone() const { return( protocols::moves::MoverOP( new DetectSymmetry( *this ) ) ); }

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, Pose const & );

private:
	inline core::Real angle_with_x_axis_proj_y( xyzVector const & v) const { return numeric::dihedral_degrees(xyzVector(v[0],1,v[2]), xyzVector(0,1,0), xyzVector(0,0,0), xyzVector(1,0,0)); }
	inline core::Real angle_with_y_axis_proj_x( xyzVector const & v) const { return numeric::dihedral_degrees(xyzVector(1,v[1],v[2]), xyzVector(1,0,0), xyzVector(0,0,0), xyzVector(0,1,0)); }

	inline core::Real angle_with_y_axis_proj_z( xyzVector const & v) const { return numeric::dihedral_degrees( xyzVector(v[0],v[1],1), xyzVector(0,0,1), xyzVector(0,0,0), xyzVector(0,1,0)); }

	inline core::Real angle_with_z_axis_proj_y( xyzVector const & v) const { return numeric::dihedral_degrees( xyzVector(v[0],1,v[2]), xyzVector(0,1,0), xyzVector(0,0,0), xyzVector(0,0,1)); }


private:
	core::Real subunit_tolerance_;
	core::Real plane_tolerance_;

};

} // symmetry
} // simple_moves
} // protocols


#endif
