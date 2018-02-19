// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/PointWaterPotential.hh
/// @brief  Statistical point water potential
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_PointWaterPotential_hh
#define INCLUDED_core_scoring_PointWaterPotential_hh

// Unit Headers
#include <core/scoring/PointWaterPotential.fwd.hh>
#include <core/scoring/DerivVectorPair.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/interpolation/spline/BicubicSpline.hh>

namespace core {
namespace scoring {

struct tableID {
	tableID( chemical::AA aa_in, std::string name1_in, std::string name2_in) :
		name1( name1_in), name2(name2_in), aa(aa_in)
	{ }

	std::string name1;
	std::string name2;
	chemical::AA aa;
};

class PointWaterPotential : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	PointWaterPotential();
	~PointWaterPotential() {}

	Real
	eval_pointwater_score(
		AA const res_aa1,
		core::conformation::Residue const &res,
		core::Vector O_pos
	) const;

	Real
	eval_pointwater_derivs(
		AA const res_aa1,
		core::conformation::Residue const &res,
		core::Vector O_pos,
		utility::vector1< DerivVectorPair > & r_atom_derivs,
		utility::vector1< DerivVectorPair > & O_derivs,
		core::Real wt
	) const;

protected:

	void read_pointwater_tables();
	void read_table_from_stream( utility::io::izstream &, ObjexxFCL::FArray2D< Real > &);
	void setup_interpolation( ObjexxFCL::FArray2D< Real > &, numeric::interpolation::spline::BicubicSpline  &);

	void dampen_longrange_energies( ObjexxFCL::FArray2D< Real > &table );

	// spline interpolated E(d,theta)
	utility::vector1< numeric::interpolation::spline::BicubicSpline > aa_potentials_;

	numeric::interpolation::spline::BicubicSpline lig_potential_don_;
	numeric::interpolation::spline::BicubicSpline lig_potential_acc_sp2_;
	numeric::interpolation::spline::BicubicSpline lig_potential_acc_sp3_;
	//numeric::interpolation::spline::BicubicSpline lig_potential_acc_ring_;  // to add

	utility::vector1< tableID > potential_ids_;
};

}
}

#endif
