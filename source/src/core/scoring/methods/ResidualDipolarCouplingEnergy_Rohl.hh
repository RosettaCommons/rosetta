// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ResidualDipolarCouplingEnergy.hh
/// @brief  RDC energy - comparing experimental RDC values to calculated values
/// @author Srivatsan Raman


#ifndef INCLUDED_core_scoring_methods_ResidualDipolarCouplingEnergy_Rohl_hh
#define INCLUDED_core_scoring_methods_ResidualDipolarCouplingEnergy_Rohl_hh

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
// AUTO-REMOVED #include <core/scoring/ResidualDipolarCoupling_Rohl.hh>

#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreType.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

//Objexx headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <ObjexxFCL/char.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/Fmath.hh>

#include <core/scoring/ResidualDipolarCoupling_Rohl.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>



// Utility headers


namespace core {
namespace scoring {
namespace methods {

///
class ResidualDipolarCouplingEnergy_Rohl : public WholeStructureEnergy  {
public:
	typedef WholeStructureEnergy  parent;

public:

	ResidualDipolarCouplingEnergy_Rohl();

	//clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

 private:

	ResidualDipolarCoupling_Rohl const & rdc_from_pose(
		pose::Pose & pose
	) const;

	Real eval_dipolar(
		pose::Pose & pose
	) const;

	void
	assemble_datamatrix(
		pose::Pose const & pose,
		utility::vector1< core::scoring::RDC_Rohl > const & All_RDC_lines,
		ObjexxFCL::FArray2D< Real > & A,
		ObjexxFCL::FArray1D< Real > & b,
		ObjexxFCL::FArray1D< Real > & weights //for alignment tensor calculation
	) const;


	void calc_ordermatrix(
		Size const & nrow,
		Size const & ORDERSIZE,
		ObjexxFCL::FArray2D< Real > & A,
		ObjexxFCL::FArray1D< Real > & b,
		ObjexxFCL::FArray1D< Real > & x,
		ObjexxFCL::FArray1D< Real > & weights, //for alignment tensor calculation
		bool & reject
	) const;

	void svdcmp(
		ObjexxFCL::FArray2D< Real > & a,
		Size const & m,
		Size const & n,
		ObjexxFCL::FArray1D< Real > & w,
		ObjexxFCL::FArray2D< Real > & v
	) const;

	Real pythag(
		Real const & a,
		Real const & b
	) const;

	void svbksb(
		ObjexxFCL::FArray2D< Real > const & u,
		ObjexxFCL::FArray1D< Real > const & w,
		ObjexxFCL::FArray2D< Real > const & v,
		Size const & m,
		Size const & n,
		ObjexxFCL::FArray1D< Real > const & b,
		ObjexxFCL::FArray1D< Real > & x
	) const;

	void calc_orderparam(
		ObjexxFCL::FArray1D< Real > x,
		ObjexxFCL::FArray2D< Real > vec,
		Real & Azz,
		Real & eta
	) const;

	Real calc_dipscore(
	ObjexxFCL::FArray2D< Real > const & A,
	ObjexxFCL::FArray1D< Real > const & x,
	ObjexxFCL::FArray1D< Real > const & b,
	utility::vector1< core::scoring::RDC_Rohl > const & All_RDC_lines,
	Size const & ORDERSIZE,
	Real const & Azz
	) const;
virtual
core::Size version() const;

};

} //methods
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
