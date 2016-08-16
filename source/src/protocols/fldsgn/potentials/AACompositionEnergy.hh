// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief amino acid composition energy
/// @author Nobuyasu Koga


#ifndef INCLUDED_protocols_fldsgn_potentials_AACompositionEnergy_HH
#define INCLUDED_protocols_fldsgn_potentials_AACompositionEnergy_HH

// Unit Headers
#include <protocols/fldsgn/potentials/AACompositionEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility headers
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {

/// @brief AACompositionEnergy
class AACompositionEnergy : public core::scoring::methods::ContextDependentOneBodyEnergy {
public:

	typedef core::scoring::methods::ContextDependentOneBodyEnergy  parent;

public: // typedef


	//typedef std::string String;
	typedef core::chemical::AA AA;
	typedef core::conformation::Residue Residue;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreTypes ScoreTypes;
	typedef core::scoring::methods::EnergyMethodOP EnergyMethodOP;


public: // constructor/destructor


	/// @brief default constructor
	AACompositionEnergy();

	/// @brief value constructor
	AACompositionEnergy( std::map< AA, std::pair< Real, Real >  > const & comp_constraint_aas );

	/// @brief copy constructor
	AACompositionEnergy( AACompositionEnergy const & src );

	/// @brief destructor
	virtual ~AACompositionEnergy();

	/// @brief clone
	virtual EnergyMethodOP clone() const;


private:

	void initialize();

public: // mutator


	void set_comp_constraint_aa( std::map< AA, std::pair< Real, Real > > const & comp_constraint_aas );


public:

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextDependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	residue_energy(
		Residue const & rsd,
		Pose const & pose,
		EnergyMap & emap
	) const;


	/// @brief scoring
	// virtual void setup_for_scoring( Pose &, ScoreFunction const & ) const;

	/// @brief scoring
	// virtual void finalize_total_energy( Pose & pose, ScoreFunction const &, EnergyMap & totals ) const;


	/// @brief DunbrackEnergy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	virtual
	Size version() const { return 1; }


private:

	/// @brief
	std::map< AA, std::pair< Real, Real > > comp_constraint_aas_;


};

} // namespace potentials
} // namespace fldsgn
} // namespace protocols

#endif // INCLUDED_protocols_fldsgn_AACompositionEnergy_HH
