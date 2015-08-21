// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/potentials/sspot/NatbiasSecondaryStructureEnergy.hh
/// @brief native biased centroid score for secondary structure
/// @author Nobuyasu Koga ( nobuyasau@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasSecondaryStructureEnergy_HH
#define INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasSecondaryStructureEnergy_HH

// Unit Headers
#include <protocols/fldsgn/potentials/sspot/NatbiasSecondaryStructureEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <protocols/fldsgn/topology/HelixPairing.fwd.hh>
#include <protocols/fldsgn/topology/HSSTriplet.fwd.hh>
#include <protocols/fldsgn/topology/StrandPairing.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <protocols/fldsgn/potentials/sspot/NatbiasHelicesSheetPotential.fwd.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.fwd.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasStrandPairPotential.fwd.hh>
#include <utility/vector1.hh>


// Utility headers


namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

/// @brief NatbiasSecondaryStructureEnergy
class NatbiasSecondaryStructureEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:


	typedef core::scoring::methods::WholeStructureEnergy  parent;


public: // typedef


	typedef std::string String;
	typedef core::Real Real;
	typedef core::Distance Distance;
	typedef core::pose::Pose Pose;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::methods::EnergyMethodOP EnergyMethodOP;

	typedef protocols::fldsgn::topology::StrandPairingSetOP StrandPairingSetOP;
	typedef protocols::fldsgn::topology::HelixPairingSetOP HelixPairingSetOP;
	typedef protocols::fldsgn::topology::HSSTripletSetOP HSSTripletSetOP;

	typedef protocols::fldsgn::potentials::sspot::NatbiasHelixPairPotentialOP NatbiasHelixPairPotentialOP;
	typedef protocols::fldsgn::potentials::sspot::NatbiasHelicesSheetPotentialOP NatbiasHelicesSheetPotentialOP;
	typedef protocols::fldsgn::potentials::sspot::NatbiasStrandPairPotentialOP NatbiasStrandPairPotentialOP;


public: // constructor/destructor


	/// @brief default constructor
	NatbiasSecondaryStructureEnergy();

	/// @brief copy constructor
	NatbiasSecondaryStructureEnergy( NatbiasSecondaryStructureEnergy const & src );

	/// @brief clone
	virtual EnergyMethodOP clone() const;


public: // mutator


	/// @brief set native secondary structure
	void native_secstruct( String const & secstruct );

	/// @brief set NatbiasStrandPairPotential
	void set_natbias_spairpot( StrandPairingSetOP const spairset );

	/// @brief set NatbiasHelixPairPotential
	void set_natbias_hpairpot( HelixPairingSetOP const hpairset );

	/// @brief set NatbiasHelicesSheetPotential
	void set_natbias_helices_sheet_pot( HSSTripletSetOP const hss3set );


public: // mutator


	/// @brief set native NatbiasStrandPairPotential
	void set_natbias_spairpot( NatbiasStrandPairPotentialOP const sspot );

	/// @brief set NatbiasHelixPairPotential
	void set_natbias_hpairpot( NatbiasHelixPairPotentialOP const hhpot );

	/// @brief set NatbiasHelicesSheetPotential
	void set_natbias_helices_sheet_pot( NatbiasHelicesSheetPotentialOP const hspot );


public: //


	/// @brief use use original secondary structure potential
	void use_nobias( bool const b ) { use_nobias_ = b; }


public:

	/// @brief scoring
	virtual void setup_for_scoring( Pose & pose, ScoreFunction const & scorefxn ) const;

	/// @brief scoring
	virtual void finalize_total_energy( Pose & pose, ScoreFunction const &, EnergyMap & totals ) const;

	/// @brief The NatbiasSecondaryStructureEnergy class requires that the EnergyGraph
	/// span 12 Angstroms between centroids.  The centroids residues build-in a
	/// 3 Angstrom radius each.
	virtual Distance atomic_interaction_cutoff() const;

	/// @brief
	virtual void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

	/// @brief
	virtual Size version() const;

private:


	/// @brief native secondary structure
	String native_secstruct_;

	/// @brief Is NatbiasStrandPairPotential to be used ?
	bool use_sspot_;

	/// @brief Is NatbiasHelixPairPotential to be used ?
	bool use_hhpot_;

	/// @brief Is NatbiasHelicesSheetPotential to be used ?
	bool use_hspot_;

	/// @brief use original secondary structure potential if this is true
	bool use_nobias_;


private: // mutable variables


	/// @brief pointer of NatbiasStrandPairPotential
	NatbiasStrandPairPotentialOP sspot_;

	/// @brief pointer of NatbiasHelixPairPotential
	NatbiasHelixPairPotentialOP hhpot_;

	/// @brief pointer of NatbiasHeliceesSheetPotential
	NatbiasHelicesSheetPotentialOP hspot_;


};


} // namespace sspot
} // namespace potentials
} // namespace fldsgn
} // namespace protocols

#endif
