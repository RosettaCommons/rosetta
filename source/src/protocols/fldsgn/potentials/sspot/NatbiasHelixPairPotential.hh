// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.cc
/// @brief header file of NatbiasHelixPairPotential
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasHelixPairPotential_hh
#define INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasHelixPairPotential_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/fldsgn/topology/HelixPairing.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

class NatbiasHelixPairPotential : public utility::pointer::ReferenceCount {
public:


	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::SS_Info2_OP SS_Info2_OP;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;
	typedef protocols::fldsgn::topology::HelixPairing HelixPairing;
	typedef protocols::fldsgn::topology::HelixPairings HelixPairings;
	typedef protocols::fldsgn::topology::HelixPairingSetOP HelixPairingSetOP;


public: // construct/destruct


	/// @brief default constructor
	NatbiasHelixPairPotential();

	/// @brief value constructor
	NatbiasHelixPairPotential( HelixPairingSetOP const hpairset );

	/// @brief copy constructor
	NatbiasHelixPairPotential( NatbiasHelixPairPotential const & src );

	/// @brief default destructor
	virtual ~NatbiasHelixPairPotential();


public:


	void set_params();


public: // mutator


	void set_hpairset( HelixPairingSetOP const hpairset );

	/// @brief set parameters for distance score between mid points of helices
	void set_params_distance_pot( Real w, Real d, Real s );

	/// @brief set  parameters for angle score of helix pair
	void set_params_angle_pot( Real w, Real d, Real s );

	/// @brief
	void set_dist_wts( Real const r ) { dist_wts_ = r; }

	/// @brief set distance parmeter
	void set_dist( Real const r ) { mid_dist_ = r; }

	/// @brief
	void set_dist_sigma2( Real const r ) { dist_sigma2_ = r; }

	/// @brief
	void set_angle_wts( Real const r ) { angle_wts_ = r; }

	/// @brief
	void set_angle( Real const r ) { cross_angle_ = r; }

	/// @brief
	void set_angle_sigma2( Real const r ) { angle_sigma2_ = r; }

	/// @brief set bend angle
	void set_bendangle( Real r ) { bend_angle_ = r; }


public: // useful functions


	void show( Pose const & pose ) const;

	/// @brief show parameters
	void show_params() const;


public: // scoring


	/// @brief score secondary structure
	void score( SS_Info2_COP const ss_info, Real & hh_score ) const;


private: // secondary structure data


	/// @brief
	Real bend_angle_;

	/// @brief
	Real dist_wts_;

	/// @brief optimal maximum length between mid points of helices
	Real mid_dist_;

	/// @brief
	Real dist_sigma2_;

	/// @brief
	Real angle_wts_;

	/// @brief optimal maximum angle of helix pair
	Real cross_angle_;

	/// @brief
	Real angle_sigma2_;

	/// @brief
	HelixPairingSetOP hpairset_;

	/// @brief
	mutable utility::vector1< Real > hh_scores_;


};

} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols


#endif
