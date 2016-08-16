// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/potentials/sspot/NatbiasStrandPairPotential.hh
/// @brief native biased centroid score for strand pairings
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasStrandPairPotential_hh
#define INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasStrandPairPotential_hh

#include <protocols/fldsgn/potentials/sspot/NatbiasStrandPairPotential.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <protocols/fldsgn/topology/BB_Pos.fwd.hh>
#include <protocols/fldsgn/topology/DimerPairing.fwd.hh>
#include <protocols/fldsgn/topology/StrandPairing.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray4D.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

/// @brief secondary structure scoring cut from classic rosetta structure.h/structure.cc
class NatbiasStrandPairPotential : public utility::pointer::ReferenceCount {
public:


	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef protocols::fldsgn::topology::BB_Pos BB_Pos;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::DimerPairings DimerPairings;
	typedef protocols::fldsgn::topology::StrandPairingSetOP StrandPairingSetOP;
	typedef protocols::fldsgn::topology::StrandPairingSetCOP StrandPairingSetCOP;

	typedef ObjexxFCL::FArray1D< int > FArray1D_int;
	typedef ObjexxFCL::FArray1D< Real > FArray1D_real;
	typedef ObjexxFCL::FArray4D< Real > FArray4D_real;


public: // construct/destruct


	/// @brief default constructor
	NatbiasStrandPairPotential();

	/// @brief default constructor
	NatbiasStrandPairPotential( StrandPairingSetOP const spairset );

	/// @brief default destructor
	virtual ~NatbiasStrandPairPotential();


public: // scoring


	/// @brief score secondary structure
	void score( Pose const & pose, SS_Info2 const & ss_info, Real & ss_score ) const;


public:


	// @brief
	void set_native_spairset( StrandPairingSetOP const spairset );


private: // methods


	/// @brief calculate sum of dot product of the co vectors of strand dimers ss1 and ss2
	/// @brief with the vector connecting the midpoints of the dimer vectors (vdist)
	/// @brief also determine return the sign of the dot products for each dimer
	/// @brief to determine which direction the CO groups point
	void
	pair_dp(
		Size const & ss1,
		Size const & ss2,
		BB_Pos const & bb_pos,
		Real & dp,
		Vector const & mid_vector,
		Size & sign1,
		Size & sign2
	) const;


private:


	// @brief
	Real calc_phithetascore( Real const phi, Real const theta ) const;

	// @brief
	Real calc_dotscore( Real const dpall ) const;

	/// @brief
	Real calc_rsigmascore( Real sig, Real dist, Size const sign1, Size const sign2 ) const;

	/// @brief
	static void rsigma_dot_initializer( FArray4D_real & rsigma_dot );


private: // initialization


	/// @brief load phi/theta bins for use in secondary structure scoring
	void load_phi_theta_bins( String const & ss_filename = "scoring/score_functions/SecondaryStructurePotential/phi.theta.36.SS.resmooth" );

	/// @brief
	void load_dotscore_bins();


private: // secondary structure data


	/// @brief
	Real strand_dist_cutoff_;
	/// @brief
	FArray4D_real phithetascore_;
	/// @brief
	FArray1D_real dotscore_;
	/// @brief
	FArray4D_real rsigma_dot_;

	/// @brief
	StrandPairingSetOP native_spairset_;


};

} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols


#endif
