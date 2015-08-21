// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_potentials_sspot_SSPairPotential_HH
#define INCLUDED_protocols_fldsgn_potentials_sspot_SSPairPotential_HH

// Unit header
#include <protocols/fldsgn/potentials/sspot/SSPairPotential.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <protocols/fldsgn/topology/BB_Pos.fwd.hh>
#include <protocols/fldsgn/topology/DimerPairing.fwd.hh>

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
class SSPairPotential : public utility::pointer::ReferenceCount {
public:


	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef core::scoring::EnergyGraph EnergyGraph;

	typedef protocols::fldsgn::topology::DimerPairing DimerPairing;
	typedef protocols::fldsgn::topology::BB_Pos BB_Pos;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::DimerPairings DimerPairings;

	typedef ObjexxFCL::FArray1D< int > FArray1D_int;
	typedef ObjexxFCL::FArray1D< Real > FArray1D_real;
	typedef ObjexxFCL::FArray4D< Real > FArray4D_real;


public: // construct/destruct


	/// @brief default constructor
	SSPairPotential();

	/// @brief default destructor
	virtual ~SSPairPotential();


public: // scoring


	/// @brief score secondary structure
	void score( Pose const & pose,
		SS_Info2 const & ss_info,
		DimerPairings & dimer_pairs,
		Real & ss_score ) const;


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
	Real calc_phithetascore( Size const strand_seqsep, Real const phi, Real const theta ) const;

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
	Size dimer_seqsep_cutoff_;
	/// @brief
	Size lowstrand_;
	/// @brief
	FArray4D_real phithetascore_;
	/// @brief
	FArray1D_real dotscore_;
	/// @brief
	FArray4D_real rsigma_dot_;


};

} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols


#endif
