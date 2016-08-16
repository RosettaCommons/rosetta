// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_potentials_sspot_HSPairPotential_HH
#define INCLUDED_protocols_fldsgn_potentials_sspot_HSPairPotential_HH

#include <protocols/fldsgn/potentials/sspot/HSPairPotential.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <protocols/fldsgn/topology/BB_Pos.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray4D.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

/// @brief secondary structure scoring cut from classic rosetta structure.h/structure.cc
class HSPairPotential : public utility::pointer::ReferenceCount {
public:


	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef core::scoring::EnergyGraph EnergyGraph;
	typedef protocols::fldsgn::topology::BB_Pos BB_Pos;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;

	typedef ObjexxFCL::FArray1D< int > FArray1D_int;
	typedef ObjexxFCL::FArray1D< Real > FArray1D_real;
	typedef ObjexxFCL::FArray4D< Real > FArray4D_real;


public: // construct/destruct


	/// @brief default constructor
	HSPairPotential();

	/// @brief default destructor
	virtual ~HSPairPotential();


public: // scoring


	/// @brief score secondary structure
	void score( Pose const & pose,
		SS_Info2 const & ss_info,
		Real & hs_score ) const;


private:


	// @brief
	Real calc_phithetascore( Size const strand_seqsep, Real const phi, Real const theta ) const;

	/// @brief
	void helix_end( Size const & pos1, BB_Pos const & bb_pos, Vector & p1, Vector & p2 ) const;


private: // initialization


	/// @brief load phi/theta bins for use in secondary structure scoring
	void load_phi_theta_bins( String const & ss_filename = "scoring/score_functions/SecondaryStructurePotential/phi.theta.36.HS.resmooth" );


private: // secondary structure data

	/// @brief
	Real dist_cutoff_;

	/// @brief
	FArray4D_real phithetascore_;


};

} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols


#endif
