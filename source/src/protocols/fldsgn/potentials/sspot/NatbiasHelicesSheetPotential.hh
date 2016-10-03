// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/potentials/sspot/NatbiasHelicesSheetPotential.hh
/// @brief centroid score for helices on sheet
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasHelicesSheetPotential_hh
#define INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasHelicesSheetPotential_hh

// Unit headers
#include <protocols/fldsgn/potentials/sspot/NatbiasHelicesSheetPotential.fwd.hh>

// Package headers
#include <protocols/fldsgn/topology/HelixPairing.fwd.hh>
#include <protocols/fldsgn/topology/HSSTriplet.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

class NatbiasHelicesSheetPotential : public utility::pointer::ReferenceCount {
public: // typedef


	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;
	typedef protocols::fldsgn::topology::HSSTripletOP HSSTripletOP;
	typedef protocols::fldsgn::topology::HSSTripletSetOP HSSTripletSetOP;
	typedef protocols::fldsgn::topology::HelixPairingSet HelixPairingSet;
	typedef protocols::fldsgn::topology::HelixPairingSetOP HelixPairingSetOP;
	typedef protocols::fldsgn::topology::HSSConstIterator HSSConstIterator;


public: // construct/destruct


	/// @brief default constructor
	NatbiasHelicesSheetPotential();

	/// @brief value constructor
	NatbiasHelicesSheetPotential( HSSTripletSetOP const hss3set );

	/// @brief value constructor
	NatbiasHelicesSheetPotential( HSSTripletSetOP const hss3set, HelixPairingSetOP const hpairset );

	/// @brief copy constructor
	NatbiasHelicesSheetPotential( NatbiasHelicesSheetPotential const & src );

	/// @brief default destructor
	virtual ~NatbiasHelicesSheetPotential();


public:


	/// @brief set parameters
	void set_params();


public: // mutator


	/// @brief set HSSTripletSet
	void
	hss_triplet_set( HSSTripletSetOP const hss3set );

	/// @brief set HelixPairingSet
	void
	hpairset( HelixPairingSetOP const hpairset );

	/// @brief set dist parameters for helix-strands interaction
	void
	set_atrdist_params_helix_strands(
		Real const hs_dist_wts,
		Real const hs_dist,
		Real const hs_dist_sigma2 );

	/// @brief set dist parameters for helix-sheet interaction

	void set_hs_atr_dist_wts( Real const r ) { hs_dist_wts_ = r; }
	void set_hs_atr_dist( Real const r ) { hs_dist_ = r; }
	void set_hs_atr_dist_sigma2( Real const r ) { hs_dist_sigma2_ = r; }
	void set_hs_angle_wts( Real const r ) { hs_angle_wts_ = r; }
	void set_hs_angle( Real const r ) { hs_angle_ = r; }
	void set_hs_angle_sigma2( Real const r ) { hs_angle_sigma2_ = r; }
	void set_hsheet_repl_dist( Real const r ) { hsheet_dist_repulsive_ = r; }
	void set_hh_angle_wts( Real const r ) { hh_align_angle_wts_ = r; }
	void set_hh_angle( Real const r ) { hh_align_angle_ = r;}
	void set_hh_angle_sigma2( Real const r ) { hh_align_angle_sigma2_ = r;}

	void
	set_repldist_params_helix_sheet( Real const hsheet_dist_repulsive );

	/// @brief set angle parameters for helix-sheet interaction
	void
	set_angle_params_helix_sheet(
		Real const hs_angle_wts,
		Real const hs_angle,
		Real const hs_angle_sigma2 );

	/// @brief set angle parameters for helices projected onto sheet
	void
	set_angle_params_helices_on_sheet(
		Real const hh_align_angle_wts,
		Real const hh_align_angle,
		Real const hh_align_angle_sigma2 );


public: // accessor


	/// @brief shows parameters for score calculation
	void show_params() const;


public: // scoring


	/// @brief calc score
	void score( SS_Info2_COP const ss_info, Real & hh_score, Real & hs_score ) const;

private:
	/// @brief Gets HSS Triplet containing helix helix_id.
	/// @details  Currently exits with error if more than one HSS triplet contain the helix
	HSSTripletOP
	get_hssop( Size const helix_id ) const;

private: // secondary structure data


	/// @brief HSSTripletSet
	HSSTripletSetOP hss3set_;

	/// @brief HelixPairingSet
	HelixPairingSetOP  hpairset_;

	/// @brief score between helix and sheet
	mutable utility::vector1< Real > hs_scores_;

	/// @brief weights for distant score between helix and sheet
	Real hs_dist_wts_;

	/// @brief maximum distance for the most favorable interaction between helix and sheet
	Real hs_dist_;

	/// @brief sigma for distance score between helix and sheet
	Real hs_dist_sigma2_;

	/// @brief repulsive distance between helix and sheet
	Real hsheet_dist_repulsive_;

	/// @brief weights for angle score between helix and sheet
	Real hs_angle_wts_;

	/// @brief maximum angle for the most favorable interaction between helix and sheet
	Real hs_angle_;

	/// @brief sigma for angle score between helix and sheet
	Real hs_angle_sigma2_;

	/// @brief score of helix pair on sheet
	mutable utility::vector1< Real > hh_scores_;

	/// @brief weight for angle score of helix pair
	Real hh_align_angle_wts_;

	/// @brief maximum angle for the most favorable interaction of helix-pair on sheet
	Real hh_align_angle_;

	/// @brief sigma for angle score of helix pair
	Real hh_align_angle_sigma2_;

};

} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols


#endif
