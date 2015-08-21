// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoringManager.hh
/// @brief  Scoring manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_SecondaryStructurePotential_hh
#define INCLUDED_core_scoring_SecondaryStructurePotential_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/SecondaryStructureWeights.hh>

#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray4D.hh>

// C++ headers
#include <string>

#include <core/scoring/SS_Info.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>


namespace core {
namespace scoring {


/// @brief secondary structure scoring cut from classic rosetta structure.h/structure.cc
class SecondaryStructurePotential : public utility::pointer::ReferenceCount {

public:

	typedef ObjexxFCL::FArray1D< Real > FArray1D_real;
	typedef ObjexxFCL::FArray2D< Real > FArray2D_real;
	typedef ObjexxFCL::FArray4D< Real > FArray4D_real;
	typedef ObjexxFCL::FArray1< Real > FArray1_real;
	typedef ObjexxFCL::FArray1A< Real > FArray1A_real;
	typedef ObjexxFCL::FArray1D< int > FArray1D_int;
	// Symmetry
	typedef conformation::symmetry::SymmetricConformation SymmetricConformation;
	typedef conformation::symmetry::SymmetryInfoCOP SymmetryInfoCOP;

public: // construct/destruct

	/// @brief default constructor
	SecondaryStructurePotential();

	/// @brief default destructor
	inline
	~SecondaryStructurePotential()
	{}


public: // scoring

	void
	setup_for_scoring( pose::Pose & pose ) const;

	/// @brief score secondary structure
	void
	score(
		pose::Pose const & pose,
		SecondaryStructureWeights const & weights,
		Real & hs_score,
		Real & ss_score,
		Real & rsigma_score,
		Real & sheet_score
	) const;


	static
	void
	rsigma_dot_initializer(
		FArray4D_real & rsigma_dot
	);

private:

	/// @brief score hspair
	void
	hspair(
		pose::Pose const & pose,
		Real & hs_score
	) const;

	/// @brief score sspair
	void
	sspair(
		pose::Pose const & pose,
		SecondaryStructureWeights const & weights,
		Real & ss_score,
		Real & rsigma_score
	) const;


	float
	hairpin_killing_score(
		pose::Pose const & pose,
		Size const & pos1,
		Size const & pos2,
		Real const & theta,
		float const & ss_score
	) const;

	/// @brief This function takes a set of dimer neighbors, and determines how
	/// many sheets there, and how many strands are in each sheet
	/// This information is then used to calculate the "poker hand" score,
	/// which reflects to probability of that distribution of strands and
	/// sheets.

	Real
	sheets_from_dimers(
		pose::Pose const & pose
	) const;


private: // methods

	/// @brief identifies secondary structure at the start of scoring and loads the Strands/Helices
	/// data
	void
	identify_ss(
		pose::Pose const & pose,
		Helices & helices,
		Strands & strands
	) const;

	/// @brief function reads in two points in sequence and returns two points in space,
	/// the endpoints of the axis through an alpha-helix
	void
	helix_end(
		int const & pos1,
		BB_Pos const & bb_pos,
		Vector & p1,
		Vector & p2
	) const;

	/// @brief calculate sum of dot product of the co vectors of strand dimers ss1 and ss2
	/// with the vector connecting the midpoints of the dimer vectors (vdist)
	/// also determine return the sign of the dot products for each dimer
	/// to determine which direction the CO groups point
	void
	pair_dp(
		int const & ss1,
		int const & ss2,
		BB_Pos const & bb_pos,
		Real & dp,
		Vector const & vdist,
		int & sign1,
		int & sign2
	) const;

	/// @brief identifies the sequence separation along the fold tree
	/// add the gap_size (default 10 ) into seqsep when there is chain break between pos1 and pos2
	int
	get_foldtree_seqsep(
		pose::Pose const & pose,
		int pos1,
		int pos2,
		int gap_size=10
	) const;


private: // static math methods

	/// @brief find the vector connecting the midpoints of two line segments
	/// defined by a1,a2 and a3,a4
	/// @param[out] dist length of v21
	/// @param[out] v21  vector connecting midpoints
	static
	void
	dist_pair(
		Vector const & a1,
		Vector const & a2,
		Vector const & a3,
		Vector const & a4,
		Real & dist,
		Vector & cen1,
		Vector & cen2,
		Vector & v21
	);

	/// @brief find the angle (theta) between two vectors v1 and v2
	/// and the angle (phi) by which v2 is rotated out of the plane
	/// defined by v1 and v21 around v1 (ie v1 is the z-axis)
	static
	void
	spherical(
		Vector const & a2,
		Vector const & a4,
		Real & phi,
		Real & theta,
		Vector const & cen1,
		Vector const & cen2,
		Vector const & v21
	);

	/// @brief find angle sigma between vectors cen1->a2 and v21
	/// @param[out] sig sigma
	static
	void
	sigma(
		Vector const & a2,
		Vector const & cen1,
		Vector const & v21,
		Real & sig
	);


private: // initialization

	/// @brief load phi/theta bins for use in secondary structure scoring
	void
	load_phi_theta_bins(
		std::string const & hs_filename = "scoring/score_functions/SecondaryStructurePotential/phi.theta.36.HS.resmooth",
		std::string const & ss_filename = "scoring/score_functions/SecondaryStructurePotential/phi.theta.36.SS.resmooth"
	);


private: // initializers for static data

	static
	void
	idsn_initializer(
		FArray1D_int & idsn
	);

	static
	void
	ids_initializer(
		FArray1D_int & ids
	);

	static
	void
	ssdist_initializer(
		FArray2D_real & ssdist
	);

	static
	void
	hs_dp_initializer(
		FArray1D_real & hs_dp
	);

	/*static
	void
	rsigma_dot_initializer(
	FArray4D_real & rsigma_dot
	);*/

	static
	void
	m_term_initializer(
		FArray1D_real & m_term
	);

	/// @brief Penalty for pairing strand dimers that are close in sequence.
	/// @details Inferred from the log ratio of pairing probabilities of strands
	/// @details in the PDB vs. in Rosetta decoys. Calculated as a function of
	/// @details strand separation.
	static
	void
	ss_penalty_initializer(
		FArray1D_real & SS_penalty
	);

private: // secondary structure data

	FArray1D_int iptsn_;
	FArray4D_real pts_;
	FArray1D_real ds_;

	FArray1D_int const idsn_;
	FArray1D_int const ids_;
	FArray2D_real const ssdist_;
	FArray1D_real const hs_dp_; // not in use?
	FArray4D_real const rsigma_dot_; // lookup for new rsigma score
	FArray1D_real const m_term_;

};

} // ns scoring
} // ns core


#endif
