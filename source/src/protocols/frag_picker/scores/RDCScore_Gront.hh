// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/RDCScore.hh
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_RDCScore_hh
#define INCLUDED_protocols_frag_picker_scores_RDCScore_hh

#include <protocols/frag_picker/scores/RDCScore.fwd.hh>

// package headers
// AUTO-REMOVED #include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.fwd.hh>
#include <protocols/frag_picker/scores/AtomBasedConstraintsScore.hh>

// mini
// AUTO-REMOVED #include <core/scoring/ResidualDipolarCoupling.hh>
// AUTO-REMOVED #include <core/scoring/func/Func.hh>
// AUTO-REMOVED #include <core/scoring/func/FuncFactory.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

//Auto Headers
#include <core/scoring/ResidualDipolarCoupling.fwd.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief Calculates a score for a fragment based on RDC experimental data
class RDCScore: public AtomBasedConstraintsScore {
public:

	/// @brief Prepare the scoring method
	RDCScore(Size, Real, bool, Size);

	/// @brief In this case caching means copying coordinates of relevant atoms from a chunk's pose
	void do_caching(VallChunkOP);

	/// @brief Erases the internal array of coordinates
	void clean_up();

	/// @brief Evaluate the score
	bool score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	utility::vector1<core::scoring::RDC> rdc_data_;
	utility::vector1<std::string> rdc_atoms_;
	utility::vector1<std::string>& rdc_atoms();

	void read_RDC_file(std::string const &, Size);
	utility::vector1<std::string> rdc_atom_names_;

	//------- RDC computing stuff ---------
	inline Real sqr(Real x) const {
		return x * x;
	}
	inline Real iprod(const Real a[3], const Real b[3]) const {
		return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
	}
	inline void mvmul(Real a[3][3], const Real src[3], Real dest[3]) const {
		dest[0] = a[0][0] * src[0] + a[0][1] * src[1] + a[0][2] * src[2];
		dest[1] = a[1][0] * src[0] + a[1][1] * src[1] + a[1][2] * src[2];
		dest[2] = a[2][0] * src[0] + a[2][1] * src[1] + a[2][2] * src[2];
	}

	void jacobi(Real a[5][5], Real d[], Real v[5][5], int *nrot);
	int m_inv_gen(Real m[5][5], int n, Real minv[5][5]);
};

} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_RDCScore_HH */

