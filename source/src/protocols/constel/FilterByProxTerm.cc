// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constel/FilterByProxTerm.cc
/// @brief implementation of class FilterByProxTerm.
/// @author Andrea Bazzoli

#include <protocols/constel/FilterByProxTerm.hh>
#include <protocols/constel/ChainTerm.hh>
#include <protocols/constel/cnl_info.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace constel {

using core::pose::Pose;
using numeric::xyzVector;
using utility::vector1;
using core::Size;
using core::Real;

vector1<ChainTerm> FilterByProxTerm::chains_;
Size FilterByProxTerm::nchains_;
Real FilterByProxTerm::max_ct_dist2_;
Size FilterByProxTerm::nres_;
std::map<char, bool> FilterByProxTerm::proxnc_;


/// @brief filter initialization.
///
/// @param[in] ps pose containing the constellations to be filtered..
/// @param[in] max_ct_dist maximum distance for a constellation to be considered
///  proximal to the termini.
/// @param[in] max_tt_dist maximum distance for chain termini to be considered
///  proximal to one another..
/// @param[in] number of residues forming either terminus of a chain.
///
void FilterByProxTerm::init(Pose const &ps, Real max_ct_dist, Real max_tt_dist, core::Size nres) {

	nres_ = nres;
	max_ct_dist2_ = max_ct_dist*max_ct_dist;

	get_chain_terms(ps, chains_);
	nchains_ = chains_.size();

	Real max_tt_dist2 = max_tt_dist*max_tt_dist;
	for ( core::Size i=1; i<=nchains_; ++i ) {
		if ( has_prox_termini(ps, chains_[i], nres_, max_tt_dist2) ) {
			proxnc_[chains_[i].get_cid()] = true;
		}
	}
}


/// @brief tells whether a constellation satisfies the filter.
///
/// @param[in] ps pose to which the constellation belongs.
/// @param[in] cnl indexes in the pose of the residues forming the
///  constellation.
///
/// @return true if 1) the constellation is sufficiently close to the termini
///  of at least one chain that it belongs to and 2) such termini are proximal
///   to one another; false otherwise.
///
/// @details in the case of proximity, the search stops at the first chain to
///  which the constellation belongs and such that an N- or C-terminal
///  residue is proximal to the constellation. The N-terminus is considered
///   first, in ascending order; the C-terminus is considered afterwards, in
///  descending order; within either terminus the search stops at the first
///   residue that is proximal to the constellation. Proximity is determined
///  based on the distance between the center of mass of the constellation
///  and the CA atom of the residue.
///
bool FilterByProxTerm::is_satisfied(Pose const &ps, vector1<core::Size> const &cnl) {

	// identify set of chains spanned by constellation
	std::map<char, bool> cspan;
	for ( core::Size i=1; i<=cnl.size(); ++i ) {
		char cid = ps.pdb_info()->chain(cnl[i]);
		cspan[cid] = true;
	}

	// compute center of mass of constellation
	xyzVector<Real> com(cnl_com(cnl, ps));

	// detect proximity between center of mass and chain termini
	typedef std::map<char, bool>::const_iterator CI;
	for ( CI i=cspan.begin(); i!=cspan.end(); ++i ) {
		char cid = i->first;
		if ( proxnc_[cid] ) {
			for ( core::Size j=1; j<=nchains_; ++j ) {
				ChainTerm const &ctm = chains_[j];
				if ( ctm.get_cid() == cid ) {

					core::Size const NPS = ctm.get_nps();
					core::Size const CPS = ctm.get_cps();
					core::Size const ALL = CPS-NPS+1;

					// try N-terminus
					core::Size const NSTA = NPS;
					core::Size const NEND = (nres_ <= ALL) ? NSTA + nres_ : CPS+1;
					for ( core::Size k=NSTA; k<NEND; ++k ) {
						if ( com.distance_squared(ps.residue(k).xyz("CA")) <= max_ct_dist2_ ) {
							return true;
						}
					}

					// try C-terminus
					core::Size const CSTA = CPS;
					core::Size const CEND = (nres_ <= ALL) ? CSTA - nres_ : NPS-1;
					for ( core::Size k=CSTA; k>CEND; --k ) {
						if ( com.distance_squared(ps.residue(k).xyz("CA")) <= max_ct_dist2_ ) {
							return true;
						}
					}

					break;
				}
			}
		}
	}

	return false;
}


/// @brief Tells whether a chain has proximal N- and C-termini.
///
/// @param[in] ps pose to which the chain belongs.
/// @param[in] chain the chain.
/// @param[in] NRES number of residues forming either terminus.
/// @param[in] DMAX2 squared maximum distance for proximity.
///
/// @return true if the chain has proximal termini; false otherwise.
///
/// @details The termini are considered to be proximal if at least one residue
///  on one terminus is within the maximum distance of a residue on the other
///  terminus. Distance is computed on C-alpha atoms.
///
bool has_prox_termini(Pose const &ps, ChainTerm const &chain, core::Size NRES,
	Real DMAX2) {

	core::Size const NSTA = chain.get_nps();
	core::Size const NEND = NSTA + NRES;
	core::Size const CSTA = chain.get_cps();
	core::Size const CEND = CSTA - NRES;

	if ( (CSTA-NSTA) <= (2*NRES - 2) ) { // overlapping terminals
		return true;
	}

	for ( core::Size j=NSTA; j<NEND; ++j ) {
		xyzVector<Real> const &nca = ps.residue(j).xyz("CA");
		for ( core::Size k=CSTA; k>CEND; --k ) {
			xyzVector<Real> const &cca = ps.residue(k).xyz("CA");
			Real d2 = nca.distance_squared(cca);
			if ( d2 <= DMAX2 ) {
				return true;
			}
		}
	}

	return false;
}

} // constel
} // protocols
