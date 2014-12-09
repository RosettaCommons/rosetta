// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief implementation of class FilterByProxTerm.
/// @author Andrea Bazzoli

#include <devel/constel/FilterByProxTerm.hh>
#include <devel/constel/ChainTerm.hh>
#include <devel/constel/cnl_info.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/vector1.hh>


namespace devel {
namespace constel {

vector1<ChainTerm> FilterByProxTerm::chains_;
Size FilterByProxTerm::nchains_;
Real FilterByProxTerm::dct_max_2_;
Size FilterByProxTerm::nres_;
std::map<char, bool> FilterByProxTerm::proxnc_;


///
/// @brief filter initialization.
///
/// @param[in] ps pose containing the constellations to be filtered..
/// @param[in] dct maximum distance for a constellation to be considered
/// 	proximal to the termini.
/// @param[in] dtt maximum distance for chain termini to be considered
/// 	proximal to one another..
/// @param[in] number of residues forming either terminus of a chain.
///
void FilterByProxTerm::init(Pose const &ps, Real dct, Real dtt, Size nres) {

	nres_ = nres;
	dct_max_2_ = dct*dct;

	get_chain_terms(ps, chains_);
	nchains_ = chains_.size();

	Real dtt2 = dtt*dtt;
	for(Size i=1; i<=nchains_; ++i)
		if(has_prox_termini(ps, chains_[i], nres_, dtt2))
			proxnc_[chains_[i].get_cid()] = true;
}


///
/// @brief tells whether a constellation satisfies the filter.
///
/// @param[in] ps pose to which the constellation belongs.
/// @param[in] cnl indexes in the pose of the residues forming the
/// 	constellation.
///
/// @return true if 1) the constellation is sufficiently close to the termini
/// 	of at least one chain that it belongs to and 2) such termini are proximal
///   to one another; false otherwise.
///
/// @details in the case of proximity, the search stops at the first chain to
/// 	which the constellation	belongs and such that an N- or C-terminal
/// 	residue is proximal to the constellation. The N-terminus is considered
///   first, in ascending order; the C-terminus is considered afterwards, in
/// 	descending order; within either terminus the search stops at the first
///   residue that is proximal to the constellation. Proximity is determined
/// 	based on the distance between the center of mass of the constellation
/// 	and the CA atom of the residue.
///
bool FilterByProxTerm::sat(Pose const &ps, vector1<Size> const &cnl) {

	// identify set of chains spanned by constellation
	std::map<char, bool> cspan;
	for(Size i=1; i<=cnl.size(); ++i) {
		char cid = ps.pdb_info()->chain(cnl[i]);
		cspan[cid] = true;
	}

	// compute center of mass of constellation
	xyzVector<Real> com(cnl_com(cnl, ps));

	// detect proximity between center of mass and chain termini
	typedef std::map<char, bool>::const_iterator CI;
	for(CI i=cspan.begin(); i!=cspan.end(); ++i) {
		char cid = i->first;
		if(proxnc_[cid])
			for(Size j=1; j<=nchains_; ++j) {
				ChainTerm const &ctm = chains_[j];
				if(ctm.get_cid() == cid) {

					Size const NPS = ctm.get_nps();
					Size const CPS = ctm.get_cps();
					Size const ALL = CPS-NPS+1;

					// try N-terminus
					Size const NSTA = NPS;
					Size const NEND = (nres_ <= ALL) ? NSTA + nres_ : CPS+1;
					for(Size k=NSTA; k<NEND; ++k)
						if(com.distance_squared(ps.residue(k).xyz("CA")) <= dct_max_2_)
							return true;

					// try C-terminus
					Size const CSTA = CPS;
					Size const CEND = (nres_ <= ALL) ? CSTA - nres_ : NPS-1;
					for(Size k=CSTA; k>CEND; --k)
						if(com.distance_squared(ps.residue(k).xyz("CA")) <= dct_max_2_)
            	return true;

					break;
				}
			}
	}

	return false;
}


///
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
/// 	on one terminus is within the maximum distance of a residue on the other
/// 	terminus. Distance is computed on C-alpha atoms.
///
bool has_prox_termini(Pose const &ps, ChainTerm const &chain, Size NRES,
	Real DMAX2) {

	Size const NSTA = chain.get_nps();
	Size const NEND = NSTA + NRES;
	Size const CSTA = chain.get_cps();
	Size const CEND = CSTA - NRES;

	if(!(ps.residue(NSTA).has("CA"))) // non-protein chain
		return false;

	if((CSTA-NSTA) <= (2*NRES - 2)) // overlapping terminals
		return true;

	for(Size j=NSTA; j<NEND; ++j) {
		xyzVector<Real> const &nca = ps.residue(j).xyz("CA");
		for(Size k=CSTA; k>CEND; --k) {
			xyzVector<Real> const &cca = ps.residue(k).xyz("CA");
			Real d2 = nca.distance_squared(cca);
			if(d2 <= DMAX2)
				return true;
		}
	}

	return false;
}

} // constel
} // devel
