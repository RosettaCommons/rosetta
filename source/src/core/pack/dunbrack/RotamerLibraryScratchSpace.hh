// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/dunbrack/RotamerLibraryScratchSpace.hh
/// @brief  Declaration of scratch space class for Dunbrack rotamer library
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_dunbrack_RotamerLibraryScratchSpace_hh
#define INCLUDED_core_pack_dunbrack_RotamerLibraryScratchSpace_hh

// Unit headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>

// Package headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {


class RotamerLibraryScratchSpace : public utility::pointer::ReferenceCount
{
public:
	static Size const AA_OMG_INDEX = 3;
	static Size const AA_PHI_INDEX = 1;
	static Size const AA_PSI_INDEX = 2;

public:

	RotamerLibraryScratchSpace();
	//RotamerLibraryScratchSpace();


	virtual ~RotamerLibraryScratchSpace();

	Real rotprob() const { return rotprob_; }
	Real negln_rotprob() const { return negln_rotprob_; }
	Size4 const & rotwell() const { return rotwell_; }
	Real4 const & chimean() const { return chimean_; }
	Real4 const & chisd() const { return chisd_; }
	Real4 const & chidev() const { return chidev_; }
	Real4 const & chidevpen() const { return chidevpen_; }
	Real5 const & drotprob_dbb() const { return drotprob_dbb_; } // Only preserved temporarily for use in the semirotameric libraries ('08, '10)
	Real5 const & dneglnrotprob_dbb() const { return dneglnrotprob_dbb_; }
	FiveReal4 const & dchimean_dbb() const { return dchimean_dbb_; }
	FiveReal4 const & dchisd_dbb() const { return dchisd_dbb_; }
	Real5 const & dchidevpen_dbb() const { return dchidevpen_dbb_; }
	Real4 const & dchidevpen_dchi() const { return dchidevpen_dchi_; }
	Real5 const & dE_dbb() const { return dE_dbb_; }
	Real5 const & dE_dbb_dev() const { return dE_dbb_dev_; }
	Real5 const & dE_dbb_semi() const { return dE_dbb_semi_; }
	Real4 const & dE_dchi() const { return dE_dchi_; }
	Real4 const & dE_dchi_dev() const { return dE_dchi_dev_; }
	Real4 const & dE_dchi_semi() const { return dE_dchi_semi_; }

	FiveReal4 const & dE_dbb_dev_perchi() const { return dE_dbb_dev_perchi_; }

	Real fa_dun_tot() const { return fa_dun_tot_; }
	Real fa_dun_rot() const { return fa_dun_rot_; }
	Real fa_dun_semi() const { return fa_dun_semi_; }
	Real fa_dun_dev() const { return fa_dun_dev_; }


	Real   & rotprob()  { return rotprob_; }
	Real   & negln_rotprob() { return negln_rotprob_; }
	Size4  & rotwell()  { return rotwell_; }
	Real4  & chimean()  { return chimean_; }
	Real4  & chisd()    { return chisd_; }
	Real4  & chidev()   { return chidev_; }
	Real4  & chidevpen()     { return chidevpen_; }
	Real5  & drotprob_dbb()  { return drotprob_dbb_; } // TEMP
	Real5  & dneglnrotprob_dbb()  { return dneglnrotprob_dbb_; }
	FiveReal4  & dchimean_dbb() { return dchimean_dbb_; }
	FiveReal4  & dchisd_dbb()   { return dchisd_dbb_; }
	Real5  & dchidevpen_dbb()   { return dchidevpen_dbb_; }
	Real4  & dchidevpen_dchi()  { return dchidevpen_dchi_; }
	Real5  & dE_dbb()  { return dE_dbb_; }
	Real5  & dE_dbb_dev()  { return dE_dbb_dev_; }
	Real5  & dE_dbb_rot()  { return dE_dbb_rot_; }
	Real5  & dE_dbb_semi()  { return dE_dbb_semi_; }
	Real4  & dE_dchi() { return dE_dchi_; }
	Real4  & dE_dchi_dev() { return dE_dchi_dev_; }
	Real4  & dE_dchi_semi() { return dE_dchi_semi_; }

	FiveReal4 & dE_dbb_dev_perchi() { return dE_dbb_dev_perchi_; }

	Real & fa_dun_tot() { return fa_dun_tot_; }
	Real & fa_dun_rot() { return fa_dun_rot_; }
	Real & fa_dun_semi() { return fa_dun_semi_; }
	Real & fa_dun_dev() { return fa_dun_dev_; }

	// Entropy correction
	Real    entropy() const { return entropy_; }
	Real5   dentropy_dbb() const { return dentropy_dbb_; }
	Real  & entropy() { return entropy_; }
	Real5 & dentropy_dbb() { return dentropy_dbb_; }

private:

	//////////////////////////////////////////////////////
	// IMPORTANT: If you add a new member variable here,
	// be sure it's zero initialized in the constructor.
	// (Not all uses of RotamerLibraryScratchSpace will be
	// handed to a Dunbrack rotamer library before use.)
	//////////////////////////////////////////////////////

	Real rotprob_;
	Real negln_rotprob_;
	Size4 rotwell_;
	Real4 chimean_;
	Real4 chisd_;
	Real4 chidev_;
	Real4 chidevpen_;
	Real5 drotprob_dbb_;
	Real5 dneglnrotprob_dbb_;
	/*Real4 dchimean_dphi_;
	Real4 dchimean_dpsi_;
	Real4 dchisd_dphi_;
	Real4 dchisd_dpsi_;*/
	FiveReal4 dchimean_dbb_;
	FiveReal4 dchisd_dbb_;
	Real5 dchidevpen_dbb_;
	Real4 dchidevpen_dchi_;
	Real5 dE_dbb_;
	Real5 dE_dbb_dev_;
	Real5 dE_dbb_rot_;
	Real5 dE_dbb_semi_;
	Real4 dE_dchi_;
	Real4 dE_dchi_dev_;
	Real4 dE_dchi_semi_;

	// fpd per-chi components of bb derivs
	//Real4 dE_dphi_dev_;
	//Real4 dE_dphi_rot_;
	//Real4 dE_dpsi_dev_;
	FiveReal4 dE_dbb_dev_perchi_;
	//Real4 dE_dpsi_rot_;

	Real fa_dun_tot_;
	Real fa_dun_rot_;
	Real fa_dun_semi_;
	Real fa_dun_dev_;

	// Entropic correction
	Real entropy_;
	Real5 dentropy_dbb_;

	// DOUG DOUG DOUG May need to make a sub class, for peptoid rotlibs
public:
	Real4 const & dchimean_domg() const { return dchimean_domg_; }
	Real4 const & dchisd_domg() const { return dchisd_domg_; }
	Real4  & dchimean_domg() { return dchimean_domg_; }
	Real4  & dchisd_domg()   { return dchisd_domg_; }

private:
	Real4 dchimean_domg_;
	Real4 dchisd_domg_;

};

}
}
}

#endif

