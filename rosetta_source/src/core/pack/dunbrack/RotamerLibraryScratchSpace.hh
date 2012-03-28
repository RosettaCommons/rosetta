// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/RotamerLibraryScratchSpace.hh
/// @brief  Declaration of scratch space class for Dunbrack rotamer library
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_dunbrack_RotamerLibraryScratchSpace_hh
#define INCLUDED_core_pack_dunbrack_RotamerLibraryScratchSpace_hh

// Unit headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>

// Package headers
// AUTO-REMOVED #include <core/scoring/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/exit.hh>

#include <utility/fixedsizearray1.hh>


namespace core {
namespace pack {
namespace dunbrack {


class RotamerLibraryScratchSpace : public utility::pointer::ReferenceCount
{
public:

	static Size const AA_PHI_INDEX = 1;
	static Size const AA_PSI_INDEX = 2;

public:

	RotamerLibraryScratchSpace();
	virtual ~RotamerLibraryScratchSpace();

	Real rotprob() const { return rotprob_; }
	Real negln_rotprob() const { return negln_rotprob_; }
	Size4 const & rotwell() const { return rotwell_; }
	Real4 const & chimean() const { return chimean_; }
	Real4 const & chisd() const { return chisd_; }
	Real4 const & chidev() const { return chidev_; }
	Real4 const & chidevpen() const { return chidevpen_; }
	Real3 const & drotprob_dbb() const { return drotprob_dbb_; } // Only preserved temporarily for use in the semirotameric libraries ('08, '10)
	Real3 const & dneglnrotprob_dbb() const { return dneglnrotprob_dbb_; }
	Real4 const & dchimean_dphi() const { return dchimean_dphi_; }
	Real4 const & dchimean_dpsi() const { return dchimean_dpsi_; }
	Real4 const & dchisd_dphi() const { return dchisd_dphi_; }
	Real4 const & dchisd_dpsi() const { return dchisd_dpsi_; }
	Real3 const & dchidevpen_dbb() const { return dchidevpen_dbb_; }
	Real4 const & dchidevpen_dchi() const { return dchidevpen_dchi_; }
	Real3 const & dE_dbb() const { return dE_dbb_; }
	Real4 const & dE_dchi() const { return dE_dchi_; }

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
	Real3  & drotprob_dbb()  { return drotprob_dbb_; } // TEMP
	Real3  & dneglnrotprob_dbb()  { return dneglnrotprob_dbb_; }
	Real4  & dchimean_dphi() { return dchimean_dphi_; }
	Real4  & dchimean_dpsi() { return dchimean_dpsi_; }
	Real4  & dchisd_dphi()   { return dchisd_dphi_; }
	Real4  & dchisd_dpsi()   { return dchisd_dpsi_; }
	Real3  & dchidevpen_dbb()   { return dchidevpen_dbb_; }
	Real4  & dchidevpen_dchi()  { return dchidevpen_dchi_; }
	Real3  & dE_dbb()  { return dE_dbb_; }
	Real4  & dE_dchi() { return dE_dchi_; }

	Real & fa_dun_tot() { return fa_dun_tot_; }
	Real & fa_dun_rot() { return fa_dun_rot_; }
	Real & fa_dun_semi() { return fa_dun_semi_; }
	Real & fa_dun_dev() { return fa_dun_dev_; }


private:

	Real rotprob_;
	Real negln_rotprob_;
	Size4 rotwell_;
	Real4 chimean_;
	Real4 chisd_;
	Real4 chidev_;
	Real4 chidevpen_;
	Real3 drotprob_dbb_;
	Real3 dneglnrotprob_dbb_;
	Real4 dchimean_dphi_;
	Real4 dchimean_dpsi_;
	Real4 dchisd_dphi_;
	Real4 dchisd_dpsi_;
	Real3 dchidevpen_dbb_;
	Real4 dchidevpen_dchi_;
	Real3 dE_dbb_;
	Real4 dE_dchi_;

	Real fa_dun_tot_;
	Real fa_dun_rot_;
	Real fa_dun_semi_;
	Real fa_dun_dev_;

};

}
}
}

#endif

