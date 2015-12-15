// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/types.hh
/// @brief  core::scoring package type declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)

#ifndef INCLUDED_core_scoring_hbonds_types_hh
#define INCLUDED_core_scoring_hbonds_types_hh

// Package headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/hbonds/HBEvalTuple.fwd.hh>

// Utility headers
#include <utility/exit.hh>

#include <core/types.hh>
#include <core/chemical/types.hh>


#include <numeric/xyzVector.fwd.hh>


#include <ObjexxFCL/FArray3D.fwd.hh>


namespace core {
namespace scoring {
namespace hbonds {

struct HBondDerivs {
	DerivVectorPair don_deriv;    // derivative vectors for the heavyatom donor for a hydrogen
	DerivVectorPair h_deriv;      // derivative vectors for the hydrogen forming a hydrogen bond
	DerivVectorPair acc_deriv;    // derivative vectors for the heavyatom acceptor for a hydrogen
	DerivVectorPair abase_deriv;  // derivative vectors for the acceptor base
	DerivVectorPair abase2_deriv; // derivative vectors for the acceptor base 2 -- for sp2 acceptors
#ifdef    SERIALIZATION
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/////////////////////////////////////////////////////////////////////////////////
/////// WARNING WARNING WARNING
///////
/////// if you modify the hbond types please update the strings name
/////// in ScoreTypeManager.cc
///////
/////// WARNING WARNING WARNING
/////////////////////////////////////////////////////////////////////////////////


enum HBondWeightType {
	hbw_NONE = 1,
	hbw_SR_BB,
	hbw_LR_BB,
	hbw_SR_BB_SC,
	hbw_LR_BB_SC,
	hbw_SC,
	hbw_MAX = hbw_SC
};

///////////////////////////////////////////////////////////////////////////////
////// WARNING WARNING WARNING
//////
////// Changing the HBAccChemType or HBDonChemType will change the
////// sort order in hbonds::hbtrie::HBAtom::operator<(...) and
////// hbonds::hbtrie::HBAtom::operator==(...). Changing the sort
////// order will __definitely__ cause trajectory changes.  So make
////// changes very carefully!
//////
//////////////////////////////////////////////////////////////////////////////


enum HBAccChemType {
	hbacc_NONE = 1,
	hbacc_PBA, // hbacc_PROTEIN_BB_AMIDE
	hbacc_CXA, // hbacc_CARBOXAMIDE
	hbacc_CXL, // hbacc_CARBOXYL
	hbacc_IMD, // hbacc_IMIDAZOL_DELTA
	hbacc_IME, // hbacc_IMIDAZOL_EPSILON
	hbacc_AHX, // hbacc_AROMATIC_HYDROXYL
	hbacc_HXL, // hbacc_HYDROXY
	hbacc_PCA_DNA, // hbacc_PHOSPHATE_CARBONYL_DNA
	hbacc_PES_DNA, // hbacc_PHOSPHATE_ESTER_DNA
	hbacc_RRI_DNA, // hbacc_RIBOSE_RING_DNA
	hbacc_PCA_RNA, // hbacc_PHOSPHATE_CARBONYL_RNA
	hbacc_PES_RNA, // hbacc_PHOSPHATE_ESTER_RNA
	hbacc_RRI_RNA, // hbacc_RIBOSE_RING_RNA
	hbacc_H2O, // hbacc_WATER
	hbacc_GENERIC_SP2BB,
	hbacc_GENERIC_SP2SC,
	hbacc_GENERIC_SP3BB,
	hbacc_GENERIC_SP3SC,
	hbacc_GENERIC_RINGBB,
	hbacc_GENERIC_RINGSC,
	hbacc_MAX = hbacc_GENERIC_RINGSC
};

enum HBDonChemType {
	hbdon_NONE = 1,
	hbdon_PBA, // hbdon_PROTEIN_BB_AMIDE
	hbdon_CXA, // hbdon_CARBOXAMIDE
	hbdon_IMD, // hbdon_IMIDAZOL_DELTA
	hbdon_IME, // hbdon_IMIDAZOL_EPSILON
	hbdon_IND, // hbdon_INDOL
	hbdon_AMO, // hbdon_AMINO
	hbdon_GDE, // hbdon_GUANIDINIUM_EPSILON
	hbdon_GDH, // hbdon_DIHYDRO_GUANIDINIUM
	hbdon_AHX, // hbdon_AROMATIC_HYDROXYL
	hbdon_HXL, // hbdon_HYDROXYL
	hbdon_H2O, // hbdon_WATER
	hbdon_GENERIC_BB,
	hbdon_GENERIC_SC,
	hbdon_MAX = hbdon_GENERIC_SC
};


enum HBEvalType{
	hbe_UNKNOWN = 0,
	hbe_NONE=1,
	hbe_dPBAaPBAsepM4helix, hbe_dPROTEIN_BB_AMIDEaPROTEIN_BB_AMIDEsepM4helix = hbe_dPBAaPBAsepM4helix,
	hbe_dPBAaPBAsepM3turn, hbe_dPROTEIN_BB_AMIDEaPROTEIN_BB_AMIDEsepM3turn = hbe_dPBAaPBAsepM3turn,
	hbe_dPBAaPBAsepM2turn, hbe_dPROTEIN_BB_AMIDEaPROTEIN_BB_AMIDEsepM2turn = hbe_dPBAaPBAsepM2turn,
	hbe_dPBAaPBAsepPM1, hbe_dPROTEIN_BB_AMIDEaPROTEIN_BB_AMIDEsepPM1 = hbe_dPBAaPBAsepPM1,
	hbe_dPBAaPBAsepP2turn, hbe_dPROTEIN_BB_AMIDEaPROTEIN_BB_AMIDEsepP2turn = hbe_dPBAaPBAsepP2turn,
	hbe_dPBAaPBAsepP3turn, hbe_dPROTEIN_BB_AMIDEaPROTEIN_BB_AMIDEsepP3turn = hbe_dPBAaPBAsepP3turn,
	hbe_dPBAaPBAsepP4helix, hbe_dPROTEIN_BB_AMIDEaPROTEIN_BB_AMIDEsepP4helix = hbe_dPBAaPBAsepP4helix,
	hbe_dPBAaPBAsepother, hbe_dPROTEIN_BB_AMIDEaPROTEIN_BB_AMIDEsepother = hbe_dPBAaPBAsepother,
	hbe_dCXAaPBAsepPM1, hbe_dCARBOXAMIDEaPROTEIN_BB_AMIDEsepPM1 = hbe_dCXAaPBAsepPM1,
	hbe_dIMDaPBAsepPM1, hbe_dIMIDAZOL_DELTAaPROTEIN_BB_AMIDEsepPM1 = hbe_dIMDaPBAsepPM1,
	hbe_dIMEaPBAsepPM1, hbe_dIMIDAZOL_EPSILONaPROTEIN_BB_AMIDEsepPM1 = hbe_dIMEaPBAsepPM1,
	hbe_dINDaPBAsepPM1, hbe_dINDOLaPROTEIN_BB_AMIDEsepPM1 = hbe_dINDaPBAsepPM1,
	hbe_dAMOaPBAsepPM1, hbe_dAMINOaPROTEIN_BB_AMIDEsepPM1 = hbe_dAMOaPBAsepPM1,
	hbe_dGDEaPBAsepPM1, hbe_dGUANIDINIUM_EPSILONaPROTEIN_BB_AMIDEsepPM1 = hbe_dGDEaPBAsepPM1,
	hbe_dGDHaPBAsepPM1, hbe_dDIHYDRO_GUANIDINIUMaPROTEIN_BB_AMIDEsepPM1 = hbe_dGDHaPBAsepPM1,
	hbe_dAHXaPBAsepPM1, hbe_dAROMATIC_HYDROXYLaPROTEIN_BB_AMIDEsepPM1 = hbe_dAHXaPBAsepPM1,
	hbe_dHXLaPBAsepPM1, hbe_dHYDROXYLaPROTEIN_BB_AMIDEsepPM1 = hbe_dHXLaPBAsepPM1,
	hbe_dCXAaPBAsepother, hbe_dCARBOXAMIDEaPROTEIN_BB_AMIDEsepother = hbe_dCXAaPBAsepother,
	hbe_dIMDaPBAsepother, hbe_dIMIDAZOL_DELTAaPROTEIN_BB_AMIDEsepother = hbe_dIMDaPBAsepother,
	hbe_dIMEaPBAsepother, hbe_dIMIDAZOL_EPSILONaPROTEIN_BB_AMIDEsepother = hbe_dIMEaPBAsepother,
	hbe_dINDaPBAsepother, hbe_dINDOLaPROTEIN_BB_AMIDEsepother = hbe_dINDaPBAsepother,
	hbe_dAMOaPBAsepother, hbe_dAMINOaPROTEIN_BB_AMIDEsepother = hbe_dAMOaPBAsepother,
	hbe_dGDEaPBAsepother, hbe_dGUANIDINIUM_EPSILONaPROTEIN_BB_AMIDEsepother = hbe_dGDEaPBAsepother,
	hbe_dGDHaPBAsepother, hbe_dDIHYDRO_GUANIDINIUMaPROTEIN_BB_AMIDEsepother = hbe_dGDHaPBAsepother,
	hbe_dAHXaPBAsepother, hbe_dAROMATIC_HYDROXYLaPROTEIN_BB_AMIDEsepother = hbe_dAHXaPBAsepother,
	hbe_dHXLaPBAsepother, hbe_dHYDROXYLaPROTEIN_BB_AMIDEsepother = hbe_dHXLaPBAsepother,
	hbe_dH2OaPBA, hbe_WATERaPROTEIN_BB_AMIDE = hbe_dH2OaPBA,
	hbe_dPBAaCXAsepPM1, hbe_dPROTEIN_BB_AMIDEaCARBOXAMIDEsepPM1 = hbe_dPBAaCXAsepPM1,
	hbe_dPBAaCXAsepother, hbe_dPROTEIN_BB_AMIDEaCARBOXAMIDEsepother = hbe_dPBAaCXAsepother,
	hbe_dCXAaCXA, hbe_dCARBOXAMIDEaCARBOXAMIDE = hbe_dCXAaCXA,
	hbe_dIMDaCXA, hbe_dIMIDAZOL_DELTAaCARBOXAMIDE = hbe_dIMDaCXA,
	hbe_dIMEaCXA, hbe_dIMIDAZOL_EPSILONaCARBOXAMIDE = hbe_dIMEaCXA,
	hbe_dINDaCXA, hbe_dINDOLaCARBOXAMIDE = hbe_dINDaCXA,
	hbe_dAMOaCXA, hbe_dAMINOaCARBOXAMIDE = hbe_dAMOaCXA,
	hbe_dGDEaCXA, hbe_dGUANIDINIUM_EPSILONaCARBOXAMIDE = hbe_dGDEaCXA,
	hbe_dGDHaCXA, hbe_dDIHYDRO_GUANIDINIUMaCARBOXAMIDE = hbe_dGDHaCXA,
	hbe_dAHXaCXA, hbe_dAROMATIC_HYDROXYLaCARBOXAMIDE = hbe_dAHXaCXA,
	hbe_dHXLaCXA, hbe_dHYDROXYLaCARBOXAMIDE = hbe_dHXLaCXA,
	hbe_dH2OaCXA, hbe_dWATERaCARBOXAMIDE = hbe_dH2OaCXA,
	hbe_dPBAaCXLsepPM1, hbe_dPROTEIN_BB_AMIDEaCARBOXYLsepPM1 = hbe_dPBAaCXLsepPM1,
	hbe_dPBAaCXLsepother, hbe_dPROTEIN_BB_AMIDEaCARBOXYLsepother = hbe_dPBAaCXLsepother,
	hbe_dCXAaCXL, hbe_dCARBOXAMIDEaCARBOXYL = hbe_dCXAaCXL,
	hbe_dIMDaCXL, hbe_dIMIDAZOL_DELTAaCARBOXYL = hbe_dIMDaCXL,
	hbe_dIMEaCXL, hbe_dIMIDAZOL_EPSILONaCARBOXYL = hbe_dIMEaCXL,
	hbe_dINDaCXL, hbe_dINDOLaCARBOXYL = hbe_dINDaCXL,
	hbe_dAMOaCXL, hbe_dAMINOaCARBOXYL = hbe_dAMOaCXL,
	hbe_dGDEaCXL, hbe_dGUANIDINIUM_EPSILONaCARBOXYL = hbe_dGDEaCXL,
	hbe_dGDHaCXL, hbe_dDIHYDRO_GUANIDINIUMaCARBOXYL = hbe_dGDHaCXL,
	hbe_dAHXaCXL, hbe_dAROMATIC_HYDROXYLaCARBOXYL = hbe_dAHXaCXL,
	hbe_dHXLaCXL, hbe_dHYDROXYLaCARBOXYL = hbe_dHXLaCXL,
	hbe_dH2OaCXL, hbe_dWATERaCARBOXYL = hbe_dH2OaCXL,
	hbe_dPBAaIMDsepPM1, hbe_dPROTEIN_BB_AMIDEaIMIDAZOL_DELTAsepPM1 = hbe_dPBAaIMDsepPM1,
	hbe_dPBAaIMDsepother, hbe_dPROTEIN_BB_AMIDEaIMIDAZOL_DELTAsepother = hbe_dPBAaIMDsepother,
	hbe_dCXAaIMD, hbe_dCARBOXAMIDEaIMIDAZOL_DELTA = hbe_dCXAaIMD,
	hbe_dIMDaIMD, hbe_dIMIDAZOL_DELTAaIMIDAZOL_DELTA = hbe_dIMDaIMD,
	hbe_dIMEaIMD, hbe_dIMIDAZOL_EPSILONaIMIDAZOL_DELTA = hbe_dIMEaIMD,
	hbe_dINDaIMD, hbe_dINDOLaIMIDAZOL_DELTA = hbe_dINDaIMD,
	hbe_dAMOaIMD, hbe_dAMINOaIMIDAZOL_DELTA = hbe_dAMOaIMD,
	hbe_dGDEaIMD, hbe_dGUANIDINIUM_EPSILONaIMIDAZOL_DELTA = hbe_dGDEaIMD,
	hbe_dGDHaIMD, hbe_dDIHYDRO_GUANIDINIUMaIMIDAZOL_DELTA = hbe_dGDHaIMD,
	hbe_dAHXaIMD, hbe_dAROMATIC_HYDROXYLaIMIDAZOL_DELTA = hbe_dAHXaIMD,
	hbe_dHXLaIMD, hbe_dHYDROXYLaIMIDAZOL_DELTA = hbe_dHXLaIMD,
	hbe_dH2OaIMD, hbe_dWATERaIMIDAZOL_DELTA = hbe_dH2OaIMD,
	hbe_dPBAaIMEsepPM1, hbe_dPROTEIN_BB_AMIDEaIMIDAZOL_EPSILONsepPM1 = hbe_dPBAaIMEsepPM1,
	hbe_dPBAaIMEsepother, hbe_dPROTEIN_BB_AMIDEaIMIDAZOL_EPSILONsepother = hbe_dPBAaIMEsepother,
	hbe_dCXAaIME, hbe_dCARBOXAMIDEaIMIDAZOL_EPSILON = hbe_dCXAaIME,
	hbe_dIMDaIME, hbe_dIMIDAZOL_DELTAaIMIDAZOL_EPSILON = hbe_dIMDaIME,
	hbe_dIMEaIME, hbe_dIMIDAZOL_EPSILONaIMIDAZOL_EPSILON = hbe_dIMEaIME,
	hbe_dINDaIME, hbe_dINDOLaIMIDAZOL_EPSILON = hbe_dINDaIME,
	hbe_dAMOaIME, hbe_dAMINOaIMIDAZOL_EPSILON = hbe_dAMOaIME,
	hbe_dGDEaIME, hbe_dGUANIDINIUM_EPSILONaIMIDAZOL_EPSILON = hbe_dGDEaIME,
	hbe_dGDHaIME, hbe_dDIHYDRO_GUANIDINIUMaIMIDAZOL_EPSILON = hbe_dGDHaIME,
	hbe_dAHXaIME, hbe_dAROMATIC_HYDROXYLaIMIDAZOL_EPSILON = hbe_dAHXaIME,
	hbe_dHXLaIME, hbe_dHYDROXYLaIMIDAZOL_EPSILON = hbe_dHXLaIME,
	hbe_dH2OaIME, hbe_dWATERaIMIDAZOL_EPSILON = hbe_dH2OaIME,
	hbe_dPBAaAHXsepPM1, hbe_dPROTEIN_BB_AMIDEaAROMATIC_HYDROXYLsepPM1 = hbe_dPBAaAHXsepPM1,
	hbe_dPBAaAHXsepother, hbe_dPROTEIN_BB_AMIDEaAROMATIC_HYDROXYLsepother = hbe_dPBAaAHXsepother,
	hbe_dCXAaAHX, hbe_dCARBOXAMIDEaAROMATIC_HYDROXYL = hbe_dCXAaAHX,
	hbe_dIMDaAHX, hbe_dIMIDAZOL_DELTAaAROMATIC_HYDROXYL = hbe_dIMDaAHX,
	hbe_dIMEaAHX, hbe_dIMIDAZOL_EPSILONaAROMATIC_HYDROXYL = hbe_dIMEaAHX,
	hbe_dINDaAHX, hbe_dINDOLaAROMATIC_HYDROXYL = hbe_dINDaAHX,
	hbe_dAMOaAHX, hbe_dAMINOaAROMATIC_HYDROXYL = hbe_dAMOaAHX,
	hbe_dGDEaAHX, hbe_dGUANIDINIUM_EPSILONaAROMATIC_HYDROXYL = hbe_dGDEaAHX,
	hbe_dGDHaAHX, hbe_dDIHYDRO_GUANIDINIUMaAROMATIC_HYDROXYL = hbe_dGDHaAHX,
	hbe_dAHXaAHX, hbe_dAROMATIC_HYDROXYLaAROMATIC_HYDROXYL = hbe_dAHXaAHX,
	hbe_dHXLaAHX, hbe_dHYDROXYLaAROMATIC_HYDROXYL = hbe_dHXLaAHX,
	hbe_dH2OaAHX, hbe_dWATERaAROMATIC_HYDROXYL = hbe_dH2OaAHX,
	hbe_dPBAaHXLsepPM1, hbe_dPROTEINS_BB_AMIDEaHYDROXYLsepPM1 = hbe_dPBAaHXLsepPM1,
	hbe_dPBAaHXLsepother, hbe_dPROTEINS_BB_AMIDEaHYDROXYLsepother = hbe_dPBAaHXLsepother,
	hbe_dCXAaHXL, hbe_dCARBOXAMIDEaHYDROXYL = hbe_dCXAaHXL,
	hbe_dIMDaHXL, hbe_dIMIDAZOL_DELTAaHYDROXYL = hbe_dIMDaHXL,
	hbe_dIMEaHXL, hbe_dIMIDAZOL_EPSILONaHYDXROXYL = hbe_dIMEaHXL,
	hbe_dINDaHXL, hbe_dINDOLaHYDROXYL = hbe_dINDaHXL,
	hbe_dAMOaHXL, hbe_dAMINOaHYDROXYL = hbe_dAMOaHXL,
	hbe_dGDEaHXL, hbe_dGUANIDINIUM_EPSILONaHYDROXYL = hbe_dGDEaHXL,
	hbe_dGDHaHXL, hbe_dDIHYDRO_GUANIDINIUMaHYDROXYL = hbe_dGDHaHXL,
	hbe_dAHXaHXL, hbe_dAROMATIC_HYDROXYLaHYDROXYL = hbe_dAHXaHXL,
	hbe_dHXLaHXL, hbe_dHYDROXYLaHYDROXYL = hbe_dHXLaHXL,
	hbe_dH2OaHXL, hbe_dWATERaHYDROXYL = hbe_dH2OaHXL,
	hbe_dPBAaPCA_DNAsepPM1, hbe_dPROTEINS_BB_AMIDEaPHOSPHATE_CARBONYL_DNAsepPM1 = hbe_dPBAaPCA_DNAsepPM1,
	hbe_dPBAaPCA_DNAsepother, hbe_dPROTEINS_BB_AMIDEaPHOSPHATE_CARBONYL_DNAsepother = hbe_dPBAaPCA_DNAsepother,
	hbe_dCXAaPCA_DNA, hbe_dCARBOXAMIDEaPHOSPHATE_CARBONYL_DNA = hbe_dCXAaPCA_DNA,
	hbe_dIMDaPCA_DNA, hbe_dIMIDAZOL_DELTAaPHOSPHATE_CARBONYL_DNA = hbe_dIMDaPCA_DNA,
	hbe_dIMEaPCA_DNA, hbe_dIMIDAZOL_EPSILONaPHOSPHATE_CARBONYL_DNA = hbe_dIMEaPCA_DNA,
	hbe_dINDaPCA_DNA, hbe_dINDOLaPHOSPHATE_CARBONYL_DNA = hbe_dINDaPCA_DNA,
	hbe_dAMOaPCA_DNA, hbe_dAMINOaPHOSPHATE_CARBONYL_DNA = hbe_dAMOaPCA_DNA,
	hbe_dGDEaPCA_DNA, hbe_dGUANIDINIUM_EPSILONaPHOSPHATE_CARBONYL_DNA = hbe_dGDEaPCA_DNA,
	hbe_dGDHaPCA_DNA, hbe_dDIHYDRO_GUANIDINIUMaPHOSPHATE_CARBONY_DNA = hbe_dGDHaPCA_DNA,
	hbe_dAHXaPCA_DNA, hbe_dAROMATIC_HYDROXYLaPHOSPHATE_CARBONYL_DNA = hbe_dAHXaPCA_DNA,
	hbe_dHXLaPCA_DNA, hbe_dHYDROXYLaPHOSPHATE_CARBONYL_DNA = hbe_dHXLaPCA_DNA,
	hbe_dH2OaPCA_DNA, hbe_dWATERaPHOSPHATE_CARBONYL_DNA = hbe_dH2OaPCA_DNA,
	hbe_dPBAaPCA_RNAsepPM1, hbe_dPROTEINS_BB_AMIDEaPHOSPHATE_CARBONYL_RNAsepPM1 = hbe_dPBAaPCA_RNAsepPM1,
	hbe_dPBAaPCA_RNAsepother, hbe_dPROTEINS_BB_AMIDEaPHOSPHATE_CARBONYL_RNAsepother = hbe_dPBAaPCA_RNAsepother,
	hbe_dCXAaPCA_RNAsepPM1, hbe_dCARBOXAMIDEaPHOSPHATE_CARBONYL_RNAsepPM1 = hbe_dCXAaPCA_RNAsepPM1,
	hbe_dCXAaPCA_RNAsepother, hbe_dCARBOXAMIDEaPHOSPHATE_CARBONYL_RNAsepother = hbe_dCXAaPCA_RNAsepother,
	hbe_dIMDaPCA_RNAsepPM1, hbe_dIMIDAZOL_DELTAaPHOSPHATE_CARBONYL_RNAsepPM1 = hbe_dIMDaPCA_RNAsepPM1,
	hbe_dIMDaPCA_RNAsepother, hbe_dIMIDAZOL_DELTAaPHOSPHATE_CARBONYL_RNAsepother = hbe_dIMDaPCA_RNAsepother,
	hbe_dIMEaPCA_RNAsepPM1, hbe_dIMIDAZOL_EPSILONaPHOSPHATE_CARBONYL_RNAsepPM1 = hbe_dIMEaPCA_RNAsepPM1,
	hbe_dIMEaPCA_RNAsepother, hbe_dIMIDAZOL_EPSILONaPHOSPHATE_CARBONYL_RNAsepother = hbe_dIMEaPCA_RNAsepother,
	hbe_dINDaPCA_RNAsepPM1, hbe_dINDOLaPHOSPHATE_CARBONYL_RNAsepPM1 = hbe_dINDaPCA_RNAsepPM1,
	hbe_dINDaPCA_RNAsepother, hbe_dINDOLaPHOSPHATE_CARBONYL_RNAsepother = hbe_dINDaPCA_RNAsepother,
	hbe_dAMOaPCA_RNAsepPM1, hbe_dAMINOaPHOSPHATE_CARBONYL_RNAsepPM1 = hbe_dAMOaPCA_RNAsepPM1,
	hbe_dAMOaPCA_RNAsepother, hbe_dAMINOaPHOSPHATE_CARBONYL_RNAsepother = hbe_dAMOaPCA_RNAsepother,
	hbe_dGDEaPCA_RNAsepPM1, hbe_dGUANIDINIUM_EPSILONaPHOSPHATE_CARBONYL_RNAsepPM1 = hbe_dGDEaPCA_RNAsepPM1,
	hbe_dGDEaPCA_RNAsepother, hbe_dGUANIDINIUM_EPSILONaPHOSPHATE_CARBONYL_RNAsepother = hbe_dGDEaPCA_RNAsepother,
	hbe_dGDHaPCA_RNAsepPM1, hbe_dDIHYDRO_GUANIDINIUMaPHOSPHATE_CARBONY_RNAsepPM1 = hbe_dGDHaPCA_RNAsepPM1,
	hbe_dGDHaPCA_RNAsepother, hbe_dDIHYDRO_GUANIDINIUMaPHOSPHATE_CARBONY_RNAsepother = hbe_dGDHaPCA_RNAsepother,
	hbe_dAHXaPCA_RNAsepPM1, hbe_dAROMATIC_HYDROXYLaPHOSPHATE_CARBONYL_RNAsepPM1 = hbe_dAHXaPCA_RNAsepPM1,
	hbe_dAHXaPCA_RNAsepother, hbe_dAROMATIC_HYDROXYLaPHOSPHATE_CARBONYL_RNAsepother = hbe_dAHXaPCA_RNAsepother,
	hbe_dHXLaPCA_RNAsepPM1, hbe_dHYDROXYLaPHOSPHATE_CARBONYL_RNAsepPM1 = hbe_dHXLaPCA_RNAsepPM1,
	hbe_dHXLaPCA_RNAsepother, hbe_dHYDROXYLaPHOSPHATE_CARBONYL_RNAsepother = hbe_dHXLaPCA_RNAsepother,
	hbe_dH2OaPCA_RNA, hbe_dWATERaPHOSPHATE_CARBONYL_RNA = hbe_dH2OaPCA_RNA,
	hbe_dPBAaPES_DNAsepPM1, hbe_dPROTEINS_BB_AMIDEaPHOSPHATE_ESTER_DNAsepPM1 = hbe_dPBAaPES_DNAsepPM1,
	hbe_dPBAaPES_DNAsepother, hbe_dPROTEINS_BB_AMIDEaPHOSPHATE_ESTER_DNAsepother = hbe_dPBAaPES_DNAsepother,
	hbe_dCXAaPES_DNA, hbe_dCARBOXAMIDEaPHOSPHATE_ESTER_DNA = hbe_dCXAaPES_DNA,
	hbe_dIMDaPES_DNA, hbe_dIMIDAZOL_DELTAaPHOSPHATE_ESTER_DNA = hbe_dIMDaPES_DNA,
	hbe_dIMEaPES_DNA, hbe_dIMIDAZOL_EPSILONaPHOSPHATE_ESTER_DNA = hbe_dIMEaPES_DNA,
	hbe_dINDaPES_DNA, hbe_dINDOLaPHOSPHATE_ESTER_DNA = hbe_dINDaPES_DNA,
	hbe_dAMOaPES_DNA, hbe_dAMINOaPHOSPHATE_ESTER_DNA = hbe_dAMOaPES_DNA,
	hbe_dGDEaPES_DNA, hbe_dGUANIDINIUM_EPSILONaPHOSPHATE_ESTER_DNA = hbe_dGDEaPES_DNA,
	hbe_dGDHaPES_DNA, hbe_dDIHYDRO_GUANIDINIUMaPHOSPHATE_ESTER_DNA = hbe_dGDHaPES_DNA,
	hbe_dAHXaPES_DNA, hbe_dAROMATIC_HYDROXYLaPHOSPHATE_ESTER_DNA = hbe_dAHXaPES_DNA,
	hbe_dHXLaPES_DNA, hbe_dHYDROXYLaPHOSPHATE_ESTER_DNA = hbe_dHXLaPES_DNA,
	hbe_dH2OaPES_DNA, hbe_dWATERaPHOSPHATE_ESTER_DNA = hbe_dH2OaPES_DNA,
	hbe_dPBAaPES_RNAsepPM1, hbe_dPROTEINS_BB_AMIDEaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dPBAaPES_RNAsepPM1,
	hbe_dPBAaPES_RNAsepother, hbe_dPROTEINS_BB_AMIDEaPHOSPHATE_ESTER_RNAsepother = hbe_dPBAaPES_RNAsepother,
	hbe_dCXAaPES_RNAsepPM1, hbe_dCARBOXAMIDEaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dCXAaPES_RNAsepPM1,
	hbe_dCXAaPES_RNAsepother, hbe_dCARBOXAMIDEaPHOSPHATE_ESTER_RNAsepother = hbe_dCXAaPES_RNAsepother,
	hbe_dIMDaPES_RNAsepPM1, hbe_dIMIDAZOL_DELTAaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dIMDaPES_RNAsepPM1,
	hbe_dIMDaPES_RNAsepother, hbe_dIMIDAZOL_DELTAaPHOSPHATE_ESTER_RNAsepother = hbe_dIMDaPES_RNAsepother,
	hbe_dIMEaPES_RNAsepPM1, hbe_dIMIDAZOL_EPSILONaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dIMEaPES_RNAsepPM1,
	hbe_dIMEaPES_RNAsepother, hbe_dIMIDAZOL_EPSILONaPHOSPHATE_ESTER_RNAsepother = hbe_dIMEaPES_RNAsepother,
	hbe_dINDaPES_RNAsepPM1, hbe_dINDOLaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dINDaPES_RNAsepPM1,
	hbe_dINDaPES_RNAsepother, hbe_dINDOLaPHOSPHATE_ESTER_RNAsepother = hbe_dINDaPES_RNAsepother,
	hbe_dAMOaPES_RNAsepPM1, hbe_dAMINOaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dAMOaPES_RNAsepPM1,
	hbe_dAMOaPES_RNAsepother, hbe_dAMINOaPHOSPHATE_ESTER_RNAsepother = hbe_dAMOaPES_RNAsepother,
	hbe_dGDEaPES_RNAsepPM1, hbe_dGUANIDINIUM_EPSILONaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dGDEaPES_RNAsepPM1,
	hbe_dGDEaPES_RNAsepother, hbe_dGUANIDINIUM_EPSILONaPHOSPHATE_ESTER_RNAsepother = hbe_dGDEaPES_RNAsepother,
	hbe_dGDHaPES_RNAsepPM1, hbe_dDIHYDRO_GUANIDINIUMaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dGDHaPES_RNAsepPM1,
	hbe_dGDHaPES_RNAsepother, hbe_dDIHYDRO_GUANIDINIUMaPHOSPHATE_ESTER_RNAsepother = hbe_dGDHaPES_RNAsepother,
	hbe_dAHXaPES_RNAsepPM1, hbe_dAROMATIC_HYDROXYLaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dAHXaPES_RNAsepPM1,
	hbe_dAHXaPES_RNAsepother, hbe_dAROMATIC_HYDROXYLaPHOSPHATE_ESTER_RNAsepother = hbe_dAHXaPES_RNAsepother,
	hbe_dHXLaPES_RNAsepPM1, hbe_dHYDROXYLaPHOSPHATE_ESTER_RNAsepPM1 = hbe_dHXLaPES_RNAsepPM1,
	hbe_dHXLaPES_RNAsepother, hbe_dHYDROXYLaPHOSPHATE_ESTER_RNAsepother = hbe_dHXLaPES_RNAsepother,
	hbe_dH2OaPES_RNA, hbe_dWATERaPHOSPHATE_ESTER_RNAsepother = hbe_dH2OaPES_RNA,
	hbe_dPBAaRRI_DNAsepPM1, hbe_dPROTEINS_BB_AMIDEaRIBOSE_RING_DNAsepPM1 = hbe_dPBAaRRI_DNAsepPM1,
	hbe_dPBAaRRI_DNAsepother, hbe_dPROTEINS_BB_AMIDEaRIBOSE_RING_DNAsepother = hbe_dPBAaRRI_DNAsepother,
	hbe_dCXAaRRI_DNA, hbe_dCARBOXAMIDEaRIBOSE_RING_DNA = hbe_dCXAaRRI_DNA,
	hbe_dIMDaRRI_DNA, hbe_dIMIDAZOL_DELTAaRIBOSE_RING_DNA = hbe_dIMDaRRI_DNA,
	hbe_dIMEaRRI_DNA, hbe_dIMIDAZOL_EPSILONaRIBOSE_RING_DNA = hbe_dIMEaRRI_DNA,
	hbe_dINDaRRI_DNA, hbe_dINDOLaRIBOSE_RING_DNA = hbe_dINDaRRI_DNA,
	hbe_dAMOaRRI_DNA, hbe_dAMINOaRIBOSE_RING_DNA = hbe_dAMOaRRI_DNA,
	hbe_dGDEaRRI_DNA, hbe_dGUANIDINIUM_EPSILONaRIBOSE_RING_DNA = hbe_dGDEaRRI_DNA,
	hbe_dGDHaRRI_DNA, hbe_dDIHYDRO_GUANIDINIUMaRIBOSE_RING_DNA = hbe_dGDHaRRI_DNA,
	hbe_dAHXaRRI_DNA, hbe_dAROMATIC_HYDROXYLaRIBOSE_RING_DNA = hbe_dAHXaRRI_DNA,
	hbe_dHXLaRRI_DNA, hbe_dHYDROXYLaRIBOSE_RING_DNA = hbe_dHXLaRRI_DNA,
	hbe_dH2OaRRI_DNA, hbe_dWATERaRIBOSE_RING_DNA = hbe_dH2OaRRI_DNA,
	hbe_dPBAaRRI_RNAsepPM1, hbe_dPROTEINS_BB_AMIDEaRIBOSE_RING_RNAsepPM1 = hbe_dPBAaRRI_RNAsepPM1,
	hbe_dPBAaRRI_RNAsepother, hbe_dPROTEINS_BB_AMIDEaRIBOSE_RING_RNAsepother = hbe_dPBAaRRI_RNAsepother,
	hbe_dCXAaRRI_RNAsepPM1, hbe_dCARBOXAMIDEaRIBOSE_RING_RNAsepPM1 = hbe_dCXAaRRI_RNAsepPM1,
	hbe_dCXAaRRI_RNAsepother, hbe_dCARBOXAMIDEaRIBOSE_RING_RNAsepother = hbe_dCXAaRRI_RNAsepother,
	hbe_dIMDaRRI_RNAsepPM1, hbe_dIMIDAZOL_DELTAaRIBOSE_RING_RNAsepPM1 = hbe_dIMDaRRI_RNAsepPM1,
	hbe_dIMDaRRI_RNAsepother, hbe_dIMIDAZOL_DELTAaRIBOSE_RING_RNAsepother = hbe_dIMDaRRI_RNAsepother,
	hbe_dIMEaRRI_RNAsepPM1, hbe_dIMIDAZOL_EPSILONaRIBOSE_RING_RNAsepPM1 = hbe_dIMEaRRI_RNAsepPM1,
	hbe_dIMEaRRI_RNAsepother, hbe_dIMIDAZOL_EPSILONaRIBOSE_RING_RNAsepother = hbe_dIMEaRRI_RNAsepother,
	hbe_dINDaRRI_RNAsepPM1, hbe_dINDOLaRIBOSE_RING_RNAsepPM1 = hbe_dINDaRRI_RNAsepPM1,
	hbe_dINDaRRI_RNAsepother, hbe_dINDOLaRIBOSE_RING_RNAsepother = hbe_dINDaRRI_RNAsepother,
	hbe_dAMOaRRI_RNAsepPM1, hbe_dAMINOaRIBOSE_RING_RNAsepPM1 = hbe_dAMOaRRI_RNAsepPM1,
	hbe_dAMOaRRI_RNAsepother, hbe_dAMINOaRIBOSE_RING_RNAsepother = hbe_dAMOaRRI_RNAsepother,
	hbe_dGDEaRRI_RNAsepPM1, hbe_dGUANIDINIUM_EPSILONaRIBOSE_RING_RNAsepPM1 = hbe_dGDEaRRI_RNAsepPM1,
	hbe_dGDEaRRI_RNAsepother, hbe_dGUANIDINIUM_EPSILONaRIBOSE_RING_RNAsepother = hbe_dGDEaRRI_RNAsepother,
	hbe_dGDHaRRI_RNAsepPM1, hbe_dDIHYDRO_GUANIDINIUMaRIBOSE_RING_RNAsepPM1 = hbe_dGDHaRRI_RNAsepPM1,
	hbe_dGDHaRRI_RNAsepother, hbe_dDIHYDRO_GUANIDINIUMaRIBOSE_RING_RNAsepother = hbe_dGDHaRRI_RNAsepother,
	hbe_dAHXaRRI_RNAsepPM1, hbe_dAROMATIC_HYDROXYLaRIBOSE_RING_RNAsepPM1 = hbe_dAHXaRRI_RNAsepPM1,
	hbe_dAHXaRRI_RNAsepother, hbe_dAROMATIC_HYDROXYLaRIBOSE_RING_RNAsepother = hbe_dAHXaRRI_RNAsepother,
	hbe_dHXLaRRI_RNAsepPM1, hbe_dHYDROXYLaRIBOSE_RING_RNAsepPM1 = hbe_dHXLaRRI_RNAsepPM1,
	hbe_dHXLaRRI_RNAsepother, hbe_dHYDROXYLaRIBOSE_RING_RNAsepother = hbe_dHXLaRRI_RNAsepother,
	hbe_dH2OaRRI_RNA, hbe_dWATERaRIBOSE_RING_RNAsepother = hbe_dH2OaRRI_RNA,
	hbe_dPBAaH2O, hbe_dPROTEIN_BB_AMIDEaWATER = hbe_dPBAaH2O,
	hbe_dCXAaH2O, hbe_dCARBOXAMIDEaWATER = hbe_dCXAaH2O,
	hbe_dIMDaH2O, hbe_dIMIDAZOL_DELTAaWATER = hbe_dIMDaH2O,
	hbe_dIMEaH2O, hbe_dIMIDAZOL_EPSILONaWATER = hbe_dIMEaH2O,
	hbe_dINDaH2O, hbe_dINDOLaWATER = hbe_dINDaH2O,
	hbe_dAMOaH2O, hbe_dAMINOaWATER = hbe_dAMOaH2O,
	hbe_dGDEaH2O, hbe_dGUANIDINIUM_EPSILONaWATER = hbe_dGDEaH2O,
	hbe_dGDHaH2O, hbe_dDIHYDRO_GUANIDINIUMaWATER = hbe_dGDHaH2O,
	hbe_dAHXaH2O, hbe_dAROMATIC_HYDROXYLaWATER = hbe_dAHXaH2O,
	hbe_dHXLaH2O, hbe_dHYDROXYLaWATER = hbe_dHXLaH2O,
	hbe_dH2OaH2O, hbe_dWATERaWATER = hbe_dH2OaH2O,
	hbe_GENERIC_SP2BB_SR,
	hbe_GENERIC_SP2BB_LR,
	hbe_GENERIC_SP3BB_SR,
	hbe_GENERIC_SP3BB_LR,
	hbe_GENERIC_RINGBB_SR,
	hbe_GENERIC_RINGBB_LR,
	hbe_GENERIC_SP2BSC_SR,
	hbe_GENERIC_SP2BSC_LR,
	hbe_GENERIC_SP3BSC_SR,
	hbe_GENERIC_SP3BSC_LR,
	hbe_GENERIC_RINGBSC_SR,
	hbe_GENERIC_RINGBSC_LR,
	hbe_GENERIC_SP2SCSC_SR,
	hbe_GENERIC_SP2SCSC_LR,
	hbe_GENERIC_SP3SCSC_SR,
	hbe_GENERIC_SP3SCSC_LR,
	hbe_GENERIC_RINGSCSC_SR,
	hbe_GENERIC_RINGSCSC_LR,
	hbe_MAX = hbe_GENERIC_RINGSCSC_LR
};

enum HBSeqSep{
	seq_sep_other = 1, // // all other sequence separation not specified
	seq_sep_M4, // // acc_rsd.seqpos() - don_rsd.seqpos() = -4
	seq_sep_M3, // // acc_rsd.seqpos() - don_rsd.seqpos() = -3
	seq_sep_M2, // // acc_rsd.seqpos() - don_rsd.seqpos() = -2
	seq_sep_PM1, // // std::abs(acc_rsd.seqpos() - don_rsd.seqpos()) = 1
	seq_sep_P2, // // acc_rsd.seqpos() - don_rsd.seqpos() = 2
	seq_sep_P3, // // acc_rsd.seqpos() - don_rsd.seqpos() = -3
	seq_sep_P4, // // acc_rsd.seqpos() - don_rsd.seqpos() = 4
	seq_sep_MAX = seq_sep_P4
};


void
HBEval_lookup_initializer( ObjexxFCL::FArray3D<HBEvalType> & hbe );

extern ObjexxFCL::FArray3D<HBEvalType> const HBEval_lookup;


///////////////////////////////////////////////////////////
//:::DEPRICATION NOTICE::::
//
// Note backbone/sidechain identification is used
//1) packing: eg since the backbone is fixed, backbone-backbone, and
// backbone-sidechain interactions can be precomputed.  This is found
// by looking at rsd.is_backbone()
//2) environmental dependence is only applied to side chains, (ask
//Tanja or Lin, because I don't know why)
//3) to take some account of variable correlation,
//sidechain-sidechain hbonds have two different angle potentials
//that depend on the AH distance so hbe_is_BB_type is used.
//
//
//Backbone/sidechain distinctions are tricky when one looks at
//modified-proteins or non-proteins.  Since these functions are only
//used in questional circumstances they are depricated and will be
//removed in subsiquent versions of the hbond-potential.
//

bool hbe_is_BB_type( HBEvalType hbe );

bool hbe_is_SC_type( HBEvalType hbe );


HBondWeightType
get_hbond_weight_type( HBEvalType const & hbe_type );


chemical::Hybridization
get_hbe_acc_hybrid( HBEvalType const & hbe );

enum HBGeoDimType {
	hbgd_NONE = 1,

	// distance from the acceptor atom to the hydrogen
	hbgd_AHdist,

	// cosine of the base-acceptor-hydrogen angle
	// the base is acceptor hybridization dependent:
	//    sp2 hybrid -> base = res.atom_base(atm_num)
	//    sp3 hybrid -> base = res.abase2(atm_num)
	//    ring hybrid-> base = (res.atom_base(atm_num) + res.abase2(atm_num))/2
	// let BAunit = unit vector from base to the acceptor
	// let AHunit = unit vector from acceptor to hydrogen
	// cosBAH = BAunit <dot> AHunit
	hbgd_cosBAH,

	// cosine of the acceptor-hydrogen-donor angle
	// let AHunit = unit vector from acceptor to the hydrogen
	// let HDunit = unit vector from hydrogen to donor
	// cosAHD = AHDunit <dot> HDunit
	hbgd_cosAHD,


	// the angle formed by the acceptor-hydrogen-donor
	// this the interior angle measured in radians.

	// In score12, the hydrogen bond score function evaluated the cosine
	// of exterior BAH and AHD angles rather than the angles
	// themselves. This was done for two reasons:
	//
	//    1) The cosine of the exterior angle is easy to evaluate, just
	//    take the dot product of the noralized bond vectors.
	//
	//    2) When projecting uniform density density over cartesian
	//    space onto the theta angle in spherical coordinates, the
	//    resulting density is not uniform. This happens because the
	//    change in the volume of the conic section per unit angle
	//    depends on the angle itself. It turns out that the
	//    distribution of the cosine of the angle is uniform. Therefore
	//    when estimating distributions, as is done for knowledge based
	//    potentials, one should normalize the distribution by computing
	//    the distrbution in "cosine" space.
	//
	// Because density estimation should be done in cosine space and the
	// cosine of the angles is easy to evaluate, the polynomials in the
	// hydrogen bond score function were defined as functions of the
	// cosine of the angles.
	//
	// A limitation of this parametrization is that the dynamic range
	// from optimal AHD angle (180 degrees) to decent AHD angle (~160
	// degrees) is compressed. This can be seen by noticing that the
	// acos(x) around zero is steep, so a small change in x results in
	// large change in acos(x). To create polynomials that have such a
	// tight distribution requires them to be relatively high degree.
	//
	// As an alternative, the AHD angle can be used directly in the
	// parametrization.
	hbgd_AHD,

	// Torsional angle about the base-acceptor bond vector
	// Not yet implemented (11/09) but coming soon...
	hbgd_chi,

	hbgd_MAX=hbgd_chi
};

enum HBDerivType {
	hbderiv_NONE = 1, // no derivative
	hbderiv_ABE_GO, // standard hbond derivative calculation
	hbderiv_ABE_GO_GEOMSOL_OCC_ACC, // geometric solvation, occluding atom ('pseudo-water') + acceptor
	hbderiv_ABE_GO_GEOMSOL_OCC_DON, // geometric solvation, occluding atom ('pseudo-water') + donor
	hbderiv_MAX = hbderiv_ABE_GO_GEOMSOL_OCC_DON
};

extern Real DUMMY_DERIV;
extern HBondDerivs DUMMY_DERIVS;  // f1/f2 vectors for four atoms
extern HBondDerivs const ZERO_DERIV2D;

Size
hb_eval_type(
	HBDonChemType don_chem_type,
	HBAccChemType acc_chem_type,
	HBSeqSep seq_sep_type
);


} // namespace hbonds
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_types_HH
