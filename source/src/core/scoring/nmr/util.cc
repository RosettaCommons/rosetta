// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/util.cc
/// @brief   function definitions for util.hh
/// @details last Modified: 05/11/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/util.hh>
#include <core/io/nmr/AtomSelection.hh>

// Project headers
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/Exceptions.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/string.functions.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <cmath>

namespace core {
namespace scoring {
namespace nmr {

static basic::Tracer TR( "core.scoring.nmr.util" );

/// @brief the type how equivalent, ambiguous spins are averaged
NMR_VALUE_AVERAGING_TYPE
convert_string_to_averaging_type(std::string const & averaging_type) {
	std::string str = ObjexxFCL::uppercased(averaging_type);
	if ( str == "SUM" ) {
		return SUM;
	} else { // str == "MEAN"
		return MEAN;
	}
}

SINGLE_NMR_VALUE_WEIGHTING
convert_string_to_weighting_scheme(std::string const & weighting_scheme) {
	std::string str = ObjexxFCL::uppercased(weighting_scheme);
	if ( str == "CONST" ) {
		return CONST;
	} else if ( str == "SIGMA" ) {
		return SIGMA;
	} else { //  str == "OBSIG"
		return OBSIG;
	}
}

RDC_NORM_TYPE
convert_string_to_normalization_type(std::string const & computation_type) {
	std::string str = ObjexxFCL::uppercased(computation_type);
	if ( str == "NH" ) {
		return NORM_TYPE_NH;
	} else if ( str == "CH" ) {
		return NORM_TYPE_CH;
	} else { // str == "NONE"
		return NORM_TYPE_NONE;
	}
}

PRE_RATE_TYPE
convert_string_to_rate_type(std::string const & rate_type) {
	std::string str = ObjexxFCL::uppercased(rate_type);
	if ( str == "R2" ) {
		return R2_PARA;
	} else if ( str == "R1" ) {
		return R1_PARA;
	} else {
		utility_exit_with_message( "ERROR: Provided string does not match any PRE rate type. Possible options are \"R2\" or \"R1\".");
	}
}

std::string
convert_rdc_type_to_string(RDC_TYPE const & type) {
	if ( type == RDC_TYPE_NH ) {
		return "NH";
	} else if ( type == RDC_TYPE_NCO ) {
		return "NCO";
	} else if ( type == RDC_TYPE_NCA ) {
		return "NCA";
	} else if ( type == RDC_TYPE_CAHA ) {
		return "CAHA";
	} else if ( type == RDC_TYPE_CAHN ) {
		return "CAHN";
	} else if ( type == RDC_TYPE_COHN ) {
		return "COHN";
	} else if ( type == RDC_TYPE_CACO ) {
		return "CACO";
	} else { // RDC_TYPE_CACB
		return "CACB";
	}
}

std::string
convert_pre_spin_type_to_string(PRE_SPIN_TYPE const & type) {
	if ( type == PRE_SPIN_TYPE_H ) {
		return "1H";
	} else if ( type == PRE_SPIN_TYPE_N ) {
		return "15N";
	} else { // PRE_SPIN_TYPE_C
		return "13C";
	}
}

/// @brief: function that handles conversion from pseudoatom identifier for degenerate protons to fullatom name
utility::vector1<id::AtomID>
lookup_pseudoprotons(
	Size const rsd,
	std::string const & atomname,
	pose::Pose const & pose
)
{
	using namespace chemical;
	using chemical::AA;
	using id::NamedAtomID;

	runtime_assert( rsd >= 1 || rsd <= pose.total_residue() );

	utility::vector1<id::AtomID> Atoms;

	AA aa_type(pose.residue_type(rsd).aa());

	if ( (atomname == "H" || atomname == "H*") && rsd == 1 && pose.is_fullatom() ) { // N-terminus
		Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1H", rsd ), pose ) );
		Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2H", rsd ), pose ) );
		Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3H", rsd ), pose ) );
	} else if ( aa_type == aa_ala || aa_type == aa_dal ) { // L-Alanine or D-Alanine
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HB", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_arg || aa_type == aa_dar ) { // L-Arginine or D-Arginine
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else if ( atomname == "HG*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG", rsd ), pose ) );
		} else if ( atomname == "HD*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HD", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HD", rsd ), pose ) );
		} else if ( atomname == "HH1*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HH1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HH1", rsd ), pose ) );
		} else if ( atomname == "HH2*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HH2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HH2", rsd ), pose ) );
		} else if ( atomname == "HH*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HH1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HH1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HH2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HH2", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_asp || aa_type == aa_cys || aa_type == aa_his
			|| aa_type == aa_ser || aa_type == aa_trp
			|| aa_type == aa_das || aa_type == aa_dcs || aa_type == aa_dhi
			|| aa_type == aa_dse || aa_type == aa_dtr ) { // L-Aspartate or L-Cysteine or L-Histidine or L-Serine or L-Tryptophane
		if ( atomname == "HB*" ) {        // D-Aspartate or D-Cysteine or D-Histidine or D-Serine or D-Tryptophane
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_asn || aa_type == aa_dan ) { // L-Asparagine or D-Asparagine
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else if ( atomname == "HD*" || atomname == "HD2*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HD2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HD2", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_glu || aa_type == aa_dgu ) { // L-Glutamate or D-Glutamate
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else if ( atomname == "HG*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_gln || aa_type == aa_dgn ) { // L-Glutamine or D-Glutamine
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else if ( atomname == "HG*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG", rsd ), pose ) );
		} else if ( atomname == "HE*" || atomname == "HE2*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HE2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HE2", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_gly ) { // Glycine
		if ( atomname == "HA*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HA", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HA", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_ile || aa_type == aa_dil ) { // L-Isoleucine or D-Isoleucine
		if ( atomname == "HG1*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG1", rsd ), pose ) );
		} else if ( atomname == "HG2*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HG2", rsd ), pose ) );
		} else if ( atomname == "HD1*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HD1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HD1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HD1", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_leu || aa_type == aa_dle ) { // L-Leucine or D-Leucine
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else if ( atomname == "HD1*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HD1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HD1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HD1", rsd ), pose ) );
		} else if ( atomname == "HD2*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HD2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HD2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HD2", rsd ), pose ) );
		} else if ( atomname == "HD*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HD1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HD1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HD1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HD2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HD2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HD2", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_lys || aa_type == aa_dly ) { // L-Lysine or D-Lysine
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else if ( atomname == "HG*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG", rsd ), pose ) );
		} else if ( atomname == "HD*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HD", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HD", rsd ), pose ) );
		} else if ( atomname == "HE*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HE", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HE", rsd ), pose ) );
		} else if ( atomname == "HZ*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HZ", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HZ", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HZ", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_met || aa_type == aa_dme ) { // L-Methionine or D-Methionine
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else if ( atomname == "HG*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG", rsd ), pose ) );
		} else if ( atomname == "HE*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HE", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HE", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HE", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_phe || aa_type == aa_tyr || aa_type == aa_dph || aa_type == aa_dty ) { // L-Phenylalanine or L-Tyrosine or D-Phenylalanine or D-Tyrosine
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else if ( atomname == "HD*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "HD1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "HD2", rsd ), pose ) );
		} else if ( atomname == "HE*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "HE1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "HE2", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_pro || aa_type == aa_dpr ) { // L-Proline or D-Proline
		if ( atomname == "HB*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HB", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HB", rsd ), pose ) );
		} else if ( atomname == "HG*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG", rsd ), pose ) );
		} else if ( atomname == "HD*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HD", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HD", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_thr || aa_type == aa_dth ) { // L-Threonine or D-Threonine
		if ( atomname == "HG2*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HG2", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else if ( aa_type == aa_val || aa_type == aa_dva ) { // L-Valine or D-Valine
		if ( atomname == "HG1*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HG1", rsd ), pose ) );
		} else if ( atomname == "HG2*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HG2", rsd ), pose ) );
		} else if ( atomname == "HG*" ) {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HG1", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "1HG2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "2HG2", rsd ), pose ) );
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( "3HG2", rsd ), pose ) );
		} else {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ) , pose ) );
		}
	} else {
		TR.Warning << "Translation of NMR pseudoatom names is only supported for canonical L/D amino acids so far. Trying to convert atom names directly into atom IDs." << std::endl;
		try {
			Atoms.push_back(named_atom_id_to_atom_id( NamedAtomID( atomname, rsd ), pose ) );
		} catch ( id::EXCN_AtomNotFound & excn ) {
			utility_exit_with_message("ERROR while translating NMR pseudoatom names into atom IDs. " + excn.msg());
		}
	}
	return Atoms;
}

/// @brief determine the RDC type from the spins atom names
RDC_TYPE
rdc_type_from_atom_names(std::pair< core::io::nmr::AtomSelection, core::io::nmr::AtomSelection > const & spinsAB) {
	if ( ( (spinsAB.first.get_atom() == "H" || spinsAB.first.get_atom() == "h") && (spinsAB.second.get_atom() == "N" || spinsAB.second.get_atom() == "n") ) ||
			( (spinsAB.first.get_atom() == "N" || spinsAB.first.get_atom() == "n") && (spinsAB.second.get_atom() == "H" || spinsAB.second.get_atom() == "h") ) ) {
		return RDC_TYPE_NH;
	} else if ( ( (spinsAB.first.get_atom() == "N" || spinsAB.first.get_atom() == "n") && (spinsAB.second.get_atom() == "C" || spinsAB.second.get_atom() == "c") ) ||
			( (spinsAB.first.get_atom() == "C" || spinsAB.first.get_atom() == "c") && (spinsAB.second.get_atom() == "N" || spinsAB.second.get_atom() == "n") ) ) {
		return RDC_TYPE_NCO;
	} else if ( ( (spinsAB.first.get_atom() == "N" || spinsAB.first.get_atom() == "n") && (spinsAB.second.get_atom() == "CA" || spinsAB.second.get_atom() == "ca") ) ||
			( (spinsAB.first.get_atom() == "CA" || spinsAB.first.get_atom() == "ca") && (spinsAB.second.get_atom() == "N" || spinsAB.second.get_atom() == "n") ) ) {
		return RDC_TYPE_NCA;
	} else if ( ( (spinsAB.first.get_atom() == "CA" || spinsAB.first.get_atom() == "ca") && (spinsAB.second.get_atom() == "HA" || spinsAB.second.get_atom() == "ha") ) ||
			( (spinsAB.first.get_atom() == "HA" || spinsAB.first.get_atom() == "ha") && (spinsAB.second.get_atom() == "CA" || spinsAB.second.get_atom() == "ca") ) ) {
		return RDC_TYPE_CAHA;
	} else if ( ( (spinsAB.first.get_atom() == "CA" || spinsAB.first.get_atom() == "ca") && (spinsAB.second.get_atom() == "HA*" || spinsAB.second.get_atom() == "ha*") ) ||
			( (spinsAB.first.get_atom() == "HA*" || spinsAB.first.get_atom() == "ha*") && (spinsAB.second.get_atom() == "CA" || spinsAB.second.get_atom() == "ca") ) ) {
		// special case of glycine 1HA, 2HA
		return RDC_TYPE_CAHA;
	} else if ( ( (spinsAB.first.get_atom() == "CA" || spinsAB.first.get_atom() == "ca") && (spinsAB.second.get_atom() == "H" || spinsAB.second.get_atom() == "h") ) ||
			( (spinsAB.first.get_atom() == "H" || spinsAB.first.get_atom() == "h") && (spinsAB.second.get_atom() == "CA" || spinsAB.second.get_atom() == "ca") ) ) {
		return RDC_TYPE_CAHN;
	} else if ( ( (spinsAB.first.get_atom() == "C" || spinsAB.first.get_atom() == "c") && (spinsAB.second.get_atom() == "H" || spinsAB.second.get_atom() == "h") ) ||
			( (spinsAB.first.get_atom() == "H" || spinsAB.first.get_atom() == "h") && (spinsAB.second.get_atom() == "C" || spinsAB.second.get_atom() == "c") ) ) {
		return RDC_TYPE_COHN;
	} else if ( ( (spinsAB.first.get_atom() == "CA" || spinsAB.first.get_atom() == "ca") && (spinsAB.second.get_atom() == "C" || spinsAB.second.get_atom() == "c") ) ||
			( (spinsAB.first.get_atom() == "C" || spinsAB.first.get_atom() == "c") && (spinsAB.second.get_atom() == "CA" || spinsAB.second.get_atom() == "ca") ) ) {
		return RDC_TYPE_CACO;
	} else if ( ( (spinsAB.first.get_atom() == "CA" || spinsAB.first.get_atom() == "ca") && (spinsAB.second.get_atom() == "CB" || spinsAB.second.get_atom() == "cb") ) ||
			( (spinsAB.first.get_atom() == "CB" || spinsAB.first.get_atom() == "cb") && (spinsAB.second.get_atom() == "CA" || spinsAB.second.get_atom() == "ca") ) ) {
		return RDC_TYPE_CACB;
	} else {
		utility_exit_with_message( "ERROR: The RDC data type could not be determined from the atom names. Possible atom names are \"H\", \"N\", \"C\", \"CA\", \"CB\", \"HA\" or \"HA*\" for glycine." );
	}
}

PRE_SPIN_TYPE
pre_spin_type_from_atom_name(core::io::nmr::AtomSelection const & atom) {
	std::string name = ObjexxFCL::uppercased(atom.get_atom());
	if ( name.find("HN") != std::string::npos ) {   // NMR specific convention to refer to backbone amides
		return PRE_SPIN_TYPE_H;
	} else if ( name.find('C') != std::string::npos ) { // all C atoms
		return PRE_SPIN_TYPE_C;
	} else if ( name.find('N') != std::string::npos && name.find("HN") == std::string::npos ) { // all N atoms
		return PRE_SPIN_TYPE_N;
	} else if ( name.find('H') != std::string::npos ) { // Now only H atoms are left
		return PRE_SPIN_TYPE_H;
	} else {
		utility_exit_with_message( "ERROR: The PRE spin type could not be inferred from the atom names. Possible PRE spin types include: \"H\", \"HN\", \"C\" and \"N\".");
	}
}

/// @brief: creates a rotation matrix given three Euler angles
///         and a convention as determined at runtime
Matrix rotation_matrix_from_euler_angles(
	Vector const & angles,
	EULER_CONVENTION convention
)
{
	Matrix rotM;
	if ( convention == ZXZ_CONVENTION ) {
		rotM = numeric::rotation_matrix_from_euler_angles_ZXZ(angles);
	} else if ( convention == ZYZ_CONVENTION ) {
		rotM = numeric::rotation_matrix_from_euler_angles_ZYZ(angles);
	} else if ( convention == ZYX_CONVENTION ) {
		rotM = numeric::rotation_matrix_from_euler_angles_ZYX(angles);
	}
	return rotM;
}

/// @brief: determines the three Euler angles given a rotation matrix
///         and the respective Euler convention as determined at runtime
Vector
euler_angles_from_rotation_matrix(
	Matrix const & rotM,
	EULER_CONVENTION convention
)
{
	Vector euler;
	if ( convention == ZXZ_CONVENTION ) {
		euler = numeric::euler_angles_from_rotation_matrix_ZXZ(rotM);
	} else if ( convention == ZYZ_CONVENTION ) {
		euler = numeric::euler_angles_from_rotation_matrix_ZYZ(rotM);
	} else if ( convention == ZYX_CONVENTION ) {
		euler = numeric::euler_angles_from_rotation_matrix_ZYX(rotM);
	}
	return euler;
}

/// @brief orders the three eigenvalues of NMRTensor and brings NMRTensor in unique tensor representation
void order_tensor_parameters(
	utility::vector1<Real> & params,
	EULER_CONVENTION convention
)
{
	// Core code has been taken from Numbat. Schmitz C. et al. (2008) J.Biomol.NMR 41(3):179-189
	Real aEigxx = fabs(params[1]);
	Real aEigyy = fabs(params[2]);
	Real aEigzz = fabs(params[3]);
	Vector angles(params[4], params[5], params[6]);

	if ( (aEigzz >= aEigyy) && (aEigyy >= aEigxx) ) {        // case 1
		;
	} else if ( (aEigzz >= aEigxx) && (aEigxx >= aEigyy) ) { // case 2
		params[6] += 90.0;
		std::swap(params[2], params[1]);
	} else if ( (aEigyy >= aEigzz) && (aEigzz >= aEigxx) ) { // case 3
		std::swap(params[2], params[3]);
		Matrix rotX90(Matrix::rows(1,0,0, 0,0,1, 0,-1,0));
		Matrix rotM = rotation_matrix_from_euler_angles(angles, convention);
		Matrix transformed_rotM = rotX90 * rotM;
		Vector updated_angles = euler_angles_from_rotation_matrix(transformed_rotM, convention);
		params[4] = updated_angles(1);
		params[5] = updated_angles(2);
		params[6] = updated_angles(3);
	} else if ( (aEigyy >= aEigxx) && (aEigxx >= aEigzz) ) { // case 4
		angles(3) += 90.0;
		std::swap(params[2], params[1]);
		std::swap(params[3], params[1]);
		Matrix rotY90(Matrix::rows(0,0,-1, 0,1,0, 1,0,0));
		Matrix rotM = rotation_matrix_from_euler_angles(angles, convention);
		Matrix transformed_rotM = rotY90 * rotM;
		Vector updated_angles = euler_angles_from_rotation_matrix(transformed_rotM, convention);
		params[4] = updated_angles(1);
		params[5] = updated_angles(2);
		params[6] = updated_angles(3);
	} else if ( (aEigxx >= aEigzz) && (aEigzz >= aEigyy) ) { // case 5
		angles(3) += 90.0;
		std::swap(params[2], params[1]);
		std::swap(params[2], params[3]);
		Matrix rotX90(Matrix::rows(1,0,0, 0,0,1, 0,-1,0));
		Matrix rotM = rotation_matrix_from_euler_angles(angles, convention);
		Matrix transformed_rotM = rotX90 * rotM;
		Vector updated_angles = euler_angles_from_rotation_matrix(transformed_rotM, convention);
		params[4] = updated_angles(1);
		params[5] = updated_angles(2);
		params[6] = updated_angles(3);
	} else if ( (aEigxx >= aEigyy) && (aEigyy >= aEigzz) ) { // case 6
		std::swap(params[3], params[1]);
		Matrix rotY90(Matrix::rows(0,0,-1, 0,1,0, 1,0,0));
		Matrix rotM = rotation_matrix_from_euler_angles(angles, convention);
		Matrix transformed_rotM = rotY90 * rotM;
		Vector updated_angles = euler_angles_from_rotation_matrix(transformed_rotM, convention);
		params[4] = updated_angles(1);
		params[5] = updated_angles(2);
		params[6] = updated_angles(3);
	}
}

/// @brief Utility function that transforms an NMR tensor from the frame of residue_from
///        into the coordinate frame of the symmetric residue_to given a symmetric pose.
///        Returns the original matrix if the pose is not symmetric or if residue_from and
///        residue_to belong to the same subunit residue_to.
Matrix
apply_tensor_transformation(
	pose::Pose & pose,
	Matrix const & Min,
	Size resid_from,
	Size resid_to
)
{
	using namespace conformation::symmetry;
	using namespace pose::symmetry;

	if ( is_symmetric(pose) ) {
		SymmetryInfoCOP syminfo = symmetry_info( pose );
		SymmetricConformation & symmconf = dynamic_cast< SymmetricConformation & >( pose.conformation() );
		if ( syminfo->get_num_components() > 1 ) {
			runtime_assert_msg(syminfo->get_component_of_residue(resid_from) == syminfo->get_component_of_residue(resid_to), "ERROR: Residues used for matrix transformation must be in the same component.");
		}
		HT Tsym_from = symmconf.get_transformation(resid_from);
		HT Tsym_to = symmconf.get_transformation(resid_to);

		Matrix Msym_from = Tsym_from.rotation_matrix();
		Matrix Msym_to = Tsym_to.rotation_matrix();
		Matrix Q = Msym_to * Msym_from.transpose();
		Matrix Mout;

		// Since we are rotating the tensor here and not a vector, we have to multiply on both sides Mout = Q * Min * Q.T
		// In contrast, for a vector we would have to calculate Vecout = Q * Vecin
		Mout = Min * Q.transpose();
		Mout = Q * Mout;

		return Mout;
	}
	return Min;
}

/// @brief Utility function that transforms a vector of metal (or spinlabel) coordinates from
///        the frame of residue_from into the coordinate frame of the symmetric residue_to given a symmetric pose.
///        Returns the original vector if the pose is not symmetric or if residue_from and
///        residue_to belong to the same subunit residue_to.
Vector
apply_vector_transformation(
	pose::Pose & pose,
	Vector const & Vin,
	Size resid_from,
	Size resid_to
)
{
	using namespace conformation::symmetry;
	using namespace pose::symmetry;

	if ( is_symmetric(pose) ) {
		SymmetryInfoCOP syminfo = symmetry_info( pose );
		SymmetricConformation & symmconf = dynamic_cast< SymmetricConformation & >( pose.conformation() );
		if ( syminfo->get_num_components() > 1 ) {
			runtime_assert_msg(syminfo->get_component_of_residue(resid_from) == syminfo->get_component_of_residue(resid_to), "ERROR: Residues used for matrix transformation must be in the same component.");
		}
		HT Tsym_from = symmconf.get_transformation(resid_from);
		HT Tsym_to = symmconf.get_transformation(resid_to);
		Vector Vout;

		Vout = Tsym_from.inverse() * Vin;
		Vout = Tsym_to * Vout;
		return Vout;
	}
	return Vin;
}

/// @brief Utility function that rotates a vector from the frame of residue_from
///        into the coordinate frame of the symmetric residue_to given a symmetric pose.
///        Returns the original vector if the pose is not symmetric or if residue_from and
///        residue_to belong to the same subunit residue_to.
Vector
apply_vector_rotation(
	pose::Pose & pose,
	Vector const & Vin,
	Size resid_from,
	Size resid_to
)
{
	using namespace conformation::symmetry;
	using namespace pose::symmetry;

	if ( is_symmetric(pose) ) {
		SymmetryInfoCOP syminfo = symmetry_info( pose );
		SymmetricConformation & symmconf = dynamic_cast< SymmetricConformation & >( pose.conformation() );
		if ( syminfo->get_num_components() > 1 ) {
			runtime_assert_msg(syminfo->get_component_of_residue(resid_from) == syminfo->get_component_of_residue(resid_to), "ERROR: Residues used for matrix transformation must be in the same component.");
		}
		HT Tsym_from = symmconf.get_transformation(resid_from);
		HT Tsym_to = symmconf.get_transformation(resid_to);

		Matrix Msym_from = Tsym_from.rotation_matrix();
		Matrix Msym_to = Tsym_to.rotation_matrix();
		Matrix Q = Msym_to * Msym_from.transpose();
		Vector Vout;

		Vout = Q * Vin;
		return Vout;
	}
	return Vin;
}

/// @brief Create a local coordinate frame from backbone coordinates.
///        The local frame is defined by a homogeneous transform object
///        which can be created from three axes and one point.
///        The center is located at CA. The z'-axis is along CA - CB.
///        The x'-axis is within the CA-CB - CA-C plane and the
///        y'-axis is perpendicular to the CB - CA - C plane.
HT
define_local_frame(
	Vector const & CA,
	Vector const & CB,
	Vector const & C
)
{
	Vector vCACB(CB - CA);
	vCACB.normalize();
	Vector vCAC(C - CA);
	Real dotprod(vCAC.dot(vCACB));
	vCAC -= dotprod * vCACB;
	vCAC.normalize();
	Vector Zaxis(vCACB);
	Vector Xaxis(vCAC);
	Vector Yaxis(vCACB.cross(vCAC));
	return HT(Xaxis, Yaxis, Zaxis, CA);
}

/// @brief Auxiliary PCS function
/// @params par: Tensor values [xM, yM, zM, Xax, Xrh]
/// rotM: Rotation matrix to transform spin coordinates in tensor frame
/// spin_coord: Spin xyz coordinates
/// scal: optional scaling factor
Real
pcs_func(
	Vec5 const & par,
	Matrix const & rotM,
	Vector const & spin_coord,
	Real const & scal
)
{
	// vector between spin and metal center
	Real x(spin_coord.x() - par[1]);
	Real y(spin_coord.y() - par[2]);
	Real z(spin_coord.z() - par[3]);

	// transformed vector after rotation
	Real x_t(rotM(1,1)*x + rotM(1,2)*y + rotM(1,3)*z);
	Real y_t(rotM(2,1)*x + rotM(2,2)*y + rotM(2,3)*z);
	Real z_t(rotM(3,1)*x + rotM(3,2)*y + rotM(3,3)*z);

	Real r2(x_t*x_t + y_t*y_t + z_t*z_t);
	Real r5(r2 * r2 * std::sqrt(r2));
	Real value_1_12_PI_r5(10000.0 / (12.0 * numeric::constants::d::pi * r5));

	Real pcs(scal * value_1_12_PI_r5 * (par[4] * (3.0 * z_t*z_t - r2) + par[5] * 1.5 * (x_t*x_t - y_t*y_t)));
	return pcs;
}

} // nmr
} // scoring
} // core
