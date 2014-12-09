// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBEvalTuple.hh
/// @brief  Tuple describing data about the donor and acceptor in a single hbond
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// unit headers
#include <core/scoring/hbonds/HBEvalTuple.hh>

// Package headers
#include <core/scoring/hbonds/hbonds_geom.hh>

// Project headers
#include <core/conformation/Residue.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray3D.hh>

namespace core {
namespace scoring {
namespace hbonds {

HBEvalTuple::HBEvalTuple(
	HBDonChemType don,
	HBAccChemType acc,
	HBSeqSep sequence_sep
) :
	don_type_( don ),
	acc_type_( acc ),
	seq_sep_( sequence_sep )
{
	update_hbevaltype();
}

HBEvalTuple::HBEvalTuple(
	int const datm,
	core::conformation::Residue const & don_rsd,
	int const aatm,
	core::conformation::Residue const & acc_rsd
)
{
	don_type_ = get_hb_don_chem_type(datm, don_rsd);
	acc_type_ = get_hb_acc_chem_type(aatm, acc_rsd);
	seq_sep_ = get_seq_sep(don_type_, acc_type_, don_rsd.polymeric_oriented_sequence_distance(acc_rsd));
	eval_type_ = HBEval_lookup(don_type_, acc_type_, seq_sep_);
}

bool
operator==(HBEvalTuple const & a, HBEvalTuple const & b){
	return(
		a.don_type_ == b.don_type_ &&
		a.acc_type_ == b.acc_type_ &&
		a.seq_sep_ == b.seq_sep_ &&
		a.eval_type_ == b.eval_type_
	);
}


void HBEvalTuple::don_type( HBDonChemType don )
{
	don_type_ = don;
	update_hbevaltype();
}

void HBEvalTuple::acc_type( HBAccChemType acc )
{
	acc_type_ = acc;
	update_hbevaltype();
}

void HBEvalTuple::sequence_sep( HBSeqSep seqsep )
{
	seq_sep_ = seqsep;
	update_hbevaltype();
}

void HBEvalTuple::update_hbevaltype()
{
		eval_type_ = HBEval_lookup(don_type_, acc_type_, seq_sep_ );
}



} // namespace hbonds
} // namespace scoring
} // namespace core
