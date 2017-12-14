// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverOutputter.cc
/// @brief an abstract base class for handling calculation outputs from PeptideDeriverFilter
/// @author orlypolo (orlymarcu@gmail.com)

#include <protocols/peptide_deriver/PeptideDeriverOutputter.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

static basic::Tracer TR( "protocols.peptide_deriver.PeptideDeriverOutputter" );


namespace protocols {
namespace peptide_deriver {

PeptideDeriverOutputter::PeptideDeriverOutputter() = default;

PeptideDeriverOutputter::~PeptideDeriverOutputter()= default;

// helper function: avoid +- signed on zero-values
core::Real PeptideDeriverOutputter::avoid_negative_zero(core::Real const value, core::Real const threshold) {
	return ((value < 0) && (-value < threshold))? 0 : value;
}


} //protocols
} //peptide_deriver






