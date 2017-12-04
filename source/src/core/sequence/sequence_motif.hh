// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/sequence/sequtil.hh
/// @brief Utilities for sequence motifs
/// @author Jared Adolf-Bryfogle


#ifndef INCLUDED_core_sequence_sequence_motif_hh
#define INCLUDED_core_sequence_sequence_motif_hh

// C/C++ headers

// Utility headers
#include <utility/vector1.fwd.hh>

// Project headers
#include <core/types.hh>

// Package headers
#include <core/sequence/Sequence.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

namespace core {
namespace sequence {

///@brief Splits the sequence motif (Ex. N[^P][STN]) into strings of individual positions
///  (which would result in [N+, ^P, STN])
/// (Any string with + in it signifies that it does not come from the [] notation and can be accounted for as individual commands.
///
utility::vector1<std::string>
split_sequence_motif( std::string const & motif );

core::Size
get_motif_length( std::string const & motif );

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Design Sequence Motif /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


///@brief Return the description of the syntax for a design sequence motif.
std::string
get_design_sequence_motif_syntax();


} // sequence
} // core

#endif
