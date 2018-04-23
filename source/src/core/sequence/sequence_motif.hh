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

/*
"  This is slightly similar to a regex, but not quite. We are not matching a sequence,\n"
"   we are designing in a motif regardless of the current sequence, anywhere in a protein.\n"
"\n"
"   - Each letter corresponds to a position. Using [ ] indicates a more complicated expression for that position.\n"
"   - An X indicates it can be anything, and that we are designing here.\n"
"   - An AA Letter, like V, indicates that that position will be designed to a V.\n"
"   - A - charactor indicates that that position stays with whatever it is currently.  We essentially skip this position.\n"
"   - An expression like: [^PAV] indicates that we will design anything except Proline, Alanine, and Valine \n"
"   - An expression like: [NTS] indicates that that position can be Asparigine, Threonine, or Serine and \n"
"      only these will be enabled during the design.\n"
"   - RESFILE commands are accepted as well. These require a % charactor in from of the whole expression.\n"
"     For example [%POLAR] would set that position to only polar design.\n"
"     This is exactly the same as a resfile line, so you can even do NC like so: \n"
"      [%EMPTY NC R2 NC T6 NC OP5]\n"
"\n"
" EXAMPLE:\n"
"  Glycosylation N-Linked motif design: N[^P][ST]\n"
"\n";
*/


///@brief Splits the sequence motif (Ex. N[^P]-[STN]A) into strings of commands for individual positions
/// Example: ["N", "^P", "-", "STN","A",] for the three positions specified
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
