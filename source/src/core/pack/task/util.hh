// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/util.hh
/// @brief Utility functions for main task classes.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_pack_task_util_hh
#define INCLUDED_core_pack_task_util_hh

#include <core/pack/task/ResfileReader.fwd.hh>
#include <utility/vector1.hh>

#include <map>

namespace core {
namespace pack {
namespace task {

///@brief Parse a command line that is residue agnostic and has no residue information.
///  IE "PIKAA  ST"  or [EMPTY NC xxx NC xxx]
utility::vector1< ResfileCommandOP >
parse_res_agnostic_commands(
	std::string const & line,
	std::map< std::string, ResfileCommandOP > const & command_map);


///@brief Take a Design Sequence Motif (like a simple version of a resfile), and parse into a list of ResfileCommands.
/// These commands can then be used to do whatever you want with, including non-cannonicals.
///
/// This is used in the CreateSequenceMotifMover and the SequenceMotifTaskOperation interfaces.
///
///@details
///
///  This is slightly similar to a regex, but not quite. We are not matching a sequence, we are designing in a motif regardless of the current sequence, anywhere in a protein.
///
///   - Each letter corresponds to a position. Using [ ] indicates a more complicated expression for that position.
///   - An X indicates it can be anything, and that we are designing here.
///   - An AA Letter, like V, indicates that that position will be designed to a V.
///   - A - charactor indicates that that position stays with whatever it is currently.  We essentially skip this position.
///   - An expression like: [^PAV] indicates that we will design anything except Proline, Alanine, and Valine
///   - An expression like: [NTS] indicates that that position can be Asparigine, Threonine, or Serine and only of these will be enabled during the design.
///   - RESFILE commands are accepted as well.
///      These require a % Charactor to start the whole line off.
///      For example [%POLAR] is totally cool, as is [%PIKAA ST] or [%EMPTY NC R2 NC T6 NC OP5] for Noncannonicals.
///
///       (Non-cannonical, DNA, RNA, etc. should all work here!!)
///
/// EXAMPLE:
///  Glycosylation N-Linked motif design: N[^P][ST]
///
utility::vector1< utility::vector1< ResfileCommandOP > >
get_resfile_commands( std::string const & motif );



} //core
} //pack
} //task


#endif //core/pack/task_util_hh

