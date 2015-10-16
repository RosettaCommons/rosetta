// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_util_io_HH
#define INCLUDED_protocols_sewing_util_io_HH

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>

//Utility
#include <utility/vector1.hh>

//Devel
#include <protocols/sewing/conformation/Model.fwd.hh>
#include <protocols/sewing/conformation/Assembly.fwd.hh>
#include <protocols/sewing/sampling/SewGraph.hh>

//C++
#include <map>
#include <string>

namespace protocols {
namespace sewing  {

typedef std::map<core::Size, utility::vector1<std::pair<bool, core::conformation::ResidueOP> > > NativeRotamersMap;

//Create a file that saves native rotamers. Used later during design
void
write_native_residue_file(
	NativeRotamersMap native_residue_map,
	std::string filename
);

NativeRotamersMap
read_native_residue_file(
	std::string filename
);

///@brief Take the given StructureScores and save them to disk
void
write_hashing_scores_to_file(
	ScoreResults const & scores,
	std::string filename
);

utility::vector1<BasisPair>
scores_to_alignments(
	ScoreResults const & scores
);

///@brief Take the given StructureScores and save them to disk
utility::vector1<BasisPair>
read_hashing_scores_from_file(
	std::string filename
);

std::string
serialize_graph_json(
	SewGraphCOP graph,
	core::Size max_nodes
);

} //sewing namespace
} //protocols namespace

#endif
