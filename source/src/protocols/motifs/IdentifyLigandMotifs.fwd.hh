// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IdentifyLigandMotifs.fwd.hh
/// @brief .fwd.hh fils for IdentifyLigandMotifs protocol. Protocol object reads in a pdb file(s) and outputs motifs for protein-ligand interactions in .pdb and .motifs format. App originally written by mdsmith, optimized and converted to protocol by Ari Ginsparg.

#ifndef INCLUDED_protocols_motifs_IdentifyLigandMotifs_fwd_hh
#define INCLUDED_protocols_motifs_IdentifyLigandMotifs_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace motifs {

class IdentifyLigandMotifs;
typedef utility::pointer::shared_ptr< IdentifyLigandMotifs > IdentifyLigandMotifsOP;
typedef utility::pointer::shared_ptr< IdentifyLigandMotifs const > IdentifyLigandMotifsCOP;


} // namespace motifs
} // namespace protocols

#endif
