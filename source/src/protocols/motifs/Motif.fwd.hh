// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Motif.fwd.hh
/// @brief forward declaration for interaction motif between residues
/// @author havranek, sthyme (sthyme@gmail.com)

#ifndef INCLUDED_protocols_motifs_Motif_fwd_hh
#define INCLUDED_protocols_motifs_Motif_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace motifs {

class Motif;
typedef utility::pointer::shared_ptr< Motif > MotifOP;
typedef utility::pointer::shared_ptr< Motif const > MotifCOP;
typedef utility::vector1< MotifOP > MotifOPs;
typedef utility::vector1< MotifCOP > MotifCOPs;

} // namespace motifs
} // namespace protocols

#endif // INCLUDED_protocols_motifs_Motif_fwd
