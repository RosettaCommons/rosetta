// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file     protocols/membrane/VisualizeEmbeddingMoverCreator.hh
/// @brief      Visualize Embedding normal and center with Virtual Residues
/// @details    Add a set of virtual residues as an additional chain to the
///    membrane pose. This tool is strictly for visualization of
///    the implicit membrane and should not be present in modeling.
///    Last Modified: 11/20/14
/// @author  JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_visualize_VisualizeEmbeddingMoverCreator_hh
#define INCLUDED_protocols_membrane_visualize_VisualizeEmbeddingMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {
namespace visualize {

/// @brief Mover Creator
class VisualizeEmbeddingMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();

};

} // visualize
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_visualize_VisualizeEmbeddingMoverCreator_hh
