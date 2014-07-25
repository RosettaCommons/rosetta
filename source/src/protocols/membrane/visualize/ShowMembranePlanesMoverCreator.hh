// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/visualize/ShowMembranePlanesMoverCreator.hh
///
/// @brief 		Add Anchor Residues for Membrane Planes
/// @details    Add 6 virtual membrane residues to the pose as an additional
///				chain to the protein. Anchor residues are attached by jump to the membrane center
///				and are pointing downstream to allow movement of the planes. This plug in
///				works directly with the PyMol Mover
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @author		Evan H. Baugh (evan.baugh.social@gmail.com)
/// @note       Last Modified: 6/27/14

#ifndef INCLUDED_protocols_membrane_visualize_ShowMembranePlanesMoverCreator_hh
#define INCLUDED_protocols_membrane_visualize_ShowMembranePlanesMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {
namespace visualize {

/// @brief Mover Creator
class ShowMembranePlanesMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();

};

} // visualize
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_visualize_ShowMembranePlanesMoverCreator_hh
