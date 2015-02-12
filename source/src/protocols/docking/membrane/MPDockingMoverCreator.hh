// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/MPDockingMoverCreator.hh
/// @brief      Dock two membrane proteins
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @author     Rebecca Alford (rfalford12@gmail.com) - Added RosettaScripts hook
/// @note       Last Modified (2/9/15)

#ifndef INCLUDED_protocols_docking_membrane_MPDockingMoverCreator_hh
#define INCLUDED_protocols_docking_membrane_MPDockingMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace docking {
namespace membrane {
    
/// @brief Mover Creator
class MPDockingMoverCreator : public protocols::moves::MoverCreator {
    
public:
    
    virtual protocols::moves::MoverOP create_mover() const;
    virtual std::string keyname() const;
    static std::string mover_name();
    
};

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPDockingMoverCreator_hh
