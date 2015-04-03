// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/motifs/MotifDnaPackerCreator.hh
/// @brief This class will create instances of Mover MotifDnaPacker for the MoverFactory
/// @author sthyme (sthyme@gmail.com) copying Andrew Leaver-Fay via code_writer.py (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_motifs_MotifDnaPackerCreator_hh
#define INCLUDED_protocols_motifs_MotifDnaPackerCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace motifs {

class MotifDnaPackerCreator : public protocols::moves::MoverCreator {
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
};

}
}

#endif

