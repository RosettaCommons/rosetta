// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/qsar/qsarMap.fwd.hh
/// @author Sam DeLuca


#ifndef INCLUDED_protocols_qsar_qsarMap_fwd_hh
#define INCLUDED_protocols_qsar_qsarMap_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace qsar {




class qsarPoint;

class qsarMap;

typedef utility::pointer::shared_ptr<qsarMap> qsarMapOP;
typedef utility::pointer::shared_ptr<qsarMap const> qsarMapCOP;

typedef utility::pointer::shared_ptr<qsarPoint> qsarPointOP;
typedef utility::pointer::shared_ptr<qsarPoint const> qsarPointCOP;

}
}


#endif
