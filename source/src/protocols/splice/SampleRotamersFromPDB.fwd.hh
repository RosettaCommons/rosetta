// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/task_operations/SampleRotamersFromPDB.hh
/// @brief  rotamer set operation forward declaration

#ifndef INCLUDED_protocols_splice_SampleRotamersFromPDB_fwd_hh
#define INCLUDED_protocols_splice_SampleRotamersFromPDB_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <core/types.hh>


namespace protocols {
namespace splice {


class SampleRotamersFromPDB_RotamerSetOperation;
class SampleRotamersFromPDB;

//Rotamer object
struct ROT {
	std::string AA;
	utility::vector1< core::Real> chi_vec;
} ;
typedef utility::vector1 < utility::vector1< ROT > >rot_matrix;

typedef utility::pointer::shared_ptr< SampleRotamersFromPDB_RotamerSetOperation > SampleRotamersFromPDB_RotamerSetOperationOP;
typedef utility::pointer::shared_ptr< SampleRotamersFromPDB_RotamerSetOperation const > SampleRotamersFromPDB_RotamerSetOperationCOP;

typedef utility::pointer::shared_ptr< SampleRotamersFromPDB > SampleRotamersFromPDBOP;
typedef utility::pointer::shared_ptr< SampleRotamersFromPDB const > SampleRotamersFromPDBCOP;

class RotLibdb;
typedef utility::pointer::shared_ptr< RotLibdb > RotLibdbOP;
typedef utility::pointer::shared_ptr< RotLibdb const > RotLibdbCOP;
}
}



#endif
