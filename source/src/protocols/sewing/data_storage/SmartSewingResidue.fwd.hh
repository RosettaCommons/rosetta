// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/data_storage/SmartSewingResidue.fwd.hh
/// @brief a minimal container for SEWING residues
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_data_storage_SmartSewingResidue_fwd_hh
#define INCLUDED_protocols_sewing_data_storage_SmartSewingResidue_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace sewing {
namespace data_storage {

class SmartSewingResidue;

typedef utility::pointer::shared_ptr< SmartSewingResidue > SmartSewingResidueOP;
typedef utility::pointer::shared_ptr< SmartSewingResidue const > SmartSewingResidueCOP;

enum PolymericType{ polymer, ligand };

} //protocols
} //sewing
} //data_storage


#endif //INCLUDED_protocols_sewing_data_storage_SmartSewingResidue_fwd_hh





