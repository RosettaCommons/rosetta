// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/hashing/ModelFileReader.fwd.hh
/// @brief a reader of SEWING model files
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_hashing_ModelFileReader_fwd_hh
#define INCLUDED_protocols_sewing_hashing_ModelFileReader_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace sewing {
namespace hashing {

class ModelFileReader;

typedef utility::pointer::shared_ptr< ModelFileReader > ModelFileReaderOP;
typedef utility::pointer::shared_ptr< ModelFileReader const > ModelFileReaderCOP;



} //protocols
} //sewing
} //hashing


#endif //INCLUDED_protocols_sewing_hashing_ModelFileReader_fwd_hh





