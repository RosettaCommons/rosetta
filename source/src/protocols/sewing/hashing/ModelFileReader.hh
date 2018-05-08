// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/hashing/ModelFileReader.hh
/// @brief a reader of SEWING model files
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_hashing_ModelFileReader_hh
#define INCLUDED_protocols_sewing_hashing_ModelFileReader_hh

#include <protocols/sewing/hashing/ModelFileReader.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/sewing/hashing/hasher_data.hh>
#include <core/conformation/Atom.hh>
#include <utility/vector1.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>

namespace protocols {
namespace sewing {
namespace hashing {

///@brief a reader of SEWING model files
class ModelFileReader : public utility::pointer::ReferenceCount {

public:

	ModelFileReader();
	ModelFileReader(ModelFileReader const & src);

	virtual ~ModelFileReader();

	ModelFileReaderOP
	clone() const;

	static SegmentVectorOP
	read_model_file(std::string filename);

private:

};


} //protocols
} //sewing
} //hashing



#endif //INCLUDED_protocols_sewing_hashing_ModelFileReader_hh





