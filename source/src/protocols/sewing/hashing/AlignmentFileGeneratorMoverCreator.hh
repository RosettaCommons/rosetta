// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/AlignmentFileGeneratorMoverCreator.hh
/// @brief Given a model file, edge file ,and one or more input structures, generates alignment files for use with AppendAssemblyMover.
/// @author guffysl (guffy@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_hashing_AlignmentFileGeneratorMoverCreator_hh
#define INCLUDED_protocols_sewing_hashing_AlignmentFileGeneratorMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace sewing {
namespace hashing {

class AlignmentFileGeneratorMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

} //hashing
} //sewing
} //protocols
#endif //INCLUDED_protocols/sewing/hashing_AlignmentFileGeneratorMoverCreator_fwd_hh



