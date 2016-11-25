// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/benchmark/MakeCanonicalHelixCreator.hh
/// @brief Creates an ideal a-helix from sequence given a range of residues
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_benchmark_MakeCanonicalHelixCreator_hh
#define INCLUDED_protocols_membrane_benchmark_MakeCanonicalHelixCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {
namespace benchmark {

class MakeCanonicalHelixCreator : public protocols::moves::MoverCreator {

public:

	// XRW TEMP  virtual protocols::moves::MoverOP create_mover() const;
	// XRW TEMP  virtual std::string keyname() const;
	// XRW TEMP  static std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} //protocols
} //membrane
} //benchmark

#endif // INCLUDED_protocols_membrane_benchmark_MakeCanonicalHelix_fwd_hh



