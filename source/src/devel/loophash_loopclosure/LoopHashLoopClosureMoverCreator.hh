// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/loophash_loopclosure/LoopHashLoopClosureMoverCreator.hh
/// @brief This class will create instances of protocols::moves::Mover LoopHashLoopClosureMoverCreator for the protocols::moves::MoverFactory
/// @author Sachko Honda (honda@apl.washington.edu)

#ifndef INCLUDED_devel_loophash_loopclosure_LoopHashLoopClosureMoverCreator_HH
#define INCLUDED_devel_loophash_loopclosure_LoopHashLoopClosureMoverCreator_HH

#include <protocols/moves/MoverCreator.hh>

namespace devel {
namespace loophash_loopclosure {

class LoopHashLoopClosureMoverCreator : public protocols::moves::MoverCreator {
public:
	// XRW TEMP  LoopHashLoopClosureMoverCreator();
	// XRW TEMP  ~LoopHashLoopClosureMoverCreator() override;
	// XRW TEMP  protocols::moves::MoverOP create_mover() const override;
	// XRW TEMP  std::string keyname() const override;
	// XRW TEMP  static std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} // loophash_loopclosure
} // devel

#endif
