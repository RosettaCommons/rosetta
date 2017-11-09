// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/legacy_sewing/sampling/LegacyAssemblyMoverCreator.hh
///
/// @brief
/// @author Tim Jacobs

#ifdef NOT_IN_SCONS_DEPRECATED

#ifndef INCLUDED_protocols_legacy_sewing_sampling_LegacyAssemblyMoverCreator_hh
#define INCLUDED_protocols_legacy_sewing_sampling_LegacyAssemblyMoverCreator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
namespace protocols {
namespace legacy_sewing  {

class LegacyAssemblyMoverCreator : public protocols::moves::MoverCreator
{
public:
	virtual protocols::moves::MoverOP create_mover() const = 0;
	virtual std::string keyname() const;
  virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
	//static  std::string mover_name();
};

}
}

#endif

#endif
