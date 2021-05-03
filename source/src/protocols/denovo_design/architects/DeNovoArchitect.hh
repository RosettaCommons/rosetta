// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/denovo_design/architects/StructureArchitect.hh
/// @brief Designs topologies
/// @author Tom Linsky (tlinsky@uw.edu)
/// @note   This is interface: it has no fields, and only
///         pure virtual methods.  No further constructors should
///         be defined.

#ifndef INCLUDED_protocols_denovo_design_architects_DeNovoArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_DeNovoArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/DeNovoArchitect.fwd.hh>
#include <protocols/denovo_design/architects/StructureArchitect.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// Boost headers
#include <boost/lexical_cast.hpp>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace architects {

/// @brief for planning ideal pieces of structures
/// @details Derived classes must still implment their own parse_tag() and type() functions.
///          This base handles the apply() virtual, though. In parse_tag(), derived classes
///          MUST be sure to set the motif list.
class DeNovoArchitect : public StructureArchitect {
public:
	typedef components::StructureDataOP StructureDataOP;
public:
	DeNovoArchitect( std::string const & id );

	~DeNovoArchitect() override;

	// pure virtual API
public:
	std::string
	type() const override = 0;

	virtual DeNovoArchitectOP
	clone() const = 0;

	virtual components::StructureDataOP
	design( core::pose::Pose const & pose, core::Real & random ) const = 0;
	/*
	static std:string
	class_name();

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	*/
	static void
	add_common_denovo_architect_attributes( utility::tag::AttributeList & attlist );
protected:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override = 0;

public:
	// constants
	static std::string const
		DATA_MAP_NAME;

public:
	components::StructureDataOP
	apply( core::pose::Pose const & pose ) const;

}; // DeNovoArchitect

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_DeNovoArchitect_hh
