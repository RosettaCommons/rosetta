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

#ifndef INCLUDED_protocols_denovo_design_architects_StructureArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_StructureArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/StructureArchitect.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.fwd.hh>

// Boost headers
#include <boost/lexical_cast.hpp>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace architects {

typedef utility::vector1< core::Size > Lengths;
typedef utility::vector1< std::string > Abego;

struct SecStructInfo {
	SecStructInfo():
		ss(), abego() {};
	std::string ss;
	Abego abego;
};

/// @brief Designs topologies
class StructureArchitect : public utility::pointer::ReferenceCount  {
public: // Creation

	StructureArchitect( std::string const & id );

	virtual
	~StructureArchitect();

public:
	/// @brief simply returns the name of this type of architect
	virtual std::string
	type() const = 0;

protected:
	/// @brief Configuration by XML
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) = 0;

public:
	void
	parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	std::string const &
	id() const;

	void
	set_id( std::string const & new_id );

private:
	/// @brief Prevent direct instantiation: No other constructors allowed.
	StructureArchitect();

private:
	/// @brief name of this architect
	std::string id_;
}; // StructureArchitect

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_StructureArchitect_hh
