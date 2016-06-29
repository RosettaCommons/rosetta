// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/movers/SetResidueAliasMover.hh
/// @brief Sets a residue alias in the StructureData
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_denovo_design_movers_SetResidueAliasMover_hh
#define INCLUDED_protocols_denovo_design_movers_SetResidueAliasMover_hh

// Unit headers
#include <protocols/denovo_design/movers/SetResidueAliasMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

///@brief Sets a residue alias in the StructureData
class SetResidueAliasMover : public protocols::moves::Mover {

public:

	SetResidueAliasMover();

	// copy constructor (not needed unless you need deep copies)
	//SetResidueAliasMover( SetResidueAliasMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~SetResidueAliasMover();

	static std::string
	class_name();

public:
	// mover virtual API
	virtual void
	apply( core::pose::Pose & pose );

	virtual void
	show( std::ostream & output = std::cout ) const;

	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//SetResidueAliasMover & operator=( SetResidueAliasMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

private:
	std::string alias_name_;
	std::string segment_name_;
	core::Size resid_;

};

std::ostream &
operator<<( std::ostream & os, SetResidueAliasMover const & mover );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_SetResidueAliasMover_hh
