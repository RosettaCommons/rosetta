// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/connection/ExtendChain.hh
/// @brief The ExtendChain Protocol
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_connection_ExtendChain_hh
#define INCLUDED_protocols_denovo_design_connection_ExtendChain_hh

// Unit headers
#include <protocols/denovo_design/connection/ExtendChain.fwd.hh>
#include <protocols/denovo_design/connection/BridgeChains.hh>

// Package headers

// Protocol headers

// Core headers

// C++ headers

namespace protocols {
namespace denovo_design {
namespace connection {

class ExtendChain : public BridgeChains {
public:
	/// @brief default constructor
	ExtendChain();

	/// @brief virtual constructor to allow derivation
	virtual ~ExtendChain();

	/// @brief Parses the ExtendChainTags
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	/// @brief Mover virtuals
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual protocols::moves::MoverOP clone() const;

	virtual void process_permutation( components::StructureData & perm ) const;

	/// @brief Performs chain connections and assumes all connection setup has been performed
	virtual void apply_connection( components::StructureData & perm ) const;

	/// @brief this mover creates a polymer bond, so this is true
	virtual bool polymer_connection() const { return true; }

	/// @brief configures based on a permutation and saves info for building
	/// @throw EXCN_Setup if no valid connection endpoints are found
	virtual void
	setup_permutation( components::StructureData & perm ) const;

};

} // connection
} // denovo_design
} // protocols

#endif
