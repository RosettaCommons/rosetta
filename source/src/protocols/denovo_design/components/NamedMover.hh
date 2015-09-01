// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/NamedMover.hh
/// @brief NamedMover functions for building structures from components
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_NamedMover_hh
#define INCLUDED_protocols_denovo_design_components_NamedMover_hh

// Unit headers
#include <protocols/denovo_design/components/NamedMover.fwd.hh>

// Project headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>

// Core headers

// Basic/Numeric/Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ Headers

namespace protocols {
namespace denovo_design {
namespace components {

/// @brief manages information about segments of residues
class NamedMover : public protocols::moves::Mover {
public:
	NamedMover();
	NamedMover( std::string const & id, std::string const & parent_id );

	virtual ~NamedMover() {}

	/// @brief setup the parameters via an xml tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	// virtual interface functions
public:
	/// @brief performs setup and applies loop building
	/// @details Default steps to apply():
	/// 1. Pulls StructureData from the pose
	/// 2. setup_permutation() stores data about this connection which is not
	///    static for every apply() call
	/// 3. apply_permutation() uses the StructureData object to build the loop
	/// 4. check() checks the built structure
	virtual void apply( core::pose::Pose & pose );

	/// @brief uses the StructureData to set up, builds structure, and updates permutation
	/// @details Steps should be:
	///  1. Use StructureData to collect information about the structure to be built
	///  2. Build the pose as instructed by the StructureData object
	///  3. Store the pose in the StructureData object
	virtual void apply_permutation( StructureData & perm ) = 0;

	/// @brief Choose and store build options in the StructureData object
	/// @details You can store build-specified information (loop length, desired
	/// abego, etc.) in the StructureData object which will be used later when
	/// apply_permutation() is called.
	/// This function SHOULD NOT engage in any pose modification
	virtual protocols::moves::MoverStatus
	setup_permutation( StructureData & perm ) const = 0;

	/// @brief checks an unbuilt permutation modified by setup_permutation to see if it is acceptable for building
	virtual bool check_permutation( StructureData const & perm ) const = 0;

	/// @brief checks the built StructureData vs. the desired params.
	/// @details should return true if it matches, false otherwise
	virtual bool check( StructureData const & perm ) const = 0;

public:
	// member functions

	/// @brief returns an identifier for this mover
	std::string const & id() const;

	/// @brief returns the name of the parent mover
	std::string const & parent_id() const;

	/// @brief sets an identifier for this mover
	void set_id( std::string const & idval );

	/// @brief sets an identifier for the parent of this mover, for nested moves
	void set_parent_id( std::string const & parent );

protected:
	/// @brief adds prefix if necessary, returns result
	std::string add_parent_prefix( std::string const & s ) const;

private:
	std::string id_;
	std::string parent_id_;
};

} // components
} // denovo_design
} // protocols

#endif
