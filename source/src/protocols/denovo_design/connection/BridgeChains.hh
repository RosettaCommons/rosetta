// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/connection/BridgeChains.hh
/// @brief The BridgeChains Protocol
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_connection_BridgeChains_hh
#define INCLUDED_protocols_denovo_design_connection_BridgeChains_hh

// Unit headers
#include <protocols/denovo_design/connection/BridgeChains.fwd.hh>
#include <protocols/denovo_design/connection/Connection.hh>

// Package headers
#include <protocols/denovo_design/components/Picker.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Protocol headers
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace denovo_design {
namespace connection {

class BridgeChains : public Connection {
public:
	/// @brief default constructor
	BridgeChains();

	/// @brief virtual constructor to allow derivation
	virtual ~BridgeChains();

	/// @brief Parses the BridgeChainsTags
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

	/// @brief Performs chain connections and assumes all connection setup has been performed
	virtual void apply_connection( components::StructureData & perm );

	/// @brief this mover creates a polymer bond, so this is true
	virtual bool polymer_connection() const { return true; }

	/// @brief configures based on a permutation and saves info for building
	virtual protocols::moves::MoverStatus
		setup_permutation( components::StructureData & perm ) const;

	// accessor/mutator
public:
	/// @brief sets the scorefunction
	inline void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn_val ) { scorefxn_ = scorefxn_val; }

	// member functions
public:

	/// @brief creates mover that does remodeling of loop residues
	/// @details result is a ready-to-call mover
	virtual protocols::moves::MoverOP
	create_remodel_mover(
			core::pose::Pose const & pose,
			protocols::loops::LoopsOP loops,
			bool const const_fold_tree,
			std::string const & complete_ss,
			StringVec const & complete_abego,
			core::Size const left,
			core::Size const right );

	/// @brief checks to ensure that both pieces being connected are fixed relative to one another
	bool segments_fixed( components::StructureData const & perm ) const;

	/// @brief using the motif list, find the desired abego for each position in the connection
	utility::vector1< std::string > abego_insert(
			StringVec const & complete_abego,
			std::string const & connection_abego,
			core::Size const left,
			core::Size const right,
			core::Size const end1,
			core::Size const start2 ) const;

	/// @brief using the motif list and input pose, find the desired aa sequence for each position in the connection. default="V"
	std::string aa_insert(
			core::pose::Pose const & pose,
			core::Size const connection_len,
			core::Size const left,
			core::Size const right,
			core::Size const end1,
			core::Size const start2 ) const;

	/// @brief using the motif list and input pose, find the desired secondary structure for each position in the connection
	std::string ss_insert(
			core::pose::Pose const & pose,
			std::string const & connection_ss,
			core::Size const left,
			core::Size const right,
			core::Size const end1,
			core::Size const start2 ) const;

	/// @brief builds the loop
	void build_loop( components::StructureData & perm );

	/// @brief get score function
	core::scoring::ScoreFunctionCOP scorefxn();

private: // options
	core::scoring::ScoreFunctionCOP scorefxn_;

private:   // other data
	/// @brief fragment picker/cache
	components::PickerOP frag_picker_;
};

} // connection
} // denovo_design
} // protocols

#endif
