// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/DeleteChainsMover.hh
/// @brief Delete a chain from a pose
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_simple_moves_DeleteChainsMover_hh
#define INCLUDED_protocols_simple_moves_DeleteChainsMover_hh

#include <protocols/simple_moves/DeleteChainsMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.hh>

// C++ Headers
namespace protocols {
namespace simple_moves {

/// @brief Delete a chain from a Pose, or a set of chains.
/// @details
/// 1) Deletes Each chain using conformation.delete_residue_range_slow
///      2) Calls detect bonds, pseudobonds, and disulfides
///      3) Makes sure PDBInfo is still intact
///      4) Sets a new reasonable FoldTree
///
class DeleteChainsMover : public moves::Mover {
public:

	/// @brief Default Constructor
	DeleteChainsMover();

	/// @brief Constructor with chain
	DeleteChainsMover( std::string const & chains );

	void
	apply( core::pose::Pose & pose ) override;

	//////////////// Accessors + Mutators //////////////////////////////////////

public:

	/// @brief Set the chain
	///  Will figure out the chainIDs at apply time.
	void
	set_chains( std::string const & chains );

	/// @brief Get the chain
	std::string
	chains() const;

	/// @brief Set the class to detect bonds after full deletion of chains?
	///    Default True.
	///    Override for speed and if you know there are no bonds between chains other than possible disulfides.
	void set_detect_bonds(bool detect_bonds);
	bool detect_bonds() const;

	/// @brief Set the class to detect psuedobonds after full deletion of chains?
	///   Default True.
	///   Override for speed and if you know there are no pseudobonds between chains other than possible disulfides.
	void set_detect_pseudobonds(bool detect_pseudobonds);
	bool detect_pseudobonds() const;


	//////////////// Mover Methods /////////////////////////////////////////////

public:

	moves::MoverOP
	clone() const override;

	moves::MoverOP
	fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;


	void set_defaults();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:


	std::string chains_;

	bool detect_bonds_;
	bool detect_pseudobonds_;

};


} // simple_moves
} // protocols

#endif


