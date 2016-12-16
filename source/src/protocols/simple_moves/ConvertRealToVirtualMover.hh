// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/simple_moves/ConvertRealToVirtualMover.hh
/// @brief Mover for switching a residue type to all virtual
/// @author Sebastian Rämisch (raemisch@scripps.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_ConvertRealToVirtualMover_hh
#define INCLUDED_protocols_simple_moves_ConvertRealToVirtualMover_hh

// Unit headers
#include <protocols/simple_moves/ConvertRealToVirtualMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace simple_moves {

///@brief
///
///       Mover for switching a residue type to all VIRTUAL.
///       A VIRTUAL Residue is one that is not scored or output.
///
///@details
///
///      If no residue_selector is set, will work on ALL residues
///
///See also:
///
///      FAToVirtualMover
///
class ConvertRealToVirtualMover : public protocols::moves::Mover {

public:

	ConvertRealToVirtualMover();
	
	/// @brief Constructor with residue selector.  If you need a particular set of residues,
	///  use the ReturnSubsetResidueSelector.
	ConvertRealToVirtualMover( core::select::residue_selector::ResidueSelectorCOP selector );
	
	// copy constructor (not needed unless you need deep copies)
	//ConvertRealToVirtualMover( ConvertRealToVirtualMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~ConvertRealToVirtualMover();

	static std::string
	mover_name();

public:
	
	/// @brief Set the residue selector.  If you need a particular set of residues,
	///  use the ReturnSubsetResidueSelector.
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );
	
	
	// mover virtual API
	virtual void
	apply( core::pose::Pose & pose );

	virtual void
	show( std::ostream & output = std::cout ) const;

	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//ConvertRealToVirtualMover & operator=( ConvertRealToVirtualMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	
private:
	/// @brief The ResidueSelector to use.
	core::select::residue_selector::ResidueSelectorCOP selector_;
}; // ConvertRealToVirtualMover class

std::ostream &
operator<<( std::ostream & os, ConvertRealToVirtualMover const & mover );

} //protocols
} //simple_moves

#endif //protocols/simple_moves_ConvertRealToVirtualMover_hh
