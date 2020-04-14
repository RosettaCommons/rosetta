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
	ConvertRealToVirtualMover( ConvertRealToVirtualMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	~ConvertRealToVirtualMover() override;

	static std::string
	mover_name();

public:

	/// @brief Set the residue selector.  If you need a particular set of residues,
	///  use the ReturnSubsetResidueSelector.
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );


	// mover virtual API
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	std::string
	get_name() const override;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//ConvertRealToVirtualMover & operator=( ConvertRealToVirtualMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public:  // Citation Management

	/// @brief Does this mover provide information about how to cite it?
	bool
	mover_provides_citation_info() const override;

	/// @brief Provide a list of authors and their e-mail addresses, as strings.
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
	provide_authorship_info_for_unpublished() const override;

	/// @brief Although this mover has no citation info, the residue selector that it calls might have some.
	utility::vector1< basic::citation_manager::CitationCollectionCOP > provide_citation_info() const override;

private:
	/// @brief The ResidueSelector to use.
	core::select::residue_selector::ResidueSelectorCOP selector_;
}; // ConvertRealToVirtualMover class

std::ostream &
operator<<( std::ostream & os, ConvertRealToVirtualMover const & mover );

} //protocols
} //simple_moves

#endif //protocols/simple_moves_ConvertRealToVirtualMover_hh
