// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SliceResidueSelector.hh
/// @brief  The SliceResidueSelector allows slicing of the returned values of other residue selectors
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_select_residue_selector_SliceResidueSelector_HH
#define INCLUDED_core_select_residue_selector_SliceResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/SliceResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

namespace slice_enums {
enum SliceMode {
	SPARSE,
	CONTIGUOUS
};

enum OutOfBoundsMode {
	ERROR,
	WARN,
	IGNORE
};
}

/// @brief The SliceResidueSelector allows slicing of the returned values of other residue selectors
class SliceResidueSelector : public ResidueSelector {
public:
	// derived from base class
	SliceResidueSelector();

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	ResidueSelectorOP clone() const override;

	SliceResidueSelector(
		ResidueSelectorCOP const & selector,
		slice_enums::SliceMode slice_mode,
		int from,
		int to,
		slice_enums::OutOfBoundsMode oob_mode = slice_enums::ERROR
	);

	SliceResidueSelector(
		ResidueSelectorCOP selector,
		slice_enums::SliceMode slice_mode,
		utility::vector1< int > const & indices,
		slice_enums::OutOfBoundsMode oob_mode = slice_enums::ERROR
	);

	ResidueSubset apply( core::pose::Pose const & pose ) const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	std::string
	get_name() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void set_residue_selector( ResidueSelectorCOP selector );

	void set_slice_mode( slice_enums::SliceMode slice_mode );

	void set_from( int from );

	void set_to( int to );

	void set_indices( utility::vector1< int > const & indices );

	void set_out_of_bounds_behavior( slice_enums::OutOfBoundsMode oob_mode );

private:

	int
	wrap_index( int index, int elements, utility::vector1< std::string > & error_messages ) const;

	void
	show_selection_logic(
		ResidueSubset const & initial,
		ResidueSubset const & final
	) const;

public: //Functions needed for the citation manager

	/// @brief Provide the citation.
	/// @returns A vector of citation collections.  This allows the residue selector to provide citations for
	/// itself and for any modules that it invokes.
	/// @details This residue selector has no citation.  It may provide citations for the residue selector that
	/// it calls, though.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	utility::vector1< basic::citation_manager::CitationCollectionCOP > provide_citation_info() const override;

	/// @brief Does this residue selector indicate that it is unpublished (and, by extension, that the author should be
	/// included in publications resulting from it)?
	/// @details Returns true (this residue selector is unpublished).
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	bool residue_selector_is_unpublished() const override;

	/// @brief Provide a list of authors and their e-mail addresses, as strings.
	/// @returns This residue selector was created by Brian Coventry.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > provide_authorship_info_for_unpublished() const override;

private: // data members
	ResidueSelectorCOP selector_;

	slice_enums::SliceMode slice_mode_;
	int from_;
	int to_;
	utility::vector1< int > indices_;
	slice_enums::OutOfBoundsMode oob_mode_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_SliceResidueSelector )
#endif // SERIALIZATION


#endif
