// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueSelector.hh
/// @brief  The ResidueSelector class identifies a subset of residues from a Pose
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_ResidueSelector_HH
#define INCLUDED_core_select_residue_selector_ResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/citation_manager/CitationCollection.fwd.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

// Utility Headers
#include <utility/VirtualBase.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

class ResidueSelector : public utility::VirtualBase {
public:

	/// @brief Constructor.
	///
	ResidueSelector();

	/// @brief Destructor.
	///
	~ResidueSelector() override;

	/// @brief Clone operator.
	/// @details All ResidueSelectors must implement a clone() operator.  This must create a copy of the object and
	/// return a ResidueSelectorOP to the original object.
	virtual ResidueSelectorOP clone() const = 0;

	/// @brief Return a ResidueSubset indicating a selection of Residues from the
	/// input Pose; the ResidueSubset is an array of booleans where a value of
	/// "true" for position i indicates that residue i is a part of the selected
	/// subset -- and a value of "false" would indicate that it is not.
	virtual
	ResidueSubset
	apply(
		core::pose::Pose const & pose
	) const = 0;

	/// @brief Initialize any data members of this instance from an input tag
	/// and a DataMap object
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	);

	virtual
	std::string
	get_name() const = 0;

public: //Functions needed for the citation manager

	/// @brief Does this residue selector provide information about how to cite it?
	/// @details Defaults to false.  Derived classes may override this to provide citation info.  If set to
	/// true, the provide_citation_info() override should also be provided.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	virtual bool residue_selector_provides_citation_info() const;

	/// @brief Provide the citation.
	/// @returns A vector of citation collections.  This allows the residue selector to provide citations for
	/// itself and for any modules that it invokes.
	/// @details The default implementation of this function provides an empty vector.  It may be
	/// overriden by residue selectors wishing to provide citation information.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	virtual utility::vector1< basic::citation_manager::CitationCollectionCOP > provide_citation_info() const;

	/// @brief Does this residue selector indicate that it is unpublished (and, by extension, that the author should be
	/// included in publications resulting from it)?
	/// @details Defaults to false.  Derived classes may override this to provide authorship info.  If set to
	/// true, the provide_authorship_info_for_unpublished() override should also be provided.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	virtual bool residue_selector_is_unpublished() const;

	/// @brief Provide a list of authors and their e-mail addresses, as strings.
	/// @returns A list of pairs of (author, e-mail address).  Empty list if not unpublished.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	virtual utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > provide_authorship_info_for_unpublished() const;

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
CEREAL_FORCE_DYNAMIC_INIT( core_select_residue_selector_ResidueSelector )
#endif // SERIALIZATION


#endif
