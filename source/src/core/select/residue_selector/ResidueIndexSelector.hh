// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueIndexSelector.hh
/// @brief  The ResidueIndexSelector selects residues using a string containing pose indices
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

#ifndef INCLUDED_core_select_residue_selector_ResidueIndexSelector_HH
#define INCLUDED_core_select_residue_selector_ResidueIndexSelector_HH

// Unit headers
#include <core/select/residue_selector/ResidueIndexSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief The ResidueIndexSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which match the given residue index. The index is read as comma-separated
/// list of either Rosetta indices (e.g. 10) or PDB numbers (e.g. 10A, residue 10 of chain A). Detection
/// and mapping from PDB to Rosetta residue numbers is done internally.
class ResidueIndexSelector : public ResidueSelector {
public:
	// derived from base class
	ResidueIndexSelector();

	ResidueIndexSelector( std::string const & index_str );

	/// @brief Convenience constructor for a single residue index
	ResidueIndexSelector( core::Size index_in );

	/// @brief Convenience constructor for a vector of indexes
	ResidueIndexSelector( utility::vector1< core::Size > const & index_in );

	/// @brief Copy constructor
	///
	ResidueIndexSelector( ResidueIndexSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	virtual ~ResidueIndexSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//unit-specific
	/**
	* @brief sets the string by which residues are selected
	*/
	void set_index( std::string const & index_str );

	/// @brief Append an additional index (in Rosetta numbering) to the list of indices.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void append_index( core::Size index_in );

	/// @brief Append additional indexes (in Rosetta numbering) to the list of indices.
	void append_index( utility::vector1< core::Size > const & index_in );

	/// @brief Is this selector set to throw an error if an out-of-range index is selected (e.g. residue
	/// 56 of a 55-residue pose)?
	inline bool error_on_out_of_bounds_index() const { return error_on_out_of_bounds_index_; }

	/// @brief Set whether this selector set to throws an error if an out-of-range index is selected (e.g. residue
	/// 56 of a 55-residue pose).
	inline void set_error_on_out_of_bounds_index( bool const setting ) { error_on_out_of_bounds_index_ = setting; }

private: // data members
	std::string index_str_;

	/// @brief If false, then there is no error if an index that is not in the pose is
	/// selected.  True by default, which means that you get an error if you try to
	/// select an index that's not in the pose.
	bool error_on_out_of_bounds_index_;
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
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_ResidueIndexSelector )
#endif // SERIALIZATION


#endif
