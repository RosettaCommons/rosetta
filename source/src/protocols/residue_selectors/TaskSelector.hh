// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/residue_selectors/TaskSelector.hh
/// @brief  The TaskSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

#ifndef INCLUDED_protocols_residue_selectors_TaskSelector_HH
#define INCLUDED_protocols_residue_selectors_TaskSelector_HH

// Unit headers
#include <protocols/residue_selectors/TaskSelector.fwd.hh>

// Package headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers

namespace protocols {
namespace residue_selectors {

/// @brief The TaskSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which are located near the given selected residues in primary sequence space
class TaskSelector : public core::select::residue_selector::ResidueSelector {
public:
	TaskSelector();
	TaskSelector(
		core::pack::task::TaskFactoryOP tf,
		bool const select_designable,
		bool const select_packable,
		bool const select_fixed );

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual core::select::residue_selector::ResidueSelectorOP
	clone() const;

	virtual ~TaskSelector();

	virtual core::select::residue_selector::ResidueSubset
	apply( core::pose::Pose const & pose ) const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & );

	virtual
	std::string
	get_name() const;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & );

	static std::string class_name();

	//unit-specific
	void set_task_factory( core::pack::task::TaskFactoryOP tf );
	void set_select_designable( bool const sel_designable );
	void set_select_packable( bool const sel_packable );
	void set_select_fixed( bool const sel_fixed );

private: // data members
	core::pack::task::TaskFactoryOP tf_;
	bool select_designable_;
	bool select_packable_;
	bool select_fixed_;
};

} //namespace residue_selectors
} //namespace protocols


#endif
