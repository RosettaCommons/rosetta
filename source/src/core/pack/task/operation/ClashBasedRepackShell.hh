// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ClashBasedRepackShell.hh
/// @brief  ClashBasedRepackShell header file.
/// @author Kale Kundert (kale@thekunderts.net)

#ifndef INCLUDED_core_pack_task_operation_ClashBasedRepackShell_hh
#define INCLUDED_core_pack_task_operation_ClashBasedRepackShell_hh

// Unit headers
#include <core/pack/task/operation/ClashBasedRepackShell.fwd.hh>

// Core headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/residue_selector/ClashBasedShellSelector.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace task {
namespace operation {

/// @brief Create a shell of residues that can repack around a smaller group of
/// residues being repacked or designed.
///
/// More specifically, this task operation disables repacking for (i.e.
/// freezes) any position that isn't part of the aforementioned group or its
/// shell.  The most common use-case is to define a set of positions that
/// you're interested in designing (i.e. with a resfile), and limit the rest of
/// the protein to repacking (e.g. NATAA).  Then you can apply the
/// ClashBasedRepackShell to limit the repacking to only the positions that
/// really interact with the design positions.  This keeps the simulations
/// efficient and prevents the final scores from being affected by random
/// changes in packing far away from the region being designed.
///
/// This task operation is a very thin wrapper around ClashBasedShellSelector.
/// You can access and configure the selector in question via the selector()
/// method.
class ClashBasedRepackShell : public TaskOperation {

public:

	/// @brief Default constructor.
	ClashBasedRepackShell();

	/// @brief Default copy constructor.
	ClashBasedRepackShell( ClashBasedRepackShell const & ) = default;

	/// @brief Default destructor.
	virtual ~ClashBasedRepackShell() = default;

	/// @brief Return a shallow copy of this object.
	virtual TaskOperationOP clone() const;

	/// @brief apply operations to PackerTask
	virtual void apply( Pose const & pose, PackerTask & task ) const;

	/// @brief Initialize from an XML tag.
	void parse_tag( utility::tag::TagCOP, basic::datacache::DataMap & );

	/// @brief Define the expected XML options.
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Return "ClashBasedRepackShell".
	static std::string keyname();

	/// @brief Get the ResidueSelector used to define the repack shell.
	core::pack::task::residue_selector::ClashBasedShellSelectorOP selector() const;

	/// @brief Set the ResidueSelector used to define the repack shell.
	void selector(core::pack::task::residue_selector::ClashBasedShellSelectorOP);

private:
	core::pack::task::residue_selector::ClashBasedShellSelectorOP shell_selector_ = nullptr;

};


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core

#endif
