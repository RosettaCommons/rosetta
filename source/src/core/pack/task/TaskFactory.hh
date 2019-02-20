// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/PackerTaskFactory.hh
/// @brief  Factory class for the creation and initialization of PackerTask objects
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- added support for PackerPalettes.

#ifndef INCLUDED_core_pack_task_TaskFactory_hh
#define INCLUDED_core_pack_task_TaskFactory_hh

// Unit Headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/types.hh>

// Package Headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/palette/PackerPalette.fwd.hh>
#include <core/pack/palette/DefaultPackerPalette.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <list>

#include <utility/vector1.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION



namespace core {
namespace pack {
namespace task {

/// @brief  Factory class for the creation and initialization of PackerTask objects
class TaskFactory : public utility::pointer::ReferenceCount
{
public:
	typedef utility::pointer::ReferenceCount parent;
	typedef operation::TaskOperation TaskOperation;
	typedef operation::TaskOperationOP TaskOperationOP;
	typedef operation::TaskOperationCOP TaskOperationCOP;
	typedef std::list< TaskOperationOP > OperationList;
	typedef OperationList::const_iterator const_iterator;
	typedef core::pack::palette::PackerPalette PackerPalette;
	typedef core::pack::palette::PackerPaletteOP PackerPaletteOP;
	typedef core::pack::palette::PackerPaletteCOP PackerPaletteCOP;
	typedef core::pack::palette::DefaultPackerPalette DefaultPackerPalette;
	typedef core::pack::palette::DefaultPackerPaletteOP DefaultPackerPaletteOP;

public:

	TaskFactory();
	TaskFactory( TaskFactory const & );
	virtual ~TaskFactory();

	TaskFactory & operator = ( TaskFactory const & );

	virtual TaskFactoryOP clone() const;

	/// @brief  Non static version.
	PackerTaskOP
	create_task_and_apply_taskoperations( pose::Pose const & pose ) const;

	/// @brief Apply each of the TaskOperations in the TaskFactory list to the PackerTask to set it up.
	/// @details Must be called AFTER PackerTask initialization.  An uninitialized PackerTask
	/// that is modified will throw an error.
	void modify_task( core::pose::Pose const & pose, PackerTaskOP task ) const;

	/// @brief clones the input task, and pushes it back into the list
	void
	push_back( TaskOperationCOP );

	/// @brief Clones the input PackerPalette, and sets it as the PackerPalette for this TaskFactory.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void set_packer_palette( PackerPaletteCOP packer_palette_in );

	/// @brief Has a PackerPalette been set for this TaskFactory?
	inline bool has_packer_palette() const { return packer_palette_ != nullptr; }

	/// @brief begin iterator of operations_
	const_iterator
	begin() const;

	/// @brief end iterator of operations_
	const_iterator
	end() const;

	/// @brief Empties the list of TaskOperations and clears the PackerPalette, replacing it with a new DefaultPackerPalette.
	///
	void
	clear();

	/// @brief return the size of the operations_ list
	core::Size size() const;

public:

	/// @brief Static construction of a task
	/// @details Returns a new PackerTask with NO TaskOperations or PackerPalette applied.
	static
	PackerTaskOP
	create_packer_task( pose::Pose const & );

	/// @brief Static construction of a task
	/// @details Returns a new PackerTask with NO TaskOperations, and with a supplied
	/// PackerPalette (used directly -- not cloned on input).
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	static
	PackerTaskOP
	create_packer_task( pose::Pose const &pose, core::pack::palette::PackerPaletteCOP palette );

private:

	void copy_operations( TaskFactory const & src );

private:

	/// @brief The PackerPalette, defining the default packer behaviour in the absence
	/// of any TaskOperations.
	/// @details Added as part of the 2016 Rosetta Chemical XRW (eXtreme Rosetta Workshop), New Brunswick, NJ.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	PackerPaletteOP packer_palette_ = nullptr;


	/// @brief The list of TaskOperations to apply.
	///
	OperationList operations_;

#ifdef SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION



};

} //namespace task
} //namespace pack
} //namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_TaskFactory )
#endif // SERIALIZATION

#endif
