// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/OptH.hh
/// @brief  run optH
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pack_task_operation_OptH_hh
#define INCLUDED_core_pack_task_operation_OptH_hh

// unit headers
#include <core/pack/task/operation/OptH.fwd.hh>

// type headers
#include <core/types.hh>

// project headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pose/Pose.fwd.hh>

// utility headers


#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {


/// @brief run optH; disallowed positions may be set to prevent optimization for those residues
class OptH : public core::pack::task::operation::TaskOperation {


private: // typedefs


	typedef core::pack::task::operation::TaskOperation Super;


public: // typedefs


	typedef core::Size Size;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;


public: // construct/destruct


	/// @brief default constructor
	OptH();


	/// @brief copy constructor
	OptH( OptH const & rval );


	/// @brief default destructor
	virtual ~OptH();


public: // virtual constructors


	/// @brief clone this object
	virtual TaskOperationOP clone() const;


public: // methods


	/// @brief apply operations to PackerTask
	virtual void apply( Pose const & , PackerTask & task ) const;


	/// @brief prevent a position from being optimized
	void disallow_resid( Size const resid );


	/// @brief init flags from command line? (default true)
	void init_from_comand_line( bool const flag );


	/// @brief include current sidechain rotamers? (default true)
	void include_current( bool const flag );


	/// @brief allow sidechain flips of HNQ? (default false)
	void flip_HNQ( bool const flag );


	/// @brief use multicool annealer? (default false)
	void use_multicool_annealer( bool const flag );


private: // data


	/// @brief prevent these positions from being optimized
	utility::vector1< Size > disallowed_resids_;

	/* the following settings are the same as those in in core::pack::optimizeH() */

	/// @brief init flags from command line? (default true)
	bool init_from_command_line_;


	/// @brief include current sidechain rotamers? (default true)
	bool include_current_;


	/// @brief allow sidechain flips of HNQ? (default false)
	bool flip_HNQ_;


	/// @brief use multicool annealer? (default false)
	bool use_multicool_annealer_;


};


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core


#endif /* INCLUDED_core_pack_task_operation_OptH_HH */
