// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_moves/ModifyVariantTypeMover.cc
/// @brief Modify variant type to residues
/// @author Alex Ford <fordas@uw.edu>

// Unit Headers
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>
#include <protocols/simple_moves/ModifyVariantTypeMoverCreator.hh>

// Package headers

// Project headers
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/rosetta_scripts/util.hh>

// tracer
#include <basic/Tracer.hh>

// Utility Headers

// C++ Headers
#include <iostream>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <boost/algorithm/string.hpp>
#include <utility/excn/Exceptions.hh>
#include <boost/foreach.hpp>


// ObjexxFCL Headers

namespace protocols {
namespace simple_moves {

static basic::Tracer TR( "protocols.simple_moves.ModifyVariantTypeMover" );


/// ModifyVariantTypeMover; based on the protocols::moves::Mover basis class
ModifyVariantTypeMover::ModifyVariantTypeMover() :
  protocols::moves::Mover("ModifyVariantType"),
  task_factory_(NULL),
  add_target_types_(),
  remove_target_types_()
  {}

// @brief apply function here
void
ModifyVariantTypeMover::apply( core::pose::Pose & pose ) 
{
	// Create map of target residues using taskoperation
  core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );

  if ( task_factory_ != 0 )
	{
    task = task_factory_->create_task_and_apply_taskoperations( pose );
		TR.Debug << "Initializing from packer task." << std::endl;
  }
	else
	{
		TR.Debug << "No packer task specified, using default task." << std::endl;
  }

	for (core::Size resi = 1; resi <= pose.n_residue(); resi++)
	{
		if( task->pack_residue(resi) )
		{
			core::chemical::ResidueTypeSet const & rsd_set(pose.residue(resi).residue_type_set());
			core::chemical::ResidueTypeCOP new_rsd_type(pose.residue(resi).type());

			BOOST_FOREACH(std::string remove_type, remove_target_types_)
			{
				new_rsd_type = rsd_set.get_residue_type_with_variant_removed( *new_rsd_type, remove_type);
			}

			BOOST_FOREACH(std::string add_type, add_target_types_)
			{
				new_rsd_type = rsd_set.get_residue_type_with_variant_added( *new_rsd_type, add_type);
			}

			core::pose::replace_pose_residue_copying_existing_coordinates( pose, resi, *new_rsd_type );
		}
	}
}

std::string
ModifyVariantTypeMover::get_name() const {
	return "ModifyVariantType";
}

moves::MoverOP
ModifyVariantTypeMover::clone() const
{
	return new ModifyVariantTypeMover( *this );
}

moves::MoverOP
ModifyVariantTypeMover::fresh_instance() const
{
	return new ModifyVariantTypeMover;
}

void ModifyVariantTypeMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/
)
{
	add_target_types_.clear();
	remove_target_types_.clear();

	std::string add_type_value = tag->getOption< std::string >( "add_type", "");
	//boost::split(add_target_types_, add_type_value, boost::is_any_of(","));
	add_target_types_ = utility::string_split<std::string>(add_type_value,',',std::string());

	std::string remove_type_value = tag->getOption< std::string >( "remove_type", "");
	//boost::split(remove_target_types_, remove_type_value, boost::is_any_of(","));
	remove_target_types_ = utility::string_split<std::string>(remove_type_value,',',std::string());

	if (add_target_types_.size() == 0 && remove_target_types_.size() == 0)
	{
		TR.Error << "Must specify add_type and/or remove_type type in ModifyVariantTypeMover." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("Must specify add_type and/or remove_type type in ModifyVariantTypeMover.");
	}

  task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data );
}

protocols::moves::MoverOP
ModifyVariantTypeMoverCreator::create_mover() const { return new ModifyVariantTypeMover; }

std::string
ModifyVariantTypeMoverCreator::keyname() const { return "ModifyVariantType"; }

} // moves
} // protocols

