// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/util.cc
/// @brief Utility functions for constraint generators
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/constraint_generator/util.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/id/SequenceMapping.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/pose/symmetry/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.constraint_generator.util" );

namespace protocols {
namespace constraint_generator {

/// @brief gets native pose from command line option, if one is specified
core::pose::PoseOP
get_native_pose()
{
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		core::pose::PoseOP ref_pose( new core::pose::Pose() );
		std::string const & native_pdb_fname = basic::options::option[ basic::options::OptionKeys::in::file::native ]();
		core::io::pdb::build_pose_from_pdb_as_is( *ref_pose, native_pdb_fname );
		return ref_pose;
	} else {
		return core::pose::PoseOP();
	}
}

core::scoring::func::FuncOP
scalar_weighted( core::scoring::func::FuncOP func, core::Real const weight )
{
	return core::scoring::func::FuncOP( new core::scoring::func::ScalarWeightedFunc( weight, func ) );
}

/// @brief creates a function from a text definition
core::scoring::func::FuncOP
create_func( std::string const & func_def )
{
	std::stringstream ss( func_def );
	std::string tag;
	ss >> tag;

	core::scoring::func::FuncOP new_func = core::scoring::constraints::ConstraintIO::get_instance()->get_func_factory().new_func( tag );
	if ( !new_func ) {
		utility_exit_with_message( "HydrogenBondConstraintGenerator could not create func from the definition " + func_def );
	}

	new_func->read_data( ss );
	return new_func;
}

/// @brief creates a sequencemapping from pose1 to pose2
core::id::SequenceMapping
generate_seqmap_from_poses(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2 )
{
	bool const same_length = ( pose1.total_residue() == pose2.total_residue() );
	bool const same_sequence = ( pose1.sequence() == pose2.sequence() );

	if ( same_length && same_sequence ) {
		return core::id::SequenceMapping::identity( pose1.total_residue() );
	} else { // !same_sequence || !same_length
		TR << "Input structure and native differ in ";
		if ( !same_length ) TR << "length and sequence ";
		else if ( !same_sequence ) TR << "sequence ";
		TR << "- aligning on PDB identity or sequence." << std::endl;
		return core::pose::sequence_map_from_pdbinfo( pose1, pose2 );
	}
}

/// @brief find number of residues -- if pose is symetric,
/// returns nubmer of symmetry-independent residues
core::Size
compute_nres_in_asymmetric_unit( core::pose::Pose const & pose )
{
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetryInfoCOP symm_info;
		core::conformation::symmetry::SymmetricConformation const & SymmConf(
			dynamic_cast< core::conformation::symmetry::SymmetricConformation const & >( pose.conformation() ) );
		symm_info = SymmConf.Symmetry_Info();
		return symm_info->num_independent_residues();
	} else {
		return pose.total_residue();
	}
}

/// @brief parses constraint generators from a tag
/// returns vector of ConstraintGeneratorCOPs
ConstraintGeneratorCOPs
parse_constraint_generators( utility::tag::TagCOP tag, basic::datacache::DataMap const & data )
{
	ConstraintGeneratorCOPs cgs;
	std::string const generators_str = tag->getOption< std::string >( "constraint_generators", "" );
	utility::vector1< std::string > const generator_strs = utility::string_split( generators_str, ',' );
	for (const auto & generator_str : generator_strs) {
		ConstraintGeneratorCOP new_cg = data.get_ptr< ConstraintGenerator const >( "ConstraintGenerators", generator_str );
		if ( !new_cg ) {
			std::stringstream msg;
			msg << "RemoveConstraints: Could not find a constraint generator named " << generator_str <<
				" in the datamap.  Make sure it is defined in an AddConstraints mover before RemoveConstraints is defined." << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
		cgs.push_back( new_cg );
	}
	return cgs;
}

} //protocols
} //constraint_generator

