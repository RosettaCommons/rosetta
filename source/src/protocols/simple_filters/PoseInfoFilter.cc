// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/PoseInfoFilter.cc
/// @brief Filter for looking at specific atom distances
/// @author Rocco Moretti (rmoretti@uw.edu)

#include <protocols/simple_filters/PoseInfoFilter.hh>
#include <protocols/simple_filters/PoseInfoFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>


namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.filters.PoseInfoFilter" );

/// @brief default ctor
PoseInfoFilter::PoseInfoFilter() :
	parent( "PoseInfo" )
{}

/// @return Print pose information and return true
bool PoseInfoFilter::apply(core::pose::Pose const & pose ) const
{
	compute( pose );
	return true;
}

core::Real
PoseInfoFilter::compute( core::pose::Pose const & pose ) const
{
	report(TR, pose);
	return 1.0;
}

/// @return Print pose information and return true
core::Real
PoseInfoFilter::report_sm( core::pose::Pose const & pose ) const
{
	compute( pose );
	return( 1 );
}

void PoseInfoFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out << "Pose Information: " << pose.size() << " residues " << std::endl;
	pose.pdb_info()->show( out );

	out << pose.fold_tree(); // Has implicit internal std::endl;
	pose.fold_tree().show( out );

	{
		using core::pose::datacache::CacheableDataType;
		out << "Cached Data: ";
		core::pose::Pose::BasicDataCache const & datacache( pose.data() );
		for ( core::Size ii(1); ii <= CacheableDataType::num_cacheable_data_types; ++ii ) {
			if ( datacache.has( ii ) ) {
				// Just output name here, as each datatype needs to be treated differently
				// (Do it below, if necessary.)
				out << CacheableDataType::get_name( static_cast<CacheableDataType::Enum>( ii ) ) << " ";
			}
		}
		out << std::endl; // To flush the Cached data list.
	}

	// Feel free to add additional pose-related information.
	// The main reason I stopped where I did was I didn't necessarily know how to access relevant others.
	out << "Constraints: ";
	pose.constraint_set()->show(out);
	pose.constraint_set()->show_definition(out,pose);
	pose.constraint_set()->show_numbers(out);
}

void PoseInfoFilter::parse_my_tag( utility::tag::TagCOP const,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/)
{
	// Right now we don't have any options to control, so don't bother doing anything.
}

protocols::filters::FilterOP
PoseInfoFilterCreator::create_filter() const { return protocols::filters::FilterOP( new PoseInfoFilter ); }

std::string
PoseInfoFilterCreator::keyname() const { return "PoseInfo"; }


} // filters
} // protocols
