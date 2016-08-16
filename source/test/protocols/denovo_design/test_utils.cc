// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/test_utils.hh
/// @brief Common test routines for permutation/component unit tests
/// @detailed
/// @author Tom Linsky

// Unit headers
#include <test/protocols/denovo_design/test_utils.hh>

// Protocol headers

// Package headers

// Core headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Basic/Numeric/Utility Headers
#include <basic/Tracer.hh>

// C++ Headers

static THREAD_LOCAL basic::Tracer utilTR( "test.protocols.denovo_design.test_utils.cc" );

namespace protocols {
namespace denovo_design {

// void checks a component in two poses to see if relative residue positions have moved
core::Size
new_res(
	core::Size const resi,
	StructureData const & orig,
	StructureData const & test )
{
	using namespace protocols::denovo_design;
	using namespace protocols::denovo_design::components;
	std::string comp = "";
	core::Size res = 0;
	for ( SegmentNameList::const_iterator r=orig.segments_begin(); r!=orig.segments_end(); ++r ) {
		if ( ( orig.segment(*r).start() <= resi ) && ( resi <= orig.segment(*r).stop() ) ) {
			comp = *r;
			res = resi - orig.segment(*r).start() + 1;
			break;
		}
	}
	assert( comp != "" );
	assert( res );
	return test.segment(comp).start() + res - 1;
}

// checks for unwanted movement between the two given components
void
check_unwanted_movement(
	std::string const & comp1,
	std::string const & comp2,
	StructureData const & orig, core::pose::Pose const & orig_pose,
	StructureData const & test, core::pose::Pose const & test_pose )
{
	Segment res_s1 = orig.segment(comp1);
	Segment res_s4 = orig.segment(comp2);
	for ( core::Size i=res_s1.start(); i<=res_s1.stop(); ++i ) {
		for ( core::Size j=res_s4.start(); j<=res_s4.stop(); ++j ) {
			utilTR.Debug << i << "__" << j << " vs " << new_res( i, orig, test ) << "__" << new_res( j, orig, test ) << std::endl;
			TS_ASSERT_DELTA(
					orig_pose.residue(i).xyz("CA").distance(orig_pose.residue(j).xyz("CA")),
					test_pose.residue(new_res(i, orig, test)).xyz("CA").distance(test_pose.residue(new_res(j, orig, test)).xyz("CA")),
					1e-4 );
		}
	}
}

void
check_movable_group(
		StructureData const & orig, core::pose::Pose const & orig_pose,
		StructureData const & test, core::pose::Pose const & test_pose,
		core::Size const group )
{
	SegmentNames const components1 = orig.segments_in_movable_group( group );
	SegmentNames const components2 = test.segments_in_movable_group( group );
	TS_ASSERT_EQUALS( components1.size(), components2.size() );
	for ( SegmentNames::const_iterator s1=components1.begin(); s1!=components1.end(); ++s1 ) {
		TS_ASSERT_DIFFERS( std::find( components2.begin(), components2.end(), *s1 ), components2.end() );
	}

	for ( SegmentNames::const_iterator s1=components1.begin(); s1!=components1.end(); ++s1 ) {
		SegmentNames::const_iterator s2 = s1;
		++s2;
		for ( ; s2!=components1.end(); ++s2 ) {
			utilTR << "Checking for unwanted movement between " << *s1 << " and " << *s2 << std::endl;
			check_unwanted_movement( *s1, *s2, orig, orig_pose, test, test_pose );
		}
	}
}

void
check_unwanted_movement(
	StructureData const & orig, core::pose::Pose const & orig_pose,
	StructureData const & test, core::pose::Pose const & test_pose )
{
	MovableGroups const & groups = orig.movable_groups();
	for ( MovableGroups::const_iterator g=groups.begin(); g!=groups.end(); ++g ) {
		check_movable_group( orig, orig_pose, test, test_pose, *g );
	}
}

void
check_sequential(
	protocols::denovo_design::components::StructureData const & perm,
	std::string const & c1,
	std::string const & c2,
	std::string const & c3 )
{
	TS_ASSERT_EQUALS( perm.segment(c1).upper()+1, perm.segment(c2).lower() );
	TS_ASSERT_EQUALS( perm.segment(c2).upper()+1, perm.segment(c3).lower() );
}

}
}
