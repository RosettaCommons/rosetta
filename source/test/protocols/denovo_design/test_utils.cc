// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/denovo_design/test_utils.hh
/// @brief Common test routines for permutation/component unit tests
/// @detailed
/// @author Tom Linsky

// Unit headers
#include <test/protocols/denovo_design/test_utils.hh>

// Protocol headers

// Package headers

// Core headers

// Basic/Numeric/Utility Headers
#include <basic/Tracer.hh>

// C++ Headers

static THREAD_LOCAL basic::Tracer utilTR( "test.protocols.denovo_design.test_utils.cc" );

namespace protocols {
namespace denovo_design {

// void checks a component in two poses to see if relative residue positions have moved
core::Size new_res(
		core::Size const resi,
		protocols::denovo_design::components::StructureData const & orig,
		protocols::denovo_design::components::StructureData const & newp )
{
	using namespace protocols::denovo_design;
	using namespace protocols::denovo_design::components;
	std::string comp = "";
	core::Size res = 0;
	for ( StringList::const_iterator r=orig.segments_begin(), end=orig.segments_end(); r != end; ++r ) {
		if ( ( orig.segment(*r).start() <= resi ) && ( resi <= orig.segment(*r).stop() ) ) {
			comp = *r;
			res = resi - orig.segment(*r).start() + 1;
			break;
		}
	}
	assert( comp != "" );
	assert( res );
	return newp.segment(comp).start() + res - 1;
}

// checks for unwanted movement between the two given components
void check_unwanted_movement(
		std::string const & comp1,
		std::string const & comp2,
		protocols::denovo_design::components::StructureData const & orig,
		protocols::denovo_design::components::StructureData const & newp )
{
	using namespace protocols::denovo_design::components;
	Segment res_s1 = orig.segment(comp1);
	Segment res_s4 = orig.segment(comp2);
	for ( core::Size i=res_s1.start(); i<=res_s1.stop(); ++i ) {
		for ( core::Size j=res_s4.start(); j<=res_s4.stop(); ++j ) {
			utilTR.Debug << i << "__" << j << " vs " << new_res( i, orig, newp ) << "__" << new_res( j, orig, newp ) << std::endl;
			TS_ASSERT_DELTA(
					orig.pose()->residue(i).xyz("CA").distance(orig.pose()->residue(j).xyz("CA")),
					newp.pose()->residue(new_res(i, orig, newp)).xyz("CA").distance(newp.pose()->residue(new_res(j, orig, newp)).xyz("CA")),
					1e-4 );
		}
	}
}

void check_movable_group(
		protocols::denovo_design::components::StructureData const & orig,
		protocols::denovo_design::components::StructureData const & perm,
		core::Size const group )
{
	utility::vector1< std::string > components1 = orig.segments_in_movable_group( group );
	utility::vector1< std::string > components2 = perm.segments_in_movable_group( group );
	TS_ASSERT( components1.size() == components2.size() );
	for ( core::Size i=1; i<=components1.size(); ++i ) {
		core::Size found = 0;
		for ( core::Size j=1; j<=components2.size(); ++j ) {
			if ( components1[i] == components2[j] ) {
				found = j;
				break;
			}
		}
		TS_ASSERT( found );
		TS_ASSERT_EQUALS( components1[i], components2[found] );
	}
	for ( core::Size i=1; i<=components1.size(); ++i ) {
		for ( core::Size j=i+1; j<=components1.size(); ++j ) {
			utilTR << "Checking for unwanted movement between " << components1[i] << " and " << components1[j] << std::endl;
			check_unwanted_movement( components1[i], components1[j], orig, perm );
		}
	}
}

void check_unwanted_movement(
		protocols::denovo_design::components::StructureData const & orig,
		protocols::denovo_design::components::StructureData const & perm )
{
	std::set< core::Size > groups = orig.movable_groups();
	for ( std::set< core::Size >::const_iterator g = groups.begin(); g != groups.end(); ++g ) {
		check_movable_group( orig, perm, *g );
	}
}

void check_sequential(
		protocols::denovo_design::components::StructureData const & perm,
		std::string const & c1,
		std::string const & c2,
		std::string const & c3 )
{
	TS_ASSERT_EQUALS( perm.segment(c1).cterm_resi()+1, perm.segment(c2).nterm_resi() );
	TS_ASSERT_EQUALS( perm.segment(c2).cterm_resi()+1, perm.segment(c3).nterm_resi() );
}

}
}
