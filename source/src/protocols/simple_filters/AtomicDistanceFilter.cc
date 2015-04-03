// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/AtomicDistanceFilter.cc
/// @brief Filter for looking at specific atom distances
/// @author Rocco Moretti (rmoretti@uw.edu)

#include <protocols/simple_filters/AtomicDistanceFilter.hh>
#include <protocols/simple_filters/AtomicDistanceFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

//parsing
#include <utility/tag/Tag.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace simple_filters {

static thread_local basic::Tracer TR( "protocols.filters.AtomicDistanceFilter" );

/// @brief default ctor
AtomicDistanceFilter::AtomicDistanceFilter() :
	parent( "AtomicDistance" )
{}

/// @brief
AtomicDistanceFilter::AtomicDistanceFilter( core::Size const res1, core::Size const res2, std::string atom_desig1, std::string atom_desig2, bool as_type1, bool as_type2, core::Real distance) :
	parent( "AtomicDistance" ),
	residue1_( res1 ),
	residue2_( res2 ),
	atomdesg1_( atom_desig1 ),
	atomdesg2_( atom_desig2 ),
	astype1_( as_type1 ),
	astype2_( as_type2 ),
	distance_( distance )
{}

/// @return Whether the atom pair is within the cutoff
bool AtomicDistanceFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	report( TR.Debug, pose );
	if( dist <= distance_ ) return true;
	return false;
}

core::Real
AtomicDistanceFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core::conformation;

	core::Real nearest_distance( 999999 );

	assert(residue1_ <= pose.n_residue() && residue2_ <= pose.n_residue());
	Residue const& res1( pose.residue( residue1_ ) );
	Residue const& res2( pose.residue( residue2_ ) );

	core::Size a1primet(1), a1end(res1.natoms());
	if ( ! astype1_ ) { // If given by name, look only at the single atom
		if ( ! res1.type().has(atomdesg1_) ) {
			TR << "WARNING! Residue "<<residue1_<<" of type "<<res1.type().name()<<" does not have atom with name "<<atomdesg1_<<std::endl;
			return nearest_distance;
		}
		a1primet = a1end = res1.atom_index(atomdesg1_);
	}
	core::Size a2primet(1), a2end(res2.natoms());
	if ( ! astype2_ ) { // If given by name, look only at the single atom
		if ( ! res2.type().has(atomdesg2_) ) {
			TR << "WARNING! Residue "<<residue2_<<" of type "<<res2.type().name()<<" does not have atom with name "<<atomdesg2_<<std::endl;
			return nearest_distance;
		}
		a2primet = a2end = res2.atom_index(atomdesg2_);
	}

	bool found1(false), found2(false);
	for ( core::Size ii(a1primet); ii <= a1end; ++ii) {
		if ( !astype1_ || res1.atom_type(ii).name() == atomdesg1_ ) {
			found1 = true;
			for ( core::Size jj(a2primet); jj <= a2end; ++jj) {
				if ( !astype2_ || res2.atom_type(jj).name() == atomdesg2_ ) {
					found2 = true;
					core::Real const dist( res1.atom(ii).xyz().distance( res2.atom(jj).xyz() ) );
					if ( dist < nearest_distance ) {
						nearest_distance = dist;
					}
				}
			}
		}
	}

	if ( ! found1 ) {
		TR << "WARNING! Residue "<<residue1_<<" of type "<<res1.type().name()<<" does not have atom with "<<(astype1_?"type ":"name ")<<atomdesg1_<<std::endl;
	}
	else if ( ! found2 ) { // elseif because the inner loop doesn't run if the outer loop doesn't trip. (if found1 is false, found2 is always false)
		TR << "WARNING! Residue "<<residue2_<<" of type "<<res2.type().name()<<" does not have atom with "<<(astype2_?"type ":"name ")<<atomdesg2_<<std::endl;
	}
	return( nearest_distance );
}

core::Real
AtomicDistanceFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	return( dist );
}

void AtomicDistanceFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	out<<"Minimal distance between residue "<<residue1_<<" atom "<<(astype1_?"type ":"name ")<<atomdesg1_<<" and residue "<<residue2_<<" atom "<<(astype2_?"type ":"name ")<<atomdesg2_<<" is "<<dist<<std::endl;
}

void AtomicDistanceFilter::parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & pose)
{
	distance_ = tag->getOption< core::Real >( "distance", 4.0 );

	std::string const res1( tag->getOption< std::string >( "residue1" ) );
	std::string const res2( tag->getOption< std::string >( "residue2" ) );
	residue1_ = core::pose::parse_resnum( res1, pose );
	residue2_ = core::pose::parse_resnum( res2, pose );

	if (residue1_ == 0) {
		TR << "Residue number "<<res1<<" not found in pose."<<std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("Residue number not found. Check xml file");
	}
	if (residue2_ == 0) {
		TR << "Residue number "<<res2<<" not found in pose."<<std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("Residue number not found. Check xml file");
	}

	if( tag->hasOption( "atomtype1" ) ) {
		if( tag->hasOption( "atomname1" ) ) {
			throw utility::excn::EXCN_RosettaScriptsOption("Can't set both atomname1 and atomtype1. Check xml file");
		}
		atomdesg1_ = tag->getOption< std::string >( "atomtype1" );
		astype1_ = true;
	}
	else {
		atomdesg1_ = tag->getOption< std::string >( "atomname1", "CB" );
		astype1_ = false;
	}

	if( tag->hasOption( "atomtype2" ) ) {
		if( tag->hasOption( "atomname2" ) ) {
			throw utility::excn::EXCN_RosettaScriptsOption("Can't set both atomname2 and atomtype2. Check xml file");
		}
		atomdesg2_ = tag->getOption< std::string >( "atomtype2" );
		astype2_ = true;
	}
	else {
		atomdesg2_ = tag->getOption< std::string >( "atomname2", "CB" );
		astype2_ = false;
	}

	TR<<"AtomicDistance filter between residue "<<residue1_<<" atom "<<(astype1_?"type ":"name ")<<atomdesg1_<<" and residue "<<residue2_<<" atom "<<(astype2_?"type ":"name ")<<atomdesg2_<<" with distance cutoff of "<<distance_<<std::endl;
}

protocols::filters::FilterOP
AtomicDistanceFilterCreator::create_filter() const { return protocols::filters::FilterOP( new AtomicDistanceFilter ); }

std::string
AtomicDistanceFilterCreator::keyname() const { return "AtomicDistance"; }


} // filters
} // protocols
