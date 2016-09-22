// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/protein_interface_design/filters/DisulfideFilter.hh
/// @brief Filters for interfaces which could form a disulfide bond between
/// docking partners.
/// @author Sarel Fleishman (sarelf@uw.edu)

#include <protocols/simple_filters/AtomicContactFilter.hh>
#include <protocols/simple_filters/AtomicContactFilterCreator.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

//parsing
#include <utility/tag/Tag.hh>
#include <core/conformation/Residue.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.filters.AtomicContactFilter" );

/// @brief default ctor
AtomicContactFilter::AtomicContactFilter() :
	parent( "AtomicContact" ),
	protocols::moves::ResId( 0u )
{}

/// @brief Constructor with a single target residue
AtomicContactFilter::AtomicContactFilter( core::Size const res1, core::Size const res2, core::Real const distance, bool const sidechain, bool const backbone, bool const protons ) :
	parent( "AtomicContact" ),
	protocols::moves::ResId( res2 ),
	residue1_( res1 ),
	distance_( distance ),
	sidechain_( sidechain ),
	backbone_( backbone ),
	protons_( protons )
{
	range1_.push_back( residue1_ );
}

/// @return Whether a disulfide bond is possible between any of the targets
bool AtomicContactFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	report( TR.Debug, pose );
	if ( dist <= distance_ ) return true;
	return false;
}

core::Real
AtomicContactFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core::conformation;

	if ( !get_resid() ) {
		TR.Error<<"ERROR: residue2 has not been defined"<<std::endl;
		runtime_assert( get_resid() );
	}
	core::Real nearest_distance( 10000 );
	Residue const res2( pose.residue( get_resid() ) );
	for ( core::Size residue1 : range1_ ) {
		Residue const res1( pose.residue( residue1 ) );

		auto atom1_begin( res1.atom_begin() ), atom1_end( res1.atom_end() ), atom2_begin( res2.atom_begin() ), atom2_end( res2.atom_end() );
		if ( sidechain_ && !backbone_ ) {
			atom1_begin = res1.sidechainAtoms_begin();
			atom2_begin = res2.sidechainAtoms_begin();
		}
		if ( !sidechain_ && backbone_ ) {
			atom1_end = res1.sidechainAtoms_begin();
			atom2_end = res2.sidechainAtoms_begin();
		}
		if ( !protons_ ) {
			atom1_end = res1.heavyAtoms_end();
			atom2_end = res2.heavyAtoms_end();
		}
		for ( auto atom1=atom1_begin; atom1!=atom1_end; ++atom1 ) {
			for ( auto atom2=atom2_begin; atom2!=atom2_end; ++atom2 ) {
				core::Real const dist( atom1->xyz().distance( atom2->xyz() ) );
				if ( dist <= nearest_distance ) nearest_distance = dist;
			}//foreach atom2
		}//foreach atom1
	}
	return( nearest_distance );
}

core::Real
AtomicContactFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	return( dist );
}

void AtomicContactFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	out<<"Minimal distance between residues "<<residue1_<<" and "<<get_resid()<<" is "<<dist<<std::endl;
}

void AtomicContactFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose)
{
	distance_ = tag->getOption< core::Real >( "distance", 4.0 );
	if ( tag->hasOption("range1") ) {
		residue1_ = 0;
		std::istringstream range_str( tag->getOption< std::string >( "range1" ) );
		core::Size num1, num2;
		range_str >> num1 >> num2;
		if ( !range_str.good() ) TR << "cannot read parameter range1" << std::endl;
		for ( core::Size i=num1; i<=num2; ++i ) {
			range1_.push_back( i );
		}
	}
	if ( range1_.size() == 0 ) {
		std::string const res1( tag->getOption< std::string >( "residue1" ) );
		residue1_ = core::pose::parse_resnum( res1, pose );
		range1_.push_back( residue1_ );
	}
	if ( tag->hasOption( "residue2" ) ) {
		std::string const res2( tag->getOption< std::string >( "residue2" ) );
		set_resid( core::pose::parse_resnum( res2, pose ) );
		modifiable( false );
	} else {
		modifiable( true );
		TR<<"AtomicContact: residue2 was not defined. A mover/filter will have to set it during the protocol\n";
	}
	sidechain_ = tag->getOption< bool >( "sidechain", 1 );
	backbone_  = tag->getOption< bool >( "backbone",  0 );
	protons_   = tag->getOption< bool >( "protons",   0 );

	TR<<"AtomicContact filter between residues "<<residue1_<<" and "<<get_resid()<<" with distance cutoff of "<<distance_<<std::endl;
}

filters::FilterOP
AtomicContactFilterCreator::create_filter() const { return filters::FilterOP( new AtomicContactFilter ); }

std::string
AtomicContactFilterCreator::keyname() const { return "AtomicContact"; }


} // filters
} // protocols
