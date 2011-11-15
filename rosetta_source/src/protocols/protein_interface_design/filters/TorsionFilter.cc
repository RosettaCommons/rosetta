// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/TorsionFilter.hh>
#include <protocols/protein_interface_design/filters/TorsionFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>

namespace protocols {
namespace protein_interface_design{
namespace filters {

static basic::Tracer TR( "protocols.protein_interface_design.filters.Torsion" );

///@brief default ctor
Torsion::Torsion() :
	parent( "Torsion" ),
	lower_( false ),
	upper_( true ),
	resnum_( 0 ),
	torsion_( "" )
{}

bool
Torsion::apply(core::pose::Pose const & pose ) const
{
	if( resnum() == 0 ){ // just print all torsions
		for( core::Size i = 1; i <= pose.total_residue(); ++i ){
			if( torsion() == "phi" || torsion() == "" )
				TR<<"Residue "<<i<<" phi "<<pose.phi( i )<<std::endl;
			if( torsion() == "psi" || torsion() == "" )
				TR<<"Residue "<<i<<" psi "<<pose.psi( i )<<std::endl;
		}
	}
	else{
			if( torsion() == "phi" || torsion() == "" )
				TR<<"Residue "<<resnum()<<" phi "<<pose.phi( resnum() )<<std::endl;
			if( torsion() == "psi" || torsion() == "" )
				TR<<"Residue "<<resnum()<<" psi "<<pose.psi( resnum() )<<std::endl;
		}

}

core::Real
Torsion::compute( core::pose::Pose const & p ) const{
	if( resnum() > 0 ){
		if( torsion() == "phi" ) return p.phi( resnum() );
		if( torsion() == "psi" ) return p.psi( resnum() );
	}
}

core::Real
Torsion::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose ) );
}

void
Torsion::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"Torsion returns "<<compute( pose )<<std::endl;
}

void
Torsion::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose )
{
	lower( tag->getOption< core::Real >( "lower", 0 ) );
	upper( tag->getOption< core::Real >( "upper", 0 ) );
	if( tag->hasOption( "resnum" ))
		resnum( protocols::rosetta_scripts::parse_resnum( tag->getOption< std::string >( "resnum" ), pose ) );
	else resnum( 0 );
	TR<<"resnum: "<<resnum()<<" lower "<<lower()<<" upper: "<<upper()<<std::endl;
}

protocols::filters::FilterOP
Torsion::fresh_instance() const{
	return new Torsion();
}

Torsion::~Torsion(){}

protocols::filters::FilterOP
Torsion::clone() const{
	return new Torsion( *this );
}

protocols::filters::FilterOP
TorsionCreator::create_filter() const { return new Torsion; }

std::string
TorsionCreator::keyname() const { return "Torsion"; }

} // filters
} // protein_interface_design
} // protocols
