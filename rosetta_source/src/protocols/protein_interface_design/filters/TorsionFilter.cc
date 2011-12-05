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
#include <sstream>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>


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
		std::stringstream s("");
		for( core::Size i = 1; i <= pose.total_residue(); ++i ){
			if( i % 5 == 0 ) s << pose.residue( i ).name1()<<pose.pdb_info()->number( i )<<pose.pdb_info()->chain( i )<<'\t';
			TR<<"Residue "<<pose.residue( i ).name1()<<pose.pdb_info()->number( i )<<pose.pdb_info()->chain( i )<<'\t';
			if( torsion() == "phi" || torsion() == "" ){
				TR<<" phi "<<pose.phi( i )<<'\t';
				s<<pose.phi( i )<<' ';
			}
			if( torsion() == "psi" || torsion() == "" ){
				TR<<" psi "<<pose.psi( i )<<'\t';
				s<<pose.psi( i )<<' ';
			}
			if( torsion() == "omega" || torsion() == "" ){
				TR<<" omega "<<pose.omega( i )<<std::endl;
				s<<pose.omega( i )<<' ';
			}
		}
		TR<<s.str()<<std::endl;
		return true;
	}
	else{
		TR<<"Residue "<<pose.residue( resnum() ).name1()<<pose.pdb_info()->number( resnum() )<<pose.pdb_info()->chain( resnum() )<<'\t';
		if( torsion() == "phi" || torsion() == "" ){
			core::Real const phi( pose.phi( resnum() ) );
			TR<<" phi "<<phi<<std::endl;
			if( torsion() == "phi" )
				return( phi>=lower() && phi<=upper() );
		}
		if( torsion() == "psi" || torsion() == "" ){
			core::Real const psi( pose.psi( resnum() ) );
			TR<<" psi "<<pose.psi( resnum() )<<std::endl;
			if( torsion() == "psi" )
				return( psi>=lower() && psi<=upper() );
		}
	}

	return false;
}

core::Real
Torsion::compute( core::pose::Pose const & p ) const{
	if( resnum() > 0 ){
		if( torsion() == "phi" ) return p.phi( resnum() );
		if( torsion() == "psi" ) return p.psi( resnum() );
		if( torsion() == "omega" ) return p.omega( resnum() );
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
	torsion( tag->getOption< std::string >( "torsion", "" ) );
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
