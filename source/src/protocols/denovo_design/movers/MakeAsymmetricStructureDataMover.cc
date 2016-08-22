// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/movers/MakeAsymmetricStructureDataMover.cc
/// @brief Converts a StructureData for a symmetric pose into an asymmetric representation
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/movers/MakeAsymmetricStructureDataMover.hh>
#include <protocols/denovo_design/movers/MakeAsymmetricStructureDataMoverCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.MakeAsymmetricStructureDataMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

MakeAsymmetricStructureDataMover::MakeAsymmetricStructureDataMover():
	protocols::moves::Mover( MakeAsymmetricStructureDataMover::class_name() )
{

}

MakeAsymmetricStructureDataMover::~MakeAsymmetricStructureDataMover(){}

void
MakeAsymmetricStructureDataMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

}

protocols::moves::MoverOP
MakeAsymmetricStructureDataMover::clone() const
{
	return protocols::moves::MoverOP( new MakeAsymmetricStructureDataMover( *this ) );
}

protocols::moves::MoverOP
MakeAsymmetricStructureDataMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new MakeAsymmetricStructureDataMover );
}

std::string
MakeAsymmetricStructureDataMover::get_name() const
{
	return MakeAsymmetricStructureDataMover::class_name();
}

std::string
MakeAsymmetricStructureDataMover::class_name()
{
	return "MakeAsymmetricStructureDataMover";
}

void
MakeAsymmetricStructureDataMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, MakeAsymmetricStructureDataMover const & mover )
{
	mover.show(os);
	return os;
}

std::string
new_id( core::Size const subunit, std::string const & orig_name )
{
	static std::string const prefix = "S";
	static std::string const delimeter = ".";

	std::stringstream new_id;
	new_id << prefix << subunit << delimeter << orig_name;
	return new_id.str();
}

void
MakeAsymmetricStructureDataMover::apply( core::pose::Pose & pose )
{
	components::StructureDataFactory const & factory = *components::StructureDataFactory::get_instance();
	runtime_assert( pose.pdb_info() );
	components::StructureData const sd = factory.create_from_remarks( pose.pdb_info()->remarks() );
	TR << "Orig sd=" << sd << std::endl;

	components::StructureData asymm_sd;

	core::Size subunit = 1;
	core::Size resid = 0;
	while ( resid + sd.pose_length() <= pose.total_residue() ) {
		for ( SegmentNameList::const_iterator s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
			components::Segment newseg = sd.segment( *s );
			newseg.set_id( new_id( subunit, *s ) );
			if ( !newseg.lower_segment().empty() ) {
				newseg.set_lower_segment( new_id( subunit, newseg.lower_segment() ) );
			}
			if ( !newseg.upper_segment().empty() ) {
				newseg.set_upper_segment( new_id( subunit, newseg.upper_segment() ) );
			}
			asymm_sd.add_segment( newseg );
		}
		resid += sd.pose_length();
		++subunit;
	}

	TR << "Asymm SD=" << asymm_sd << std::endl;
	components::StructureDataFactory::get_instance()->save_into_pose( pose, asymm_sd );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
MakeAsymmetricStructureDataMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new MakeAsymmetricStructureDataMover );
}

std::string
MakeAsymmetricStructureDataMoverCreator::keyname() const
{
	return MakeAsymmetricStructureDataMover::class_name();
}

} //protocols
} //denovo_design
} //movers

