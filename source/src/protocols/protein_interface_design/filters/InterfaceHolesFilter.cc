// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/InterfaceHolesFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/InterfaceHolesFilter.hh>
#include <protocols/protein_interface_design/filters/InterfaceHolesFilterCreator.hh>
#include <protocols/filters/Filter.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/scoring/Interface.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/scoring/Energies.hh>

#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {


static thread_local basic::Tracer TR( "protocols.protein_interface_design.filters.InterfaceHolesFilter" );
core::Real
InterfaceHolesFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core;
	using namespace core::scoring;
	pose::Pose copy_pose = pose;
	ScoreFunctionOP sfxn( new ScoreFunction );
	sfxn->set_weight( holes_decoy, 1.0 );
	core::Real const weight( (*sfxn)[ ScoreType( holes_decoy ) ] );
	protocols::scoring::Interface iface( rb_jump_ );
	iface.calculate( copy_pose );

	(*sfxn)( copy_pose );
	Real const bound_holes( copy_pose.energies().total_energies()[ ScoreType( holes_decoy ) ]);
	TR.Debug << "Bound holes: " << bound_holes << std::endl;

	protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( copy_pose, rb_jump_ ) );
	translate->step_size( 1000.0 );
	translate->apply( copy_pose );
	(*sfxn)( copy_pose );
	Real const unbound_holes( copy_pose.energies().total_energies()[ ScoreType( holes_decoy ) ]);
	TR.Debug << "Unbound holes: " << unbound_holes << std::endl;

	Real delta_holes = (bound_holes - unbound_holes) * weight; // more negative as bound state gets better

	// Normalize for number of heavy atoms in the pose. Otherwise larger interfaces will always be larger
	//Brings this in line with scores reported by Will's holes app.
	// 8-19-09: this now appears to be done internally by the scoretype. No need to do it myself
	//core::Size nheavyatoms=0;
	//for( core::Size i=1; i<= copy_pose.total_residue(); ++i ) {
	// nheavyatoms += copy_pose.residue(i).nheavyatoms();
	// }
	//delta_holes /= nheavyatoms;
	return delta_holes;
}

InterfaceHolesFilter::~InterfaceHolesFilter() {}

bool
InterfaceHolesFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const iface_holes( compute( pose ));
	TR << "delta holes score: " << iface_holes ;
	if ( iface_holes <= threshold_ ) {
		TR<<" passing."<<std::endl;
		return( true );
	} else TR<<" failing." << std::endl;
	return( false );
}

void
InterfaceHolesFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const iface_holes( compute( pose ));
	out<<"InterfaceHoles score across jump # "<< rb_jump_ << ": " << iface_holes<<'\n';
}

core::Real
InterfaceHolesFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const iface_holes( compute( pose ));
	return( (core::Real) iface_holes );
}

void
InterfaceHolesFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	threshold_ = tag->getOption<core::Size>( "threshold", 200 );
	rb_jump_ = tag->getOption<core::Size>( "jump", 1 );

	TR<<"Using jump number " << rb_jump_ << " with threshold "<< threshold_ <<std::endl;
}

protocols::filters::FilterOP
InterfaceHolesFilterCreator::create_filter() const { return protocols::filters::FilterOP( new InterfaceHolesFilter ); }

std::string
InterfaceHolesFilterCreator::keyname() const { return "InterfaceHoles"; }


} // filters
} // protein_interface_design
} // devel


