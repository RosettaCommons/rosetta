// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/InterfaceSasaFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/InterfaceSasaFilter.hh>
#include <protocols/simple_filters/InterfaceSasaFilterCreator.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <basic/MetricValue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.hh>
// Jacob
#include <utility/string_util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

// Project Headers
#include <core/types.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer interface_sasa_filter_tracer( "protocols.simple_filters.InterfaceSasaFilter" );

protocols::filters::FilterOP
InterfaceSasaFilterCreator::create_filter() const { return new InterfaceSasaFilter; }

std::string
InterfaceSasaFilterCreator::keyname() const { return "Sasa"; }


InterfaceSasaFilter::InterfaceSasaFilter() :
	Filter( "Sasa" ),
	lower_threshold_( 0.0 ),
	hydrophobic_( false ),
	polar_( false ),
	jump_( 1 ),
	sym_dof_names_( "" ),
	upper_threshold_(100000000000.0)
{}

InterfaceSasaFilter::InterfaceSasaFilter( core::Real const lower_threshold, bool const hydrophobic/*=false*/, bool const polar/*=false*/, core::Real upper_threshold, std::string sym_dof_names ) :
	Filter( "Sasa" ),
	lower_threshold_( lower_threshold ),
	hydrophobic_( hydrophobic ),
	polar_( polar ),
	upper_threshold_(upper_threshold),
	sym_dof_names_(sym_dof_names)
{}

InterfaceSasaFilter::~InterfaceSasaFilter(){}

filters::FilterOP
InterfaceSasaFilter::clone() const{
	return new InterfaceSasaFilter( *this );
}

filters::FilterOP
InterfaceSasaFilter::fresh_instance() const{
	return new InterfaceSasaFilter;
}

void
InterfaceSasaFilter::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap &, filters::Filters_map const &,moves::Movers_map const &, core::pose::Pose const & )
{
	lower_threshold_ = tag->getOption<core::Real>( "threshold", 800 );
	upper_threshold_ = tag->getOption<core::Real>( "upper_threshold", 1000000);
	jump( tag->getOption< core::Size >( "jump", 1 ));
	sym_dof_names( tag->getOption< std::string >( "sym_dof_names" , "" ) );
	hydrophobic_ = tag->getOption<bool>( "hydrophobic", false );
	polar_ = tag->getOption<bool>( "polar", false );
	runtime_assert( !hydrophobic_ || !polar_ );
	if( jump() != 1 && ( polar_ || hydrophobic_ ) )
		utility_exit_with_message( "ERROR: presently, only total sasa is supported across a jump other than 1. Remove polar and hydrophobic flags and try again." );

	interface_sasa_filter_tracer<<"SasaFilter with lower threshold of "<<lower_threshold_<<" Ang^2 and jump "<<jump()<<'\n';
	if( hydrophobic_ )
		interface_sasa_filter_tracer<<"Only reporting hydrophobic sasa\n";
	if( polar_ )
		interface_sasa_filter_tracer<<"Only reporting polar sasa\n";
	interface_sasa_filter_tracer.flush();
}

bool
InterfaceSasaFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ) );

	interface_sasa_filter_tracer<<"sasa is "<<sasa<<". ";
	if( sasa >= lower_threshold_ && sasa <= upper_threshold_ ){
		interface_sasa_filter_tracer<<"passing." <<std::endl;
		return true;
	}
	else {
		interface_sasa_filter_tracer<<"failing."<<std::endl;
		return false;
	}
}

void
InterfaceSasaFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ));
	out<<"Sasa= "<< sasa<<'\n';
}

core::Real
InterfaceSasaFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ));
	return( sasa );
}

void
InterfaceSasaFilter::jump( core::Size const jump )
{
	jump_ = jump;
}

core::Size
InterfaceSasaFilter::jump() const
{
	return jump_;
}

void
InterfaceSasaFilter::sym_dof_names( std::string const sym_dof_names )
{
	sym_dof_names_ = sym_dof_names;
}

std::string
InterfaceSasaFilter::sym_dof_names() const
{
	return sym_dof_names_;
}

core::Real
InterfaceSasaFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::pose::metrics;
	using basic::MetricValue;
	using namespace core;
	using namespace protocols::moves;

	core::pose::Pose split_pose( pose );
	// JBB 120819
	int sym_aware_jump_id = 0;
	if ( sym_dof_names() != "" ) {
		utility::vector1<std::string> sym_dof_name_list;
		sym_dof_name_list = utility::string_split( sym_dof_names() , ',' );
		for (Size i = 1; i <= sym_dof_name_list.size(); i++) {
			sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( split_pose, sym_dof_name_list[i] );
			protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( split_pose, sym_aware_jump_id ) );
			translate->step_size( 1000.0 );
			translate->apply( split_pose );
		}
	} else {
		sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( split_pose, jump() );
		protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( split_pose, sym_aware_jump_id ) );
		translate->step_size( 1000.0 );
		translate->apply( split_pose );
	}
	//split_pose.dump_pdb("split_pose.pdb");
	// JBB 120819

	runtime_assert( !hydrophobic_ || !polar_ );
	if( !hydrophobic_ && !polar_ ){
		MetricValue< core::Real > mv_sasa;

		pose.metric( "sasa", "total_sasa", mv_sasa);
		core::Real const bound_sasa( mv_sasa.value() );
		split_pose.metric( "sasa", "total_sasa", mv_sasa );
		core::Real const unbound_sasa( mv_sasa.value() );
		if( core::pose::symmetry::is_symmetric( pose )) {
			core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
			core::Real const buried_sasa( (unbound_sasa - bound_sasa) /(sym_info->subunits()));
			interface_sasa_filter_tracer << "subunits: " << sym_info->subunits() << std::endl;
			return( buried_sasa );
		} else {
			core::Real const buried_sasa(unbound_sasa - bound_sasa);
			return( buried_sasa );
		}
	}
	else{
		MetricValue< id::AtomID_Map< core::Real > > atom_sasa;
		pose.metric( "sasa_interface", "delta_atom_sasa", atom_sasa );
		core::Real polar_sasa( 0.0 ), hydrophobic_sasa( 0.0 );
		for( core::Size pos(1); pos<=pose.total_residue(); ++pos ){
			for( core::Size atomi( 1 ); atomi <= atom_sasa.value().n_atom( pos ); ++atomi ){
				core::Real const atomi_delta_sasa( atom_sasa.value()( pos, atomi ) );
				core::conformation::Residue const pos_rsd( pose.residue( pos ) );
				core::chemical::AtomType const atom_type( pos_rsd.atom_type( atomi ) );
				bool const is_polar( atom_type.is_donor() || atom_type.is_acceptor() || atom_type.is_polar_hydrogen() );
				if( is_polar ) polar_sasa += atomi_delta_sasa;
				else hydrophobic_sasa += atomi_delta_sasa;
			}
		}
		if( hydrophobic_ ) return hydrophobic_sasa;
		if( polar_ ) return polar_sasa;
	}
	utility_exit_with_message( "Execution should never have reached this point." );
	return( 0 );
}

}
}
