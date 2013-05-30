// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Jacob Bale (balej@uw.edu)
#include <devel/matdes/AverageInterfaceEnergyFilter.hh>
#include <devel/matdes/AverageInterfaceEnergyFilterCreator.hh>

#include <protocols/simple_filters/EnergyPerResidueFilter.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/elscripts/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <utility/string_util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pose/symmetry/util.hh>
#include <ObjexxFCL/format.hh>


namespace devel {
namespace matdes {

static basic::Tracer TR( "devel.matdes.AverageInterfaceEnergyFilter" );

///@brief default ctor
AverageInterfaceEnergyFilter::AverageInterfaceEnergyFilter() :
	parent( "AverageInterfaceEnergy" ),
	task_factory_( NULL ),
	scorefxn_( NULL ),
	threshold_( 100000 ),
	bb_bb_( false ),
	sym_dof_names_( "" ),
	jump_( 0 ),
	return_total_( false )
{}

core::pack::task::TaskFactoryOP
AverageInterfaceEnergyFilter::task_factory() const{
	return task_factory_;
}

void
AverageInterfaceEnergyFilter::task_factory( core::pack::task::TaskFactoryOP task_factory ){
	task_factory_ = task_factory;
}

core::scoring::ScoreFunctionOP
AverageInterfaceEnergyFilter::scorefxn() const{
	return scorefxn_;
}

void
AverageInterfaceEnergyFilter::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::scoring::ScoreType
AverageInterfaceEnergyFilter::score_type() const{
	return score_type_;
}

void
AverageInterfaceEnergyFilter::score_type( core::scoring::ScoreType const st ){
	score_type_ = st;
}

core::Real
AverageInterfaceEnergyFilter::threshold() const{
	return threshold_;
}

void
AverageInterfaceEnergyFilter::threshold( core::Real const thresh ){
	threshold_ = thresh;
}

bool
AverageInterfaceEnergyFilter::bb_bb() const{
	return bb_bb_;
}

void
AverageInterfaceEnergyFilter::bb_bb( bool const bb ){
	bb_bb_ = bb;
}

bool AverageInterfaceEnergyFilter::unbound() const { return unbound_; }
void AverageInterfaceEnergyFilter::unbound( bool const unbound ) { unbound_ = unbound; }
std::string AverageInterfaceEnergyFilter::sym_dof_names() const { return sym_dof_names_; }
void AverageInterfaceEnergyFilter::sym_dof_names( std::string const sym_dofs ) { sym_dof_names_ = sym_dofs; }
core::Size AverageInterfaceEnergyFilter::jump() const { return jump_; }
void AverageInterfaceEnergyFilter::jump( core::Size const jump ) { jump_ = jump; }
bool AverageInterfaceEnergyFilter::return_total() const { return return_total_; }
void AverageInterfaceEnergyFilter::return_total( bool const total ) { return_total_ = total; }
std::string AverageInterfaceEnergyFilter::score_type_name() const { return score_type_name_; }
void AverageInterfaceEnergyFilter::score_type_name( std::string const name ) { score_type_name_ = name; }

bool
AverageInterfaceEnergyFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real aie(compute( pose ));
	if( aie <= threshold_) {
		TR<<"passing."<<std::endl;
		return true;
	} 
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

core::Real
AverageInterfaceEnergyFilter::compute( core::pose::Pose const & pose ) const{
	runtime_assert( task_factory() );
	core::pack::task::PackerTaskCOP packer_task( task_factory()->create_task_and_apply_taskoperations( pose ) );
	core::Size total_residue;
	if(core::pose::symmetry::is_symmetric( pose )) { 
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		total_residue = symm_info->num_independent_residues();
	} else {
		total_residue = pose.total_residue(); 
	}
	std::string interface_pos("interface positions considered: ");
	core::Real tie=0;
	core::Real interface_energy = 0;
	core::Size interface_size=0;
	core::pose::Pose p = pose;

	// If unbound energies are being evaluated, create the unbound state
	if ( unbound() ) {
	  int sym_aware_jump_id = 0;
	  if ( sym_dof_names() != "" ) {
	    utility::vector1<std::string> sym_dof_name_list;
	    sym_dof_name_list = utility::string_split( sym_dof_names() , ',' );
	    for (Size i = 1; i <= sym_dof_name_list.size(); i++) {
	      sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name_list[i] );
	      protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( p, sym_aware_jump_id ) );
	      translate->step_size( 1000.0 );
	      translate->apply( p );
	    }
  	} else if ( jump() != 0 ) {
	    sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( pose, jump() );
	    protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( p, sym_aware_jump_id ) );
	    translate->step_size( 1000.0 );
		  translate->apply( p );
  	} else {
			utility_exit_with_message( "If unbound is set to true, you must provide either sym_dof_names or a jump in order to create the unbound pose!" );
		}
	}

	scorefxn()->score( p );
  core::scoring::EnergyMap em;
	for( core::Size resi=1; resi<=total_residue; ++resi ){
		if( packer_task->being_packed( resi ) ) {
			interface_size += 1;
			interface_pos.append(ObjexxFCL::string_of(resi) + "+");
			interface_energy += p.energies().residue_total_energy( resi );
    	em += p.energies().residue_total_energies( resi );
		}
	}
	if ( score_type_name() == "total_score" ) {
		tie = interface_energy;
	} else {
	  em *= scorefxn()->weights();
		tie = em[score_type()];
	}
	core::Real aie = tie / interface_size;
	TR.Debug <<interface_pos<<std::endl;
	if ( return_total() )
		return ( tie );
	return( aie );
}

core::Real
AverageInterfaceEnergyFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real aie(compute( pose ));
	return( aie );
}

void
AverageInterfaceEnergyFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"AverageInterfaceEnergyFilter returns "<<compute( pose )<<std::endl;
}

void
AverageInterfaceEnergyFilter::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & )
{
	TR << "AverageInterfaceEnergyFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	scorefxn(protocols::rosetta_scripts::parse_score_function(tag, data));
	score_type(core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) ) );
	score_type_name( tag->getOption<std::string>( "score_type", "total_score" ) );
	threshold( tag->getOption< core::Real >( "cutoff", 100000 ) );
	bb_bb( tag->getOption< bool >( "bb_bb", false ) );
	unbound( tag->getOption< bool >( "unbound", false ) );
	sym_dof_names( tag->getOption< std::string >( "sym_dof_names", "" ) );
	jump( tag->getOption< core::Size >( "jump", 0 ) );
	return_total( tag->getOption< bool >( "return_total", false ) );
	TR<<"with options scoretype: "<<score_type()<<" and cutoff: "<<threshold()<<std::endl;
}
void AverageInterfaceEnergyFilter::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & score_fxns,
				utility::lua::LuaObject const & tasks ) {
	TR << "AverageInterfaceEnergyFilter"<<std::endl;
	task_factory( protocols::elscripts::parse_taskdef( def["tasks"], tasks ));
	if( def["scorefxn"] ) {
		scorefxn( protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns ) );
	} else {
		scorefxn( score_fxns["score12"].to<core::scoring::ScoreFunctionSP>()->clone()  );
	}
	score_type( core::scoring::score_type_from_name( def["score_type"] ? def["score_type"].to<std::string>() : "total_score" ) );
	threshold( def["cutoff"] ? def["cutoff"].to<core::Real>() : 100000 );
	bb_bb( def["bb_bb"] ? def["bb_bb"].to<bool>() : false );
	TR<<"with options scoretype: "<<score_type()<<" and cutoff: "<<threshold()<<std::endl;
}

protocols::filters::FilterOP
AverageInterfaceEnergyFilter::fresh_instance() const{
	return new AverageInterfaceEnergyFilter();
}

AverageInterfaceEnergyFilter::~AverageInterfaceEnergyFilter(){}

protocols::filters::FilterOP
AverageInterfaceEnergyFilter::clone() const{
	return new AverageInterfaceEnergyFilter( *this );
}

protocols::filters::FilterOP
AverageInterfaceEnergyFilterCreator::create_filter() const { return new AverageInterfaceEnergyFilter; }

std::string
AverageInterfaceEnergyFilterCreator::keyname() const { return "AverageInterfaceEnergy"; }

} // matdes
} // devel
