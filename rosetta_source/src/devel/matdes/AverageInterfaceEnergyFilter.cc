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
	bb_bb_( false )
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
	core::Size interface_size=0;
	protocols::simple_filters::EnergyPerResidueFilter epr;
	epr.scorefxn(scorefxn());
	epr.score_type(score_type());
	epr.threshold(threshold());
	epr.bb_bb(bb_bb());
	for( core::Size resi=1; resi<=total_residue; ++resi ){
		if( packer_task->being_packed( resi ) ) {
			interface_size += 1;
			epr.resnum(resi);
			tie += epr.compute( pose ); 
			interface_pos.append(ObjexxFCL::string_of(resi) + "+");   
		}
	}
	core::Real aie = tie / interface_size;
	TR.Debug <<interface_pos<<std::endl;
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
	std::string const scorefxn_name( tag->getOption< std::string >( "scorefxn", "score12" ) );
	scorefxn( data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ) );
	score_type(core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) ) );
	threshold( tag->getOption< core::Real >( "cutoff", 100000 ) );
	bb_bb( tag->getOption< bool >( "bb_bb", false ) );
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
