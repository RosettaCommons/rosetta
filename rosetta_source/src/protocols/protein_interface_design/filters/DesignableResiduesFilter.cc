// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/DesignableResiduesFilter.hh>
#include <protocols/protein_interface_design/filters/DesignableResiduesFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design{
namespace filters {

static basic::Tracer TR( "protocols.protein_interface_design.filters.DesignableResiduesFilter" );

///@brief default ctor
DesignableResiduesFilter::DesignableResiduesFilter() :
	parent( "DesignableResidues" ),
	task_factory_( NULL ),
	packable_( false ),
	designable_( true )
{}

bool
DesignableResiduesFilter::packable() const{
	return packable_;
}

bool
DesignableResiduesFilter::designable() const{
	return designable_;
}

void
DesignableResiduesFilter::designable( bool const d ){
	designable_ = d;
}

void
DesignableResiduesFilter::packable( bool const p ){
	packable_ = p;
}

core::pack::task::TaskFactoryOP
DesignableResiduesFilter::task_factory() const
{
	return task_factory_;
}

void
DesignableResiduesFilter::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

bool
DesignableResiduesFilter::apply(core::pose::Pose const & pose ) const
{
	compute( pose );
	return( true );
}

core::Real
DesignableResiduesFilter::compute( core::pose::Pose const & pose ) const{
	runtime_assert( task_factory() );
	runtime_assert( packable() || designable() );
	core::pack::task::PackerTaskCOP packer_task( task_factory()->create_task_and_apply_taskoperations( pose ) );
	if( designable() ){
		TR<<"Designable residues:"<<std::endl;
		for( core::Size resi=1; resi<=pose.total_residue(); ++resi ){
			if( packer_task->being_designed( resi ) )
				TR<<pose.residue( resi ).name3()<<" "<< pose.pdb_info()->number( resi )<<std::endl;
		}
	}
	if( packable() ){
		TR<<"Repackable residues:"<<std::endl;
		for( core::Size resi=1; resi<=pose.total_residue(); ++resi ){
			if( packer_task->being_packed( resi ) )
				TR<<pose.residue( resi ).name3()<<" "<<pose.pdb_info()->number( resi )<<std::endl;
		}
	}
	return( 0.0 );
}

core::Real
DesignableResiduesFilter::report_sm( core::pose::Pose const & ) const
{
	return( 0.0 );
}

void
DesignableResiduesFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"DesignableResiduesFilter returns "<<compute( pose )<<std::endl;
}

void
DesignableResiduesFilter::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & )
{
	TR << "DesignableResiduesFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	packable( tag->getOption< bool >( "packable", false ) );
	designable( tag->getOption< bool >( "designable", true ) );
	runtime_assert( designable() || packable() );
	TR<<"with options designable: "<<designable()<<" and repackable "<<packable()<<std::endl;
}

protocols::filters::FilterOP
DesignableResiduesFilter::fresh_instance() const{
	return new DesignableResiduesFilter();
}

DesignableResiduesFilter::~DesignableResiduesFilter(){}

protocols::filters::FilterOP
DesignableResiduesFilter::clone() const{
	return new DesignableResiduesFilter( *this );
}

protocols::filters::FilterOP
DesignableResiduesFilterCreator::create_filter() const { return new DesignableResiduesFilter; }

std::string
DesignableResiduesFilterCreator::keyname() const { return "DesignableResidues"; }

} // filters
} // protein_interface_design
} // protocols
