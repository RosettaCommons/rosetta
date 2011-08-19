// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/SequenceRecoveryFilter.hh>
#include <protocols/protein_interface_design/filters/SequenceRecoveryFilterCreator.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/protein_interface_design/dock_design_filters.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/pdb/pose_io.hh>
#include <protocols/protein_interface_design/design_utils.hh>

//Auto Headers

namespace protocols {
namespace protein_interface_design{
namespace filters {

static basic::Tracer TR( "protocols.protein_interface_design.filters.SequenceRecoveryFilter" );

///@brief default ctor
SequenceRecoveryFilter::SequenceRecoveryFilter() :
	parent( "SequenceRecovery" ),
	task_factory_( NULL ),
	reference_pose_( NULL ),
	rate_threshold_( 0.0 )
{}

core::pack::task::TaskFactoryOP
SequenceRecoveryFilter::task_factory() const
{
	return task_factory_;
}

void
SequenceRecoveryFilter::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

core::Real
SequenceRecoveryFilter::rate_threshold() const
{
	return( rate_threshold_ );
}

void
SequenceRecoveryFilter::rate_threshold( core::Real const rate )
{
	runtime_assert( rate_threshold_ >= 0 && rate_threshold_ <= 1.0 );
	rate_threshold_ = rate;
}

core::pose::PoseCOP
SequenceRecoveryFilter::reference_pose() const
{
	return reference_pose_;
}

void
SequenceRecoveryFilter::reference_pose( core::pose::PoseCOP pose )
{
	reference_pose_ = pose;
}

void
SequenceRecoveryFilter::reference_pose( core::pose::Pose const & pose )
{
	reference_pose_ = new core::pose::Pose( pose );
}

bool
SequenceRecoveryFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const recovery_rate( compute( pose ) );
	TR<<"Sequence recovery rate evaluates to "<<recovery_rate<<". ";
	if( recovery_rate <= rate_threshold_ ){
		TR<<"Failing."<<std::endl;
		return false;
	}
	TR<<"Success."<<std::endl;
	return true;
}

core::Real
SequenceRecoveryFilter::compute( core::pose::Pose const & pose ) const{
	runtime_assert( task_factory() );
	runtime_assert( reference_pose() );
	if( reference_pose()->total_residue() != pose.total_residue() )
		utility_exit_with_message( "Reference pose and current pose have a different number of residues" );
	core::pack::task::PackerTaskOP packer_task( task_factory_->create_task_and_apply_taskoperations( pose ) );
	core::Size designable_count( 0 );
	for( core::Size resi=1; resi<=pose.total_residue(); ++resi )
		if( packer_task->being_designed( resi ) ) ++designable_count;

	if( !designable_count )
		utility_exit_with_message( "No designable residues identified in pose. Are you sure you have set the correct task operations?" );

	using namespace core::scoring;
  protocols::protein_interface_design::ReportSequenceDifferences rsd( ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH ) );
  rsd.calculate( *reference_pose(), pose );
  std::map< core::Size, std::string > const res_names1( rsd.res_name1() );
  core::Size const mutated( res_names1.size() );
  core::Real const rate( 1.0 - (core::Real) mutated / designable_count );
  TR<<"Your design mover mutated "<<mutated<<" positions out of "<<designable_count<<" designable positions. Sequence recovery is: "<<rate<<std::endl;
	return( rate );
}

core::Real
SequenceRecoveryFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose ) );
}

void
SequenceRecoveryFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"SequenceRecoveryFilter returns "<<compute( pose )<<std::endl;
}

void
SequenceRecoveryFilter::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & pose )
{
	TR << "SequenceRecoveryFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	rate_threshold( tag->getOption< core::Real >( "rate_threshold", 0.0 ) );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if( option[ in::file::native ].user() ){
		std::string const reference_pdb = option[ in::file::native ]();
		core::pose::PoseOP temp_pose( new core::pose::Pose );
		core::import_pose::pose_from_pdb( *temp_pose, reference_pdb );
		reference_pose( temp_pose );
		TR<<"Using native pdb "<<reference_pdb<<" as reference.";
	}
	else{
		TR<<"Using starting pdb as reference. You could use -in::file::native to specify a different pdb for reference";
		reference_pose( pose );
	}
	TR<<std::endl;
}

protocols::filters::FilterOP
SequenceRecoveryFilter::fresh_instance() const{
	return new SequenceRecoveryFilter();
}

SequenceRecoveryFilter::~SequenceRecoveryFilter(){}


protocols::filters::FilterOP
SequenceRecoveryFilter::clone() const{
	return new SequenceRecoveryFilter( *this );
}

protocols::filters::FilterOP
SequenceRecoveryFilterCreator::create_filter() const { return new SequenceRecoveryFilter; }

std::string
SequenceRecoveryFilterCreator::keyname() const { return "SequenceRecovery"; }


} // filters
} // protein_interface_design
} // protocols
