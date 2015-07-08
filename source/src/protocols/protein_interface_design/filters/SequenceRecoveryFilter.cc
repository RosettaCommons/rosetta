// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/elscripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <ObjexxFCL/format.hh>
#include <protocols/protein_interface_design/design_utils.hh>
#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace protein_interface_design{
namespace filters {

static thread_local basic::Tracer TR( "protocols.protein_interface_design.filters.SequenceRecoveryFilter" );

/// @brief default ctor
SequenceRecoveryFilter::SequenceRecoveryFilter() :
	parent( "SequenceRecovery" ),
	task_factory_( /* NULL */ ),
	reference_pose_( /* NULL */ ),
	rate_threshold_( 0.0 ),
	mutation_threshold_( 100 ),
	mutations_( 0 ),
	verbose_( 0 )
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
	reference_pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) );
}

core::Size
SequenceRecoveryFilter::mutation_threshold() const
{
	return( mutation_threshold_ );
}

void
SequenceRecoveryFilter::mutation_threshold( core::Size const mut )
{
	mutation_threshold_ = mut;
}

bool
SequenceRecoveryFilter::mutations() const
{
	return( mutations_ );
}

void
SequenceRecoveryFilter::mutations( bool const muts )
{
	mutations_ = muts;
}

bool
SequenceRecoveryFilter::verbose() const
{
	return( verbose_ );
}

void
SequenceRecoveryFilter::verbose( bool const verb )
{
	verbose_ = verb;
}

bool
SequenceRecoveryFilter::write2pdb() const
{
	return( write2pdb_ );
}

void
SequenceRecoveryFilter::write2pdb( bool const write )
{
	write2pdb_ = write;
}

bool
SequenceRecoveryFilter::apply(core::pose::Pose const & pose ) const
{
	if ( mutations_ ) {
		core::Size const num_mutations( (core::Size) compute( pose, false ) );
		TR<<"The designed pose possesses "<<num_mutations<<" compared to the reference pose. ";
		if( num_mutations <= mutation_threshold_ ){
			TR<<"Success."<<std::endl;
			return true;
		} else {
			TR<<"Failing."<<std::endl;
			return false;
		}		
	} else {
		core::Real const recovery_rate( compute( pose, false ) );
		TR<<"Sequence recovery rate evaluates to "<<recovery_rate<<". ";
		if( recovery_rate <= rate_threshold_ ){
			TR<<"Failing."<<std::endl;
			return false;
		} else {
		TR<<"Success."<<std::endl;
			return true;
		}
	}
}
void SequenceRecoveryFilter::apply( core::io::serialization::PipeMap & pmap ) {
	if( pmap["reference"]->size() == 1 )
		reference_pose( *(*pmap["reference"])[0] );
	Filter::apply( pmap );
}

core::Real
SequenceRecoveryFilter::compute( core::pose::Pose const & pose, bool const & write ) const{
	runtime_assert( task_factory() != 0 );
	runtime_assert( reference_pose() != 0 );
	core::Size total_residue_ref;
	core::pose::Pose asym_ref_pose;
	if(core::pose::symmetry::is_symmetric( *reference_pose() )) { 
		core::pose::symmetry::extract_asymmetric_unit( *reference_pose(), asym_ref_pose);
  	for (core::Size i = 1; i <= asym_ref_pose.total_residue(); ++i) {
    	if (asym_ref_pose.residue_type(i).name() == "VRT") {
				asym_ref_pose.conformation().delete_residue_slow(asym_ref_pose.total_residue());
			}
		}
		total_residue_ref = asym_ref_pose.total_residue(); 
	} else {
		total_residue_ref = reference_pose()->total_residue(); 
		asym_ref_pose = *reference_pose();
	} 
	core::Size total_residue;
	core::pose::Pose asym_pose;
	if (core::pose::symmetry::is_symmetric( pose )) { 
		core::pose::symmetry::extract_asymmetric_unit(pose, asym_pose);
  	for (core::Size i = 1; i <= asym_pose.total_residue(); ++i) {
    	if (asym_pose.residue_type(i).name() == "VRT") {
				asym_pose.conformation().delete_residue_slow(asym_pose.total_residue());
			}
		}
		total_residue = asym_pose.total_residue(); 
	} else {
		total_residue = pose.total_residue(); 
		asym_pose = pose;
	}
	if( total_residue_ref != total_residue )
		utility_exit_with_message( "Reference pose and current pose have a different number of residues" );
	core::pack::task::PackerTaskOP packer_task( task_factory_->create_task_and_apply_taskoperations( pose ) );
	core::Size designable_count( 0 );
	core::Size packable_count( 0 );
	for( core::Size resi=1; resi<=total_residue; ++resi ) {
		if( packer_task->being_designed( resi ) ) {
			designable_count++;
		}
		if( packer_task->being_packed( resi ) ) {
			packable_count++;
		}
	}
	if( !designable_count ) {
		TR<<"Warning: No designable residues identified in pose. Are you sure you have set the correct task operations?"<<std::endl;
		if( !packable_count ) {
			utility_exit_with_message("No designable or packable residues identified in pose. Are you sure you have set the correct task operations?" );
		}
	}
	using namespace core::scoring;
  protocols::protein_interface_design::ReportSequenceDifferences rsd( get_score_function() );
  rsd.calculate( asym_ref_pose, asym_pose );
  std::map< core::Size, std::string > const res_names1( rsd.res_name1() );
  std::map< core::Size, std::string > const res_names2( rsd.res_name2() );
  core::Size const mutated( res_names1.size() );
  // AMW: cppcheck notices that if there are packable residues but no designable residues
  // we divide by zero here.
  core::Real const rate( 1.0 - (core::Real) mutated / designable_count );
  TR<<"Your design mover mutated "<<mutated<<" positions out of "<<designable_count<<" designable positions. Sequence recovery is: "<<rate<<std::endl;
	if ( verbose_ ) {
		rsd.report( TR );
		TR.flush();
	}	
	if ( write ) {
		write_to_pdb( res_names1, res_names2 );
	}
	if ( mutations_ ) {
		return( (core::Real) mutated );
	} else {
		return( rate );
	}
}

/// @brief Add each mutation to the output pdb if desired
void
SequenceRecoveryFilter::write_to_pdb( std::map< core::Size, std::string > const & res_names1, std::map< core::Size, std::string > const & res_names2 ) const {

  protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
  std::string user_name = this->get_user_defined_name();
	std::map< Size, std::string >::const_iterator it_name1 = res_names1.begin();
	std::map< Size, std::string >::const_iterator it_name2 = res_names2.begin();
	while( it_name1 != res_names1.end() ) {
		std::string output_string = "SequenceRecoveryFilter " + user_name + ": " + it_name2->second + ObjexxFCL::string_of( it_name1->first ) + it_name1->second;
		job->add_string( output_string );
		++it_name1; ++it_name2;
 	}

}

core::Real
SequenceRecoveryFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose, write2pdb() ) );
}

void
SequenceRecoveryFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"SequenceRecoveryFilter returns "<<compute( pose, false )<<std::endl;
}

void
SequenceRecoveryFilter::parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & pose )
{
	TR << "SequenceRecoveryFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	rate_threshold( tag->getOption< core::Real >( "rate_threshold", 0.0 ) );
	mutation_threshold( tag->getOption< core::Size >( "mutation_threshold", 100 ) );
	mutations( tag->getOption< bool >( "report_mutations", 0 ) );
	verbose( tag->getOption< bool >( "verbose", 0 ) );
	write2pdb( tag->getOption< bool >( "write2pdb", 0 ) );

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

void SequenceRecoveryFilter::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & ,
				utility::lua::LuaObject const & tasks ) {
	TR << "SequenceRecoveryFilter"<<std::endl;
	task_factory( protocols::elscripts::parse_taskdef( def["tasks"], tasks ));
	rate_threshold( def["rate_threshold"] ? def["rate_threshold"].to<core::Real>() : 0.0 );
	mutation_threshold( def["mutation_threshold"] ? def["mutation_threshold"].to<core::Size>() : 100 );
	mutations( def["report_mutations"] ? def["report_mutations"].to<bool>() : false );
	verbose( def["verbose"] ? def["verbose"].to<bool>() : false );
	write2pdb( def["write2pdb"] ? def["write2pdb"].to<bool>() : false );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if( option[ in::file::native ].user() ){
		std::string const reference_pdb = option[ in::file::native ]();
		core::pose::PoseOP temp_pose( new core::pose::Pose );
		core::import_pose::pose_from_pdb( *temp_pose, reference_pdb );
		reference_pose( temp_pose );
		TR<<"Using native pdb "<<reference_pdb<<" as reference.";
	}
	TR<<std::endl;
}
protocols::filters::FilterOP
SequenceRecoveryFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new SequenceRecoveryFilter() );
}

SequenceRecoveryFilter::~SequenceRecoveryFilter(){}


protocols::filters::FilterOP
SequenceRecoveryFilter::clone() const{
	return protocols::filters::FilterOP( new SequenceRecoveryFilter( *this ) );
}

protocols::filters::FilterOP
SequenceRecoveryFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SequenceRecoveryFilter ); }

std::string
SequenceRecoveryFilterCreator::keyname() const { return "SequenceRecovery"; }


} // filters
} // protein_interface_design
} // protocols
