// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/matdes/TaskAwareAlaScan.cc
/// @brief
/// @author Neil King (neilking@uw.edu)
// Project Headers

#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <map>
#include <numeric/random/random.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <devel/matdes/TaskAwareAlaScan.hh>
#include <devel/matdes/TaskAwareAlaScanCreator.hh>


static basic::Tracer TR( "devel.matdes.TaskAwareAlaScan" );

namespace devel {
namespace matdes {

// @brief default constructor
TaskAwareAlaScan::TaskAwareAlaScan():
  task_factory_( NULL ),
  jump_( 1 ),
  repeats_( 3 ),
  scorefxn_( NULL ),
	repack_( 1 ),
	report_diffs_( 0 )
{}

// @brief constructor with arguments
TaskAwareAlaScan::TaskAwareAlaScan(
	core::pack::task::TaskFactoryOP task_factory,
	core::Size jump,
	core::Size const repeats,
	core::scoring::ScoreFunctionCOP /*scorefxn*/,
	bool repack,
	bool report_diffs
):
		Filter( "TaskAwareAlaScan" ),
		task_factory_( task_factory ),
		jump_( jump ),
		repeats_( repeats ),
		repack_( repack ),
		report_diffs_( report_diffs )
{}

// @brief copy constructor
TaskAwareAlaScan::TaskAwareAlaScan( TaskAwareAlaScan const & rval ):
		Filter( rval ),
		task_factory_( rval.task_factory_ ),
		jump_( rval.jump_ ),
		repeats_( rval.repeats_ ),
		scorefxn_( rval.scorefxn_ ),
		repack_( rval.repack_ ),
		report_diffs_( rval.report_diffs_ )
{}

protocols::filters::FilterOP TaskAwareAlaScan::clone() const { return new TaskAwareAlaScan( *this ); }
protocols::filters::FilterOP TaskAwareAlaScan::fresh_instance() const { return new TaskAwareAlaScan(); }

// @brief setters
void TaskAwareAlaScan::task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
void TaskAwareAlaScan::jump( core::Size const j ) { jump_ = j; }
void TaskAwareAlaScan::repeats( core::Size const r ) { repeats_ = r; }
void TaskAwareAlaScan::scorefxn( core::scoring::ScoreFunctionOP const scorefxn ) { scorefxn_ = scorefxn; }
void TaskAwareAlaScan::repack( bool const repack ) { repack_ = repack; }
void TaskAwareAlaScan::report_diffs( bool const report_diffs ) { report_diffs_ = report_diffs; }
void TaskAwareAlaScan::write2pdb( bool const write ) { write2pdb_ = write; }

// @brief getters
core::pack::task::TaskFactoryOP TaskAwareAlaScan::task_factory() const { return task_factory_; }
core::Size TaskAwareAlaScan::jump() const { return jump_; }
core::Size TaskAwareAlaScan::repeats() const { return repeats_; }
bool TaskAwareAlaScan::repack() const { return repack_; }
bool TaskAwareAlaScan::report_diffs() const { return report_diffs_; }
bool TaskAwareAlaScan::write2pdb() const { return write2pdb_; }

// @brief Dummy Filter apply function
bool TaskAwareAlaScan::apply( core::pose::Pose const & ) const { return true; }

// @brief parse xml
void
TaskAwareAlaScan::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/
){
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	jump( tag->getOption< Size >( "jump", 1 ) );
	repeats( tag->getOption< core::Size >( "repeats", 1 ) );
	std::string const scorefxn_name( tag->getOption< std::string >( "scorefxn", "score12" ));
	repack( tag->getOption< bool >( "repack", 1 ) );
	report_diffs( tag->getOption< bool >("report_diffs", 1) );
	write2pdb( tag->getOption< bool >("write2pdb", 0) );
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data, scorefxn_name );
  std::string unparsed_exempt_identities = tag->getOption< std::string >( "exempt_identities" );
  if( unparsed_exempt_identities != "" ){
    utility::vector1< std::string > const ids( utility::string_split( unparsed_exempt_identities , ',' ) );
    exempt_identities_.clear();
    foreach( std::string const id, ids ){
        exempt_identities_.insert( id );
    }
  }
}

// @brief Calculate the ddG for an alanine mutation at the specified position
core::Real
TaskAwareAlaScan::ddG_for_single_residue( core::pose::Pose const & const_pose, core::Size const resi ) const
{
  if( !const_pose.residue( resi ).is_protein() ){
    TR<<"WARNING: Non-protein residue "<< resi<<" was requested for ala-scan. Returning 0"<<std::endl;
    return 0.0;
  }
  core::Size const rb_jump( jump_ );
  core::pose::Pose pose( const_pose );

	// First, mutate the residue in question to alanine
  utility::vector1< bool > allowed_aas;
  allowed_aas.assign( core::chemical::num_canonical_aas, false );
  std::string mut_name = "ALA";
  if ( exempt_identities_.find( pose.residue( resi ).name3() ) != exempt_identities_.end() ) {
    allowed_aas[ pose.residue( resi ).aa() ] = true;
		mut_name = pose.residue( resi ).name3();
	} else {
    allowed_aas[ core::chemical::aa_ala ] = true;
	}
  using namespace core::pack::task;
  PackerTaskOP task = TaskFactory::create_packer_task( pose );
  task->initialize_from_command_line().or_include_current( true );
  for( core::Size resj=1; resj<=pose.total_residue(); ++resj ){
    if( resi == resj )
      task->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
    else
      task->nonconst_residue_task( resj ).prevent_repacking();
  }
	bool symmetric = 0;
	if (core::pose::symmetry::is_symmetric(pose)) {
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task);
		symmetric = 1;
	}
  core::pack::pack_rotamers( pose, *scorefxn_, task );

	// Create DdgFilter and ScoreTypeFilter
  protocols::simple_filters::DdgFilter ddg_filter( 10000, scorefxn_, rb_jump, 1, symmetric );
  if( repack() )
    TR << "Energy calculations are carried out with repacking in the bound and unbound states" << std::endl;
  else
    TR << "Energy calculations are carried out without repacking in the bound and unbound states" << std::endl;
  ddg_filter.repack( repack() );
  protocols::simple_filters::ScoreTypeFilter const energy_filter( scorefxn_, core::scoring::total_score, 0 );

	// Calculate either ddG or total score and return value
  core::Real accumulate_ddg = 0;
  for( core::Size r=1; r<=repeats_; ++r )
    accumulate_ddg += (rb_jump==0 ? energy_filter.compute( pose ) : ddg_filter.compute( pose ) );
  core::Real const mut_ddg( accumulate_ddg / repeats_ );

	TR << protocols::jd2::current_output_name() << " ala scan ddG for mutation " << const_pose.residue( resi ).name3() << resi << mut_name << " = " << mut_ddg << std::endl;

  TR.flush();
  return( mut_ddg );
}

// @brief calculate and report the per-residue ddGs
void
TaskAwareAlaScan::report( std::ostream & out, core::pose::Pose const & const_pose ) const
{

	core::Size const rb_jump( jump_ );
	core::pose::Pose pose( const_pose );
	bool symmetric = 0;
  if (core::pose::symmetry::is_symmetric(pose)) { symmetric = 1; }

	protocols::simple_filters::DdgFilter const ddg_filter( 10000, scorefxn_, rb_jump, 1, symmetric );
  protocols::simple_filters::ScoreTypeFilter const energy_filter( scorefxn_, core::scoring::total_score, 0 );

	// Get the baseline energy/ddG of the input structure
	core::Real accumulate_ddg( 0 );
	for( core::Size r=1; r<=repeats_; ++r )
		accumulate_ddg += (rb_jump==0 ? energy_filter.compute( const_pose ) : ddg_filter.compute( const_pose ) );
	core::Real const wt_ddg( accumulate_ddg / repeats_ );

	// Apply TaskOperations from xml
  core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
  if ( task_factory_ != 0 ) {
    task = task_factory_->create_task_and_apply_taskoperations( pose );
  } else {
    TR << "Warning: You have not provided any TaskOperations. A default will be used." << std::endl;
  }
	if (symmetric) {
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task);
	}

	// For packable residues, calculate the ddG/score upon mutation to alanine
  for (core::Size resi = 1; resi <= pose.n_residue(); resi++) {
		if( !pose.residue( resi ).is_protein() ) continue;
		if( task->pack_residue( resi ) ) {
			core::Real const mut_ddg( ddG_for_single_residue( const_pose, resi ) );
			core::Real const diff_ddg( mut_ddg - wt_ddg );

			core::pose::PDBInfoCOP pose_info( const_pose.pdb_info() );
			char const chain( pose_info->chain( resi ) );
			std::string const res_type( const_pose.residue( resi ).name3() );
			core::Real output_ddG = ( report_diffs() == 1 ) ? diff_ddg : mut_ddg;
			// Output to pdb file
			if ( write2pdb() ) { write_to_pdb( resi, res_type, output_ddG ); }
			// Output to tracer
			out<<" "<<res_type<<" "<<resi<<" "<<chain<<" : "<< ObjexxFCL::fmt::F (9,4,output_ddG)<<std::endl;
		}
	}
	out<<std::endl;
}

void TaskAwareAlaScan::write_to_pdb( core::Size const & residue, std::string const & residue_name, core::Real const & ddG ) const
{
  
  protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
  std::string filter_name = this->name();
  std::string user_name = this->get_user_defined_name();
	std::string mut_name = "ALA";
	if ( exempt_identities_.find( residue_name ) != exempt_identities_.end() ) { mut_name = residue_name; }
  std::string output_string = filter_name + " " + user_name + ": " + residue_name + ObjexxFCL::string_of(residue) + mut_name + " = " + ObjexxFCL::fmt::F (9,4,ddG);
  job->add_string(output_string);

}

protocols::filters::FilterOP
TaskAwareAlaScanCreator::create_filter() const { return new TaskAwareAlaScan; }

std::string
TaskAwareAlaScanCreator::keyname() const { return "TaskAwareAlaScan"; }

}
}
