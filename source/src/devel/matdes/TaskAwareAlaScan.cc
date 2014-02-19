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
#include <core/chemical/ResidueConnection.hh>
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
#include <core/pack/task/operation/TaskOperation.hh>
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
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/elscripts/util.hh>
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

/// @brief default constructor
TaskAwareAlaScan::TaskAwareAlaScan():
		task_factory_( NULL ),
		ddG_task_factory_( NULL ),
		use_custom_task_( false ),
		jump_( 1 ),
		sym_dof_name_( "" ),
		repeats_( 3 ),
		scorefxn_( NULL ),
		repack_( 1 ),
		report_diffs_( 0 )
{}

/// @brief constructor with arguments
TaskAwareAlaScan::TaskAwareAlaScan(
	core::pack::task::TaskFactoryOP task_factory,
	core::pack::task::TaskFactoryOP ddG_task_factory,
	bool use_custom_task,
	core::Size jump,
	std::string sym_dof_name,
	core::Size const repeats,
	core::scoring::ScoreFunctionCOP /*scorefxn*/,
	bool repack,
	bool report_diffs
):
		Filter( "TaskAwareAlaScan" ),
		task_factory_( task_factory ),
		ddG_task_factory_( ddG_task_factory ),
		use_custom_task_( use_custom_task ),
		jump_( jump ),
		sym_dof_name_( sym_dof_name ),
		repeats_( repeats ),
		repack_( repack ),
		report_diffs_( report_diffs )
{}

/// @brief copy constructor
TaskAwareAlaScan::TaskAwareAlaScan( TaskAwareAlaScan const & rval ):
		Filter( rval ),
		task_factory_( rval.task_factory_ ),
		ddG_task_factory_( rval.ddG_task_factory_ ),
		use_custom_task_( rval.use_custom_task_ ),
		jump_( rval.jump_ ),
		sym_dof_name_( rval.sym_dof_name_ ),
		repeats_( rval.repeats_ ),
		scorefxn_( rval.scorefxn_ ),
		repack_( rval.repack_ ),
		report_diffs_( rval.report_diffs_ )
{}

/// @brief destructor
TaskAwareAlaScan::~TaskAwareAlaScan() {}

protocols::filters::FilterOP TaskAwareAlaScan::clone() const { return new TaskAwareAlaScan( *this ); }
protocols::filters::FilterOP TaskAwareAlaScan::fresh_instance() const { return new TaskAwareAlaScan(); }

// setters
void TaskAwareAlaScan::task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
void TaskAwareAlaScan::ddG_task_factory( core::pack::task::TaskFactoryOP ddG_task_factory ) { ddG_task_factory_ = ddG_task_factory; }
void TaskAwareAlaScan::use_custom_task( bool const uct ) { use_custom_task_ = uct; }
void TaskAwareAlaScan::jump( core::Size const j ) { jump_ = j; }
void TaskAwareAlaScan::sym_dof_name( std::string const j ) { sym_dof_name_ = j; }
void TaskAwareAlaScan::repeats( core::Size const r ) { repeats_ = r; }
void TaskAwareAlaScan::scorefxn( core::scoring::ScoreFunctionOP const scorefxn ) { scorefxn_ = scorefxn; }
void TaskAwareAlaScan::repack( bool const repack ) { repack_ = repack; }
void TaskAwareAlaScan::report_diffs( bool const report_diffs ) { report_diffs_ = report_diffs; }
void TaskAwareAlaScan::write2pdb( bool const write ) { write2pdb_ = write; }

// getters
core::pack::task::TaskFactoryOP TaskAwareAlaScan::task_factory() const { return task_factory_; }
core::pack::task::TaskFactoryOP TaskAwareAlaScan::ddG_task_factory() const { return ddG_task_factory_; }
bool TaskAwareAlaScan::use_custom_task() const { return use_custom_task_; }
core::Size TaskAwareAlaScan::jump() const { return jump_; }
std::string TaskAwareAlaScan::sym_dof_name() const { return sym_dof_name_; }
core::Size TaskAwareAlaScan::repeats() const { return repeats_; }
bool TaskAwareAlaScan::repack() const { return repack_; }
bool TaskAwareAlaScan::report_diffs() const { return report_diffs_; }
bool TaskAwareAlaScan::write2pdb() const { return write2pdb_; }

/// @brief Dummy Filter apply function
bool TaskAwareAlaScan::apply( core::pose::Pose const & ) const { return true; }

/// @brief parse xml
void
TaskAwareAlaScan::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/
)
{
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	jump( tag->getOption< Size >( "jump", 1 ) );
	sym_dof_name( tag->getOption< std::string >( "sym_dof_name", "" ) );
	repeats( tag->getOption< core::Size >( "repeats", 1 ) );
	repack( tag->getOption< bool >( "repack", 1 ) );
	report_diffs( tag->getOption< bool >("report_diffs", 1) );
	write2pdb( tag->getOption< bool >("write2pdb", 0) );
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	// Handle exempt identities
  std::string unparsed_exempt_identities = tag->getOption< std::string >( "exempt_identities", "" );
  if( unparsed_exempt_identities != "" ) {
    utility::vector1< std::string > const ids( utility::string_split( unparsed_exempt_identities , ',' ) );
    exempt_identities_.clear();
    foreach( std::string const id, ids ) {
        exempt_identities_.insert( id );
    }
  }
	// Handle task_operations for ddG repacking
	use_custom_task( tag->hasOption("ddG_task_operations") );
	if ( use_custom_task() ) {
  	parse_ddG_task_operations( tag, data );
	}
}
void TaskAwareAlaScan::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & score_fxns,
				utility::lua::LuaObject const & tasks ) {
	task_factory( protocols::elscripts::parse_taskdef( def["tasks"], tasks ));
	jump( def["jump"] ? def["jump"].to<core::Size>() : 1 );
	repeats( def["repeats"] ? def["repeats"].to<core::Size>() : 1 );
	if( def["scorefxn"] ) {
		scorefxn_ = protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns );
	} else {
		scorefxn_ = score_fxns["score12"].to<core::scoring::ScoreFunctionSP>()->clone();
	}
	repack( def["repack"] ? def["repack"].to<bool>() : 1 );
	report_diffs( def["report_diffs"] ? def["report_diffs"].to<bool>() : 1 );
	write2pdb( def["write2pdb"] ? def["write2pdb"].to<bool>() : 0 );
	if( def["exempt_identites"] ) {
    exempt_identities_.clear();
		for (utility::lua::LuaIterator i=def["exempt_identities"].begin(), end; i != end; ++i) {
			exempt_identities_.insert( (*i).to<std::string>() );
		}
	}
	// i dont think you need all this use_custom_task stuff now....
	use_custom_task( def["ddG_task_operations"] );
	if( use_custom_task() ) {
		ddG_task_factory( protocols::elscripts::parse_taskdef( def["ddG_tasks"], tasks ));
		// that was easy!
	}
}

/// @brief Calculate the ddG for an alanine mutation at the specified position
core::Real
TaskAwareAlaScan::ddG_for_single_residue( core::pose::Pose const & const_pose, core::Size const resi ) const
{
	if( !const_pose.residue( resi ).is_protein() ){
		TR<<"WARNING: Non-protein residue "<< resi<<" was requested for ala-scan. Returning 0"<<std::endl;
		return 0.0;
	}
	core::pose::Pose pose( const_pose );
	core::Size sym_aware_jump_id = 0;
	if ( sym_dof_name() != "" ) {
		sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name() );
	} else {
		sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( pose, jump() );
	}

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

	if (core::pose::symmetry::is_symmetric(pose)) {
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task);
	}
	core::pack::pack_rotamers( pose, *scorefxn_, task );

	// Create ddG mover, calculate, and return value
	protocols::simple_moves::ddG ddG( scorefxn_, sym_aware_jump_id );
	if ( use_custom_task() ) {
		ddG.use_custom_task( use_custom_task() );
		ddG.task_factory( ddG_task_factory() );
	}
	core::Real average( 0.0 );
	for( core::Size i=1; i<=repeats(); ++i ) {
		ddG.calculate( pose );
		average += ddG.sum_ddG();
		ddG.report_ddG( TR );
	}
	core::Real const mut_ddg( average / (core::Real)repeats() );

	TR << protocols::jd2::current_output_name() << " ala scan ddG for mutation " <<
			const_pose.residue( resi ).name3() << resi << mut_name << " = " << mut_ddg << std::endl;

	TR.flush();
	return( mut_ddg );
}

/// @brief calculate and report the per-residue ddGs
void
TaskAwareAlaScan::report( std::ostream & out, core::pose::Pose const & const_pose ) const
{

	// Set up, get data
	core::pose::Pose pose( const_pose );
	bool symmetric = 0;
  if ( core::pose::symmetry::is_symmetric(pose) ) symmetric = 1;
	core::Size sym_aware_jump_id = 0;
	if ( sym_dof_name() != "" ) {
		sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name() );
	} else {
		sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( pose, jump() );
	}

	// Calculate the ddG of the input pose
	protocols::simple_moves::ddG ddG( scorefxn_, sym_aware_jump_id /*, symmetric*/ ); //ddG autodetects symmetry now
	if ( use_custom_task() ) {
		ddG.use_custom_task( use_custom_task() );
		ddG.task_factory( ddG_task_factory() );
	}
	core::Real average( 0.0 );
	for( core::Size i=1; i<=repeats(); ++i ) {
		ddG.calculate( pose );
		average += ddG.sum_ddG();
		ddG.report_ddG( TR );
	}
	core::Real const wt_ddg( average / (core::Real)repeats() );

	// Apply TaskOperations from xml to pick which residues you'll be alanine scanning
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
			out<<" "<<res_type<<" "<<resi<<" "<<chain<<" : "<< ObjexxFCL::format::F (9,4,output_ddG)<<std::endl;
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
  std::string output_string = filter_name + " " + user_name + ": " + residue_name + ObjexxFCL::string_of(residue) + mut_name + " = " + ObjexxFCL::format::F (9,4,ddG);
  job->add_string(output_string);

}

/// @details Adapted from protocols/rosetta_scripts/util.cc
/// Since RosettaScripts can only return one TaskFactory from one set of task_operations,
/// a separate function is needed to provide an auxiliary TaskFactory that can be
/// passed along to the ddG mover to control repacking during ddG calculations.
/// @note Do NOT use STMStoredTasks generated by the StoreTaskMover as ddG_task_operations.
/// These will prevent mutation of the residue in question to alanine.
void
TaskAwareAlaScan::parse_ddG_task_operations( utility::tag::TagCOP const tag, basic::datacache::DataMap const & data )
{

  core::pack::task::TaskFactoryOP new_task_factory = new core::pack::task::TaskFactory;
  std::string const t_o_val( tag->getOption<std::string>("ddG_task_operations") );
  typedef utility::vector1< std::string > StringVec;
  StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
  TR<<"Passing the following task operations to ddG mover from "<<tag->getName()<<" called "<<tag->getOption<std::string>( "name", "no_name" )<<":\n";
  for ( StringVec::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
        t_o_key != end; ++t_o_key ) {
    if ( data.has( "task_operations", *t_o_key ) ) {
      new_task_factory->push_back( data.get< core::pack::task::operation::TaskOperation * >( "task_operations", *t_o_key ) );
      TR<<*t_o_key<<' ';
    } else {
      utility_exit_with_message("TaskOperation " + *t_o_key + " not found in basic::datacache::DataMap.");
    }
  }
  TR<<std::endl;
	ddG_task_factory( new_task_factory );
}

protocols::filters::FilterOP
TaskAwareAlaScanCreator::create_filter() const { return new TaskAwareAlaScan; }

std::string
TaskAwareAlaScanCreator::keyname() const { return "TaskAwareAlaScan"; }

}
}
