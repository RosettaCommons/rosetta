// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/FilterScan.hh>
#include <protocols/protein_interface_design/filters/FilterScanCreator.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/pose/PDBInfo.hh>
#include <fstream>
#include <utility/file/FileName.hh>
#include <iostream>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
// Auto-header: duplicate removed #include <protocols/rosetta_scripts/util.hh>
#include <protocols/jd2/util.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/jd2/JobDistributor.hh>

namespace protocols {
namespace protein_interface_design{
namespace filters {

static basic::Tracer TR( "protocols.protein_interface_design.filters.FilterScanFilter" );
static basic::Tracer TR_residue_scan( "ResidueScan" );

///@brief default ctor
FilterScanFilter::FilterScanFilter() :
	parent( "FilterScan" ),
	task_factory_( NULL ),
	triage_filter_( NULL ),
	filter_( NULL ),
	resfile_general_property_( "nataa" ),
	relax_mover_( NULL ),
	scorefxn_( NULL ),
	delta_( false ),
	unbound_( false ),
	report_all_( false ),
	jump_( 0 ),
	dump_pdb_( false )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string temp_resfile_name( protocols::jd2::current_output_name() );
	temp_resfile_name = temp_resfile_name + ".resfile";
	resfile_name( temp_resfile_name );
}

bool
FilterScanFilter::delta() const{
	return delta_;
}

void
FilterScanFilter::delta( bool const d ){
	delta_ = d;
}

bool
FilterScanFilter::unbound() const{
	return unbound_;
}

void
FilterScanFilter::unbound( bool const u ){
	unbound_ = u;
}

bool
FilterScanFilter::report_all() const{
	return report_all_;
}

void
FilterScanFilter::report_all( bool const ra ){
	report_all_ = ra;
}

core::Size
FilterScanFilter::jump() const{
	return jump_;
}

void
FilterScanFilter::jump( core::Size const j ){
	jump_ = j;
}

protocols::moves::MoverOP
FilterScanFilter::relax_mover() const{
	return relax_mover_;
}

void
FilterScanFilter::relax_mover( protocols::moves::MoverOP mover ){
	relax_mover_ = mover;
}

protocols::filters::FilterOP
FilterScanFilter::triage_filter() const{
	return triage_filter_;
}

void
FilterScanFilter::triage_filter( protocols::filters::FilterOP filter ){
	triage_filter_ = filter;
}
protocols::filters::FilterOP
FilterScanFilter::filter() const{
	return filter_;
}

void
FilterScanFilter::filter( protocols::filters::FilterOP filter ){
	filter_ = filter;
}

std::string
FilterScanFilter::resfile_general_property() const{
	return resfile_general_property_;
}

void
FilterScanFilter::resfile_general_property( std::string const resfile_general_property ){
	resfile_general_property_ = resfile_general_property;
}

std::string
FilterScanFilter::resfile_name() const{
	return resfile_name_;
}

void
FilterScanFilter::resfile_name( std::string const resfile_name ){
	resfile_name_ = resfile_name;
}

core::pack::task::TaskFactoryOP
FilterScanFilter::task_factory() const
{
	return task_factory_;
}

void
FilterScanFilter::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

void
FilterScanFilter::unbind( core::pose::Pose & pose ) const{
	if( !unbound() ) return;
	protocols::moves::RigidBodyTransMover rbtm( pose, jump() );
	rbtm.step_size( 10000.0 );
	rbtm.apply( pose );
}

bool
FilterScanFilter::apply(core::pose::Pose const & p ) const
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	core::pose::Pose pose_orig( p );

	PackerTaskCOP task = task_factory()->create_task_and_apply_taskoperations( pose_orig );
	utility::vector1< core::Size > being_designed;
	being_designed.clear();

	for( core::Size resi = 1; resi <= pose_orig.total_residue(); ++resi ){
		if( task->residue_task( resi ).being_designed() && pose_orig.residue(resi).is_protein() )
			being_designed.push_back( resi );
	}
	if( being_designed.empty() ) {
		TR.Warning << "WARNING: No residues are listed as designable." << std::endl;
		return true;
	}
	std::map< core::Size, utility::vector1< AA > > residue_id_map;
	std::map< std::pair< core::Size, AA >, std::pair< core::Real, bool > > residue_id_val_map; // position, aa identity : report filter value, triage filter accepted?
	residue_id_map.clear(); residue_id_val_map.clear();
	unbind( pose_orig );
	core::pose::Pose pose( pose_orig );
	//compute baseline
	TR<<"Computing baseline filter value\n";
	PackerTaskOP repack = task_factory()->create_task_and_apply_taskoperations( pose );
	repack->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
  core::pack::pack_rotamers( pose, *scorefxn(), repack );
	pose_orig = pose;
	relax_mover()->apply( pose );
	core::Real const baseline( filter()->report_sm( pose ) );
	TR<<"Starting pose's filter value (baseline) is: "<<baseline<<std::endl;
	pose = pose_orig;

	foreach( core::Size const resi, being_designed ){
  	 typedef std::list< ResidueTypeCAP > ResidueTypeCAPList;
		 ResidueTypeCAPList const & allowed( task->residue_task( resi ).allowed_residue_types() );
		 utility::vector1< AA > allow_temp;
		 allow_temp.clear();
		 foreach( ResidueTypeCAP const t, allowed ){
		 	allow_temp.push_back( t->aa() );
		 }
		foreach( AA const target_aa, allow_temp ){
			 pose = pose_orig;
    	 utility::vector1< bool > allowed_aas;
    	 allowed_aas.clear();
    	 allowed_aas.assign( num_canonical_aas, false );
    	 allowed_aas[ target_aa ] = true;
			 core::pack::task::TaskFactoryOP mut_res = new core::pack::task::TaskFactory( *task_factory() );
			 protocols::toolbox::task_operations::DesignAroundOperationOP dao = new protocols::toolbox::task_operations::DesignAroundOperation;///restrict repacking to 8.0A around target res to save time
			 dao->design_shell( 8.0 );
			 dao->include_residue( resi );
			 dao->repack_on( false );
			 mut_res->push_back( dao );
    	 PackerTaskOP mutate_residue = mut_res->create_task_and_apply_taskoperations( pose );
    	 mutate_residue->initialize_from_command_line().or_include_current( true );
    	 for( core::Size resj = 1; resj <= pose.total_residue(); ++resj ){
    	   if( resi != resj )
    	     mutate_residue->nonconst_residue_task( resj ).restrict_to_repacking();
    	   else
    	     mutate_residue->nonconst_residue_task( resj ).restrict_absent_canonical_aas( allowed_aas );
    	 }
			 pose = pose_orig;//unbind if necessary
    	 TR<<"Mutating residue "<<pose.residue( resi ).name3()<<resi<<" to ";
    	 core::pack::pack_rotamers( pose, *scorefxn(), mutate_residue );
    	 TR<<pose.residue( resi ).name3()<<". Now relaxing..."<<std::endl;
    	 relax_mover()->apply( pose );
			 bool const triage_filter_pass( triage_filter()->apply( pose ) );
			 if( !triage_filter_pass ){
				 TR<<"Triage filter fails"<<std::endl;
				 if(report_all_) {
					 residue_id_val_map[ std::pair< core::Size, AA >( resi, target_aa ) ] = std::pair< core::Real, bool >(filter()->report_sm( pose ), false);
				 }
				 continue;
			 }
			 TR<<"Triage filter succeeds"<<std::endl;
			 residue_id_val_map[ std::pair< core::Size, AA >( resi, target_aa ) ] = std::pair< core::Real, bool >(filter()->report_sm( pose ), true);
    	 residue_id_map[ resi ].push_back( target_aa );
			 if( dump_pdb() ){
			 	using namespace protocols::jd2;
			 	JobOP job( JobDistributor::get_instance()->current_job() );
			 	std::stringstream fname;
				fname << job->input_tag() << pose_orig.residue( resi ).name3() << resi << pose.residue( resi ).name3();
			 	TR<<"Saving pose "<<fname.str();
				pose.dump_scored_pdb( fname.str(), *scorefxn() );
			 }
			 TR.flush();
		}//foreach target_aa
	}//foreach resi
	if( resfile_name() != "" ){
		std::ofstream resfile;
		resfile.open( resfile_name().c_str(), std::ios::out );
		resfile << resfile_general_property()<<"\nstart\n";
		for( std::map< core::Size, utility::vector1< AA > >::const_iterator pair = residue_id_map.begin(); pair != residue_id_map.end(); ++pair ){
			resfile << pose.pdb_info()->number( pair->first )<<'\t'<<pose.pdb_info()->chain( pair->first )<<"\tPIKAA\t";
			foreach( AA const aa, pair->second )
				resfile<<oneletter_code_from_aa( aa );
			resfile<<'\n';
		}
		resfile.close();
	}//fi resfile_name()
	for( std::map< std::pair< core::Size, AA >, std::pair< core::Real, bool > >::const_iterator pair = residue_id_val_map.begin(); pair != residue_id_val_map.end(); ++pair ){
		core::conformation::Residue const native_res( pose.conformation().residue( pair->first.first ) );
		TR_residue_scan<<resfile_name()<<'\t'
									 <<pair->first.first<<'\t'
									 <<oneletter_code_from_aa( pair->first.second )<<'\t'
									 <<( delta() ? pair->second.first - baseline : pair->second.first )
									 <<(pair->second.second?"":"\tTRIAGED")<<std::endl;
	}
	TR.flush();
	return true;
}

core::Real
FilterScanFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( 1 );
}

void
FilterScanFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
}

void
FilterScanFilter::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose )
{
	TR << "FilterScanFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	std::string const triage_filter_name( tag->getOption< std::string >( "triage_filter", "true_filter" ) );
	protocols::filters::Filters_map::const_iterator triage_filter_it( filters.find( triage_filter_name ) );
	if( triage_filter_it == filters.end() )
		utility_exit_with_message( "Triage filter "+triage_filter_name+" not found" );

	triage_filter( triage_filter_it->second );

	std::string const filter_name( tag->getOption< std::string >( "filter" ) );
	protocols::filters::Filters_map::const_iterator filter_it( filters.find( filter_name ) );
	if( filter_it == filters.end() )
		utility_exit_with_message( "Filter "+filter_name+" not found" );
	filter( filter_it->second );
	std::string const relax_mover_name( tag->getOption< std::string >( "relax_mover", "null" ) );
	protocols::moves::Movers_map::const_iterator mover_it( movers.find( relax_mover_name ) );
	if( mover_it == movers.end() )
		utility_exit_with_message( "Relax mover "+relax_mover_name+" not found" );
	relax_mover( mover_it->second );

	delta( tag->getOption< bool >( "delta", false ) );
	report_all( tag->getOption< bool >( "report_all", false ) );
//	if( delta() ){
//		unbound( tag->getOption< bool >( "unbound", false ) );
//		if( unbound() )
//			jump( tag->getOption< core::Size >( "jump", 1 ) );

//		core::pose::Pose p( pose );
//		unbind( p );
//		core::scoring::ScoreFunctionOP score12 = data.get< core::scoring::ScoreFunction *>( "scorefxns", "score12" );
//		(*score12)(p);

//		relax_mover()->apply( p );
//		baseline( filter()->report_sm( p ) );
//	}

	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	resfile_name( tag->getOption< std::string >( "resfile_name",resfile_name() ) );
	resfile_general_property( tag->getOption< std::string >( "resfile_general_property", "nataa" ) );
	dump_pdb( tag->getOption< bool >( "dump_pdb", false ) );
	TR<<"with options resfile_name: "<<resfile_name()<<" resfile_general_property "<<resfile_general_property()<<" unbound "<<unbound()<<" jump "<<jump()<<" delta "<<delta()<<" filter "<<filter_name<<" dump_pdb "<<dump_pdb()<<std::endl;
}

core::scoring::ScoreFunctionOP
FilterScanFilter::scorefxn() const{
	return scorefxn_;
}

void
FilterScanFilter::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

protocols::filters::FilterOP
FilterScanFilter::fresh_instance() const{
	return new FilterScanFilter();
}

FilterScanFilter::~FilterScanFilter(){}

protocols::filters::FilterOP
FilterScanFilter::clone() const{
	return new FilterScanFilter( *this );
}

protocols::filters::FilterOP
FilterScanFilterCreator::create_filter() const { return new FilterScanFilter; }

std::string
FilterScanFilterCreator::keyname() const { return "FilterScan"; }

void
FilterScanFilter::dump_pdb( bool const d ){
	dump_pdb_ = d;
}

bool
FilterScanFilter::dump_pdb() const{
	return dump_pdb_;
}

} // filters
} // protein_interface_design
} // protocols
