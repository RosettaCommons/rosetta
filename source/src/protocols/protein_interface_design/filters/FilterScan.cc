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
// AUTO-REMOVED #include <utility/file/FileName.hh>
#include <iostream>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/util.hh>
#include <boost/foreach.hpp>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.hh>
#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_filters/DeltaFilter.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <basic/options/keys/OptionKeys.hh>

namespace protocols {
namespace protein_interface_design{
namespace filters {

static thread_local basic::Tracer TR( "protocols.protein_interface_design.filters.FilterScanFilter" );
static thread_local basic::Tracer TR_residue_scan( "ResidueScan" );

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
	dump_pdb_( false ),
	rtmin_( false )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string temp_resfile_name( protocols::jd2::current_output_name() );
	temp_resfile_name = temp_resfile_name + ".resfile";
	resfile_name( temp_resfile_name );
	delta_filters_.clear();
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

std::string
FilterScanFilter::score_log_file() const{
	return score_log_file_;
}

void
FilterScanFilter::resfile_name( std::string const resfile_name ){
	resfile_name_ = resfile_name;
}

void
FilterScanFilter::score_log_file( std::string const score_log_file ){
	score_log_file_ = score_log_file;
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
	protocols::rigid::RigidBodyTransMover rbtm( pose, jump() );
	rbtm.step_size( 10000.0 );
	rbtm.apply( pose );
}

///@brief introduces a single-point subsitution and then performs the repack, rtmin, and relax moves that are requested.
void
FilterScanFilter::single_substitution( core::pose::Pose & pose, core::Size const resi, core::chemical::AA const target_aa ) const{
	using namespace core::chemical;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
  using namespace protocols::toolbox::task_operations;

  utility::vector1< bool > allowed_aas;
  allowed_aas.clear();
  allowed_aas.assign( num_canonical_aas, false );
  allowed_aas[ target_aa ] = true;
  TaskFactoryOP mut_res = new TaskFactory( *task_factory() );
  DesignAroundOperationOP dao = new DesignAroundOperation;///restrict repacking to 8.0A around target res to save time
  dao->design_shell( 0.1 );
  dao->include_residue( resi );
  mut_res->push_back( dao );
  PackerTaskOP mutate_residue = mut_res->create_task_and_apply_taskoperations( pose );
  mutate_residue->initialize_from_command_line().or_include_current( true );
  for( core::Size resj = 1; resj <= pose.total_residue(); ++resj ){
    if( resi != resj )
      mutate_residue->nonconst_residue_task( resj ).restrict_to_repacking();
    else
      mutate_residue->nonconst_residue_task( resj ).restrict_absent_canonical_aas( allowed_aas );
  }
  TR<<"Mutating residue "<<pose.residue( resi ).name3()<<resi<<" to ";
  protocols::simple_moves::PackRotamersMoverOP pack;
  if( core::pose::symmetry::is_symmetric( pose ) ) {
		pack =  new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn(), mutate_residue );
	} else{
		pack = new protocols::simple_moves::PackRotamersMover( scorefxn(), mutate_residue );
  }
  pack->apply( pose );
  if( rtmin() ) {
		// definition/allocation of RTmin mover must flag dependant, as some scoreterms are incompatable with RTmin initilization
		if( core::pose::symmetry::is_symmetric( pose ) ) {
			protocols::simple_moves::symmetry::SymRotamerTrialsMover rt( scorefxn(), *mutate_residue );
			rt.apply( pose );
		} else{
			protocols::simple_moves::RotamerTrialsMinMoverOP rtmin;
			rtmin = new protocols::simple_moves::RotamerTrialsMinMover( scorefxn(), *mutate_residue );
			rtmin->apply(pose);
		}
	}
  TR<<pose.residue( resi ).name3()<<". Now relaxing..."<<std::endl;
	if ( relax_mover() ) {
		relax_mover()->apply( pose );
	}
}

utility::vector1< protocols::simple_filters::DeltaFilterOP > FilterScanFilter::delta_filters() const { return delta_filters_; }
void FilterScanFilter::delta_filters( utility::vector1< protocols::simple_filters::DeltaFilterOP > const d ){ delta_filters_ = d; }

bool
FilterScanFilter::apply(core::pose::Pose const & p ) const
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	core::pose::Pose pose( p );

	PackerTaskCOP task = task_factory()->create_task_and_apply_taskoperations( pose );
	utility::vector1< core::Size > being_designed;
	being_designed.clear();

	for( core::Size resi = 1; resi <= pose.total_residue(); ++resi ){
		if( task->residue_task( resi ).being_designed() && pose.residue(resi).is_protein() )
			being_designed.push_back( resi );
	}
	if( being_designed.empty() ) {
		TR.Warning << "WARNING: No residues are listed as designable." << std::endl;
		return true;
	}
	std::map< core::Size, utility::vector1< AA > > residue_id_map;
	std::map< std::pair< core::Size, AA >, std::pair< core::Real, bool > > residue_id_val_map; // position, aa identity : report filter value, triage filter accepted?
	residue_id_map.clear(); residue_id_val_map.clear();
	unbind( pose );
	core::pose::Pose const pose_orig( pose );// const to ensure that nothing silly happens along the way...
	BOOST_FOREACH( core::Size const resi, being_designed ){
		pose = pose_orig;
		///compute baseline
// SJF 29Aug13 in the following we try to extract a stable baseline value for the substitution. The repacking steps create large levels of noise where mutations to self appear show high ddG differences. To prevent that, we cheat Rosetta by asking it to replace the residue to self three times in a row and take the baseline after the third computation. If you use ddG, it is recommended to use many repeats (>=5) to lower noise levels further. I think that if mutations to self show noise levels <=0.5R.e.u. the protocol is acceptable.
		single_substitution( pose, resi, pose.residue( resi ).aa() ); /// mutates to self. This simply activates packing/rtmin/relax at the site. By doing this on a per residue basis we ensure that the baseline is computed in exactly the same way as the mutations
		single_substitution( pose, resi, pose.residue( resi ).aa() );
		core::pose::Pose const pose_ref( pose ); // 10Jul14 Adi debugging: save pose after two substitutions to self but compute baseline relative to 3 substitutions
		single_substitution( pose, resi, pose.residue( resi ).aa() );
//		pose.dump_scored_pdb( "at_baseline.pdb", *scorefxn() );
		BOOST_FOREACH( protocols::simple_filters::DeltaFilterOP const delta_filter, delta_filters_ ){
			std::string const fname( delta_filter->get_user_defined_name() );
			core::Real const fbaseline( delta_filter->filter()->report_sm( pose ) );
			delta_filter->baseline( fbaseline );
			TR<<"Computed baseline at position "<<resi<<" with filter "<<fname<<" is "<<fbaseline<<std::endl;
		}
    typedef std::list< ResidueTypeCOP > ResidueTypeCOPList;
    ResidueTypeCOPList const & allowed( task->residue_task( resi ).allowed_residue_types() );
    utility::vector1< AA > allow_temp;
    allow_temp.clear();
    BOOST_FOREACH( ResidueTypeCOP const t, allowed ){
    	allow_temp.push_back( t->aa() );
    }
//		core::pose::Pose const pose_ref( pose ); //10Jul14 Adi debugging: original line for pose_ref
		BOOST_FOREACH( AA const target_aa, allow_temp ){
			 pose = pose_ref; //29Aug13 previously was pose_orig; here
			 single_substitution( pose, resi, target_aa );
//				pose.dump_scored_pdb( "after_mut.pdb", *scorefxn() );
			 bool triage_filter_pass( false );
			 if( delta_filters_.size() > 0 ){
				 triage_filter_pass=true;
				 BOOST_FOREACH( protocols::simple_filters::DeltaFilterCOP const delta_filter, delta_filters_ ){
					 triage_filter_pass = delta_filter->apply( pose );
					 if( delta_filters_.size() == 1 )
							residue_id_val_map[ std::pair< core::Size, AA >( resi, target_aa ) ] = std::pair< core::Real, bool >(delta_filter()->report_sm( pose ), triage_filter_pass);

					 if( !triage_filter_pass )
						 break;
				 }
				 if( triage_filter_pass )
    	 		 residue_id_map[ resi ].push_back( target_aa );
			 }
			 else
			 	triage_filter_pass = triage_filter()->apply( pose );
			 if( !triage_filter_pass ){
				 TR<<"Triage filter fails"<<std::endl;
				 if(report_all_ && !delta()) {
					 residue_id_val_map[ std::pair< core::Size, AA >( resi, target_aa ) ] = std::pair< core::Real, bool >(filter()->report_sm( pose ), false);
				 }
				 continue;
			 }
			 TR<<"Triage filter succeeds"<<std::endl;
			 if( !delta() ){
			 		residue_id_val_map[ std::pair< core::Size, AA >( resi, target_aa ) ] = std::pair< core::Real, bool >(filter()->report_sm( pose ), true);
    	 		residue_id_map[ resi ].push_back( target_aa );
			 }
			 if( dump_pdb() ){
			 	using namespace protocols::jd2;
			 	JobOP job( JobDistributor::get_instance()->current_job() );
			 	std::stringstream fname;
				fname << job->input_tag() << pose_orig.residue( resi ).name3() << resi << pose.residue( resi ).name3()<<".pdb";
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
			BOOST_FOREACH( AA const aa, pair->second )
				resfile<<oneletter_code_from_aa( aa );
			resfile<<'\n';
		}
		resfile.close();
	} //fi resfile_name()

	if ( score_log_file() != "" ) {
		std::ofstream scorefile;
		scorefile.open( score_log_file().c_str(), std::ios::out );

		using namespace ObjexxFCL::format;
		for( std::map< std::pair< core::Size, AA >, std::pair< core::Real, bool > >::const_iterator pair = residue_id_val_map.begin(); pair != residue_id_val_map.end(); ++pair ){
			core::conformation::Residue const native_res( pose.conformation().residue( pair->first.first ) );
			scorefile
			<< pair->first.first << '\t'
			<< p.residue( pair->first.first ).name1() <<'\t'
			<< oneletter_code_from_aa( pair->first.second )<<'\t'
			<< F(9,6, pair->second.first) <<std::endl;
		}
		scorefile.close();
	} // fi score_log_file()

	for( std::map< std::pair< core::Size, AA >, std::pair< core::Real, bool > >::const_iterator pair = residue_id_val_map.begin(); pair != residue_id_val_map.end(); ++pair ){
		core::conformation::Residue const native_res( pose.conformation().residue( pair->first.first ) );
		TR_residue_scan<<resfile_name()<<'\t'
									 << pair->first.first<<'\t'
									 << oneletter_code_from_aa( pair->first.second )<<'\t'
		               <<pair->second.first
                   <<(pair->second.second?"":"\tTRIAGED")<<std::endl;
	}
	TR.flush();
	return true;
}

core::Real
FilterScanFilter::report_sm( core::pose::Pose const & ) const
{
	return( 1 );
}

void
FilterScanFilter::report( std::ostream &, core::pose::Pose const & ) const
{
}

void
FilterScanFilter::parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & )
{
	TR << "FilterScanFilter"<<std::endl;
	runtime_assert( tag->hasOption( "filter" ) || tag->hasOption( "delta_filters" ));
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	std::string const triage_filter_name( tag->getOption< std::string >( "triage_filter", "true_filter" ) );
	protocols::filters::Filters_map::const_iterator triage_filter_it( filters.find( triage_filter_name ) );

	//These #ifdefs are a terrible hack to work around a compiler bug in mpicxx. sorry
#ifdef USEMPI
	if( triage_filter_it == filters.end() )
		utility_exit_with_message( "Triage filter "+triage_filter_name+" not found" );
#endif
#ifndef USEMPI
	if( triage_filter_it == filters.end() )
		throw utility::excn::EXCN_RosettaScriptsOption( "Triage filter "+triage_filter_name+" not found" );
#endif

	triage_filter( triage_filter_it->second );

	std::string const filter_name( tag->getOption< std::string >( "filter", "true_filter" ) );
	protocols::filters::Filters_map::const_iterator filter_it( filters.find( filter_name ) );

#ifdef USEMPI
	if( filter_it == filters.end() )
		utility_exit_with_message( "Filter "+filter_name+" not found" );
#endif
#ifndef USEMPI
	if( filter_it == filters.end() )
		throw utility::excn::EXCN_RosettaScriptsOption( "Filter "+filter_name+" not found" );
#endif

	filter( filter_it->second );
	std::string const relax_mover_name( tag->getOption< std::string >( "relax_mover", "null" ) );
	protocols::moves::Movers_map::const_iterator mover_it( movers.find( relax_mover_name ) );

#ifdef USEMPI
	if( mover_it == movers.end() )
		utility_exit_with_message( "Relax mover "+relax_mover_name+" not found" );
#endif
#ifndef USEMPI
	if( mover_it == movers.end() )
		throw utility::excn::EXCN_RosettaScriptsOption( "Relax mover "+relax_mover_name+" not found" );
#endif

	relax_mover( mover_it->second );

	delta( tag->getOption< bool >( "delta", false ) );
	report_all( tag->getOption< bool >( "report_all", false ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	resfile_name( tag->getOption< std::string >( "resfile_name",resfile_name() ) );
	resfile_general_property( tag->getOption< std::string >( "resfile_general_property", "nataa" ) );
	rtmin( tag->getOption< bool >( "rtmin", false ) );
	score_log_file( tag->getOption< std::string >( "score_log_file",score_log_file() ) );
	dump_pdb( tag->getOption< bool >( "dump_pdb", false ) );


	utility::vector1< std::string > delta_filter_names;
	delta_filter_names.clear();
	if( tag->hasOption( "delta_filters" ) ){
		delta_filter_names = utility::string_split( tag->getOption< std::string >( "delta_filters" ), ',' );
		TR<<"Using delta filters: ";
		BOOST_FOREACH( std::string const fname, delta_filter_names ){
			delta_filters_.push_back( dynamic_cast< protocols::simple_filters::DeltaFilter * >( protocols::rosetta_scripts::parse_filter( fname, filters )() ) );
			TR<<fname<<",";
		}
		TR<<std::endl;
	}
	TR<<"with options resfile_name: "<<resfile_name()<<" resfile_general_property "<<resfile_general_property()<<" unbound "<<unbound()<<" jump "<<jump()<<" delta "<<delta()<<" filter "<<filter_name<<" dump_pdb "<<dump_pdb()<<" rtmin "<<rtmin()<<std::endl;
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
