// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/simple_filters/RelativePoseFilter.hh>
#include <protocols/simple_filters/RelativePoseFilterCreator.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <boost/foreach.hpp>
#include <utility/string_util.hh>
#include <protocols/simple_moves/DumpPdb.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmDataFactory.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/pose/util.hh>

namespace protocols {
namespace simple_filters {

static thread_local basic::Tracer TR( "protocols.simple_filters.RelativePoseFilter" );

///@brief default ctor
RelativePoseFilter::RelativePoseFilter() :
	parent( "RelativePose" ),
	filter_( /* NULL */ ),
	relax_mover_( /* NULL */ ),
	dump_pose_fname_( "" ),
	pose_( /* NULL */ ),
	scorefxn_( /* NULL */ ),
	packing_shell_( 8.0 ),
	thread_( true ),
	baseline_( true ),
	baseline_val_( -9999 ),
	unbound_( false ),
	copy_stretch_( false ),
	symmetry_definition_(""),
	rtmin_( false )
{
	alignment_.clear();
	copy_comments_.clear();
}

void
RelativePoseFilter::thread( bool const t ){
	thread_ = t;
}

bool
RelativePoseFilter::thread() const{
	return thread_;
}

void
RelativePoseFilter::packing_shell( core::Real const s ){
	packing_shell_ = s;
}

core::Real
RelativePoseFilter::packing_shell() const{
	return packing_shell_;
}

std::string
RelativePoseFilter::dump_pose_fname() const{
	return dump_pose_fname_;
}

void
RelativePoseFilter::dump_pose_fname( std::string const s ){
	dump_pose_fname_ = s;
}

core::pose::PoseOP
RelativePoseFilter::pose() const{
	return pose_;
}

void
RelativePoseFilter::pose( core::pose::PoseOP pose ){
	pose_= pose;
}

void
RelativePoseFilter::filter( protocols::filters::FilterOP filter ){
	filter_ = filter;
}

protocols::filters::FilterOP
RelativePoseFilter::filter() const{
	return filter_;
}

bool
RelativePoseFilter::apply(core::pose::Pose const & p ) const
{
	core::pose::PoseCOP threaded_pose( thread_seq( p ) );
	return( filter()->apply( *threaded_pose ) );
}

core::pose::PoseOP
RelativePoseFilter::thread_seq( core::pose::Pose const & p ) const{
	using namespace core::chemical;
	using namespace protocols::toolbox::task_operations;

	core::pose::PoseOP copy_pose( new core::pose::Pose( *pose() ) ); // don't let the pose drift
	// if(symmetry_definition()!=""){
	// 	if(core::pose::symmetry::is_symmetric(*copy_pose)) {
	// 		core::pose::symmetry::extract_asymmetric_unit(*copy_pose,*copy_pose);
	// 	}
	// 	core::pose::symmetry::make_symmetric_pose(*pose_,*symmdata_);
	// }
	if( unbound() ){
	  	protocols::rigid::RigidBodyTransMover rbtm( *copy_pose, 1 );
	  	rbtm.step_size( 10000.0 );
  		rbtm.apply( *copy_pose );
	}
	if( copy_comments().size() ){//copy relevant comments to copy_pose
		BOOST_FOREACH( std::string const key, copy_comments() ){
			std::string val;
			bool const success = core::pose::get_comment( p, key, val );
			if( success ){
				core::pose::add_comment( *copy_pose, key, val );
				TR<<"Added comment key/val: "<<key<<'/'<<val<<" to the relative pose"<<std::endl;
			}
			else
				TR<<"Failed to find comment key/val: "<<key<<'/'<<val<<" in reference pose"<<std::endl;
		}
	}//fi copy_comments.size()

	if( copy_stretch() ){ // just copy the aligned stretch, and then go straight to relax. No repacking
		copy_pose->copy_segment( alignment_.size()/*how many residues*/, p/*src*/, alignment_.begin()->first/*start on target*/, alignment_.begin()->second/*start on src*/ );
		copy_pose->conformation().detect_disulfides();
	}
	else{ // no copy_stretch. Repack etc. carefully
		DesignAroundOperationOP dao( new DesignAroundOperation );
		dao->design_shell( packing_shell() );
		std::vector< core::Size > diffs;
		diffs.clear();
		TR<<"differences at positions: ";
		for( std::map< core::Size, core::Size >::const_iterator aln=alignment_.begin(); aln!=alignment_.end(); ++aln ){
			char const res1_name(pose()->conformation().residue( aln->first ).name1());
			char const res2_name(p.conformation().residue( aln->second ).name1());
			if( res1_name != res2_name ) {
				diffs.push_back( aln->first );
				TR<<res1_name<<aln->first<<res2_name<<", ";
			}
		}
		TR<<std::endl;
		if( baseline() )
			TR<<"baseline: "<<baseline_val()<<std::endl;
		BOOST_FOREACH( core::Size const d, diffs )
			dao->include_residue( d );
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		TaskFactoryOP tf( new TaskFactory );
		tf->push_back( dao );
		tf->push_back( TaskOperationCOP( new IncludeCurrent ) );
		tf->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
		core::pack::task::PackerTaskOP pack = tf->create_task_and_apply_taskoperations( *pose() );
		TR<<"prevent repacking (p); restrict to repacking (r), design (d): ";
		for( core::Size i = 1; i<=pose()->total_residue(); ++i ){
			if( !pack->nonconst_residue_task( i ).being_designed() ){ // prevent repacking on all non-designable residues
				pack->nonconst_residue_task( i ).prevent_repacking();
				TR<<i<<"(p), ";
			}
			else if( std::find( diffs.begin(), diffs.end(), i ) == diffs.end() ){
				pack->nonconst_residue_task( i ).restrict_to_repacking();
				TR<<i<<"(r), ";
			}
			else{//design!
				utility::vector1< bool > allowed_aas( num_canonical_aas, false );
				if( thread() ){
					//allow the aa at the corresonding position in the alignment with the disk pose p
					//we can't use the alignment_ map's [] operator because it can add new elements and this is a const function!
					allowed_aas[ p.residue( alignment_.find( i )->second ).aa() ] = true;
				} else
					//or allow the original aa in the current pose
					allowed_aas[ pose()->residue( i ).aa() ] = true;
				pack->nonconst_residue_task( i ).restrict_absent_canonical_aas( allowed_aas );
				TR<<i<<"(d), ";
			}
		}//for i=1->total_residue
		TR<<std::endl;
		using namespace protocols::simple_moves;
		PackRotamersMoverOP prm;
		RotamerTrialsMinMoverOP rtmin;
		if( core::pose::symmetry::is_symmetric( *copy_pose ) )
			prm = PackRotamersMoverOP( new symmetry::SymPackRotamersMover( scorefxn(), pack ) );
		else
			prm = PackRotamersMoverOP( new PackRotamersMover( scorefxn(), pack ) );
		prm->apply( *copy_pose );
		if( rtmin.get() ){
			rtmin = RotamerTrialsMinMoverOP( new RotamerTrialsMinMover( scorefxn(), *pack ) );
			rtmin->apply( *copy_pose );
			TR<<"finished rtmin"<<std::endl;
		}
	}/// end else no copy_stretch
	relax_mover()->apply( *copy_pose );
	return( copy_pose );
}


core::Real
RelativePoseFilter::compute( core::pose::Pose const & p ) const{
	core::pose::PoseOP threaded_pose( thread_seq( p ) );
	core::Real const filter_val( filter()->report_sm( *threaded_pose ) );
	if( dump_pose_fname() != "" ){
		protocols::simple_moves::DumpPdb dump( dump_pose_fname() );
		dump.set_scorefxn( scorefxn() );
		dump.apply( *threaded_pose );
	}
	TR<<"filter "<<filter_name()<<" reports value pos: "<<filter_val<<", "<<std::endl;
	if( baseline() ){
		TR<<"filter "<<filter_name()<<" reports value, baseline: "<<filter_val<<", "<<baseline_val()<<std::endl;
		return( filter_val - baseline_val() );
	}
	else{
		TR<<"filter "<<filter_name()<<" val: "<<filter_val<<std::endl;
		return( filter_val );
	}
}

core::Real
RelativePoseFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

void
RelativePoseFilter::report( std::ostream & out, core::pose::Pose const & ) const
{
	out<<"Dummy: All reports are dummy for this filter. Only report_sm is defined"<<std::endl;
}

void
RelativePoseFilter::pdb_name( std::string const pdb_name ){
	pose( core::import_pose::pose_from_pdb( pdb_name, false /*read foldtree*/ ) );
}


/// alignment is expecting X1:Y1,X2:Y2,X3:Y3... where X is the protein on disk (target) and Y is the active structure (starting structure). When no alignment is given it is implied that the poses are trivially aligned 1..nres
/// Feb2012 added option to align entire chains: A:B,D:C. Notice that no testing is made to ensure correct lengths etc., simply aligns from the start to end of the chains sequentially.
void
RelativePoseFilter::parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & p )
{
	using namespace protocols::rosetta_scripts;

	TR << "RelativePoseFilter"<<std::endl;
	std::string pose_fname ("");
	bool use_native ( tag->getOption< bool >( "use_native_pdb", false ));
	if (use_native)	pose( core::import_pose::pose_from_pdb(basic::options::option[ basic::options::OptionKeys::in::file::native ], false));

	if ( tag->hasOption( "pdb_name" ) ) {
			pose_fname = tag->getOption< std::string >( "pdb_name" );
			pdb_name( pose_fname );
	}
	if( tag->hasOption( "symmetry_definition" ) ){
		symmetry_definition_ = tag->getOption< std::string >( "symmetry_definition", "" );
		if(core::pose::symmetry::is_symmetric(*pose_)) {
			core::pose::symmetry::extract_asymmetric_unit(*pose_,*pose_);
		}
		symmdata_ = core::conformation::symmetry::SymmDataOP( new core::conformation::symmetry::SymmData( pose_->n_residue(), pose_->num_jump() ) );
		symmdata_->read_symmetry_data_from_file(symmetry_definition_);
		core::pose::symmetry::make_symmetric_pose(*pose_,*symmdata_);
	}
	relax_mover( parse_mover( tag->getOption< std::string >( "relax_mover", "null" ), movers ) );
	filter( parse_filter( tag->getOption< std::string >( "filter" ), filters ) );
	filter_name( tag->getOption< std::string> ( "name", "" ) );
	baseline( tag->getOption< bool >( "baseline", 1 ));
	if( baseline() ){
		relax_mover()->apply( *pose() );
		baseline_val( filter()->report_sm( *pose() ) );
		TR<<"The baseline value for the pose read from disk is: "<<baseline_val()<<std::endl;
	}
	else
		TR<<"Baseline turned off. Is that intended?"<<std::endl;
	dump_pose_fname( tag->getOption< std::string >( "dump_pose", "" ) );
	if( tag->hasOption( "alignment" ) ){
		utility::vector1< std::string > const residue_pairs( utility::string_split( tag->getOption< std::string >( "alignment", "" ), ',' ) );
		BOOST_FOREACH( std::string const residue_pair, residue_pairs ){
			utility::vector1< std::string > const residues( utility::string_split( residue_pair, ':' ) );
			runtime_assert( residues.size() == 2 );
			char const residues1_cstr( residues[ 1 ].c_str()[ 0 ] ), residues2_cstr( residues[ 2 ].c_str()[ 0 ] ); // these may hold the chain designators
			if( residues[ 1 ].length() == 1 &&
				( residues1_cstr <= 'Z' && residues2_cstr >= 'A' ) &&
				( residues[ 2 ].length() == 1 &&
				( residues2_cstr <= 'Z' && residues2_cstr >= 'A' ) ) ){ // are we aligning two chains to one another?
					core::pose::PDBInfoCOP pdbinfo1( pose()->pdb_info() ), pdbinfo2( p.pdb_info() );
					core::Size pose_res( 1 ), p_res( 1 );
					for(; pose_res <= pose()->total_residue(); ++pose_res )// find chain1 start
						if( pdbinfo1->chain( pose_res ) == residues1_cstr ) break;
					for(; p_res <= p.total_residue(); ++p_res ) // find chain2 start
						if( pdbinfo2->chain( p_res ) == residues2_cstr ) break;
					for( core::Size index = 0; ; ++index ){/// push aligned residues
						alignment_[ pose_res ] = p_res;
						if( pose_res == pose()->total_residue() || p_res == p.total_residue() )//end of chains -> stop aligning
							break;
						pose_res++; p_res++;
						if( pdbinfo1->chain( pose_res ) != residues1_cstr ||
							pdbinfo2->chain( p_res ) != residues2_cstr ) /// end of aligned chains
							break;
					}
				}//fi aligning two chains
			else /// aligning individual residues
				alignment_[ core::pose::parse_resnum( residues[ 1 ], *pose() ) ] = core::pose::parse_resnum( residues[ 2 ], p );
		}
	}
	else{
//		runtime_assert( pose()->total_residue() == p.total_residue() || core::pose::symmetry::is_symmetric( p ) );
		for( core::Size i=1; i<=p.total_residue(); ++i )
			alignment_[ i ] = i;
	}
	scorefxn( parse_score_function( tag, data ) );
	packing_shell( tag->getOption( "packing_shell", packing_shell() ));
	rtmin( tag->getOption< bool >("rtmin", 0 ) );
	runtime_assert( packing_shell() >= 0 );
	thread( tag->getOption< bool >( "thread", thread() ) );
	unbound( tag->getOption< bool >( "unbound", false ) );
	copy_stretch( tag->getOption< bool >( "copy_stretch", false ) );
	if( tag->hasOption( "copy_comments" ) )
		copy_comments( utility::string_split( tag->getOption< std::string >( "copy_comments","" ), ',' ) );
	TR<<"with pdb: "<<pose_fname<<" dumping fname "<<dump_pose_fname()<<" thread: "<<thread()<<" unbound "<<unbound()<<" copy_stretch: "<<copy_stretch()<<" rtmin "<<rtmin()<<" and packing_shell: "<<packing_shell();
	if(symmetry_definition_!="") TR<<" symmetry: "<<symmetry_definition_<<std::endl;
	TR<<std::endl;
}

protocols::filters::FilterOP
RelativePoseFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new RelativePoseFilter() );
}

RelativePoseFilter::~RelativePoseFilter(){}

protocols::filters::FilterOP
RelativePoseFilter::clone() const{
	return protocols::filters::FilterOP( new RelativePoseFilter( *this ) );
}

protocols::filters::FilterOP
RelativePoseFilterCreator::create_filter() const { return protocols::filters::FilterOP( new RelativePoseFilter ); }

std::string
RelativePoseFilterCreator::keyname() const { return "RelativePose"; }

protocols::moves::MoverOP
RelativePoseFilter::relax_mover() const{
	return relax_mover_;
}

void
RelativePoseFilter::relax_mover( protocols::moves::MoverOP const m ){
	relax_mover_ = m;
}

void
RelativePoseFilter::scorefxn( core::scoring::ScoreFunctionOP const scorefxn ){
	scorefxn_ = scorefxn;
}

core::scoring::ScoreFunctionOP
RelativePoseFilter::scorefxn() const{
	return( scorefxn_ );
}

bool
RelativePoseFilter::baseline() const{
	return baseline_;
}

void
RelativePoseFilter::baseline( bool const b ){
	baseline_ = b;
}

core::Real
RelativePoseFilter::baseline_val() const{
	runtime_assert( baseline() );
	return baseline_val_;
}

void
RelativePoseFilter::baseline_val( core::Real const b ){
	runtime_assert( baseline() );
	baseline_val_ = b;
}

void
RelativePoseFilter::symmetry_definition( std::string const s ){
	symmetry_definition_ = s;
}

std::string
RelativePoseFilter::symmetry_definition() const{
	return symmetry_definition_;
}

void
RelativePoseFilter::filter_name( std::string const s ){
	filter_name_ = s;
}

std::string
RelativePoseFilter::filter_name()const{ return filter_name_; }

bool
RelativePoseFilter::rtmin() const{ return rtmin_; }

void
RelativePoseFilter::rtmin( bool const b ){ rtmin_ = b; }

} // simple_filters
} // protocols
