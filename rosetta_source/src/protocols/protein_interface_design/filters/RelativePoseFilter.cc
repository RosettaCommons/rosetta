// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/RelativePoseFilter.hh>
#include <protocols/protein_interface_design/filters/RelativePoseFilterCreator.hh>
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
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <utility/string_util.hh>
#include <protocols/protein_interface_design/movers/DumpPdb.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design{
namespace filters {

static basic::Tracer TR( "protocols.protein_interface_design.filters.RelativePoseFilter" );

///@brief default ctor
RelativePoseFilter::RelativePoseFilter() :
	parent( "RelativePose" ),
	filter_( NULL ),
	relax_mover_( NULL ),
	dump_pose_fname_( "" ),
	pose_( NULL ),
	scorefxn_( NULL ),
	packing_shell_( 8.0 ),
	thread_( true ),
	baseline_( true ),
	baseline_val_( -9999 ),
	unbound_( false ),
	copy_stretch_( false )
{
	alignment_.clear();
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
	TR<<"filter's value: "<<filter()->report_sm( *threaded_pose )<<std::endl;
	return( filter()->apply( *threaded_pose ) );
}

core::pose::PoseOP
RelativePoseFilter::thread_seq( core::pose::Pose const & p ) const{
	using namespace core::chemical;
	using namespace protocols::toolbox::task_operations;

	core::pose::PoseOP copy_pose( new core::pose::Pose( *pose() ) ); // don't let the pose drift
	if( unbound() ){
  	protocols::rigid::RigidBodyTransMover rbtm( *copy_pose, 1 );
  	rbtm.step_size( 10000.0 );
  	rbtm.apply( *copy_pose );
	}
	if( copy_stretch() ) // just copy the aligned stretch, and then go straight to relax. No repacking
		copy_pose->copy_segment( alignment_.size()/*how many residues*/, p/*src*/, alignment_.begin()->first/*start on target*/, alignment_.begin()->second/*start on src*/ );
	else{ // no copy_stretch. Repack etc. carefully
		DesignAroundOperationOP dao = new DesignAroundOperation;
		dao->design_shell( packing_shell() );
		std::vector< core::Size > diffs;
		diffs.clear();
		for( std::map< core::Size, core::Size >::const_iterator aln=alignment_.begin(); aln!=alignment_.end(); ++aln )
			if( pose()->conformation().residue( aln->first ).aa() != p.conformation().residue( aln->second ).aa() ) diffs.push_back( aln->first );

		TR<<"baseline: "<<baseline_val()<<std::endl;
		TR<<"differences at positions: ";
		foreach( core::Size const d, diffs ){
			dao->include_residue( d );
			TR<<d<<", ";
		}
		TR<<std::endl;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		TaskFactoryOP tf = new TaskFactory;
		tf->push_back( dao );
		tf->push_back( new IncludeCurrent );
		tf->push_back( new InitializeFromCommandline );
		core::pack::task::PackerTaskOP pack = tf->create_task_and_apply_taskoperations( *pose() );
		for( core::Size i = 1; i<=pose()->total_residue(); ++i ){
			if( !pack->nonconst_residue_task( i ).being_designed() ) // prevent repacking on all non-designable residues
				pack->nonconst_residue_task( i ).prevent_repacking();
			else if( std::find( diffs.begin(), diffs.end(), i ) == diffs.end() )
				pack->nonconst_residue_task( i ).restrict_to_repacking();
			else{//design!
				utility::vector1< bool > allowed_aas( num_canonical_aas, false );
				if( thread() )
					allowed_aas[ p.residue( i ).aa() ] = true;
				else
					allowed_aas[ pose()->residue( i ).aa() ] = true;
				pack->nonconst_residue_task( i ).restrict_absent_canonical_aas( allowed_aas );
			}
		}//for i=1->total_residue
		using namespace protocols::simple_moves;
		PackRotamersMoverOP prm;
		if( core::pose::symmetry::is_symmetric( *copy_pose ) )
			prm = new symmetry::SymPackRotamersMover( scorefxn(), pack );
		else
			prm = new PackRotamersMover( scorefxn(), pack );
		prm->apply( *copy_pose );
	}/// end else no copy_stretch
	relax_mover()->apply( *copy_pose );
	return( copy_pose );
}


core::Real
RelativePoseFilter::compute( core::pose::Pose const & p ) const{
	core::pose::PoseOP threaded_pose( thread_seq( p ) );
	core::Real const filter_val( filter()->report_sm( *threaded_pose ) );
	if( dump_pose_fname() != "" ){
		protocols::protein_interface_design::movers::DumpPdb dump( dump_pose_fname() );
		dump.set_scorefxn( scorefxn() );
		dump.apply( *threaded_pose );
	}
	if( baseline() )
		return( filter_val - baseline_val() );
	else
		return( filter_val );
}

core::Real
RelativePoseFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

void
RelativePoseFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"Dummy: All reports are dummy for this filter. Only report_sm is defined"<<std::endl;
}

/// alignment is expecting X1:Y1,X2:Y2,X3:Y3... where X is the protein on disk (target) and Y is the active structure (starting structure). When no alignment is given it is implied that the poses are trivially aligned 1..nres
/// Feb2012 added option to align entire chains: A:B,D:C. Notice that no testing is made to ensure correct lengths etc., simply aligns from the start to end of the chains sequentially.
void
RelativePoseFilter::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & p )
{
	using namespace protocols::rosetta_scripts;

	TR << "RelativePoseFilter"<<std::endl;
	std::string const pose_fname( tag->getOption< std::string >( "pdb_name" ) );
	pose( core::import_pose::pose_from_pdb( pose_fname, false /*read foldtree*/ ) );
	relax_mover( parse_mover( tag->getOption< std::string >( "relax_mover", "null" ), movers ) );
	filter( parse_filter( tag->getOption< std::string >( "filter" ), filters ) );
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
		foreach( std::string const residue_pair, residue_pairs ){
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
						pose_res++; p_res++;
						if( pdbinfo1->chain( pose_res ) != residues1_cstr ||
							pdbinfo2->chain( p_res ) != residues2_cstr ) /// end of aligned chains
							break;
					}
				}//fi aligning two chains
			else /// aligning individual residues
				alignment_[ parse_resnum( residues[ 1 ], *pose() ) ] = parse_resnum( residues[ 2 ], p );
		}
	}
	else{
//		runtime_assert( pose()->total_residue() == p.total_residue() || core::pose::symmetry::is_symmetric( p ) );
		for( core::Size i=1; i<=p.total_residue(); ++i )
			alignment_[ i ] = i;
	}
	scorefxn( parse_score_function( tag, data ) );
	packing_shell( tag->getOption( "packing_shell", packing_shell() ));
	runtime_assert( packing_shell() >= 0 );
	thread( tag->getOption< bool >( "thread", thread() ) );
	unbound( tag->getOption< bool >( "unbound", false ) );
	copy_stretch( tag->getOption< bool >( "copy_stretch", false ) );
	TR<<"with pdb: "<<pose_fname<<" dumping fname "<<dump_pose_fname()<<" thread: "<<thread()<<" unbound "<<unbound()<<" copy_stretch: "<<copy_stretch()<<" and packing_shell: "<<packing_shell()<<std::endl;
}

protocols::filters::FilterOP
RelativePoseFilter::fresh_instance() const{
	return new RelativePoseFilter();
}

RelativePoseFilter::~RelativePoseFilter(){}

protocols::filters::FilterOP
RelativePoseFilter::clone() const{
	return new RelativePoseFilter( *this );
}

protocols::filters::FilterOP
RelativePoseFilterCreator::create_filter() const { return new RelativePoseFilter; }

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

} // filters
} // protein_interface_design
} // protocols
