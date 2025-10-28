// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_ddg/AlaScan.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
// Project Headers

#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_ddg/DdgFilter.hh>
#include <protocols/simple_ddg/AlaScan.hh>
#include <protocols/simple_ddg/AlaScanCreator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <core/pack/task/ResidueLevelTask.hh> // AUTO IWYU For ResidueLevelTask


namespace protocols {
namespace simple_ddg {

static basic::Tracer TR( "protocols.simple_filters.AlaScan" );




void
AlaScan::scorefxn( core::scoring::ScoreFunctionOP scorefxn )
{
	scorefxn_ = scorefxn;
}

AlaScan::AlaScan():
	Filter( "AlaScan" ),
	chain1_( false ),
	chain2_( true ),
	repeats_( 1 ),
	distance_threshold_( 8.0 ),
	jump_( 1 ),
	symmetry_( false ),
	repack_( true )
{}

AlaScan::AlaScan( bool const chain1, bool const chain2, core::Size const repeats, core::Real const dist, core::scoring::ScoreFunctionCOP scorefxn, core::Size const jump=1, bool const symmetry=false ) : Filter( "AlaScan" ),
	chain1_( chain1 ),
	chain2_( chain2 ),
	repeats_( repeats ),
	distance_threshold_( dist ),
	jump_( jump ),
	symmetry_( symmetry ),
	repack_( true )
{
	scorefxn_ = scorefxn->clone();
}

AlaScan::~AlaScan() = default;

bool
AlaScan::repack() const
{
	return repack_;
}

void
AlaScan::repack( bool const repack )
{
	repack_ = repack;
}

void
AlaScan::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	distance_threshold_ = tag->getOption<core::Real>( "interface_distance_cutoff", distance_threshold_ );
	chain1_ = tag->getOption< bool >( "partner1", chain1_ );
	chain2_ = tag->getOption< bool >( "partner2", chain2_ );
	jump_ = tag->getOption< core::Size >( "jump", jump_ );
	runtime_assert( chain1_ || chain2_ );
	repeats_ = tag->getOption< core::Size >( "repeats", repeats_ );
	symmetry_ = tag->getOption< bool >( "symmetry", symmetry_ );
	repack( tag->getOption< bool >( "repack", repack_ ) );

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	if ( symmetry_ ) {
		TR<<"Symmetric AlaScan with distance threshold of "<<distance_threshold_<<" Ang "<<". jump="<<jump_<<" partner1="<<chain1_<<", partner2="<<chain2_<<" using "<<repeats_<<" repeats."<<std::endl;
	} else {
		TR<<"AlaScan with distance threshold of "<<distance_threshold_<<" Ang "<<". jump="<<jump_<<" partner1="<<chain1_<<", partner2="<<chain2_<<" using "<<repeats_<<" repeats repack "<<repack()<<std::endl;
	}
}


core::Real
AlaScan::ddG_for_single_residue( core::pose::Pose const & const_pose, core::Size const resi ) const
{
	if ( !const_pose.residue( resi ).is_protein() ) {
		TR.Warning << "Non-protein residue "<< resi<<" was requested for ala-scan. Returning 0"<<std::endl;
		return 0.0;
	}
	core::Size const rb_jump( jump_ );
	core::pose::Pose pose( const_pose );

	simple_ddg::DdgFilter ddg_filter( 100/*ddg_threshold*/, scorefxn_, rb_jump, 1 /*repeats*/ );
	if ( repack() ) {
		TR<<"Energy calculations are carried out with repacking in the bound and unbound states (ddG)\n";
	} else {
		TR<<"Energy calculations are carried out without repackign in the bound and unbound states (dG)\n";
	}
	ddg_filter.repack( repack() );

	utility::vector1< bool > allowed_aas;
	allowed_aas.assign( core::chemical::num_canonical_aas, false );
	allowed_aas[ core::chemical::aa_ala ] = true;
	using namespace core::pack::task;

	PackerTaskOP task = TaskFactory::create_packer_task( pose );
	task->initialize_from_command_line().or_include_current( true );
	for ( core::Size resj=1; resj<=pose.size(); ++resj ) {
		if ( resi == resj ) {
			task->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
		} else {
			task->nonconst_residue_task( resj ).prevent_repacking();
		}
	}
	core::pack::pack_rotamers( pose, *scorefxn_, task );
	core::Real accumulate_ddg = 0;

	for ( core::Size r=1; r<=repeats_; ++r ) {
		accumulate_ddg += (rb_jump==0 ? ( *scorefxn_ )( pose ) : ddg_filter.compute( pose ) );
	}
	core::Real const mut_ddg( accumulate_ddg / repeats_ );

	TR.flush();
	return( mut_ddg );
}

void
AlaScan::report( std::ostream & out, core::pose::Pose const & const_pose ) const
{
	if ( symmetry_ ) {
		report_symmetry( out, const_pose );
		return;
	}

	core::Size const rb_jump( jump_ );
	core::pose::Pose pose( const_pose );

	core::kinematics::FoldTree const fold_tree = pose.conformation().fold_tree();

	core::Size upstream_jump_res, downstream_jump_res;
	upstream_jump_res = fold_tree.upstream_jump_residue( jump_ );
	downstream_jump_res = fold_tree.downstream_jump_residue( jump_ );

	core::Size const chain_begin( chain1_ ? 1 : downstream_jump_res );
	core::Size const chain_end  ( chain2_ ? pose.size() : upstream_jump_res );

	protocols::scoring::Interface interface_obj;
	interface_obj.jump( rb_jump == 0 ? 1 : rb_jump ); // 0 plays badly with interface obj.
	pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
	interface_obj.distance( distance_threshold_ );
	interface_obj.calculate( pose );

	simple_ddg::DdgFilter const ddg_filter( 100/*ddg_threshold*/, scorefxn_, rb_jump, 1 /*repeats*/ );

	core::Real accumulate_ddg( 0 );
	for ( core::Size r=1; r<=repeats_; ++r ) {
		accumulate_ddg += (rb_jump==0 ? ( *scorefxn_ )( pose ) : ddg_filter.compute( const_pose ) );
	}

	core::Real const wt_ddg( accumulate_ddg / repeats_ );
	for ( core::Size resi = chain_begin; resi <= chain_end; ++resi ) {
		if ( !pose.residue( resi ).is_protein() ) continue;
		if ( interface_obj.is_interface( resi ) ) {
			core::Real const mut_ddg( ddG_for_single_residue( const_pose, resi ) );
			core::Real const diff_ddg( mut_ddg - wt_ddg );

			core::pose::PDBInfoCOP pose_info( const_pose.pdb_info() );
			std::string chain( pose_info->chain( resi ) );
			int const number( pose_info->number( resi ) );
			std::string const res_type( const_pose.residue( resi ).name3() );
			out<<" "<<res_type<<" "<<number<<" "<<chain<<" : "<< ObjexxFCL::format::F (9,4,diff_ddg)<<'\n';
		}
	}
	out<<std::endl;
}

void
AlaScan::report_symmetry( std::ostream & out, core::pose::Pose const & const_pose ) const
{
	core::pose::Pose pose( const_pose );

	debug_assert( core::pose::symmetry::is_symmetric( pose ));
	auto & symm_conf (
		dynamic_cast<core::conformation::symmetry::SymmetricConformation & > ( pose.conformation()) );

	protocols::scoring::Interface interface_obj(1);
	pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
	interface_obj.distance( distance_threshold_ );
	interface_obj.calculate( pose );

	simple_ddg::DdgFilter const ddg( 100/*ddg_threshold*/, scorefxn_, 1, 1 /*repeats*/ /*, true */ ); //DdfFilter autodetects symmetry from input now
	core::Real accumulate_ddg( 0 );
	for ( core::Size r=1; r<=repeats_; ++r ) {
		accumulate_ddg += ddg.compute( const_pose );
	}
	core::Real const wt_ddg( accumulate_ddg / repeats_ );

	//core::Real const wt_ddg( ddg.compute( const_pose ) );
	utility::vector1< bool > allowed_aas;
	allowed_aas.assign( core::chemical::num_canonical_aas, false );
	allowed_aas[ core::chemical::aa_ala ] = true;
	for ( core::Size resi = 1; resi <= pose.size(); ++resi ) {
		if ( !symm_conf.Symmetry_Info()->bb_is_independent(resi) ) continue;
		if ( !pose.residue( resi ).is_protein() ) continue;
		if ( interface_obj.is_interface( resi ) ) {
			using namespace core::pack::task;

			PackerTaskOP task = TaskFactory::create_packer_task( pose );
			task->initialize_from_command_line().or_include_current( true );
			for ( core::Size resj=1; resj<=pose.size(); ++resj ) {
				if ( !pose.residue( resi ).is_protein() ) continue;
				if ( resi == resj ) {
					task->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
				} else {
					task->nonconst_residue_task( resj ).prevent_repacking();
				}
			}
			core::pack::pack_rotamers( pose, *scorefxn_, task );
			accumulate_ddg = 0;
			for ( core::Size r=1; r<=repeats_; ++r ) accumulate_ddg += ddg.compute( pose );
			core::Real const mut_ddg( accumulate_ddg / repeats_ );
			//core::Real const mut_ddg( ddg.compute( const_pose ) );
			core::Real const diff_ddg( mut_ddg - wt_ddg );

			core::pose::PDBInfoCOP pose_info( const_pose.pdb_info() );
			//char const chain( pose_info->chain( resi ) );
			//core::Size const number( pose_info->number( resi ) );
			std::string const res_type( const_pose.residue( resi ).name3() );
			//out<<" "<<res_type<<" "<<number<<" "<<chain<<" : "<< F (9,4,diff_ddg)<<'\n';
			out<<" "<<res_type<<" "<< resi <<" : "<< ObjexxFCL::format::F (9,4,diff_ddg)<<'\n';
			pose=const_pose;
		}
	}
}

std::string AlaScan::name() const {
	return class_name();
}

std::string AlaScan::class_name() {
	return "AlaScan";
}

void AlaScan::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute( "interface_distance_cutoff", xsct_real, "Distance from chain 1 to define chain 2's interface, and vice versa" )
		+ XMLSchemaAttribute( "partner1", xsct_rosetta_bool, "Perform alanine mutations on chain 1" )
		+ XMLSchemaAttribute( "partner2", xsct_rosetta_bool, "Perform alanine mutations on chain 2" )
		+ XMLSchemaAttribute( "jump", xsct_non_negative_integer, "An inter-chain jump over which the calculation is performed, numbered sequentially from 1" )
		+ XMLSchemaAttribute( "repeats", xsct_non_negative_integer, "Number of repeats over which to average" )
		+ XMLSchemaAttribute( "symmetry", xsct_rosetta_bool, "Assume a symmetric pose" )
		+ XMLSchemaAttribute( "repack", xsct_rosetta_bool, "Repack in a radius around the mutated residue" );

	rosetta_scripts::attributes_for_parse_score_function( attlist );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Perform alanine scanning over a particular jump-defined interface, with a certain number of repeats with which to average results. Repacking is suggested.", attlist );
}

std::string AlaScanFilterCreator::keyname() const {
	return AlaScan::class_name();
}

protocols::filters::FilterOP
AlaScanFilterCreator::create_filter() const {
	return utility::pointer::make_shared< AlaScan >();
}

void AlaScanFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AlaScan::provide_xml_schema( xsd );
}



}
}
