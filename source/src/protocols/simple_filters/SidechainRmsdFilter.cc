// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SidechainRmsdFilter.cc
/// @brief A filter based on automorphic sidechain RMSD
/// @author Noah Ollikainen

#include <protocols/simple_filters/SidechainRmsdFilter.hh>
#include <protocols/simple_filters/SidechainRmsdFilterCreator.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/rms_util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer sidechain_rmsd_filter_tracer( "protocols.simple_filters.SidechainRmsdFilter" );

protocols::filters::FilterOP
SidechainRmsdFilterCreator::create_filter() const { return new SidechainRmsdFilter; }

std::string
SidechainRmsdFilterCreator::keyname() const { return "SidechainRmsd"; }

SidechainRmsdFilter::SidechainRmsdFilter() : filters::Filter( "SidechainRmsd"  ) {}

SidechainRmsdFilter::SidechainRmsdFilter( core::Size const res1, core::Size const res2, core::Real const rmsd_threshold ) :
		Filter( "SidechainRmsd" ), res1_( res1 ), res2_( res2 ), rmsd_threshold_( rmsd_threshold ) {}

SidechainRmsdFilter::~SidechainRmsdFilter(){}

void
SidechainRmsdFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	res1_ = core::pose::get_resnum( tag, pose, "res1_" );
	res2_ = core::pose::get_resnum( tag, pose, "res2_" );
	rmsd_threshold_ = tag->getOption<core::Real>("threshold", 1.0);
	include_backbone_ = tag->getOption<bool>("include_backbone", false);
	
	if( tag->hasOption("reference_name") ){
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data_map);
	}
	else{
		reference_pose_ = new core::pose::Pose( pose );
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() )
			core::import_pose::pose_from_pdb( *reference_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ] );
	}

	//residue_distance_filter_tracer<<"ResidueDistanceFilter with distance threshold of "<<distance_threshold_<<" between residues "<<res1_<<" and "<<res2_<<std::endl;
}

bool
SidechainRmsdFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ) );

	sidechain_rmsd_filter_tracer<<"Sidechain RMSD of residue "<<res1_<<" is "<<rmsd<<std::endl;
	return( rmsd<=rmsd_threshold_ );
}

void
SidechainRmsdFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ) );

	out<<"Sidechain RMSD of residue "<<res1_<<" is "<<rmsd<<'\n';
}

core::Real
SidechainRmsdFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ) );

	return( rmsd );
}
core::Real
SidechainRmsdFilter::compute( core::pose::Pose const & pose ) const {
	core::conformation::Residue const res_res1( pose.conformation().residue( res1_ ) );
	core::conformation::Residue const res_res2( reference_pose_->conformation().residue( res2_ ) );
	core::Real rmsd (0.0);
	
	// make sure we're comparing the same amino acid type
	runtime_assert( res_res1.aa() == res_res2.aa() );
	
	if( include_backbone_ ){
		rmsd = core::scoring::automorphic_rmsd( res_res1, res_res2, false /*superimpose*/ );
	} else {
		core::chemical::ResidueTypeSet const & res1_set( res_res1.residue_type_set() );
		core::chemical::ResidueType const & working_res1_type(
			res1_set.get_residue_type_with_variant_added( res_res1.type(), "VIRTUAL_BB" ) );
		core::conformation::ResidueOP working_res1 = 
			core::conformation::ResidueFactory::create_residue( working_res1_type );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms(
			res_res1, *working_res1, pose.conformation() );

		core::chemical::ResidueTypeSet const & res2_set( res_res2.residue_type_set() );
		core::chemical::ResidueType const & working_res2_type(
			res2_set.get_residue_type_with_variant_added( res_res2.type(), "VIRTUAL_BB" ) );
		core::conformation::ResidueOP working_res2 = 
			core::conformation::ResidueFactory::create_residue( working_res2_type );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms(
			res_res2, *working_res2, reference_pose_->conformation() );

		rmsd = core::scoring::automorphic_rmsd(*working_res1, *working_res2, false /*superimpose*/);
	}

	return( rmsd );
}

filters::FilterOP SidechainRmsdFilter::clone() const {
	return new SidechainRmsdFilter( *this );
}

filters::FilterOP SidechainRmsdFilter::fresh_instance() const{
	return new SidechainRmsdFilter();
}

}
}
