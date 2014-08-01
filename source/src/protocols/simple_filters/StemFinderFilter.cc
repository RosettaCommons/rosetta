// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/StemFinderFilter.cc
/// @brief
/// @author Ravit Netzer & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/StemFinderFilter.hh>
#include <protocols/simple_filters/StemFinderFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <utility/vector1.hh>
#include <string>
#include <boost/foreach.hpp>
#include <utility/string_util.hh>

namespace protocols{
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.StemFinder" );

protocols::filters::FilterOP
StemFinderFilterCreator::create_filter() const { return new StemFinder; }

std::string
StemFinderFilterCreator::keyname() const { return "StemFinder"; }

//default ctor
StemFinder::StemFinder() :
protocols::filters::Filter( "StemFinder" ),
from_res_( 0 ),
to_res_( 0 ),
rmsd_( 0.7 ),
stems_on_sse_( false ),
stems_are_neighbors_( true ),
neighbor_distance_( 4.0 ),
neighbor_separation_( 10 )
{
	filenames_.clear();
}

StemFinder::~StemFinder() {}

void
StemFinder::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	from_res( tag->getOption< core::Size >( "from_res", 1 ) );
	to_res( tag->getOption< core::Size >( "to_res", pose.total_residue() ) );
	rmsd( tag->getOption< core::Real >( "rmsd", 0.7 ) );
	stems_on_sse( tag->getOption< bool >( "stems_on_sse", false ) );
	stems_are_neighbors( tag->getOption< bool >( "stems_are_neighbors", true ) );
	if( stems_are_neighbors() ){
		neighbor_distance( tag->getOption< core::Real >( "neighbor_distance", 4.0 ) );
		neighbor_separation( tag->getOption< core::Size >( "neighbor_separation", 10 ) );
	}
	filenames_ = utility::string_split( tag->getOption< std::string > ( "filenames" ), ',' );

	TR<<"StemFinder with options: from_res: "<<from_res()<<" to_res: "<<to_res()<<" rmsd: "<<rmsd()<<" stems_on_sse: "<<stems_on_sse();
	BOOST_FOREACH( std::string const f, filenames_ )
		TR<<f<<' ';
	TR<<std::endl;
}

utility::vector1< core::Size >
positions_in_secstruct( core::pose::Pose const & pose ){
  core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced(); // switch to simplified H E L notation
  std::string const sec_struct( dssp.get_dssp_secstruct() );

	utility::vector1< core::Size > positions;
	positions.clear();
	for( core::Size i = 0; i < sec_struct.length(); ++i ){
		if( sec_struct[ i ] == 'H' || sec_struct[ i ] == 'E' )
			positions.push_back( i + 1 );
	}
	return positions;
}

std::string
dssp( core::pose::Pose const & pose ){
  core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced(); // switch to simplified H E L notation
  std::string const sec_struct( dssp.get_dssp_secstruct() );
	return( sec_struct );
}

utility::vector1< core::pose::PoseOP >
load_poses( utility::vector1< std::string > const filenames ){
	utility::vector1< core::pose::PoseOP > poses;
	poses.clear();
	TR<<"Loading "<<filenames.size()<<" poses from disk: ";
	BOOST_FOREACH( std::string const f, filenames ){
		core::pose::PoseOP new_pose( new core::pose::Pose );
	  core::import_pose::pose_from_pdb( *new_pose, f );
		poses.push_back( new_pose );
		TR<<f<<' ';
	}
	TR<<std::endl;
	return poses;
}

core::Real
atom_distance( core::pose::Pose const & p1, core::Size const r1, std::string const a1,
							 core::pose::Pose const & p2, core::Size const r2, std::string const a2 ){
	return( p1.conformation().residue( r1 ).xyz( a1 ).distance( p2.conformation().residue( r2 ).xyz( a2 ) ) );
}

core::Real
res_res_min_distance( core::pose::Pose const &p1, core::Size const r1,
											core::pose::Pose const &p2, core::Size const r2 ){
	core::conformation::Residue const res1( p1.conformation().residue( r1 ) );
	core::conformation::Residue const res2( p2.conformation().residue( r2 ) );

	core::Size const natoms1( res1.natoms() ), natoms2( res2.natoms() );
	core::Real min_dist( 1000000.0 );
	for( core::Size i1 = 1; i1 <= natoms1; ++i1 ){
		for( core::Size i2 = 1; i2 <= natoms2; ++i2 ){
			core::Real const dist( res1.atom( i1 ).xyz().distance( res2.atom( i2 ).xyz() ) );
			if( dist <= min_dist )
				min_dist = dist;
		}
	}
	return min_dist;
}

bool
StemFinder::apply( core::pose::Pose const & pose ) const {
	using utility::vector1;
	using core::Size;
	using core::Real;
	using namespace core::pose;
	vector1< PoseOP > poses( load_poses( filenames_ ) );

	core::conformation::Conformation const conf( pose.conformation() );
	std::string const template_dssp( dssp( pose ) );
//	for( core::Size i = 1; i <= pose.total_residue(); ++i ){
//		TR<<pose.conformation().residue( i ).name1()<<i<<' '<<template_dssp[ i - 1 ]<<'\n';
//	}

	vector1< std::string > poses_dssp;
	if( stems_on_sse() ){ // dssp calculations are time consuming
		poses_dssp.clear();
		BOOST_FOREACH( core::pose::PoseOP p, poses )
			poses_dssp.push_back( dssp( *p ) );
	}
	vector1< Real > distances( pose.total_residue(), 0.0 );
	for( Size pos = 1; pos <= pose.total_residue(); ++pos ){
		if( pos <= from_res() || pos >= to_res() || template_dssp[ pos - 1 ] == 'L' ){ // position is not in secondary structure element
			distances[ pos ] = 99999999999.9;
			continue;
		}
		for( Size pose_idx = 1; pose_idx <= poses.size(); ++pose_idx ){
			Size position_on_target( 0 );
			Real dist( 0.0 );
			protocols::rosetta_scripts::find_nearest_res( *poses[ pose_idx ], pose, pos, position_on_target, dist );
			bool sec_struct_agreement( true );
			if( stems_on_sse() )
				sec_struct_agreement = poses_dssp[ pose_idx ][ pos - 1 ] == template_dssp[ pos - 1 ];
			if( position_on_target > 0  && sec_struct_agreement )
				distances[ pos ] += dist;
			else{
				distances[ pos ] += 9999999.9; // this position is taken out of consideration
				break;
			}
			if( distances[ pos ] >= rmsd() * poses.size() ){ // past the rmsd limit for sure. No point computing any more
				distances[ pos ] += 99999999.9;
				break;
			}
		}
		if( distances[ pos ] <= 999999.9 )
			distances[ pos ] = (Real) distances[pos ] / (Real) poses.size();
	}
	TR<<"Candidate positions on template:\n";
	if( stems_are_neighbors() ){//find positions that are at most a cutoff distance from other positions
		vector1< Size > neighbor_idxs;
		neighbor_idxs.clear();
		for( Size dist_idx = 1; dist_idx < distances.size() - neighbor_separation() + 1; ++dist_idx ){
			if( distances[ dist_idx ] <= 999999.9 ){
				for( Size j = dist_idx + neighbor_separation() - 1; j <= distances.size(); ++j ){
					if( std::find( neighbor_idxs.begin(), neighbor_idxs.end(), dist_idx ) != neighbor_idxs.end() &&
							std::find( neighbor_idxs.begin(), neighbor_idxs.end(), j )        != neighbor_idxs.end() )
						continue;
					if( distances[ j ] <= 99999.9 ){
						Real const res_res_dist = res_res_min_distance( pose, dist_idx, pose, j );
						if( res_res_dist <= neighbor_distance() ){
							neighbor_idxs.push_back( dist_idx );
							neighbor_idxs.push_back( j );
						}
					}
				}
			}
		}
		std::sort( neighbor_idxs.begin(), neighbor_idxs.end() );
		vector1< Size >::iterator it = std::unique( neighbor_idxs.begin(), neighbor_idxs.end() );
		neighbor_idxs.resize( std::distance( neighbor_idxs.begin(), it ) );
		BOOST_FOREACH( Size const n, neighbor_idxs )
			TR<<conf.residue( n ).name1()<<n<<' '<<distances[ n ]<<'\n';
		TR<<std::endl;
		if( neighbor_idxs.size() == 0 ){
			TR<<"No pairs found"<<std::endl;
			return false;
		}
		//Find nearby pairs: A->B that cover large segments
		std::map< Size/*start stem*/, Size/*end stem*/ > stem_pairs;
		stem_pairs.clear();
		for( Size i = 1; i <= neighbor_idxs.size() - 1; ++i ){
			for( Size j = i + 1; j <= neighbor_idxs.size(); ++j ){
				if( neighbor_idxs[ i ] + neighbor_separation() > neighbor_idxs[ j ] )
					continue;
				Real const res_res_dist = res_res_min_distance( pose, neighbor_idxs[ i ], pose, neighbor_idxs[ j ] );
				if( res_res_dist <= neighbor_distance() ){
					stem_pairs[ neighbor_idxs[ i ] ] = neighbor_idxs[ j ];
					TR<<"Pair: "<<conf.residue( neighbor_idxs[ i ] ).name1()<<neighbor_idxs[ i ]<<' '<<conf.residue( neighbor_idxs[ j ] ).name1()<<neighbor_idxs[ j ]<<' '<<res_res_dist<<'\n';
				}
			}
		}
		TR<<std::endl;
		//Find contiguous triplets that cover large segments: A->B->C
		TR<<"Connected pairs:\n";
		for( std::map< Size, Size >::const_iterator mit = stem_pairs.begin(); mit != stem_pairs.end(); ++mit ){
			std::map< Size, Size >::const_iterator third_partner( stem_pairs.find( mit->second ) );
			if( third_partner != stem_pairs.end() )
					TR<<"Triplet: "<<conf.residue( mit->first ).name1()<<mit->first<<' '<<conf.residue( mit->second ).name1()<<mit->second<<' '<<conf.residue( third_partner->second ).name1()<<third_partner->second<<'\n';
		}
		TR<<std::endl;
	}
	else {
		for( Size dist_idx = 1; dist_idx <= distances.size(); ++dist_idx ){
			if( distances[ dist_idx ] <= 999999.9 )
				TR<<conf.residue( dist_idx ).name1()<<dist_idx<<' '<<distances[ dist_idx ]<<'\n';
		}
		TR<<std::endl;
	}
	return true;
}

void
StemFinder::report( std::ostream &, core::pose::Pose const & ) const {
}

core::Real
StemFinder::report_sm( core::pose::Pose const & ) const {
	return( 0.0 );
}

}
}
