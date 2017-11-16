// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <utility/string_util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.StemFinder" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP StemFinderFilterCreator::create_filter() const { return protocols::filters::FilterOP( new StemFinder ); }

// XRW TEMP std::string
// XRW TEMP StemFinderFilterCreator::keyname() const { return "StemFinder"; }

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

StemFinder::~StemFinder() = default;

void
StemFinder::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	from_res( tag->getOption< core::Size >( "from_res", 1 ) );
	to_res( tag->getOption< core::Size >( "to_res", pose.size() ) );
	rmsd( tag->getOption< core::Real >( "rmsd", 0.7 ) );
	stems_on_sse( tag->getOption< bool >( "stems_on_sse", false ) );
	stems_are_neighbors( tag->getOption< bool >( "stems_are_neighbors", true ) );
	if ( stems_are_neighbors() ) {
		neighbor_distance( tag->getOption< core::Real >( "neighbor_distance", 4.0 ) );
		neighbor_separation( tag->getOption< core::Size >( "neighbor_separation", 10 ) );
	}
	filenames_ = utility::string_split( tag->getOption< std::string > ( "filenames" ), ',' );

	TR<<"StemFinder with options: from_res: "<<from_res()<<" to_res: "<<to_res()<<" rmsd: "<<rmsd()<<" stems_on_sse: "<<stems_on_sse();
	for ( std::string const & f : filenames_ ) {
		TR<<f<<' ';
	}
	TR<<std::endl;
}

std::string
dssp( core::pose::Pose const & pose ){
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced(); // switch to simplified H E L notation
	std::string const sec_struct( dssp.get_dssp_secstruct() );
	return( sec_struct );
}

utility::vector1< core::Size >
positions_in_secstruct( core::pose::Pose const & pose ){
	std::string const sec_struct( dssp( pose ) );

	utility::vector1< core::Size > positions;
	positions.clear();
	for ( core::Size i = 0; i < sec_struct.length(); ++i ) {
		if ( sec_struct[ i ] == 'H' || sec_struct[ i ] == 'E' ) {
			positions.push_back( i + 1 );
		}
	}
	return positions;
}

utility::vector1< core::pose::PoseOP >
load_poses( utility::vector1< std::string > const & filenames ){
	utility::vector1< core::pose::PoseOP > poses;
	poses.clear();
	TR<<"Loading "<<filenames.size()<<" poses from disk: ";
	for ( std::string const & f : filenames ) {
		core::pose::PoseOP new_pose( new core::pose::Pose );
		core::import_pose::pose_from_file( *new_pose, f , core::import_pose::PDB_file);
		poses.push_back( new_pose );
		TR<<f<<' ';
	}
	TR<<std::endl;
	return poses;
}

core::Real
atom_distance( core::pose::Pose const & p1, core::Size const r1, std::string const & a1,
	core::pose::Pose const & p2, core::Size const r2, std::string const & a2 ){
	return( p1.conformation().residue( r1 ).xyz( a1 ).distance( p2.conformation().residue( r2 ).xyz( a2 ) ) );
}

core::Real
res_res_min_distance( core::pose::Pose const &p1, core::Size const r1,
	core::pose::Pose const &p2, core::Size const r2 ){
	core::conformation::Residue const res1( p1.conformation().residue( r1 ) );
	core::conformation::Residue const res2( p2.conformation().residue( r2 ) );

	core::Size const natoms1( res1.natoms() ), natoms2( res2.natoms() );
	core::Real min_dist( 1000000.0 );
	for ( core::Size i1 = 1; i1 <= natoms1; ++i1 ) {
		for ( core::Size i2 = 1; i2 <= natoms2; ++i2 ) {
			core::Real const dist( res1.atom( i1 ).xyz().distance( res2.atom( i2 ).xyz() ) );
			if ( dist <= min_dist ) {
				min_dist = dist;
			}
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

	core::conformation::Conformation const & conf( pose.conformation() );
	std::string const template_dssp( dssp( pose ) );
	// for( core::Size i = 1; i <= pose.size(); ++i ){
	//  TR<<pose.conformation().residue( i ).name1()<<i<<' '<<template_dssp[ i - 1 ]<<'\n';
	// }

	vector1< std::string > poses_dssp;
	if ( stems_on_sse() ) { // dssp calculations are time consuming
		poses_dssp.clear();
		for ( core::pose::PoseOP p : poses ) {
			poses_dssp.push_back( dssp( *p ) );
		}
	}
	vector1< Real > distances( pose.size(), 0.0 );
	for ( Size pos = 1; pos <= pose.size(); ++pos ) {
		if ( pos <= from_res() || pos >= to_res() || template_dssp[ pos - 1 ] == 'L' ) { // position is not in secondary structure element
			distances[ pos ] = 99999999999.9;
			continue;
		}
		for ( Size pose_idx = 1; pose_idx <= poses.size(); ++pose_idx ) {
			Size position_on_target( 0 );
			Real dist( 0.0 );
			protocols::rosetta_scripts::find_nearest_res( *poses[ pose_idx ], pose, pos, position_on_target, dist );
			bool sec_struct_agreement( true );
			if ( stems_on_sse() ) {
				sec_struct_agreement = poses_dssp[ pose_idx ][ pos - 1 ] == template_dssp[ pos - 1 ];
			}
			if ( position_on_target > 0  && sec_struct_agreement ) {
				distances[ pos ] += dist;
			} else {
				distances[ pos ] += 9999999.9; // this position is taken out of consideration
				break;
			}
			if ( distances[ pos ] >= rmsd() * poses.size() ) { // past the rmsd limit for sure. No point computing any more
				distances[ pos ] += 99999999.9;
				break;
			}
		}
		if ( distances[ pos ] <= 999999.9 ) {
			distances[ pos ] = (Real) distances[pos ] / (Real) poses.size();
		}
	}
	TR<<"Candidate positions on template:\n";
	if ( stems_are_neighbors() ) { //find positions that are at most a cutoff distance from other positions
		vector1< Size > neighbor_idxs;
		neighbor_idxs.clear();
		for ( Size dist_idx = 1; dist_idx < distances.size() - neighbor_separation() + 1; ++dist_idx ) {
			if ( distances[ dist_idx ] <= 999999.9 ) {
				for ( Size j = dist_idx + neighbor_separation() - 1; j <= distances.size(); ++j ) {
					if ( std::find( neighbor_idxs.begin(), neighbor_idxs.end(), dist_idx ) != neighbor_idxs.end() &&
							std::find( neighbor_idxs.begin(), neighbor_idxs.end(), j )        != neighbor_idxs.end() ) {
						continue;
					}
					if ( distances[ j ] <= 99999.9 ) {
						Real const res_res_dist = res_res_min_distance( pose, dist_idx, pose, j );
						if ( res_res_dist <= neighbor_distance() ) {
							neighbor_idxs.push_back( dist_idx );
							neighbor_idxs.push_back( j );
						}
					}
				}
			}
		}
		std::sort( neighbor_idxs.begin(), neighbor_idxs.end() );
		auto it = std::unique( neighbor_idxs.begin(), neighbor_idxs.end() );
		neighbor_idxs.resize( std::distance( neighbor_idxs.begin(), it ) );
		for ( Size const n : neighbor_idxs ) {
			TR<<conf.residue( n ).name1()<<n<<' '<<distances[ n ]<<'\n';
		}
		TR<<std::endl;
		if ( neighbor_idxs.size() == 0 ) {
			TR<<"No pairs found"<<std::endl;
			return false;
		}
		//Find nearby pairs: A->B that cover large segments
		std::map< Size/*start stem*/, Size/*end stem*/ > stem_pairs;
		stem_pairs.clear();
		for ( Size i = 1; i <= neighbor_idxs.size() - 1; ++i ) {
			for ( Size j = i + 1; j <= neighbor_idxs.size(); ++j ) {
				if ( neighbor_idxs[ i ] + neighbor_separation() > neighbor_idxs[ j ] ) {
					continue;
				}
				Real const res_res_dist = res_res_min_distance( pose, neighbor_idxs[ i ], pose, neighbor_idxs[ j ] );
				if ( res_res_dist <= neighbor_distance() ) {
					stem_pairs[ neighbor_idxs[ i ] ] = neighbor_idxs[ j ];
					TR<<"Pair: "<<conf.residue( neighbor_idxs[ i ] ).name1()<<neighbor_idxs[ i ]<<' '<<conf.residue( neighbor_idxs[ j ] ).name1()<<neighbor_idxs[ j ]<<' '<<res_res_dist<<'\n';
				}
			}
		}
		TR<<std::endl;
		//Find contiguous triplets that cover large segments: A->B->C
		TR<<"Connected pairs:\n";
		for ( auto const & stem_pair : stem_pairs ) {
			auto const third_partner( stem_pairs.find( stem_pair.second ) );
			if ( third_partner != stem_pairs.end() ) {
				TR<<"Triplet: "<<conf.residue( stem_pair.first ).name1()<<stem_pair.first<<' '<<conf.residue( stem_pair.second ).name1()<<stem_pair.second<<' '<<conf.residue( third_partner->second ).name1()<<third_partner->second<<'\n';
			}
		}
		TR<<std::endl;
	} else {
		for ( Size dist_idx = 1; dist_idx <= distances.size(); ++dist_idx ) {
			if ( distances[ dist_idx ] <= 999999.9 ) {
				TR<<conf.residue( dist_idx ).name1()<<dist_idx<<' '<<distances[ dist_idx ]<<'\n';
			}
		}
		TR<<std::endl;
	}
	return true;
}

void
StemFinder::add_filename( std::string const & s ) {
	filenames_.push_back( s );
}

void
StemFinder::report( std::ostream &, core::pose::Pose const & ) const {
}

core::Real
StemFinder::report_sm( core::pose::Pose const & ) const {
	return( 0.0 );
}

std::string StemFinder::name() const {
	return class_name();
}

std::string StemFinder::class_name() {
	return "StemFinder";
}

void StemFinder::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("from_res", xsct_non_negative_integer, "starting template position (in rosetta numbering) in which to search for stems", "1")
		+ XMLSchemaAttribute("to_res", xsct_non_negative_integer, "ending template positions (in rosetta numbering) in which to search for stems")
		+ XMLSchemaAttribute::attribute_w_default("rmsd", xsct_real, "cutoff for the average rmsd between a given position in the template and all of the closest positions in the homologs.", "0.7")
		+ XMLSchemaAttribute::attribute_w_default("stems_on_sse", xsct_rosetta_bool, "demand that in each of the homologs the candidate stems are on 2ary structural elements", "false")
		+ XMLSchemaAttribute::attribute_w_default("stems_are_neighbors", xsct_rosetta_bool, "should we eliminate stems that are farther than neighbor_distance from one another?", "true")
		+ XMLSchemaAttribute::attribute_w_default("neighbor_distance", xsct_real, "minimal atomic distance between any pair of atoms on each of the residues", "4.0")
		+ XMLSchemaAttribute::attribute_w_default("neighbor_seperation", xsct_non_negative_integer, "minimal aa separation between candidate stem sites", "10")
		+ XMLSchemaAttribute::required_attribute("filenames", xs_string, "PDB structures that are well aligned to the template");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Compare a set of homologous but structurally heterogeneous PDBs to a template PDB and find structurally highly conserved sites that can serve as stems for splicing segments", attlist );
}

std::string StemFinderFilterCreator::keyname() const {
	return StemFinder::class_name();
}

protocols::filters::FilterOP
StemFinderFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new StemFinder );
}

void StemFinderFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StemFinder::provide_xml_schema( xsd );
}


}
}
