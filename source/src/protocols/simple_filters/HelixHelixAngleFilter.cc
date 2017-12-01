// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/HelixHelixAngleFilter.cc
/// @brief  HelixHelixAngleFilter computes wither the crossing angle or distance between two TMs
/// @details HelixHelixAngleFilter computes 1 of 3 parameters: 1. crossing angle 2. atomic distance 3. vector distance
/// 1. the crossing angle at the closest point between two vectors representing the two helices
/// 2. the shortest distance between two atoms on either TM
/// 3. the shortest distance (point of approach) between the vectors representing the two helices
/// @author Jonathan Weinstein


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
#include <core/kinematics/FoldTree.hh>
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
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/HelixHelixAngleFilter.hh>
#include <protocols/simple_filters/HelixHelixAngleFilterCreator.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <fstream>

#include <core/scoring/dssp/Dssp.hh>
#include <numeric/conversions.hh>
#include <basic/svd/SVD_Solver.hh>
#include <numeric/PCA.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.HelixHelixAngleFilter" );

protocols::filters::FilterOP
HelixHelixAngleFilterCreator::create_filter() const { return protocols::filters::FilterOP( new HelixHelixAngleFilter ); }

std::string
HelixHelixAngleFilterCreator::keyname() const { return "HelixHelixAngle"; }

HelixHelixAngleFilter::~HelixHelixAngleFilter(){}

void
HelixHelixAngleFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	start_helix_1_ = tag->getOption< core::Size >( "start_helix_1", 0 );
	start_helix_2_ = tag->getOption< core::Size >( "start_helix_2", 0  );
	end_helix_1_ = tag->getOption< core::Size >( "end_helix_1", 0  );
	end_helix_2_ = tag->getOption< core::Size >( "end_helix_2", 0  );
	angle_min_ = tag->getOption< core::Real >( "angle_min", 40.0 );
	angle_max_ = tag->getOption< core::Real >( "angle_max", 100.0 );
	dist_min_ = tag->getOption< core::Real >( "dist_min", 0.0 );
	dist_max_ = tag->getOption< core::Real >( "dist_max", 5.0 );
	angle_or_dist_ = tag->getOption< std::string >( "angle_or_dist", "angle" );
	dist_by_atom_ = tag->getOption< bool >( "dist_by_atom", true );

	if ( tag->hasOption( "helix_num_1" ) ) {
		helix_num_1_ = tag->getOption< core::Size > ( "helix_num_1", 1 );
		helix_num_2_ = tag->getOption< core::Size > ( "helix_num_2", 2 );
		core::scoring::dssp::Dssp dssp( pose );
		by_helices_ = true;
		std::string dssp_str( dssp.get_dssp_secstruct() );
		utility::vector1 < std::pair< core::Size, core::Size > > helices;
		helices.clear();
		bool open_helix(false);
		core::Size start( 0 );
		core::Size end( 0 );
		int chain( pose.chain( 1 ) );
		int helix_chain( 0 );
		for ( core::Size i = 0; i < dssp_str.size(); ++i ) {
			if ( helix_chain != 0 && helix_chain != (int) pose.chain( i + 1 ) && open_helix ) {
				open_helix = false;
				end = i + 1;
				helices.push_back( std::pair< core::Size, core::Size>( start, end ) );
			}
			if ( dssp_str[ i ] != 'H' and ! open_helix ) {
				chain = pose.chain( i + 1 );
				continue;
			}
			if ( dssp_str[ i ] != 'H' && open_helix ) {
				end = i + 1;
				helices.push_back( std::pair< core::Size, core::Size>( start, end ) );
				open_helix = false;
				chain = pose.chain( i + 1 );
				continue;
			}
			if ( dssp_str[ i ] == 'H' && open_helix ) {
				chain = pose.chain( i + 1 );
				continue;
			}
			if ( dssp_str[ i ] == 'H' && ! open_helix ) {
				start = i;
				helix_chain = pose.chain( i + 1 );
				open_helix = true;
				continue;
			}
		}
		if ( helices.size() == 0 ) {
			utility_exit_with_message("DSSP has no helices, defne actual residue indices using start_helix_1, end_helix_1 etc.");
		}
		TR << chain << " " << helices << std::endl; // necessary for compilation.
		start_helix_1_ = helices[ helix_num_1_ ].first;
		end_helix_1_ = helices[ helix_num_1_ ].second;
		start_helix_2_ = helices[ helix_num_2_ ].first;
		end_helix_2_ = helices[ helix_num_2_ ].second;
	} else {
		by_helices_ = false;
	}
	if ( start_helix_1_ >= end_helix_1_ || start_helix_2_ >= end_helix_2_ ) utility_exit_with_message( "illegal helices" );
	if ( angle_min_ >= angle_max_ ) utility_exit_with_message( "illegal min and max angles" );
	if ( dist_min_ >= dist_max_ ) utility_exit_with_message( "illegal min and max distances" );
}

bool
HelixHelixAngleFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real angle = compute( pose );
	bool status( false );
	if ( angle_or_dist_ == "angle" ) {
		status = ( angle >= angle_min_ && angle <= angle_max_ );
	} else {
		status = ( angle >= dist_min_ && angle <= dist_max_ );
	}
	return status;
}

void
HelixHelixAngleFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	core::Real angle = compute( pose );
	if ( angle_or_dist_ == "angle" ) {
		out << "the angle between the helices is " << angle << std::endl;
	} else {
		out << "the distance between the helices is " << angle << std::endl;
	}
}

core::Real
HelixHelixAngleFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real angle = compute( pose );
	return( angle );
}

core::Real
HelixHelixAngleFilter::compute( core::pose::Pose const & pose ) const
{
	if ( angle_or_dist_ == "dist" && dist_by_atom_ ) {
		core::Real dist = calc_shortest_dist_by_atoms( pose );
		return( dist );
	}
	TR << "helices defined " << start_helix_1_ << "->" << end_helix_1_ << " and " << start_helix_2_ << "->" << end_helix_2_ << std::endl;
	utility::vector1< numeric::xyzVector < core::Real > > l1( find_helix_vector( pose, start_helix_1_, end_helix_1_ ) );
	utility::vector1< numeric::xyzVector < core::Real > > l2( find_helix_vector( pose, start_helix_2_, end_helix_2_ ) );

	std::pair<numeric::xyzVector<core::Real>, numeric::xyzVector<core::Real>> closest_pnts( find_closest_pnts( l1, l2 ) );

	if ( angle_or_dist_ == "dist" ) {
		core::Real distance( closest_pnts.first.distance( closest_pnts.second ) );
		return( distance );
	} else {

		// calculate the crossing angle at the point of approach
		core::Real hoist_angle_degrees = numeric::angle_degrees( l1[1], l1[2], l2[1] );
		core::Real packing_angle_degrees = numeric::angle_degrees( l1[1], l1[2], l2[1], l2[2] );
		core::Real meridian_angle_degrees = numeric::dihedral_degrees( l1[1], l1[2], l2[1], l2[2] );

		core::Real cross_angle = numeric::dihedral_degrees( l1[1], closest_pnts.first, closest_pnts.second, l2[1] );

		TR << "hoist_angle " << hoist_angle_degrees << std::endl;
		TR << "packi_angle " << packing_angle_degrees << std::endl;
		TR << "merid_angle " << meridian_angle_degrees << std::endl;
		TR << "cross " << cross_angle << std::endl;
		return( cross_angle );
	}
}

///
// @brief iterate all atoms in both TMs, return the minimal distance between them
core::Real
HelixHelixAngleFilter::calc_shortest_dist_by_atoms( core::pose::Pose const & pose ) const
{
	core::Real dist( 100000000 );
	core::Real temp_dist( 0.0 );
	core::Size sel_i( 0 ), sel_j( 0 ), sel_i_atm( 0 ), sel_j_atm( 0 );
	for ( core::Size i=start_helix_1_; i<=end_helix_1_; ++i ) {
		for ( core::Size j=start_helix_2_; j<=end_helix_2_; ++j ) {
			for ( core::Size i_atm=1; i_atm <= pose.residue( i ).nheavyatoms(); ++i_atm ) {
				if ( pose.residue( i ).atom_is_backbone( i_atm ) ) {
					for ( core::Size j_atm=1; j_atm<=pose.residue( j ).nheavyatoms(); ++j_atm ) {
						if ( pose.residue( j ).atom_is_backbone( j_atm ) ) {
							temp_dist = pose.residue( i ).xyz( i_atm ).distance( pose.residue( j ).xyz( j_atm ) );
							if ( temp_dist < dist ) {
								dist = temp_dist;
								sel_i = i;
								sel_j = j;
								sel_i_atm = i_atm;
								sel_j_atm = j_atm;
							}
						}
					}
				}
			}
		}
	}
	TR << "the distance between atoms:" << std::endl;
	TR << "res " << sel_i << " atom " << pose.residue( sel_i ).atom_name( sel_i_atm ) << std::endl;
	TR << "res " << sel_j << " atom " << pose.residue( sel_j ).atom_name( sel_j_atm ) << std::endl;
	return( dist );
}

// @brief find point of approach (closest points) on the vectors, and return the distance
std::pair< numeric::xyzVector< core::Real >, numeric::xyzVector< core::Real > >
HelixHelixAngleFilter::find_closest_pnts( utility::vector1< numeric::xyzVector < core::Real > > l1,
	utility::vector1< numeric::xyzVector < core::Real > > l2 ) const {
	std::pair< numeric::xyzVector<core::Real>, numeric::xyzVector<core::Real>> pnts;

	numeric::xyzVector< core::Real > u( l1[2] - l1[1] );
	numeric::xyzVector< core::Real > v( l2[2] - l2[1] );
	numeric::xyzVector< core::Real > w( l1[1] - l2[1] );

	core::Real a( u.dot( u ) );
	core::Real b( u.dot( v ) );
	core::Real c( v.dot( v ) );
	core::Real d( u.dot( w ) );
	core::Real e( v.dot( w ) );
	core::Real D( ( a * c ) - ( b * b ) );

	core::Real sc( 0.0 );
	core::Real tc( 0.0 );

	if ( D < 0.000001 ) {
		sc = 0.0;
		tc = ( b > c ? d / b : e / c );
	} else {
		sc = ( ( b * e ) - ( c * d ) ) / D;
		tc = ( ( a * e ) - ( b * d ) ) / D;
	}

	numeric::xyzVector< core::Real > v1( l1[1] + sc * ( l1[2] - l1[1] ) );
	numeric::xyzVector< core::Real > v2( l2[1] + tc * ( l2[2] - l2[1] ) );
	TR << "v1 " << v1.x() << " " << v1.y() << " " << v1.z() << std::endl;
	TR << "v2 " << v2.x() << " " << v2.y() << " " << v2.z() << std::endl;

	pnts.first = v1;
	pnts.second = v2;
	return( pnts );
}


// @brief claculate a vector representing the TM (a 3D linear regression)
utility::vector1< numeric::xyzVector < core::Real > >
HelixHelixAngleFilter::find_helix_vector( core::pose::Pose const & pose, core::Size start, core::Size end ) const
{
	utility::vector1< numeric::xyzVector <core::Real > > vec;
	vec.empty();

	utility::vector1< numeric::xyzVector< core::Real> > bb_coords; bb_coords.clear();

	for ( core::Size i=start; i <= end; i++ ) {
		core::conformation::Residue rsd( pose.residue( i ) );
		for ( core::Size j=1; j <= rsd.nheavyatoms(); ++j ) {
			if ( rsd.atom_is_backbone( j ) ) {
				core::Vector const & bbvec( rsd.xyz( j ) );
				bb_coords.push_back( bbvec );
			}
		}
	}
	numeric::xyzVector<core::Real> com = center_of_mass( bb_coords );
	numeric::xyzVector<core::Real> first_principal_component = numeric::first_principal_component( bb_coords );

	numeric::xyzVector<core::Real> com_principal_component = first_principal_component += com;

	//p0 is point on principal component vector closest to N-term
	numeric::xyzVector<core::Real> p0 = numeric::closest_point_on_line( com, com_principal_component, bb_coords[1]);

	//p1 is point on principal component vector closest to C-term
	numeric::xyzVector<core::Real> p1 = numeric::closest_point_on_line( com, com_principal_component, bb_coords.back());
	//TR << "com " << com.x() << " " << com.y() << " " << com.z() << std::endl;
	//TR << "p0 " << p0.x() << " " << p0.y() << " " << p0.z() << std::endl;
	//TR << "p1 " << p1.x() << " " << p1.y() << " " << p1.z() << std::endl;
	vec.push_back( p0 );
	vec.push_back( p1 );
	return( vec );
}

// return a pair of residues, which are the contact point between the helices
std::pair< core::Size, core::Size >
HelixHelixAngleFilter::find_closest_res( core::pose::Pose const & pose, core::Size s1, core::Size e1, core::Size s2, core::Size e2  ) const
{
	std::pair< core::Size, core::Size > result( 0, 0 );
	core::Real distance(0);
	core::Real min_dist(10000);
	for ( core::Size res_i=s1; res_i <= e1; res_i++ ) {
		for ( core::Size res_j=s2; res_j <= e2; res_j++ ) {
			// residue_i is from 1st span, residue_j is from 2nd span
			distance = pose.residue(res_i).atom( pose.residue(res_i).nbr_atom() ).xyz().distance_squared( pose.residue(res_j).atom( pose.residue(res_j).nbr_atom() ).xyz() );

			if ( distance < min_dist ) {
				min_dist = distance;
				result.first = res_i;
				result.second = res_j;
			}

		}
	}
	TR << "found the closest residues to be " << result.first << " and " << result.second << std::endl;
	return( result );
}

std::string HelixHelixAngleFilter::name() const {
	return class_name();
}

std::string HelixHelixAngleFilter::class_name() {
	return "HelixHelixAngle";
}

void HelixHelixAngleFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "start_helix_1" , xsct_non_negative_integer , "residue where helix 1 starts" , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "start_helix_2" , xsct_non_negative_integer , "residue where helix 2 starts" , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "end_helix_1" , xsct_non_negative_integer , "residue where helix 1 ends" , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "end_helix_2" , xsct_non_negative_integer , "residue where helix 2 ends" , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "min_angle" , xsct_real , "minimal angle to pass" , "40" )
		+ XMLSchemaAttribute::attribute_w_default( "max_angle" , xsct_real , "maximal angle to pass" , "100" )
		+ XMLSchemaAttribute::attribute_w_default( "dist_min" , xsct_real , "minimal distance to pass" , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "dist_max" , xsct_real , "maximal distance to pass" , "5" )
		+ XMLSchemaAttribute::attribute_w_default( "angle_or_dist" , xs_string , "calculate angle or distance" , "angle" )
		+ XMLSchemaAttribute::attribute_w_default( "dist_by_atom" , xsct_rosetta_bool , "whether to calcualte distance between the closeset BB atoms, or the helix vectors", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "helix_num_2" , xsct_non_negative_integer , "choose helices by DSSP helix number" , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "helix_num_1" , xsct_non_negative_integer , "choose helices by DSSP helix number" , "2" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "calculate the distance or angle between helices", attlist );
}

void HelixHelixAngleFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HelixHelixAngleFilter::provide_xml_schema( xsd );
}

}
}
