// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/ParametricSheet.cc
/// @brief Parametric generation of beta sheets
/// @detailed
/// @author Tom Linsky


//Unit Headers
#include <devel/denovo_design/ParametricSheet.hh>

//Project Headers

//Core Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>

//Protocol Headers
#include <protocols/idealize/IdealizeMover.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>

//Basic Headers
#include <basic/Tracer.hh>

//Utility Headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

//C++ Headers
#include <math.h>
#include <fstream>

#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("devel.denovo_design.ParametricSheet");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace devel {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////


///  ---------------------------------------------------------------------------------
///  ParametricSheet main code:
///  ---------------------------------------------------------------------------------

/// @brief initialize static const data
core::Size const ParametricSheet::default_strand_length_ = 6;
int const ParametricSheet::default_register_shift_ = 0;
std::string const ParametricSheet::default_orientation_ = "P";

/// @brief default constructor
ParametricSheet::ParametricSheet() :
	utility::pointer::ReferenceCount(),
	strand_dist_( 4.5 ),
	c_dist_( 3.4 ),
	twist_( 0.0 ),
	c_coil_( 0.0 ),
	h_coil_( 0.0 ),
	ztwist_( 0.0 ),
	// 0.19 was determined to be optimal for flat, parallel sheets
	strand_rotation_( 0.19 )
{
	ca_coords_.clear();
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
ParametricSheet::~ParametricSheet() {}

void
ParametricSheet::parse_tags( utility::vector1< utility::tag::TagCOP > const & tags )
{
	// get strand info from subtags
	strand_name_to_index_.clear();
	strand_data_.clear();
	for ( core::Size i=1; i<=tags.size(); ++i ) {
		std::string const type = tags[i]->getName(); // strand is the only allowed subtag
		if ( type == "Strand" ) {
			if ( tags[i]->hasOption( "name" ) ) {
				std::string const name( tags[i]->getOption< std::string >( "name" ) );
				strand_name_to_index_[ name ] = i;
				strand_data_.push_back( StrandData( name,
							tags[i]->getOption< core::Size >( "length", default_strand_length_ ),
							tags[i]->getOption< core::Size >( "register_shift", default_register_shift_ ),
							tags[i]->getOption< std::string >( "orientation", default_orientation_ ) ) );
			} else {
				utility_exit_with_message( "Strand name option is not specified, and must be." );
			}
			TR << "tag type is " << type << std::endl;
		} else {
			utility_exit_with_message( type + " is not a valid subtag for BuildSheet" );
		}
	}
	TR << " Strand info: " << std::endl;
	for ( core::Size i=1; i<=strand_data_.size(); ++i ) {
		TR << "Length: " << strand_data_[i].length << " Reg. Shift: " << strand_data_[i].register_shift << " Orient: " << strand_data_[i].orientation << std::endl;
	}
}

/// @brief setup the parameters based by parsing an input pose which contains ONLY the sheet
void
ParametricSheet::from_pose( core::pose::Pose const & pose )
{
	// read in strand parameters
	utility::vector1< core::pose::PoseOP > strands( pose.split_by_chain() );
	char chain_name( 'A' - (char)1 );
	strand_name_to_index_.clear();
	strand_data_.clear();
	for ( core::Size i=1; i<=strands.size(); ++i ) {
		std::string const name( boost::lexical_cast< std::string >( chain_name + i ) );
		strand_name_to_index_[ name ] = i;
		strand_data_.push_back( StrandData( name,
																				strands[i]->total_residue(),
																				0,
																				"P" ) );
	}

	// initialize ca_coords from the pose
	init_ca_coords();
	for ( core::Size i=2; i<=strands.size()+1; ++i ) {
		for ( core::Size j=3; j<=strands[i-1]->total_residue()+2; ++j ) {
			ca_coords_[i][j] = strands[i-1]->residue(j-2).xyz( "CA" );
		}
		if ( i - 2 >= 1 ) {
			create_ca_point( i, i-1, i-2, 2, 3, 4 );
			create_ca_point( i, i-1, i-2, 1, 2, 3 );
			create_ca_point( i, i-1, i-2, ca_coords_[i].size()-1, ca_coords_[i].size()-2, ca_coords_[i].size()-3 );
			create_ca_point( i, i-1, i-2, ca_coords_[i].size(), ca_coords_[i].size()-1, ca_coords_[i].size()-2 );
		} else if ( i + 2 <= ca_coords_.size() ) {
			create_ca_point( i, i+1, i+2, 2, 3, 4 );
			create_ca_point( i, i+1, i+2, 1, 2, 3 );
			create_ca_point( i, i+1, i+2, ca_coords_[i].size()-1, ca_coords_[i].size()-2, ca_coords_[i].size()-3 );
			create_ca_point( i, i+1, i+2, ca_coords_[i].size(), ca_coords_[i].size()-1, ca_coords_[i].size()-2 );
		}
	}
	// create dummy positions
	generate_strand( 1, 2, 3 );
	generate_strand( strands.size() + 2, strands.size() + 1, strands.size() );
	print_vector_of_vectors( ca_coords_ );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Basic Accessors
core::Size ParametricSheet::num_strands() const { return strand_data_.size(); }

core::Size ParametricSheet::register_shift( core::Size const strand ) const
{
	assert( strand );
	assert( strand <= strand_data_.size() );
	return strand_data_[ strand ].register_shift;
}

core::Size ParametricSheet::strand_len( core::Size const strand ) const
{
	assert( strand );
	assert( strand <= strand_data_.size() );
	return strand_data_[ strand ].length;
}

std::string const & ParametricSheet::orientation( core::Size const strand ) const
{
	assert( strand );
	assert( strand <= strand_data_.size() );
	return strand_data_[ strand ].orientation;
}

core::Size ParametricSheet::ca_coords_size() const { return ca_coords_.size(); }

core::Size ParametricSheet::ca_coords_size( core::Size const strand ) const
{
	assert( strand );
	assert( strand <= ca_coords_.size() );
	return ca_coords_[ strand ].size();
}

core::Vector const & ParametricSheet::ca_coords( core::Size const strand, core::Size const resi ) const
{
	return ca_coords_[strand][resi];
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Basic Mutators
void ParametricSheet::c_dist( core::Real const cdist_val ) { c_dist_ = cdist_val; }
void ParametricSheet::strand_dist( core::Real const sdist_val ) { strand_dist_ = sdist_val; }
void ParametricSheet::twist( core::Real const twist_val ) { twist_ = twist_val; }
void ParametricSheet::c_coil( core::Real const c_coil_val ) { c_coil_ = c_coil_val; }
void ParametricSheet::h_coil( core::Real const h_coil_val ) { h_coil_ = h_coil_val; }
void ParametricSheet::strand_rotation( core::Real const strand_rotation_val ) { strand_rotation_ = strand_rotation_val; }

void
ParametricSheet::strand_data( utility::vector1< StrandData > const & strand_data_val )
{
	strand_data_.clear();
	strand_data_ = strand_data_val;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper Functions
void
print_vector( numeric::xyzVector< core::Real > const & vec, std::string const & name = "" )
{
	TR << name << " [ " << vec.x() << ", " << vec.y() << ", " << vec.z() << " ] " << std::endl;
}

void
print_vector_of_vectors( utility::vector1< utility::vector1< core::Vector > > const & vec ) {
	for ( core::Size i=1; i<=vec.size(); ++i ) {
		for ( core::Size j=1; j<=vec[i].size(); ++j ) {
			TR << i << "," << j << " " << vec[i][j].x() << " " << vec[i][j].y() << " " << vec[i][j].z() << std::endl;
		}
	}
}

void
print_matrix( numeric::xyzMatrix< core::Real > const & mat )
{
	for ( core::Size i=1; i<=3; ++i ) {
		print_vector( mat.row( i ) );
	}
}

/// @brief given two previous h points, determine the h direction
core::Vector
ParametricSheet::get_h_axis( core::Vector const & p1,
														 core::Vector const & p2,
														 core::Vector const & /*c_axis*/,
														 core::Size const /*resi*/ ) const
{
	core::Vector h_point_dir( p1 - p2 );
	//return rotate_about_axis( c_axis, up_down( resi )*strand_rotation_ )*h_point_dir;
	return h_point_dir;
}

/// @brief creates a vector for a new residue in the sheet based on the previous c and h vectors
void
ParametricSheet::new_strand( core::Size const strand_num,
		core::Size const prev_strand,
		core::Size const prev2_strand,
		core::Size const center_resi )
{
	runtime_assert( strand_num );
	runtime_assert( strand_num <= ca_coords_.size() );

	std::fill( ca_coords_[strand_num].begin(), ca_coords_[strand_num].end(), core::Vector( -999.9, -999.9, -999.9 ) );

	int multiplier( 1 );
	if ( prev_strand > strand_num ) {
		multiplier = -1;
	}

	// setup h axis reference vector
	core::Vector href1( 0.0, 0.0, 0.0 );
	core::Vector href2( 0.0, -strand_dist_*multiplier*cos(strand_rotation_), -strand_dist_*multiplier*sin(strand_rotation_) );
	bool const prev_exists( check_strand_and_resi_values( prev_strand, center_resi ) );
	bool const prev2_exists( check_strand_and_resi_values( prev2_strand, center_resi ) );
	if ( prev_exists && prev2_exists ) {
		href2 = ca_coords_[prev2_strand][center_resi];
		href1 = ca_coords_[prev_strand][center_resi];
		TR << "Both exist" << std::endl;
	} else if ( prev_exists ) {
		href1 = ca_coords_[prev_strand][center_resi];
		href2 = href1 + href2;
		TR << "Prev exists" << std::endl;
	}
	// setup c axis reference vector
	core::Vector c_vec( 3.8021, 0.0, 0.0 );

	print_vector( href1, "href1" );
	print_vector( href2, "href2" );
	core::Vector h_vec( get_h_axis( href1, href2, c_vec, 1 ) );
	c_vec.normalize();
	h_vec.normalize();

	// setup orthogonal axis reference vector
	core::Vector cross( c_vec.cross( h_vec ).normalize() );

	// the new vector will have c,h coordinates -- this vector will convert that into xyz
	numeric::xyzMatrix< core::Real > const rotate( numeric::xyzMatrix< core::Real >::rows( c_vec.x(), h_vec.x(), cross.x(),
				c_vec.y(), h_vec.y(), cross.y(),
				c_vec.z(), h_vec.z(), cross.z() ) );

	// create the initial strand point and rotate it into place
	core::Vector const new_point_ch( 0.0, strand_dist_*cos(h_coil_), strand_dist_*sin(h_coil_) );
	print_vector( new_point_ch, "new_point_ch" );
	ca_coords_[strand_num][center_resi] = rotate*new_point_ch + href1;
	TR << "Added point: " << std::endl;
	print_vector( ca_coords_[strand_num][center_resi] );
}

void
ParametricSheet::test_matrix( utility::vector1< utility::vector1< core::Vector > > const & vec ) const
{
	for ( core::Size i=1; i<=vec.size(); ++i ) {
		for ( core::Size j=2; j<=vec[i].size(); ++j ) {
			// should be 3.8 angstroms from predecessor
			if ( vec[i][j].x() > -999 && vec[i][j-1].x() > -999 ) {
				core::Real const dist( vec[i][j].distance( vec[i][j-1] ) );
				//TR << i << "," << j << " dist = " << dist << std::endl;
				if ( ( dist < 3.81 ) || ( dist > 3.79 ) ) dump_ca_coords( "ca_coords_failed.pdb" );
				runtime_assert( dist < 3.81 );
				runtime_assert( dist > 3.79 );
			}
		}
	}
}

/// @brief convert a point in (c-direction, h-direction, z-direction local coordinates to xyz coordinates
/// c_point1 is assumed to be 0,0,0 in chz coordinates
core::Vector
ParametricSheet::chz_to_xyz( core::Vector const & chz_point,
		core::Vector const & c_point1,
		core::Vector const & c_point2,
		core::Vector const & h_point1,
		core::Vector const & h_point2,
		core::Size const resi ) const
{
	core::Vector c_vec( c_point1 - c_point2 );
	core::Vector h_vec( get_h_axis( h_point1, h_point2, c_vec, resi ) );
	h_vec.normalize();
	c_vec.normalize();
	core::Vector const cross( c_vec.cross( h_vec ).normalize() );

	TR << "CHZ vectors" << std::endl;
	print_vector( c_vec, "c_vec" );
	print_vector( h_vec, "h_vec" );
	print_vector( cross, "cross" );

	// the new vector will have c,h coordinates -- this vector will convert that into xyz
	numeric::xyzMatrix< core::Real > const rotate( numeric::xyzMatrix< core::Real >::rows( c_vec.x(), h_vec.x(), cross.x(),
				c_vec.y(), h_vec.y(), cross.y(),
				c_vec.z(), h_vec.z(), cross.z() ) );
	core::Vector const retval( rotate*chz_point + c_point1 );
	TR << " chz->xyz: [ " << chz_point.x() << ", " << chz_point.y() << ", " << chz_point.z() << " ] --> [ " << retval.x() << ", " << retval.y() << ", " << retval.z() << " ]" << std::endl;
	return retval;
}

numeric::xyzMatrix< core::Real >
ParametricSheet::rotate_about_axis( core::Vector const & twist_axis, core::Real const theta ) const
{
	// rotate about the given vector
	core::Vector axis( twist_axis );
	axis.normalize();
	core::Real const xx = axis.x()*axis.x()*(1-cos(theta)) + cos(theta);
	core::Real const xy = axis.x()*axis.y()*(1-cos(theta)) - axis.z()*sin(theta);
	core::Real const xz = axis.x()*axis.z()*(1-cos(theta)) + axis.y()*sin(theta);
	core::Real const yx = axis.y()*axis.x()*(1-cos(theta)) + axis.z()*sin(theta);
	core::Real const yy = axis.y()*axis.y()*(1-cos(theta)) + cos(theta);
	core::Real const yz = axis.y()*axis.z()*(1-cos(theta)) - axis.x()*sin(theta);
	core::Real const zx = axis.z()*axis.x()*(1-cos(theta)) - axis.y()*sin(theta);
	core::Real const zy = axis.z()*axis.y()*(1-cos(theta)) + axis.x()*sin(theta);
	core::Real const zz = axis.z()*axis.z()*(1-cos(theta)) + cos(theta);
	numeric::xyzMatrix< core::Real > const rotate( numeric::xyzMatrix< core::Real >::rows( xx, xy, xz,
				yx, yy, yz,
				zx, zy, zz ) );
	return rotate;
}

bool
ParametricSheet::check_strand_and_resi_values( core::Size const strand, core::Size const resi ) const
{
	if ( strand < 1 ) return false;
	if ( strand > ca_coords_.size() ) return false;
	if ( ! ca_coords_[strand].size() ) return false;
	if ( resi < 1 ) return false;
	if ( resi > ca_coords_[strand].size() ) return false;
	// actually check the value
	if ( ca_coords_[strand][resi] < -999.0 ) return false;
	return true;
}

/// @brief returns 1 if the C=O for this residue points toward strand i+1, -1 if it points toward strand i-1
int
ParametricSheet::up_down( core::Size const resi ) const
{
	if ( resi % 2 ) {
		return -1;
	} else {
		return 1;
	}
}

/// @brief creates a vector for a new residue in the reference strand based on the previous c and h vectors
/// assumption: points_[i][j-1] has been generated
void
ParametricSheet::create_ca_point( core::Size const strand,
		core::Size const prev_strand,
		core::Size const prev2_strand,
		core::Size const resi,
		core::Size const prev_resi,
		core::Size const prev2_resi )
{
	runtime_assert( resi >= 1 );
	runtime_assert( strand >= 1 );
	runtime_assert( strand <= ca_coords_.size() );
	runtime_assert( resi <= ca_coords_[strand].size() );

	core::Real const r( 3.8021 );

	int c_multiplier( 1 );
	if ( prev_resi > resi ) {
		c_multiplier = -1;
	}

	// check to see what exists
	bool const prev_strand_exists( check_strand_and_resi_values( prev_strand, resi ) );
	bool const prev2_strand_exists( check_strand_and_resi_values( prev2_strand, resi ) );
	core::Vector h_point1;
	core::Vector h_point2;
	if ( prev_strand_exists  && prev2_strand_exists ) {
		h_point1 = ca_coords_[prev_strand][resi];
		h_point2 = ca_coords_[prev2_strand][resi];
	} else if ( prev_strand_exists ) {
		h_point1 = ca_coords_[prev_strand][resi];
		h_point2 = h_point1 - core::Vector( 0.0, strand_dist_, 0.0 );
	} else {
		h_point1 = core::Vector( 0.0, 0.0, 0.0 );
		h_point2 = h_point1 - core::Vector( 0.0, strand_dist_, 0.0 );
	}
	print_vector( h_point1, "H point1" );
	print_vector( h_point2, "H point2" );

	core::Vector c_point1;
	core::Vector c_point2;
	// find c-vector
	bool const prev_resi_exists( check_strand_and_resi_values( strand, prev_resi ) );
	bool const prev2_resi_exists( check_strand_and_resi_values( strand, prev2_resi ) );
	if ( prev_resi_exists && prev2_resi_exists ) {
		c_point1 = ca_coords_[strand][prev_resi];
		c_point2 = ca_coords_[strand][prev2_resi];
	} else if ( prev_resi_exists ) {
		c_point1 = ca_coords_[strand][prev_resi];
		c_point2 = c_point1 - core::Vector( r*c_multiplier, 0.0, 0.0 );
	} else {
		c_point1 = ( h_point1 - h_point2 ) + h_point1 - core::Vector( r*c_multiplier, 0.0, 0.0 );
		c_point2 = c_point1 - core::Vector( r*c_multiplier, 0.0, 0.0 );
	}

	print_vector( c_point1, "C point1" );
	print_vector( c_point2, "C point2" );

	// find reference point = strand-1, resi-1
	core::Vector reference_point;
	if ( check_strand_and_resi_values( prev_strand, prev_resi ) ) {
		reference_point = ca_coords_[prev_strand][prev_resi];
	} else {
		reference_point = c_point1 - ( h_point1 - h_point2 );
	}
	print_vector( reference_point, "Reference point" );

	int up_down_mult( up_down( resi ) );
	up_down_mult *= c_multiplier;
	//core::Real const theta( up_down_mult*2*acos(c_dist_/r) + c_coil_ );
	core::Real const theta( 2*acos(c_dist_/r) + c_coil_*up_down_mult );
	// phi is rotation angle about c direction
	core::Real phi( 0.0 );
	if ( up_down_mult < 0 ) phi += numeric::constants::f::pi;
	TR << "theta=" << theta << ", acos=" << 2*acos(c_dist_/r) << ", c_coil=" << c_coil_ << ", up_down=" << up_down_mult << ", phi=" << phi << std::endl;
	core::Vector const new_point_chz( r*cos(theta), r*sin(theta)*sin(phi), r*sin(theta)*cos(phi) );
	TR << "strand_rotation=" << strand_rotation_ << std::endl;
	print_vector( new_point_chz, "new point_chz" );
	core::Vector const new_point_rotated( chz_to_xyz( new_point_chz, c_point1, c_point2, h_point1, h_point2, resi ) );
	print_vector( new_point_rotated, "rotated point" );

	// correct the h-direction previous strand point in the case where the previous strand isn't made yet
	if ( ! prev_strand_exists ) {
		h_point1 = new_point_rotated - ( h_point1 - h_point2 );
	}
	core::Vector const new_point_twisted( twist( new_point_rotated, c_point1, h_point1, reference_point ) );
	ca_coords_[strand][resi] = new_point_twisted;
	//print_vector_of_vectors( ca_coords_ );
	test_matrix( ca_coords_ );
}

void
ParametricSheet::generate_strand( core::Size const strand, core::Size const prev_strand, core::Size const prev2_strand )
{
	// TODO : fix this shit. Generation of sheets should start in the center and move outward
	// starting on the end of strands is bad because twist doesn't get properly accounted for
	// the strand-1 is to correct 
	core::Size const center_resi( ( ca_coords_[ strand ].size() / 2) + 2 );
	//core::Size const center_resi( strand_len_+4 );
	//core::Size const center_resi(1);
	new_strand( strand, prev_strand, prev2_strand, center_resi );
	for ( core::Size j=center_resi-1; j>=1; --j ) {
		create_ca_point( strand, prev_strand, prev2_strand, j, j+1, j+2 );
	}
	for ( core::Size j=center_resi+1; j<=ca_coords_[ strand ].size(); ++j ) {
		create_ca_point( strand, prev_strand, prev2_strand, j, j-1, j-2 );
	}
}

/// @brief Take the previous point in the c direction, the previous point in the h direction and a reference point and twist about the axis of the two points with the reference point as 0,0,0
core::Vector
ParametricSheet::twist( core::Vector const & xyz_point,
		core::Vector const & c_prev,
		core::Vector const & h_prev,
		core::Vector const & reference ) const
{
	print_vector( xyz_point, "xyz_point" );
	print_vector( c_prev, "c_prev" );
	print_vector( h_prev, "h_prev" );
	print_vector( reference, "reference" );

	TR << "p->c dist=" << xyz_point.distance( c_prev ) << ", h->ref dist=" << h_prev.distance( reference ) << std::endl;
	TR << "p->h dist=" << xyz_point.distance( h_prev ) << ", c->ref dist=" << c_prev.distance( reference ) << std::endl;
	core::Vector c_point( c_prev - c_prev );
	core::Vector h_point( h_prev - c_prev );
	core::Vector ref( reference - c_prev );
	core::Vector new_point( c_point + (h_point-ref) );

	print_vector( c_point, "c_point" );
	print_vector( h_point, "h_point" );
	print_vector( ref, "ref" );
	print_vector( new_point, "new_point" );

	// rotate the point about the c_point -> h_point axis if the new point is in the upper right quadrant
	numeric::xyzMatrix< core::Real > const twist_mat( rotate_about_axis( h_point, twist_ ) );
	print_matrix( twist_mat );

	core::Vector const retval( twist_mat*new_point + c_prev );
	TR << " twist: [ " << xyz_point.x() << ", " << xyz_point.y() << ", " << xyz_point.z() << " ] --> [ " << retval.x() << ", " << retval.y() << ", " << retval.z() << " ]" << std::endl;

	TR << "p->c dist=" << retval.distance( c_prev ) << ", h->ref dist=" << h_prev.distance( reference ) << std::endl;
	TR << "p->h dist=" << retval.distance( h_prev ) << ", c->ref dist=" << c_prev.distance( reference ) << std::endl;

	return retval;
}

/// @brief Ca points are generated in a rectangular grid. This function computes the maximum required height of that grid in the C-direction based on strand lengths and register shifts
core::Size
ParametricSheet::max_strand_len() const
{
	core::Size max_strand_len( 0 );
	for ( core::Size i=1; i<=strand_data_.size(); ++i ) {
		core::Size const cur_value( strand_data_[i].length + strand_data_[i].register_shift );
		if ( cur_value > max_strand_len ) max_strand_len = cur_value;
	}
	return max_strand_len;
}

/// @brief initializes the set of ca coords
void
ParametricSheet::init_ca_coords()
{
	core::Size const strands_to_build( num_strands() + 2 );
	ca_coords_.clear();
	for ( core::Size i=1; i<=strands_to_build; ++i ) {
		ca_coords_.push_back( utility::vector1< core::Vector >( max_strand_len()+4, core::Vector( -9999.9, -9999.9, -9999.9 ) ) );
	}
}

/// @brief generates parametric coordinates that basically define the axes of the sheet
/// sets the vector of vectors of points that make up each strand
void
ParametricSheet::generate_residue_positions()
{
	core::Size const center_low( ( num_strands() / 2 ) + 1 );
	core::Size center_high( ( num_strands() / 2 ) + 2 );

	// initialize points, ca_coords data
	init_ca_coords();

	// strands should be build from the center outward
	for ( core::Size i=center_low; i>=1; --i ) {
		generate_strand( i, i+1, i+2 );
	}
	for ( core::Size i=center_high; i<=ca_coords_.size(); ++i ) {
		generate_strand( i, i-1, i-2 );
	}
}

/// @brief generates a pdb file consisting of Ca positions
void
ParametricSheet::dump_ca_coords( std::string const & pdb_name ) const
{
	std::ofstream outfile( pdb_name.c_str() );
	// set proper floating point precisions
	outfile << std::fixed;
	core::Size atomid( 1 );
	for ( core::Size i=1; i<=ca_coords_.size(); ++i ) {
		for ( core::Size j=1; j<=ca_coords_[i].size(); ++j ) {
			if ( ( ca_coords_[i][j].x() > -999.0 ) &&
					( ca_coords_[i][j].y() > -999.0 ) &&
					( ca_coords_[i][j].z() > -999.0 ) ) {
				// 1-6 record name
				outfile << "HETATM";
				// 7-11 atom number
				outfile << std::setfill(' ') << std::setw(5) << atomid;
				// 13-16 atom name
				outfile << " C" << std::setfill('0') << std::setw(2) << j;
				// 17 alternate location
				outfile << " ";
				// 18-20 residue name, 22 chain
				outfile << "STR A";
				// 23-26 reside number
				outfile << std::setfill(' ') << std::setw( 4 ) << i;
				// 27-30 nothing
				outfile << "    ";
				// 31-38 x, 39-46 y, 47-54 z
				outfile << std::setfill(' ') << std::setw( 8 ) << std::setprecision( 3 ) << ca_coords_[i][j].x();
				outfile << std::setfill(' ') << std::setw( 8 ) << std::setprecision( 3 ) << ca_coords_[i][j].y();
				outfile << std::setfill(' ') << std::setw( 8 ) << std::setprecision( 3 ) << ca_coords_[i][j].z() << std::endl;;
				++atomid;
			}
		}
	}
	outfile.close();
}

/// @brief compute the angle between two vectors using atan2 -- both are assumed to come from 0,0,0
core::Real
vector_angle( core::Vector const & v1, core::Vector const & v2 )
{
	core::Vector const cross( v1.cross( v2 ) );
	core::Real const dot( v1.dot( v2 ) );
	return atan2( cross.length(), dot );
}

/// @brief helper function to calculate the mean
core::Real mean( utility::vector1< core::Real > const & vec )
{
	core::Real avg( 0.0 );
	for ( core::Size i=1; i<=vec.size(); ++i ) {
		avg += vec[i];
	}
	return avg / vec.size();
}

/// @brief helper function to calculate std dev given mean and vector
core::Real std_dev( utility::vector1< core::Real > const & vec, core::Real const mean )
{
	core::Real sum2( 0.0 );
	for ( core::Size i=1; i<=vec.size(); ++i ) {
		sum2 += ( vec[i] - mean )*( vec[i] - mean );
	}
	sum2 /= vec.size();
	return sqrt( sum2 );
}

/// @brief calculates the mean and std. dev of ca distance along the c axis
std::pair< core::Real, core::Real >
ParametricSheet::calc_cdist() const
{
	utility::vector1< core::Real > dists;
	// limit ourselves to the actual points included in the sheet
	for ( core::Size i=2; i<=ca_coords_.size()-1; ++i ) {
		for ( core::Size j=3; j<=ca_coords_[i].size()-2; ++j ) {
				dists.push_back( calc_cdist( i, j, j-1, j+1 ) );
		}
	}
	core::Real avg( mean( dists ) );
	core::Real dev( std_dev( dists, avg ) );
	return std::pair< core::Real, core::Real >( avg, dev );
}

/// @brief calculates the ca distance along the c axis between the given strand and residue numbers
core::Real
ParametricSheet::calc_cdist( core::Size const strand, core::Size const resi,
		core::Size const resi_prev, core::Size const resi_next ) const
{
	// get c-axis estimate from resi_prev, resi_next
	core::Vector c_axis( ca_coords_[strand][resi_next] - ca_coords_[strand][resi_prev] );
	c_axis.normalize();
	core::Vector const v_next( ca_coords_[strand][resi_next] - ca_coords_[strand][resi] );
	core::Vector const v_prev( ca_coords_[strand][resi_prev] - ca_coords_[strand][resi] );
	core::Real const dist_next( v_next.dot( c_axis ) );
	core::Real const dist_prev( v_prev.dot( c_axis ) );
	TR << "distances = " << dist_next << " and " << dist_prev << std::endl;
	return ( ( dist_next - dist_prev ) / 2.0 );
}

/// @brief calculates mean and std dev of twist
std::pair< core::Real, core::Real >
ParametricSheet::calc_twist() const
{
	utility::vector1< core::Real > twists;
	for ( core::Size i=2; i<=ca_coords_.size()-1; ++i ) {
		for ( core::Size j=3; j<=ca_coords_[i].size()-2; ++j ) {
			core::Real const t( calc_twist( ca_coords_[i][j], ca_coords_[i][j-1], ca_coords_[i-1][j], ca_coords_[i-1][j-1] ) );
			twists.push_back( t );
			TR << "Twist=" << t << std::endl;
		}
	}
	core::Real const avg( mean( twists ) );
	return std::pair< core::Real, core::Real >( avg, std_dev( twists, avg ) );
}

/// @brief given four points, calculate the twist
core::Real
ParametricSheet::calc_twist( core::Vector const & point, core::Vector const & c_prev, core::Vector const & h_prev, core::Vector const & ref ) const
{
	core::Vector ref_vec( ( h_prev - ref ) + ( c_prev - ref ) );
	core::Vector test_vec( ( point - h_prev ) + ( point - c_prev ) );
	return vector_angle( test_vec, ref_vec );
}

/// @brief calculates sheet geometry
ParametricSheet::SheetGeometry
ParametricSheet::calc_geometry() const
{
	SheetGeometry geom;

	// h distances
	utility::vector1< core::Real > dists;
	for ( core::Size i=2; i<=ca_coords_.size()-1; ++i ) {
		for ( core::Size j=3; j<=ca_coords_[i].size()-2; ++j ) {
			dists.push_back( ca_coords_[i][j].distance( ca_coords_[i-1][j] ) );
		}
	}

	geom.hdist_mean = mean( dists );
	geom.hdist_dev = std_dev( dists, geom.hdist_mean );

	// twist
	std::pair< core::Real, core::Real > const twist_geom( calc_twist() );
	geom.twist_mean = twist_geom.first;
	geom.twist_dev = twist_geom.second;

	// c-distance
	std::pair< core::Real, core::Real > const cdist_geom( calc_cdist() );
	geom.cdist_mean = cdist_geom.first;
	geom.cdist_dev = cdist_geom.second;
	return geom;
}

/// @brief prints out sheet geometry
void
ParametricSheet::SheetGeometry::show( std::ostream & out ) const
{
	out << "CDist: " << cdist_mean << " +/- " << cdist_dev << std::endl;
	out << "HDist: " << hdist_mean << " +/- " << hdist_dev << std::endl;
	out << "Twist: " << twist_mean << " +/- " << twist_dev << std::endl;
}

} // namespace denovo_design

} // namespace devel
