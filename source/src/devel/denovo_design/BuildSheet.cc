// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/BuildSheet.cc
/// @brief Parametric generation of beta sheets
/// @detailed
/// @author Tom Linsky


//Unit Headers
#include <devel/denovo_design/BuildSheet.hh>
#include <devel/denovo_design/BuildSheetCreator.hh>

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

static basic::Tracer TR("devel.denovo_design.BuildSheet");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace devel {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////


using namespace ObjexxFCL;

std::string
BuildSheetCreator::keyname() const
{
	return BuildSheetCreator::mover_name();
}

protocols::moves::MoverOP
BuildSheetCreator::create_mover() const {
	return new BuildSheet();
}

std::string
BuildSheetCreator::mover_name()
{
	return "BuildSheet";
}

///  ---------------------------------------------------------------------------------
///  BuildSheet main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
BuildSheet::BuildSheet() :
	protocols::moves::Mover( "BuildSheet" ),
	minimize_( false ),
	phi1_( -125.0 ),
	phi2_( -125.0 ),
	psi1_( 125.0 ),
	psi2_( 125.0 ),
	num_cached_dihedrals_( 30 ),
	scorefxn_( NULL )
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
BuildSheet::~BuildSheet() {}


/// Return a copy of ourselves
protocols::moves::MoverOP
BuildSheet::clone() const
{
	return new BuildSheet(*this);
}

std::string
BuildSheet::get_name() const {
	return BuildSheetCreator::mover_name();
}

void
BuildSheet::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & )
{
	if ( tag->hasOption( "c_dist" ) ) sheet_.c_dist( tag->getOption< core::Real >( "c_dist" ) );
	if ( tag->hasOption( "strand_dist" ) ) sheet_.strand_dist( tag->getOption< core::Real >( "strand_dist" ) );
	if ( tag->hasOption( "twist" ) ) sheet_.twist( tag->getOption< core::Real >( "twist" ) );
	if ( tag->hasOption( "c_coil" ) ) sheet_.c_coil( tag->getOption< core::Real >( "c_coil" ) );
	if ( tag->hasOption( "h_coil" ) ) sheet_.h_coil( tag->getOption< core::Real >( "h_coil" ) );
	if ( tag->hasOption( "strand_rotation" ) ) sheet_.strand_rotation( tag->getOption< core::Real >( "strand_rotation" ) );
	minimize_ = tag->getOption< bool >( "minimize", minimize_ );

	phi1_ = tag->getOption< core::Real >( "phi1", phi1_ );
	phi2_ = tag->getOption< core::Real >( "phi2", phi2_ );
	psi1_ = tag->getOption< core::Real >( "psi1", psi1_ );
	psi2_ = tag->getOption< core::Real >( "psi2", psi2_ );
	if ( tag->hasOption( "scorefxn" ) ) {
		scorefxn_ = data.get< core::scoring::ScoreFunction * >( "scorefxns", tag->getOption< std::string >( "scorefxn" ) );
	}
	sheet_.parse_tags( tag->getTags() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Basic Accessors


////////////////////////////////////////////////////////////////////////////////////////////////////
void BuildSheet::apply( core::pose::Pose & pose )
{
	core::pose::Pose newpose;
	newpose.pdb_info( new core::pose::PDBInfo( pose, true ) );
	sheet_.generate_residue_positions();
	sheet_.dump_ca_coords( "ca_coords.pdb" );
	ParametricSheet::SheetGeometry const geom( sheet_.calc_geometry() );
	geom.show( TR ); TR.flush();
	//sheet_.from_pose( pose );
	//sheet_.dump_ca_coords( "pose_ca_coords.pdb" );
	//ParametricSheet::SheetGeometry const geom_pose( sheet_.calc_geometry() );
	//geom_pose.show( TR );
	//utility_exit();
	read_dihedral_map();
	for ( core::Size i=1; i<=dihedral_map_.size(); ++i ) {
		for ( core::Size j=1; j<=dihedral_map_[i].size(); ++j ) {
			for ( core::Size k=1; k<=dihedral_map_[i][j].size(); ++k ) {
				if ( dihedral_map_[i][j][k].set ) {
					TR << "YART " << i << " " << j << " " << k << " " << dihedral_map_[i][j][k].phi << "  " << dihedral_map_[i][j][k].psi << " " << dihedral_map_[i][j][k].h_angle << std::endl;
				}
			}
		}
	}
	build_sheet( newpose );
	pose = newpose;
	pose.dump_pdb( "ideal_sheet.pdb" );
}

void
print_vectorBS( numeric::xyzVector< core::Real > const & vec, std::string const & name = "" )
{
	TR << name << " [ " << vec.x() << ", " << vec.y() << ", " << vec.z() << " ] " << std::endl;
}

/// @brief given a vector, return spherical coordinate for theta for that vector
core::Real calc_spherical_theta( core::Vector const & vec )
{
	return acos( vec.z() / vec.length() );
}

/// @brief given a vector, return spherical coordinate for theta for that vector
core::Real calc_spherical_phi( core::Vector const & vec )
{
	return atan2( vec.y(), vec.x() );
}

/// @brief given a set of four Ca-points on strand i ( j-2, j-1, j, j+1 ) estimate the appropriate dihedrals
BuildSheet::AngleData const &
BuildSheet::lookup_dihedrals( core::Size const i, core::Size const j ) const
{
	runtime_assert( j-2 > 0 );
	runtime_assert( i > 0 );
	runtime_assert( i <= sheet_.ca_coords_size() );
	runtime_assert( j+1 <= sheet_.ca_coords_size( i ) );
	numeric::xyzMatrix< core::Real > const rotate( align_matrix( sheet_.ca_coords( i, j-2 ),
																															 sheet_.ca_coords( i, j-1 ),
																															 sheet_.ca_coords( i, j ) ) );
	core::Vector const next_ca_vec_xyz( sheet_.ca_coords( i, j+1 ) - sheet_.ca_coords( i, j ) );
	core::Vector const next_ca_vec( rotate*next_ca_vec_xyz );
	TR << "looking up dihedral; i=" << i << ", j=" << j << std::endl;
	print_vectorBS( next_ca_vec_xyz, "next ca vec pre-rotation" );
	print_vectorBS( next_ca_vec, "next ca vec" );
	print_vectorBS( rotate * ( sheet_.ca_coords( i, j-1 ) - sheet_.ca_coords( i, j ) ), "i-1" );
	core::Real const theta( calc_spherical_theta( next_ca_vec ) );
	core::Real const phi( calc_spherical_phi( next_ca_vec ) );
	//core::Real const h_angle( sheet_.calc_h_angle( i, j ) );
	core::Real const h_angle( 0.0 );
	TR << "theta=" << theta << ", phi=" << phi << ", h_angle=" << h_angle << "; i1=" << index_from_theta( theta ) << ", i2=" << index_from_phi( phi ) << ", i3=" << index_from_hangle( h_angle ) << std::endl;
	return dihedral_map_[index_from_theta(theta)][index_from_phi(phi)][index_from_hangle(h_angle)];
}

/// @brief physically builds a strand based on the provided ca coordinates
void
BuildSheet::build_sheet( core::pose::Pose & pose ) const
{
	core::chemical::ResidueTypeSetCAP restype_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );
	std::string const res_name( "VAL" );
	core::pose::Pose dummy_pose( build_ideal_strand( restype_set, res_name, 3 ) );
	runtime_assert( dummy_pose.total_residue() >= 3 );

	// translate dummy_pose so that CA of residue 2 is at the origin
	dummy_pose.apply_transform_Rx_plus_v( numeric::xyzMatrix< core::Real >::identity(), -dummy_pose.residue(2).xyz( "CA" ) );
	dummy_pose.dump_pdb( "ideal_strand.pdb" );

	// keep track of which residues go with which strands
	std::map< std::pair< core::Size, core::Size >, core::Size > strand_res_to_pose;
	for ( core::Size i=2; i<=sheet_.ca_coords_size()-1; ++i ) {
		bool parallel( true );
		core::Size const min_value( 3 + sheet_.register_shift( i-1 ) );
		core::Size const max_value( sheet_.strand_len( i-1 ) + 2 + sheet_.register_shift( i-1 ) );
		utility::vector1< core::Size > j_values;
		if ( sheet_.orientation( i-1 ) == "P" || sheet_.orientation( i-1 ) == "p" ) {
			for ( core::Size j=min_value; j<=max_value; ++j ) {
				j_values.push_back( j );
			}
		} else {
			for ( core::Size j=max_value; j>=min_value; --j ) {
				j_values.push_back( j );
			}
			parallel = false;
		}
		for ( core::Size jindex=1; jindex<=j_values.size(); ++jindex ) {
			core::Size const j( j_values[jindex] );
			// set dummy pose phi-psi angles based on Ca positions
			AngleData const & dihedrals( lookup_dihedrals( i, j ) );
			// dihedral data should be set
			if ( !dihedrals.set ) {
				pose.dump_pdb( "ideal_sheet_incomplete.pdb" );
				utility_exit_with_message( "Could not build a sheet with the given parameters" );
			}
			dummy_pose.set_phi( 2, dihedrals.phi );
			dummy_pose.set_psi( 2, dihedrals.psi );
			TR << "set phi to " << dihedrals.phi << " and psi to " << dihedrals.psi << std::endl;

			// axes will be c-1, c+1, cross
			core::Vector const c2( dummy_pose.residue(2).xyz( "CA" ) );
			core::Vector c1( dummy_pose.residue(1).xyz( "CA" ) - c2 );
			core::Vector c3( dummy_pose.residue(3).xyz( "CA" ) - c2 );
			core::Vector cross( c1.cross( c3 ) );
			c1.normalize();
			c3.normalize();
			cross.normalize();

			// this transformation matrix converts vectors of proper length to the three bb atoms
			numeric::xyzMatrix< core::Real > const T1( numeric::xyzMatrix< core::Real >::rows( c1.x(), c3.x(), cross.x(),
						c1.y(), c3.y(), cross.y(),
						c1.z(), c3.z(), cross.z() ) );

			numeric::xyzMatrix< core::Real > const T1_inv( numeric::inverse( T1 ) );

			// Center the coordinates about the Ca residue of the residue to be superimposed.
			core::Vector const ideal_this( sheet_.ca_coords( i, j ) );
			core::Vector ideal_prev2, ideal_prev, ideal_next;
			if ( parallel ) {
				ideal_prev2 = sheet_.ca_coords( i, j-2 ) - ideal_this;
				ideal_prev = sheet_.ca_coords( i, j-1 ) - ideal_this;
				ideal_next = sheet_.ca_coords( i, j+1 ) - ideal_this;
			} else {
				ideal_prev2 = sheet_.ca_coords( i, j+2 ) - ideal_this;
				ideal_prev = sheet_.ca_coords( i, j+1 ) - ideal_this;
				ideal_next = sheet_.ca_coords( i, j-1 ) - ideal_this;
			}
			ideal_prev.normalize();
			ideal_next.normalize();
			core::Vector ideal_prevstrand( ideal_prev.cross( ideal_next ) );
			ideal_prevstrand.normalize();
			// this transformation matrix convers unit vectors to the three "ideal" ca atoms
			numeric::xyzMatrix< core::Real > const T_ideal( numeric::xyzMatrix< core::Real >::rows( ideal_prev.x(), ideal_next.x(), ideal_prevstrand.x(),
						ideal_prev.y(), ideal_next.y(), ideal_prevstrand.y(),
						ideal_prev.z(), ideal_next.z(), ideal_prevstrand.z() ) );
			core::conformation::ResidueOP new_res( new core::conformation::Residue( dummy_pose.residue(2) ) );
			for ( core::Size k=1; k<=new_res->natoms(); ++k ) {
				numeric::xyzVector< core::Real > const & cur_pos( new_res->atom(k).xyz() );
				TR << " strand " << i << " res " << j << " atom " << k << " old,new xyz=" << cur_pos.x() << " " << cur_pos.y() << " " << cur_pos.z() << " ";
				numeric::xyzVector< core::Real > const transformed_dummy( T1_inv * cur_pos );
				print_vectorBS( transformed_dummy, "Transformed_dummy" );
				core::Vector const transformed_two( ( T_ideal * transformed_dummy ) + ideal_this );
				new_res->atom(k).xyz( transformed_two );
			}
			if ( jindex == 1 ) {
				pose.append_residue_by_jump( *new_res, pose.total_residue(), "", "", true );
			} else {
				pose.conformation().safely_append_polymer_residue_after_seqpos( *new_res, pose.total_residue(), false );
			}
			strand_res_to_pose[ std::pair< core::Size, core::Size >( i, j ) ] = pose.total_residue();
			if ( minimize_ ) {
				utility::vector1< core::Size > strands;
				utility::vector1< core::Size > resis;
				strands.push_back( i );
				strands.push_back( i );
				strands.push_back( i-1 );
				resis.push_back( j );
				utility::vector1< core::Size > pose_residues;
				pose_residues.push_back( pose.total_residue() );
				// if this isn't the first residue
				if ( jindex > 1 ) {
					pose_residues.push_back( strand_res_to_pose[ std::pair< core::Size, core::Size >( i, j_values[jindex-1] ) ] );
					resis.push_back( j_values[jindex-1] );
				} else {
					pose_residues.push_back( 0 );
					resis.push_back( 0 );
				}
				// if this isn't the first strand
				if ( i > 2 ) {
					pose_residues.push_back( strand_res_to_pose[ std::pair< core::Size, core::Size >( i-1, j ) ] );
				} else {
				pose_residues.push_back( 0 );
			}
			resis.push_back( j );
			// minimize to try to get hydrogen bonds perfect
			if ( i >= 2 ) {
				minimize_residues( pose, pose_residues, strands, resis );
				pose.dump_pdb("min" + boost::lexical_cast< std::string >( i ) + "-" + boost::lexical_cast< std::string >( j ) + ".pdb");
			}
			//minimize_residues( pose, residues, i, j );
			//minimize_residues( pose, residues, i, j );
			}
		}
	}
}

/// @brief minimizes residue resid with cartesian minimization, including strong Ca constraints
/// assumes that the Ca positions are correct
/// the passed vector contains residue ids for the current residue, previous residue in the C-direction, and previous residue in the H-direction, IN THAT ORDER
/// the "strands" vector contains ca_coords_ indices for the same three residues, in the same order
/// the "resis" vector contains ca_coords_ residue indices for the same three residues, in the same order
void
BuildSheet::minimize_residues( core::pose::Pose & pose,
		utility::vector1< core::Size > const & pose_residues,
		utility::vector1< core::Size > const & strands,
		utility::vector1< core::Size > const & resis ) const
{
	TR << "Residues are: " << pose_residues[1] << " " << pose_residues[2] << " " << pose_residues[3] << std::endl;
	runtime_assert( pose_residues.size() == 3 );
	runtime_assert( strands.size() == 3 );
	runtime_assert( resis.size() == 3 );
	runtime_assert( pose_residues[1] > 0 );
	runtime_assert( pose_residues[1] <= pose.total_residue() );

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
	mm->set_bb( false );
	mm->set_chi( false );
	mm->set_jump( false );
	mm->set( core::id::D, true );
	mm->set( core::id::THETA, true );
	core::Size const anchor_ca_index( pose.residue_type(1).atom_index("CA") );

	// functions for coordinate constraints
	core::scoring::func::HarmonicFuncOP weak_func( new core::scoring::func::HarmonicFunc(0.0, 2.0 ) );
	core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc(0.0, 0.1 ) );

	// first element is the residue being minimized
	// second element is the previous residue in the c direction
	// third element is the previous residue in the h direction

	// enable moveable residues and add constraints
	for ( core::Size i=1; i<=pose_residues.size(); ++i ) {
		if ( pose_residues[i] ) {
			runtime_assert( pose_residues[i] <= pose.total_residue() );
			runtime_assert( strands[i] > 0 );
			runtime_assert( strands[i] <= sheet_.ca_coords_size() );
			runtime_assert( resis[i] > 0 );
			runtime_assert( resis[i] <= sheet_.ca_coords_size( strands[i] ) );
			mm->set_bb( pose_residues[i], true );

			// add all-atom constraints
			for ( core::Size j=1; j<=pose.residue(pose_residues[i]).natoms(); ++j ) {
				core::scoring::constraints::ConstraintOP cst;
				cst = new core::scoring::constraints::CoordinateConstraint( core::id::AtomID(j,pose_residues[i]),
						core::id::AtomID(anchor_ca_index,1),
						pose.residue(pose_residues[i]).atom(j).xyz(),
						weak_func );
				pose.add_constraint( cst );
			}

			// add ca constraints
			if ( pose.residue_type( pose_residues[i] ).has("CA") ) {
				core::Size const ca_index( pose.residue_type( pose_residues[i] ).atom_index("CA") );
				core::scoring::constraints::ConstraintOP cst;
				cst = new core::scoring::constraints::CoordinateConstraint( core::id::AtomID( ca_index, pose_residues[i] ),
																																		core::id::AtomID( anchor_ca_index, 1 ),
																																		sheet_.ca_coords( strands[i], resis[i] ),
																																		func );
				pose.add_constraint( cst );
			}
		}
	}
	mm->show( TR ); TR.flush();

	// get a default centroid scorefunction
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::getScoreFunction( "fa_standard" ) );
	scorefxn->set_weight( core::scoring::cart_bonded, 0.5 );
	scorefxn->set_weight( core::scoring::coordinate_constraint, 1.0 );
	protocols::simple_moves::MinMover minimize;
	minimize.cartesian( true );
	minimize.movemap( mm );
	minimize.score_function( *scorefxn_ );
	minimize.min_type( "lbfgs_armijo_nonmonotone" );
	minimize.tolerance( 0.00001 );
	scorefxn->show( TR, pose ); TR.flush();
	minimize.apply( pose );
	scorefxn->show( TR, pose ); TR.flush();

	pose.remove_constraints();
}

/// @brief given a dihedral theta angle, return a proper index in the AngleData vector
core::Size
BuildSheet::index_from_theta( core::Real const theta ) const
{
	// can be from 0 <= theta <= pi
	return index_from_angle( theta, numeric::constants::f::pi );
}

/// @brief given a dihedral phi angle, return a proper index in the AngleData vector
core::Size
BuildSheet::index_from_phi( core::Real const phi ) const
{
	// can be from 0 <= phi <= 2*pi
	return index_from_angle( phi, numeric::constants::f::pi*2 );
}

/// @brief given a dihedral phi angle, return a proper index in the AngleData vector
core::Size
BuildSheet::index_from_hangle( core::Real const h_angle ) const
{
	return index_from_angle( h_angle, numeric::constants::f::pi );
}

/// @brief given a dihedral angle and periodicity, return an integer that refers to the proper index in the AngleData vector
core::Size
BuildSheet::index_from_angle( core::Real angle, core::Real const max_value ) const
{
	// make angle periodic from 0-max_value
	while ( angle < 0 ) {
		angle += max_value;
	}
	while ( angle > max_value ) {
		angle -= max_value;
	}
	runtime_assert( angle <= max_value );
	runtime_assert( angle >= 0.0 );
	core::Size const val( angle/max_value*num_cached_dihedrals_ );
	runtime_assert( val+1 > 0 );
	runtime_assert( val+1 <= num_cached_dihedrals_ );
	return val+1;
}

/// @brief calculates an angle which represents the middle of the given hash index
core::Real
BuildSheet::angle_from_index( core::Size const index, core::Real const max_value ) const
{
	core::Real const val( (core::Real)index * max_value / (core::Real)num_cached_dihedrals_ );
	return val;
}

/// @brief calculates the angle between the two vectors
/// @brief creates a rotation matrix that places p1, p2, and p_center in the XY plane with p_center at the origin and p2 on the x axis
numeric::xyzMatrix< core::Real >
align_matrix( core::Vector const & p1,
		core::Vector const & p2,
		core::Vector const & p_center )
{
	// define axes for coordinate system
	// 1. p2 will lie at -r, 0, 0
	// 2. h axis is the cross of vectors from p_center to p2 and p4
	// 3. all points will lie on xy plane -- y = the cross of h axis and p2
	// 4. point p1 will lie in -y space, not +y
	core::Vector p1_vec( p1 - p_center );
	core::Vector p2_vec( p2 - p_center );
	print_vectorBS( p2_vec, "p2_vec" );
	core::Real const theta( atan2( p2_vec.y(), p2_vec.x() ) );
	// rotate about z axis
	TR << "theta = " << theta << std::endl;
	numeric::xyzMatrix< core::Real > rotate_z( numeric::xyzMatrix< core::Real >::rows( cos(theta), sin(theta), 0.0,
				-sin(theta), cos(theta), 0.0,
				0.0, 0.0, 1.0 ) );
	core::Vector const trans_point( rotate_z * p2_vec );
	print_vectorBS( trans_point, "trans_point" );
	core::Real const angle2( numeric::constants::f::pi + atan2( trans_point.z(), trans_point.x() ) );
	TR << "phi = " << angle2 << std::endl;
	numeric::xyzMatrix< core::Real > rotate_y( numeric::xyzMatrix< core::Real >::rows( cos(angle2), 0.0, sin(angle2),
				0.0, 1.0, 0.0,
				-sin(angle2), 0.0, cos(angle2) ) );

	core::Vector new_p1_vec( rotate_y * rotate_z * p1_vec );
	print_vectorBS( new_p1_vec, "new_p1_vec" );
	// after these rotations, p2 will lie along the x axis
	// now we rotate about the x axis to place p1 on the xy plane
	core::Real const angle3( atan2( new_p1_vec.z(), new_p1_vec.y() ) );
	numeric::xyzMatrix< core::Real > rotate_x( numeric::xyzMatrix< core::Real >::rows( 1.0, 0.0, 0.0,
				0.0, cos(angle3), sin(angle3),
				0.0, -sin(angle3), cos(angle3) ) );
	print_vectorBS( rotate_x * new_p1_vec, "fully_rotated_p1" );
	print_vectorBS( rotate_x * rotate_y * rotate_z * p2_vec, "fully_rotated_p2" );
	return rotate_x*rotate_y*rotate_z;
}

/// @brief aligns residues 1-3 of the pose to the X-Y plane, with CA3 as 0,0,0
void
BuildSheet::align_pose( core::pose::Pose & pose ) const
{
	runtime_assert( pose.total_residue() >= 4 );
	// translate Ca-3 to the origin
	core::Vector const c3_pos( pose.residue(3).xyz( "CA" ) );
	pose.apply_transform_Rx_plus_v( numeric::xyzMatrix< core::Real >::identity(), -c3_pos );

	// rotate such that Ca-1 and Ca-2 lie on the X axis
	// Ca-1 should be at (-r, 0, 0) in XYZ and Ca-2 should be at (0, 0, 0) in XYZ
	// Ca-0 should lie in the XZ plane
	numeric::xyzMatrix< core::Real > const rotate( align_matrix( pose.residue(1).xyz( "CA" ),
				pose.residue(2).xyz( "CA" ),
				pose.residue(3).xyz( "CA" ) ) );
	pose.apply_transform_Rx_plus_v( rotate, core::Vector( 0.0, 0.0, 0.0 ) );

	print_vectorBS( pose.residue(1).xyz( "CA" ), "C-2 AFter rotate" );
	print_vectorBS( pose.residue(2).xyz( "CA" ), "C-1 After rotate" );
	print_vectorBS( pose.residue(3).xyz( "CA" ), "C-0 After rotate" );
	print_vectorBS( pose.residue(4).xyz( "CA" ), "C1  After rotate" );
}

/// @brief initializes and rests the dihedral map
void
BuildSheet::init_dihedral_map()
{
	dihedral_map_.clear();
	for ( core::Size i=1; i<=num_cached_dihedrals_; ++i ) {
		dihedral_map_.push_back( utility::vector1< utility::vector1< AngleData > >() );
		for ( core::Size j=1; j<=num_cached_dihedrals_; ++j ) {
			utility::vector1< AngleData > dihedrals( num_cached_dihedrals_ );
			dihedral_map_[i].push_back( dihedrals );
		}
	}
}

/// @brief reads a map of Ca coordinates vs. phi/psi angles from a db file
/// if indexed is true, the first two columns in the file refer to pre-hashed indices of the dihedral_map_
/// otherwise, they are directly spherical coordinates
void
BuildSheet::read_dihedral_map()
{
	init_dihedral_map();
	std::ifstream infile( "../dihedral_db" );
	// if file can't be found, create a dihedral map
	if ( !infile.good() ) {
		build_dihedral_map();
		return;
	}
	core::Size num_bins;
	infile >> num_bins;
	if ( num_bins ) {
		num_cached_dihedrals_ = num_bins;
	}

	while ( !infile.eof() ) {
		core::Size i, j, h;
		core::Real phi, psi, angle_diff;
		infile >> i >> j >> h >>  phi >> psi >> angle_diff;
		runtime_assert( phi >= -180.0 );
		runtime_assert( phi <= 180.0 );
		runtime_assert( psi >= -180.0 );
		runtime_assert( psi <= 180.0 );
		if ( infile.good() ) {
			runtime_assert( num_bins );
			runtime_assert( i <= num_cached_dihedrals_ );
			runtime_assert( i > 0.9999 );
			runtime_assert( j <= num_cached_dihedrals_ );
			runtime_assert( j > 0.9999 );
			runtime_assert( h <= num_cached_dihedrals_ );
			runtime_assert( h > 0.9999 );
			AngleData adata( phi, psi, angle_diff );
			runtime_assert( ! dihedral_map_[i][j][h].set );
			dihedral_map_[i][j][h] = adata;
			TR << "Read [ " << i << " " << j << " " << h << " ] -> " << phi << " " << psi << " " << angle_diff << std::endl;
		}
	}
}

/// @brief builds a map of Ca coordinates vs. phi/psi angles
void
BuildSheet::build_dihedral_map()
{
	init_dihedral_map();
	core::chemical::ResidueTypeSetCAP restype_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" ) );
	core::pose::Pose dummy( build_ideal_strand( restype_set, "VAL", 4 ) );
	// put pose in the proper orientation
	align_pose( dummy );
	print_vectorBS( dummy.residue(1).xyz( "CA" ), "C-2 CA POS" );
	core::Real const r_ref( dummy.residue(2).xyz( "CA" ).distance( dummy.residue(3).xyz( "CA" ) ) );
	core::Size const nstep( 1440 ); // 0.25 degree increments
	for ( core::Size psi2_step = 1; psi2_step <= nstep; ++psi2_step ) {
		for ( core::Size phi2_step = 1; phi2_step <= nstep; ++phi2_step ) {
			core::Vector const c1( dummy.residue(1).xyz("CA") );
			core::Vector const c2( dummy.residue(2).xyz("CA") );
			core::Vector const c3( dummy.residue(3).xyz("CA") );
			core::Real const phi_value( (phi2_step-1)*360.0/nstep - 180.0 );
			core::Real const psi_value( (psi2_step-1)*360.0/nstep - 180.0 );
			dummy.set_psi( 3, psi_value );
			dummy.set_phi( 3, phi_value );
			TR << "[ phi, psi ] = " << phi_value << ", " << psi_value << std::endl;
			//dummy.dump_pdb( "test" + boost::lexical_cast< std::string >( psi2_step*phi2_step ) + ".pdb" );
			core::Vector const ca( dummy.residue(4).xyz( "CA" ) );
			print_vectorBS( dummy.residue(1).xyz( "CA" ), "C-2 CA POS" );
			print_vectorBS( ca, "[ x, y, z ] = " );
			core::Real const r( ca.length() );
			core::Real const a( calc_spherical_theta( ca ) );
			core::Real const b( calc_spherical_phi( ca ) );
			TR << "[ r, a, b ] = [ " << r << ", " << a << ", " << b << " ]" << std::endl;
			if ( ( r > r_ref+0.0001 ) || ( r < r_ref-0.0001 ) ) {
				utility_exit_with_message( "distance between Ca atoms is invalid: " + boost::lexical_cast< std::string >( r ) + " != " + boost::lexical_cast< std::string >( r_ref ) );
			}
			core::Size const i( index_from_theta(a) );
			core::Size const j( index_from_phi(b) );
			// try taking the result with the oxygen atom parallel to the h-direction
			core::Vector const cprev_axis( dummy.residue(2).xyz("CA") - dummy.residue(3).xyz("CA") );
			core::Vector const cnext_axis( dummy.residue(4).xyz("CA") - dummy.residue(3).xyz("CA") );
			core::Vector h_axis( cprev_axis.cross( cnext_axis ) );
			h_axis.normalize();
			// h axis should be co-directional with the o-c vector
			core::Vector const co( dummy.residue(3).xyz("O") - dummy.residue(3).xyz("C") );
			// get angle of CO with h axis
			core::Real const h_angle( vector_angle( h_axis, co ) );
			print_vectorBS( co, "C-O Vector" );
			print_vectorBS( h_axis, "H axis" );
			for ( core::Size h=1; h<=dihedral_map_[i][j].size(); ++h ) {
				core::Real const angle_diff( std::abs( h_angle - angle_from_index( h, numeric::constants::f::pi ) ) );
				AngleData adata( phi_value, psi_value, angle_diff );
				if ( dihedral_map_[i][j][h].set ) {
					TR << "(phi_value, psi_value)=" << phi_value << ", " << psi_value << " value for " << i << ", " << j << ", " << h << " is already set to " << dihedral_map_[i][j][h].phi << ", " << dihedral_map_[i][j][h].psi << "..." << "h_angle " << angle_diff << " saved angle is " << dihedral_map_[i][j][h].h_angle << " ";
					if ( angle_diff < dihedral_map_[i][j][h].h_angle ) {
						TR << "OVERWRITING";
						dihedral_map_[i][j][h] = adata;
					}
					TR << std::endl;
				} else {
					dihedral_map_[i][j][h] = adata;
				}
			}
			TR << "XAWK " << a << " " << b << " " << phi_value << " " << psi_value << " " << h_angle << std::endl;
		}
	}
}

/// @brief builds a small idealized strand fragment to be used for parameterized stuff
core::pose::Pose
BuildSheet::build_ideal_strand( core::chemical::ResidueTypeSetCAP restype_set,
		std::string const & res_name,
		core::Size const len ) const {
	core::pose::Pose pose;
	for ( core::Size i=1; i<=len; ++i ) {
		core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue( restype_set->name_map( res_name ) );
		if ( pose.total_residue() == 0 ) {
			pose.append_residue_by_jump( *new_res, pose.total_residue() );
		} else {
			pose.conformation().safely_append_polymer_residue_after_seqpos( *new_res, pose.total_residue(), true );
		}
	}

	// set phi, psi to idealized values
	for ( core::Size i=1; i<=len; ++i ) {
		pose.set_omega( i, 180.0 );
		if ( sheet_.up_down(i) == -1 ) {
			pose.set_phi( i, phi1_ );
			pose.set_psi( i, psi1_ );
		} else {
			pose.set_phi( i, phi2_ );
			pose.set_psi( i, psi2_ );
		}
	}
	return pose;
}

} // namespace denovo_design

} // namespace devel
