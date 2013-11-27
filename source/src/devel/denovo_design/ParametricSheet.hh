// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/ParametricSheet.hh
/// @brief Parametric representation of a beta sheet
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_devel_denovo_design_ParametricSheet_hh
#define INCLUDED_devel_denovo_design_ParametricSheet_hh

// Unit headers

// Protocol headers

// Package headers
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <basic/datacache/DataMap.hh>

//// C++ headers
#include <string>

// Utility Headers
#include <numeric/xyzMatrix.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace devel {
namespace denovo_design {

class ParametricSheet : public utility::pointer::ReferenceCount {
public:
	struct StrandData {
		StrandData( std::string const name_val,
								core::Size const len,
								core::Size const reg_shift,
								std::string const orient )
			: name( name_val),
				length( len ),
				register_shift( reg_shift ),
				orientation( orient )
		{}
		std::string name;
		core::Size length;
		core::Size register_shift;
		std::string orientation;
	};

	struct SheetGeometry {
		SheetGeometry() :
			cdist_mean( 0.0 ),
			cdist_dev( 0.0 ),
			hdist_mean( 0.0 ),
			hdist_dev( 0.0 ),
			twist_mean( 0.0 ),
			twist_dev( 0.0 )
		{}
		void show( std::ostream & out ) const;
		core::Real cdist_mean;
		core::Real cdist_dev;
		core::Real hdist_mean;
		core::Real hdist_dev;
		core::Real twist_mean;
		core::Real twist_dev;
	};

	/// @brief default constructor
	ParametricSheet();

  /// @brief virtual constructor to allow derivation
	virtual ~ParametricSheet();

	/// @brief setup the parameters via a set of xml tags
	void parse_tags( utility::vector1< utility::tag::TagCOP > const & tags );

	/// @brief setup the parameters based by parsing an input pose
	void from_pose( core::pose::Pose const & pose );

	// Accessors
public:
	/// @brief returns the number of strands in this sheet
	core::Size num_strands() const;

	/// @brief returns the length of strand i
	core::Size strand_len( core::Size const strand ) const;

	/// @brief register shift of strand i
	core::Size register_shift( core::Size const strand ) const;

	/// @brief orientation of strand i
	std::string const & orientation( core::Size const strand ) const;

	/// @brief returns number of rows in the matrix of ca coords
	core::Size ca_coords_size() const;

	/// @brief returns number of ca atoms in strand i
	core::Size ca_coords_size( core::Size const strand ) const;

	/// @brief returns the coordinates for strand strand, residue resi
	core::Vector const & ca_coords( core::Size const strand, core::Size const resi ) const;

	// Mutators
public:
	void c_dist( core::Real const cdist_val );
	void strand_dist( core::Real const sdist_val );
	void twist( core::Real const twist_val );
	void c_coil( core::Real const c_coil_val );
	void h_coil( core::Real const h_coil_val );
	void strand_rotation( core::Real const strand_rotation_val );
	void strand_data( utility::vector1< StrandData > const & strand_data_val );

	/// @brief generates parametric coordinates that basically define the axes of the sheet
	/// returns a vector of vectors of points that make up each strand
	void
	generate_residue_positions();

	/// @brief generates a strand given the strand number we wish to generate and the previous strand
	void
	generate_strand( core::Size const strand, core::Size const prev_strand, core::Size const prev2_strand );

	/// @brief creates a vector for a new residue in the sheet based on the previous c and h vectors
	void
	new_strand( core::Size const strand_num,
							core::Size const prev_strand,
							core::Size const prev2_strand,
							core::Size const center_resi );

	/// @brief generates a transformation matrix to rotate about an arbitrary axis which runs through 0,0,0
	numeric::xyzMatrix< core::Real >
	rotate_about_axis( core::Vector const & twist_axis, core::Real const theta ) const;

	/// @brief generates a Ca coordinate based on a position on the beta sheet ribbon
	void
	create_ca_point( core::Size const strand, core::Size const prev_strand, core::Size const prev2_strand,
									 core::Size const resi, core::Size const prev_resi, core::Size const prev2_resi );

	/// @brief builds a small idealized strand fragment to be used for parameterized stuff
	core::pose::Pose
	build_ideal_strand( core::chemical::ResidueTypeSetCAP restype_set,
											std::string const & res_name,
											core::Size const len ) const;

	/// @brief convert a point in (c-direction, h-direction, z-direction local coordinates to xyz coordinates
	/// c_point1 is assumed to be 0,0,0 in chz coordinates
	core::Vector
	chz_to_xyz( core::Vector const & chz_point,
							core::Vector const & c_point1,
							core::Vector const & c_point2,
							core::Vector const & h_point1,
							core::Vector const & h_point2,
							core::Size const resi ) const;

	/// @brief Take the previous point in the c direction, the previous point in the h direction and a reference point and twist about the axis of the two points with the reference point as 0,0,0
	core::Vector
	twist( core::Vector const & xyz_point,
				 core::Vector const & c_prev,
				 core::Vector const & h_prev,
				 core::Vector const & reference ) const;

	/// @brief initializes the set of ca coords
	void
	init_ca_coords();

	/// @brief dumps ca coordinates and puts them into a pdb file
	void
	dump_ca_coords( std::string const & pdb_name ) const;

	bool
	check_strand_value( core::Size const strand ) const;

	bool
	check_strand_and_resi_values( core::Size const strand, core::Size const resi ) const;

	/// @brief returns 1 if the C=O for this residue points toward strand i+1, -1 if it points toward strand i-1
	int
	up_down( core::Size const resi ) const;

	/// @brief given two previous h points, determine the h direction
	core::Vector
	get_h_axis( core::Vector const & p1,
							core::Vector const & p2,
							core::Vector const & c_axis,
							core::Size const resi ) const;

	/// @brief Ca points are generated in a rectangular grid. This function computes the maximum required height of that grid in the C-direction based on strand lengths and register shifts
	core::Size
	max_strand_len() const;

	/// @brief calculates the mean and std. dev of ca distance along the c axis
	std::pair< core::Real, core::Real >
	calc_cdist() const;

	/// @brief calculates the ca distance along the c axis between the given strand and residue numbers
	core::Real
	calc_cdist( core::Size const strand1, core::Size const resi1,
			core::Size const strand2, core::Size const resi2 ) const;

	/// @brief calculates mean and std dev of twist
	std::pair< core::Real, core::Real >
	calc_twist() const;

	/// @brief given four points, calculate the twist
	core::Real
	calc_twist( core::Vector const & point, core::Vector const & c_prev, core::Vector const & h_prev, core::Vector const & ref ) const;

	/// @brief calculates distance in the h-direction
	SheetGeometry
	calc_geometry() const;

private: // private functions
	void
	test_matrix( utility::vector1< utility::vector1< core::Vector > > const & vec ) const;

private: // options/parameters
	core::Real strand_dist_;
	core::Real c_dist_;
	core::Real twist_;
	core::Real c_coil_;
	core::Real h_coil_;
	core::Real ztwist_;
	utility::vector1< StrandData > strand_data_;
	/// @brief rotation of strands about the C axis
	core::Real strand_rotation_;

private:   // other data
	utility::vector1< utility::vector1< core::Vector > > ca_coords_;
	std::map< std::string, core::Size > strand_name_to_index_;

private: // static const data
	static core::Size const default_strand_length_;
	static int const default_register_shift_;
	static std::string const default_orientation_;
};

// helper functions

void print_vector( numeric::xyzVector< core::Real > const & vec, std::string const & name );

void
print_vector_of_vectors( utility::vector1< utility::vector1< core::Vector > > const & vec );

/// @brief compute the angle between two vectors using atan2 -- both are assumed to come from 0,0,0
core::Real
vector_angle( core::Vector const & v1, core::Vector const & v2 );

} // denovo_design
} // devel

#endif
