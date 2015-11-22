// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/select/util/SelectResiduesByLayer.hh
/// @brief Select residues depending on layer: core, boundary, and surface
/// and the layer are defined by accessible sufrace area of each residue
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )
/// @author Gabe Rocklin (sidechain neighbour selection)
/// @author Vikram K. Mulligan (vmullig@uw.edu -- moving this class to core and refactoring for noncanonicals)

#ifndef INCLUDED_core_select_util_SelectResiduesByLayer_hh
#define INCLUDED_core_select_util_SelectResiduesByLayer_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/AA.hh>
#include <core/select/util/SelectResiduesByLayer.fwd.hh>

#include <map>
#include <utility/vector1.hh>

#include <utility/pointer/ReferenceCount.hh>


namespace core {
namespace select {
namespace util {

class SelectResiduesByLayer : public utility::pointer::ReferenceCount {
public:


	typedef core::Size Size;
	typedef core::Real Real;
	typedef std::string String;
	typedef core::chemical::AA AA;
	typedef core::pose::Pose Pose;


public: // constructor/destructor

	/// @brief Default constructor
	///
	SelectResiduesByLayer();

	/// @brief Copy constructor
	///
	SelectResiduesByLayer( SelectResiduesByLayer const &src );

	/// @brief Value constructor
	///
	SelectResiduesByLayer( bool const pick_core, bool const pick_boundary, bool const pick_surface );

	/// @brief Value constructor
	/// @detail Selected layer can be given by string, for example, core_surface or core_boundary_surface
	SelectResiduesByLayer( String const &pick );

	/// @brief Destructor
	///
	virtual ~SelectResiduesByLayer();

	/// @brief Clone operator.
	/// @details Construct a copy of this object and return an owning pointer to the copy.
	SelectResiduesByLayerOP clone() const;

public:


	void initialize( Real const &burial, Real const &surface );


public:  // mutator


	/// @brief
	inline void
	set_design_layer( bool const pick_core, bool const pick_boundary, bool const pick_surface )
	{
		pick_core_ = pick_core;
		pick_boundary_ = pick_boundary;
		pick_surface_ = pick_surface;
	}


	/// @brief set pore radius for colculating asa
	inline void pore_radius( Real const ps )
	{
		pore_radius_ = ps;
	}

	/// @brief accessible surface for evaluating residues are in surface or not
	void sasa_surface( Real const &r, String const &ss="" );

	/// @brief accessible surface for evaluating residues are in core or not
	void sasa_core( Real const &r, String const &ss="" );

	/// @brief accessible surface for evaluating residues are in core or not
	void make_rasmol_format_file( bool const b )
	{
		make_rasmol_format_file_ = b;
	}

	inline void use_sidechain_neighbors( bool const b )
	{
		use_sidechain_neighbors_ = b;
	}

	/// @brief
	void exclude_aatypes_for_selection( utility::vector1< AA > const & aas );

	/// @brief
	void restrict_aatypes_for_selection( utility::vector1< AA > const & aas );


public:  // accessor

	/// @brief accessbile surface are of each residue
	Real rsd_sasa( Size const i ) const;

	/// @brief return defined layer for each residue
	String layer( Size const i ) const;

	/// @brief selected residues on boundary
	utility::vector1< Size > const & selected_boundary_residues() const;

	/// @brief selected residues in core
	utility::vector1< Size > const & selected_core_residues() const;

	/// @brief selected residues on surface
	utility::vector1< Size > const & selected_surface_residues() const;

	bool use_sidechain_neighbors() const;

	/// @brief Midpoint of the distance falloff sigmoid.
	/// @details Defaults to 9.0.  Only used by the sidchain_neighbors code.
	inline core::Real const & dist_midpoint() const { return dist_midpoint_; }

	/// @brief Set the midpoint of the distance falloff sigmoid.
	///
	inline void set_dist_midpoint( core::Real const &val ) { dist_midpoint_ = val; return; }

	/// @brief Factor by which number of residue neighbors is divided.
	/// @details Defaults to 1.0.  Only used by the sidchain_neighbors code.
	inline core::Real const & rsd_neighbor_denominator() const { return rsd_neighbor_denominator_; }

	/// @brief Set the value by which the number of residue neighbors is divided.
	///
	inline void set_rsd_neighbor_denominator( core::Real const &val ) { rsd_neighbor_denominator_ = val; return; }

	/// @brief Get the shift factor for the angle term.
	/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
	/// vector.  The shift factor (default 0.5) is then added, and the resulting value is raised to the angle_exponent and multiplied by the distance factor.
	inline core::Real const & angle_shift_factor() const { return angle_shift_factor_; }

	/// @brief Set the shift factor for the angle term.
	/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
	/// vector.  The shift factor (default 0.5) is then added, and the resulting value is raised to the angle_exponent and multiplied by the distance factor.
	inline void set_angle_shift_factor( core::Real const &val ) { angle_shift_factor_ = val; return; }

	/// @brief Get the exponent for the angle term, which affects how close other atoms have to be to the CA-CB line to be counted fully.
	/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
	/// vector.  The shift factor is then added, and the resulting value is raised to the angle_exponent (default 2.0) and multiplied by the distance factor.
	inline core::Real const & angle_exponent() const { return angle_exponent_; }

	/// @brief Set the exponent for the angle term, which affects how close other atoms have to be to the CA-CB line to be counted fully.
	/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
	/// vector.  The shift factor is then added, and the resulting value is raised to the angle_exponent (default 2.0) and multiplied by the distance factor.
	inline void set_angle_exponent ( core::Real const &val ) { angle_exponent_ = val; return; }

	/// @brief Get the exponent for the distance term, which affects how sharp the falloff is with distance.
	/// @details The distance term is: 1/(1+exp(n*(d - m))), where d is the distance, n is the exponent set by this term,
	/// and m is the midpoint of the falloff.  The n value sets the sharpness.  Defaults to 1.0.
	inline core::Real const & dist_exponent() const { return dist_exponent_; }

	/// @brief Set the exponent for the distance term, which affects how sharp the falloff is with distance.
	/// @details The distance term is: 1/(1+exp(n*(d - m))), where d is the distance, n is the exponent set by this term,
	/// and m is the midpoint of the falloff.  The n value sets the sharpness.  Defaults to 1.0.
	inline void set_dist_exponent ( core::Real const &val ) { dist_exponent_ = val; return; }

public: // main operation


	/// @brief apply
	utility::vector1< Size > const
	compute( Pose const & pose, String const &secstruct="", bool const skip_dssp=false );


private: // helper functions

	/// @brief
	utility::vector1< Real > const
	calc_rsd_sasa( Pose const & pose ) const;

	utility::vector1< Real > const
	calc_sc_neighbors( Pose const & pose ) const;

protected:


	/// @brief design core ?
	bool pick_core_;

	/// @brief design boundary ?
	bool pick_boundary_;

	/// @brief design surface ?
	bool pick_surface_;

	/// @brief pore radius for calculating asa( accessible surface area )
	Real pore_radius_;

	/// @brief asa pore radius for calculating asa
	std::map< char, Real > burial_;

	/// @brief pore radius for calculating asa
	std::map< char, Real > surface_;

	/// @brief amino acid types excluded for selection
	utility::vector1< AA > excluded_aatypes_for_selection_;

	/// @brief amino acid types restricted for selection
	utility::vector1< AA > restricted_aatypes_for_selection_;

	/// @brief selected residues in core
	utility::vector1< Size > selected_core_residues_;

	/// @brief selected residues at boundary
	utility::vector1< Size > selected_boundary_residues_;

	/// @brief selected residues on surface
	utility::vector1< Size > selected_surface_residues_;

	/// @brief create rasmol format script for selected residues; red: surface, blue: core, green: boundary
	bool make_rasmol_format_file_;

	bool use_sidechain_neighbors_;

	/// @brief output in rasmol format


	/// @brief asa for each residue
	utility::vector1< Real > rsd_sasa_;

	/// @brief
	utility::vector1< String > rsd_layer_;

	/// @brief Midpoint for the distance falloff sigmoid.  Defaults to 9.0.
	/// @details Only used by the sidechain_neighbors code.
	core::Real dist_midpoint_;

	/// @brief Divisor for the number of residue neighbours.  Defaults to 1.0.
	/// @details Only used by the sidechain_neighbors code.
	core::Real rsd_neighbor_denominator_;

	/// @brief Shift factor for the angle term.
	/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
	/// vector.  The shift factor (default 0.5) is then added, and the resulting value is raised to the angle_exponent and multiplied by the distance factor.
	core::Real angle_shift_factor_;

	/// @brief Exponent for the angle term, which affects how close other atoms have to be to the CA-CB line to be counted fully.
	/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
	/// vector.  The shift factor is then added, and the resulting value is raised to the angle_exponent (default 2.0) and multiplied by the distance factor.
	core::Real angle_exponent_;

	/// @brief Exponent for the distance term, which affects how sharp the falloff is with distance.
	/// @details The distance term is: 1/(1+exp(n*(d - m))), where d is the distance, n is the exponent set by this term,
	/// and m is the midpoint of the falloff.  The n value sets the sharpness.  Defaults to 1.0.
	core::Real dist_exponent_;

};


} // util
} // select
} // core

#endif
