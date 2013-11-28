// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/BuildSheet.hh
/// @brief The BuildSheet Protocol
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_devel_denovo_design_BuildSheet_hh
#define INCLUDED_devel_denovo_design_BuildSheet_hh

// Unit headers
#include <devel/denovo_design/BuildSheet.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <devel/denovo_design/ParametricSheet.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <basic/datacache/DataMap.hh>

//// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>


namespace devel {
namespace denovo_design {

class BuildSheet : public protocols::moves::Mover {
public:
	struct AngleData {
		AngleData() : psi( 0.0 ), phi( 0.0 ), h_angle( 0.0 ), set( false ) {}
		AngleData( core::Real const phi_val, core::Real const psi_val, core::Real const h_angle_val )
			: psi( psi_val ), phi( phi_val ), h_angle( h_angle_val ), set( true ) {}
		core::Real psi;
		core::Real phi;
		core::Real h_angle;
		bool set;
	};

	struct StrandData {
		StrandData( std::string const name_val,
				core::Size const len,
				int const reg_shift,
				std::string const orient )
			: name( name_val),
			length( len ),
			register_shift( reg_shift ),
			orientation( orient )
		{}
		std::string name;
		core::Size length;
		int register_shift;
		std::string orientation;
	};

	/// @brief default constructor
	BuildSheet();

  /// @brief virtual constructor to allow derivation
	virtual ~BuildSheet();

  /// @brief Parses the BuildSheetTags
	void parse_my_tag(
	  utility::tag::TagCOP const tag,
	  basic::datacache::DataMap & data,
	  protocols::filters::Filters_map const &,
	  protocols::moves::Movers_map const &,
	  core::pose::Pose const &
	);

  /// @brief Return the name of this mover.
  virtual std::string get_name() const;

  /// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::moves::MoverOP clone() const;

  /// @brief Apply the BuildSheet. Overloaded apply function from mover base class.
	virtual void apply( core::pose::Pose & pose );

	/// @brief physically builds a strand based on the provided ca coordinates
	void
	build_sheet( core::pose::Pose & pose ) const;

	/// @brief builds a small idealized strand fragment to be used for parameterized stuff
	core::pose::Pose
	build_ideal_strand( core::chemical::ResidueTypeSetCAP restype_set,
											std::string const & res_name,
											core::Size const len ) const;

	/// @brief minimizes residue resid with cartesian minimization, including strong Ca constraints
	/// assumes that the Ca positions are correct
	/// the "residues" vector contains pose residue ids for the current residue, previous residue in the C-direction, and previous residue in the H-direction, IN THAT ORDER
	/// the "strands" vector contains ca_coords_ indices for the same three residues, in the same order
	/// the "resis" vector contains ca_coords_ residue indices for the same three residues, in the same order
	void
	minimize_residues( core::pose::Pose & pose, utility::vector1< core::Size > const & residues,
								utility::vector1< core::Size > const & strands,
								utility::vector1< core::Size > const & resis ) const;

	/// @brief initializes and rests the dihedral map
	void init_dihedral_map();

	/// @brief reads a map of Ca coordinates vs. phi/psi angles from a db file
	void read_dihedral_map();

	/// @brief builds a map of Ca coordinates vs. phi/psi angles
	void build_dihedral_map();

	/// @brief given a dihedral theta angle, return a proper index in the AngleData vector
	core::Size
	index_from_theta( core::Real const theta ) const;

	/// @brief given a dihedral phi angle, return a proper index in the AngleData vector
	core::Size
	index_from_phi( core::Real const phi ) const;

	/// @brief given a dihedral phi angle, return a proper index in the AngleData vector
	core::Size
	index_from_hangle( core::Real const h_angle ) const;

	/// @brief given a dihedral angle, return an integer that refers to the proper index in the AngleData vector
	core::Size
	index_from_angle( core::Real angle, core::Real const max_value ) const;

	/// @brief calculates an angle which represents the middle of the given hash index
	core::Real
	angle_from_index( core::Size const index, core::Real const max_value ) const;

	/// @brief aligns residues 1-3 of the pose to the X-Y plane, with CA3 as 0,0,0
	void
	align_pose( core::pose::Pose & pose ) const;

	/// @brief given a set of four Ca-points on strand i ( j-2, j-1, j, j+1 ) estimate the appropriate dihedrals
	AngleData const &
	lookup_dihedrals( core::Size const i, core::Size const j ) const;

private: // options/parameters
	/// @brief if true, residues will be minimized as they are added to the sheet
	bool minimize_;

	core::Real phi1_;
	core::Real phi2_;
	core::Real psi1_;
	core::Real psi2_;
	core::Size num_cached_dihedrals_;

private:   // other data
	ParametricSheet sheet_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< utility::vector1< utility::vector1< AngleData > > > dihedral_map_;
};

// helper functions

/// @brief given a vector, return spherical coordinate for theta for that vector
core::Real calc_spherical_theta( core::Vector const & vec );

/// @brief given a vector, return spherical coordinate for theta for that vector
core::Real calc_spherical_phi( core::Vector const & vec );

/// @brief creates a rotation matrix that places p1, p2, and p_center in the XY plane with p_center at the origin and p2 on the x axis
numeric::xyzMatrix< core::Real >
align_matrix( core::Vector const & p1,
							core::Vector const & p2,
							core::Vector const & p_center );

} // denovo_design
} // devel

#endif
