// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/GeometrySecMatchRPE.hh
/// @brief
/// @author Florian Richter, floric@u.washington.edu, june 09


#ifndef INCLUDED_protocols_match_downstream_GeometrySecMatchRPE_hh
#define INCLUDED_protocols_match_downstream_GeometrySecMatchRPE_hh

// Unit headers
#include <protocols/match/downstream/GeometrySecMatchRPE.fwd.hh>

// Package headers
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

#include <core/chemical/ResidueType.fwd.hh>
#include <utility/vector1_bool.hh>


namespace protocols {
namespace match {
namespace downstream {


/// @brief helper class for GeometrySec
/// abstract base class for distance, angle, and dihedral derived classes
class AtomGeometrySecMatchRPE : public SecMatchResiduePairEvaluator
{

public:

	typedef std::pair< core::Size, core::Size > SizePair;

	typedef core::Real Real;

	AtomGeometrySecMatchRPE( protocols::toolbox::match_enzdes_util::GeomSampleInfo const & gsi );

	~AtomGeometrySecMatchRPE();

	virtual
	bool
	evaluate_residues(
		core::conformation::Residue const & candidate_res,
		core::conformation::Residue const & target_res
	) const = 0;

	virtual
	bool
	require_all_target_residue_atom_coordinates() const;

	virtual
	bool
	require_target_atom_coordinate( Size target_atom_id ) const;


	/// @brief determines if the passed in value is between lowval and highval
	bool
	check_value(
		core::Real value
	) const;

	utility::vector1< SizePair > const &
	at_inds() const {
		return at_inds_; }

	void
	add_at_ind(
		core::Size which_cst_res,
		core::Size atom_ind_in_res
	);

	core::Real
	lowval() const {
		return lowval_;
	}

	core::Real
	highval() const {
		return highval_;
	}

	virtual
	std::string
	print(
		core::chemical::ResidueTypeCOP candidate_restype,
		core::chemical::ResidueTypeCOP target_restype
	) const = 0;


protected:

	void
	clear_at_inds();

	void
	set_lowval( core::Real lowval );

	void
	set_highval( core::Real highval );

private:

	Real lowval_;
	Real highval_;
	utility::vector1< SizePair > at_inds_;

};

/// @brief RPE to figure out if two atoms are within a given distance
/// atoms need to be set through the parent class add_at_ind function
class AtomDistanceSecMatchRPE : public AtomGeometrySecMatchRPE
{

public:

	/// @brief the constructor for this rpe sets the lowval and highval to
	/// the squared values of the initial values, that should make distance
	/// evaluation faster
	AtomDistanceSecMatchRPE( protocols::toolbox::match_enzdes_util::GeomSampleInfo const & gsi );

	virtual
	bool
	evaluate_residues(
		core::conformation::Residue const & candidate_res,
		core::conformation::Residue const & target_res
	) const;

	virtual
	bool
	require_candidate_residue_atoms_to_lie_near_target_atom( Size target_atom_id ) const;

	virtual
	utility::vector1< Size >
	candidate_res_atoms_reqd_near_target_atom(
		Size target_atom_id
	) const;

	virtual
	Real
	max_separation_dist_to_target_atom( Size target_atom_id ) const;

	virtual
	std::string
	print(
		core::chemical::ResidueTypeCOP candidate_restype,
		core::chemical::ResidueTypeCOP target_restype
	) const;

};

/// @brief RPE to figure out if three atoms are within a given angle
/// atoms need to be set through the parent class add_at_ind function
class AtomAngleSecMatchRPE : public AtomGeometrySecMatchRPE
{

public:

	AtomAngleSecMatchRPE( protocols::toolbox::match_enzdes_util::GeomSampleInfo const & gsi );

	virtual
	bool
	evaluate_residues(
		core::conformation::Residue const & candidate_res,
		core::conformation::Residue const & target_res
	) const;

	virtual
	std::string
	print(
		core::chemical::ResidueTypeCOP candidate_restype,
		core::chemical::ResidueTypeCOP target_restype
	) const;

};


/// @brief RPE to figure out if four atoms are within a given dihedral angle
/// atoms need to be set through the parent class add_at_ind function
/// also checks whether a dihedral is periodic, i.e. multiple minima
class AtomDihedralSecMatchRPE : public AtomGeometrySecMatchRPE
{

public:

	AtomDihedralSecMatchRPE( protocols::toolbox::match_enzdes_util::GeomSampleInfo const & gsi );

	virtual
	bool
	evaluate_residues(
		core::conformation::Residue const & candidate_res,
		core::conformation::Residue const & target_res
	) const;

	virtual
	std::string
	print(
		core::chemical::ResidueTypeCOP candidate_restype,
		core::chemical::ResidueTypeCOP target_restype
	) const;

private:

	bool check_periodicity_;

	Real periodicity_;
	Real offset_;
};


/// @brief holds a list of AtomGeometrySecMatchRPEs, that get evaluated in sequence
/// when an instance of this class is asked to evaluate two residues.
class GeometrySecMatchRPE : public SecMatchResiduePairEvaluator
{

public:

	/// @brief convenience constructor from mcfi
	/// the downstream_inds and upstream_inds vector must
	/// contain atoms D1(U1), D2(U2), and D3(U3), respectively,
	/// in that order
	GeometrySecMatchRPE(
		protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo const & mcfi,
		utility::vector1< core::Size > const & downstream_inds,
		utility::vector1< core::Size > const & upstream_inds
	);

	/// @brief empty constructor
	GeometrySecMatchRPE(){}

	/// @brief performs a logical AND for all of the AtomGeometry evaluators.
	virtual
	bool
	evaluate_residues(
		core::conformation::Residue const & candidate_res,
		core::conformation::Residue const & target_res
	) const;

	void
	add_atomgeom_evaluator(
		AtomGeometrySecMatchRPECOP evaluator
	);

	virtual
	bool
	require_all_target_residue_atom_coordinates() const;

	virtual
	bool
	require_target_atom_coordinate( Size target_atom_id ) const;

	virtual
	bool
	require_candidate_residue_atoms_to_lie_near_target_atom( Size target_atom_id ) const;

	virtual
	utility::vector1< Size >
	candidate_res_atoms_reqd_near_target_atom(
		Size target_atom_id
	) const;

	virtual
	Real
	max_separation_dist_to_target_atom( Size target_atom_id ) const;

	utility::vector1< AtomGeometrySecMatchRPECOP > const &
	atom_geom_rpes() const {
		return atom_geom_rpes_;
	}

private:

	utility::vector1< AtomGeometrySecMatchRPECOP > atom_geom_rpes_;

};

}
}
}

#endif
