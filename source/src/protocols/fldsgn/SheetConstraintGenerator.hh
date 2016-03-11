// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/fldsgn/SheetConstraintGenerator.hh
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky (tlinsky@uw.edu), Nov 2012


#ifndef INCLUDED_protocols_fldsgn_SheetConstraintGenerator_hh
#define INCLUDED_protocols_fldsgn_SheetConstraintGenerator_hh

// Unit Header
#include <protocols/fldsgn/SheetConstraintGenerator.fwd.hh>

// Package Header
#include <protocols/fldsgn/topology/StrandPairing.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

// Proeject Header
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/types.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>


#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {

struct ResiduePair : public std::pair< core::Size, core::Size > {
	ResiduePair( core::Size res1, core::Size const res2, char const orient ) :
		std::pair< core::Size, core::Size >( res1, res2 ), orientation( orient ) {}
	char orientation;

	friend std::ostream & operator<<( std::ostream & os, ResiduePair const & rp ) {
		os << "{ " << rp.first << " " << rp.second << " " << rp.orientation << " }";
		return os;
	}
private:
	ResiduePair() {};
};
typedef utility::vector1< ResiduePair > ResiduePairs;

class SheetConstraintGenerator : public protocols::forge::remodel::RemodelConstraintGenerator {
public:

	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef protocols::jd2::parser::BluePrintOP BluePrintOP;


public:
	SheetConstraintGenerator();

	virtual ~SheetConstraintGenerator();

	virtual void
	parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	virtual std::string
	get_name() const;

	virtual protocols::moves::MoverOP
	fresh_instance() const;

	virtual protocols::moves::MoverOP
	clone() const;

	virtual core::scoring::constraints::ConstraintCOPs
	generate_constraints( Pose const & pose );

	/// @brief sets the secondary structure to be used for constraint generation
	void set_secstruct( std::string const & ss );

	/// @brief sets the strand pair string to be used for constraint generation (e.g. "2-3.P.1" )
	///        means strands 2 and 3 are paired in a parallel orientation with register shift of 1
	void set_strand_pairs( std::string const & spairs );

	/// @brief if set, a square-well, "flat bottom" function will be used for the constraints.
	///        Otherwise, harmonic constraints will be used.
	void set_flat_bottom_constraints( bool const flat_csts );

	/// @brief set the weight of the sheet constraints
	void set_weight( Real const coef );

	/// @brief set the maximum Ca-Ca distance between paired residues
	void set_distance( Real const dist );

	/// @brief sets the tolerance for the func to all atom pair constraints
	void set_distance_tolerance( Real const dist_tol );

	/// @brief set the flat-bottom tolerance for the backbone angle between strands for each pair
	/// This is N1-C1-C2 and N2-C2-C1 for parallel sheets, and N1-C1-N2/N2-C2-N1 for antiparallel.
	void set_angle_tolerance( Real const angle_tolerance );

	/// @brief set the flat-bottom tolerance for the Cb1-Ca1-Ca2-Cb2 dihedral angle (0 = optimal)
	void set_cacb_dihedral_tolerance( Real const dihedral_tolerance );

	/// @brief set the flat-bottom tolerance for the backbone dihedrals (0=optimal)
	/// Dihedral 1 = O1-N1-C1-C2, Dihedral 2 = O2-N2-C2-C1
	void set_bb_dihedral_tolerance( Real const dihedral_tolerance );

	/// @brief sets whether we should constrain distance only, and not generate dihedral and angle constraints
	void set_constrain_ca_ca_dist( bool const constrain_dist_only );

	/// @brief sets whether we should constrain cb_ca_ca_cb dihedrals
	void set_constrain_bb_cacb_dihedral( bool const constrain_dihedral );

	/// @brief sets whether we should constrain bb dihedrals
	void set_constrain_bb_dihedral( bool const constrain_dihedral );

	/// @brief sets whether we should constrain bb angle
	void set_constrain_bb_angle( bool const constrain_angle );

	/// @brief initialize from a blueprint
	void initialize_from_blueprint( std::string const & blueprint_file );
	void initialize_from_blueprint( protocols::jd2::parser::BluePrintCOP bp );

	std::pair< std::string, std::string > get_secstruct_and_strandpairings( Pose const & pose ) const;

	/// @brief computes and returns a vector of residues that are paired in the sheet
	ResiduePairs
	compute_residue_pairs( topology::StrandPairings const & spairs ) const;

private:
	// func generation
	core::scoring::func::FuncOP weighted_func( core::scoring::func::FuncOP func ) const;
	core::scoring::func::FuncOP create_caca_atom_pair_func( core::Real const ideal_dist ) const;
	core::scoring::func::FuncOP create_bb_angle_func( core::Real const ideal_angle ) const;
	core::scoring::func::FuncOP create_bb_dihedral_func( core::Real const ideal_dihedral ) const;
	core::scoring::func::FuncOP create_cacb_dihedral_func( core::Real const ideal_dihedral ) const;

	// constraint generation
	core::scoring::constraints::ConstraintOP create_bb_dihedral_constraint(
		core::pose::Pose const & pose,
		core::Size const res1,
		core::Size const res2,
		core::scoring::func::FuncOP func1,
		core::scoring::func::FuncOP func2 ) const;

	core::scoring::constraints::ConstraintOP
	create_ca_ca_atom_pair_constraint(
		core::id::AtomID const & atom1,
		core::id::AtomID const & atom2,
		core::scoring::func::FuncOP func ) const;

	core::scoring::constraints::ConstraintOP
	create_bb_angle_constraint(
		core::id::AtomID const & atom1,
		core::id::AtomID const & atom2,
		core::id::AtomID const & atom3,
		core::scoring::func::FuncOP func ) const;

	core::scoring::constraints::ConstraintOP
	create_bb_cacb_dihedral_constraint(
		core::id::AtomID const & atom1,
		core::id::AtomID const & atom2,
		core::id::AtomID const & atom3,
		core::id::AtomID const & atom4,
		core::scoring::func::FuncOP func ) const;

private:

	Real weight_;
	Real dist_;
	Real dist_tolerance_;
	Real angle_tolerance_;
	Real cacb_dihedral_tolerance_;
	Real bb_dihedral_tolerance_;
	bool constrain_ca_ca_dist_;
	bool constrain_bb_cacb_dihedral_;
	bool constrain_bb_dihedral_;
	bool constrain_bb_angle_;
	bool flat_bottom_constraints_;

	BluePrintOP blueprint_;
	std::string secstruct_;
	std::string spairs_;
	std::string sheet_name_;
}; //class SheetConstraintGenerator

} //namespace fldsgn
} //namespace protocols


#endif // INCLUDED_protocols_fldsgn_SheetConstraintGenerator_HH
