// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/constraints/SheetConstraintsRCG.hh
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky (tlinsky@uw.edu), Nov 2012


#ifndef INCLUDED_protocols_fldsgn_SheetConstraintsRCG_hh
#define INCLUDED_protocols_fldsgn_SheetConstraintsRCG_hh

// Unit Header
#include <protocols/fldsgn/SheetConstraintsRCG.fwd.hh>

// Package Header
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

// Proeject Header
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>


#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {

class SheetConstraintsRCG : public protocols::forge::remodel::RemodelConstraintGenerator {
public:

	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef protocols::jd2::parser::BluePrintOP BluePrintOP;


public:
	SheetConstraintsRCG();

	SheetConstraintsRCG( SheetConstraintsRCG const & rval );

	SheetConstraintsRCG( BluePrintOP const & blue );

	SheetConstraintsRCG( BluePrintOP const & blue, Real const coef );

	SheetConstraintsRCG( BluePrintOP const & blue, Real const coef, Real const dist );

	virtual ~SheetConstraintsRCG();

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

	virtual
	void generate_remodel_constraints( Pose const & pose );

	/// @brief sets teh blueprint file used for determining proper sheet pairing
	/// This function will create the blueprint object for you from the file and use it
	void set_blueprint( std::string const & blueprint_file );

	/// @brief sets the blueprint object used for determining proper sheet pairing
	void set_blueprint( BluePrintOP const & blue );

	/// @brief set the weight of the sheet constraints
	void set_weight( Real const coef );

	/// @brief set the maximum Ca-Ca distance between paired residues
	void set_distance( Real const dist );

	/// @brief set the flat-bottom tolerance for the backbone angle between strands for each pair
	/// This is N1-C1-C2 and N2-C2-C1 for parallel sheets, and N1-C1-N2/N2-C2-N1 for antiparallel.
	void set_angle_tolerance( Real const angle_tolerance );

	/// @brief set the flat-bottom tolerance for the Cb1-Ca1-Ca2-Cb2 dihedral angle (0 = optimal)
	void set_cacb_dihedral_tolerance( Real const dihedral_tolerance );

	/// @brief set the flat-bottom tolerance for the backbone dihedrals (0=optimal)
	/// Dihedral 1 = O1-N1-C1-C2, Dihedral 2 = O2-N2-C2-C1
	void set_bb_dihedral_tolerance( Real const dihedral_tolerance );

	/// @brief sets whether we should constrain distance only, and not generate dihedral and angle constraints
	void set_constrain_dist_only( bool const constrain_dist_only );
private:

	Real weight_;
	Real dist_;
	Real angle_tolerance_;
	Real cacb_dihedral_tolerance_;
	Real bb_dihedral_tolerance_;
	bool constrain_dist_only_;
	BluePrintOP blueprint_;
}; //class SheetConstraintsRCG


} //namespace fldsgn
} //namespace protocols


#endif // INCLUDED_protocols_fldsgn_SheetConstraintsRCG_HH
