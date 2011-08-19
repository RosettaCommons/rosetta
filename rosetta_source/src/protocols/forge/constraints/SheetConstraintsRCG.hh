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


#ifndef INCLUDED_protocols_forge_constraints_SheetConstraintsRCG_hh
#define INCLUDED_protocols_forge_constraints_SheetConstraintsRCG_hh

// Unit Header
#include <protocols/forge/constraints/SheetConstraintsRCG.fwd.hh>

// Package Header
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

// Proeject Header
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/fldsgn/BluePrint.fwd.hh>

#include <string>

namespace protocols{
namespace forge{
namespace constraints{

class SheetConstraintsRCG : public protocols::forge::remodel::RemodelConstraintGenerator {
public:

	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef protocols::fldsgn::BluePrintOP BluePrintOP;


public:

	SheetConstraintsRCG( BluePrintOP const & blue );

	SheetConstraintsRCG( BluePrintOP const & blue, Real const coef );

	SheetConstraintsRCG( BluePrintOP const & blue, Real const coef, Real const dist );

	virtual ~SheetConstraintsRCG();

	virtual
	void generate_remodel_constraints( Pose const & pose );

	void set_blueprint( BluePrintOP const & blue );

	void set_weight( Real const coef );

	void set_distance( Real const dist );

private:

	BluePrintOP blueprint_;
	Real coef_;
	Real dist_;


}; //class SheetConstraintsRCG


} //namespace remodel
} //namespace forge
} //namespace protocols




#endif // INCLUDED_protocols_forge_remodel_SheetConstraintsRCG_HH
