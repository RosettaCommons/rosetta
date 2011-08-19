// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelConstraintGenerator.hh
///
/// @brief  abstract base class that generates constraints during forge loop remodelling
/// @author Florian Richter, floric@u.washington.edu, april 2009




#ifndef INCLUDED_protocols_forge_remodel_RemodelConstraintGenerator_hh
#define INCLUDED_protocols_forge_remodel_RemodelConstraintGenerator_hh

//unit headers
#include <protocols/forge/components/VarLengthBuild.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.fwd.hh>
// AUTO-REMOVED #include <protocols/forge/build/BuildInstruction.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh> // WIN32 INCLUDE
#endif
#include <core/id/SequenceMapping.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
#include <utility/vector1_bool.hh>


namespace protocols {
namespace forge{
namespace remodel{


/// @brief pure virtual base class
class RemodelConstraintGenerator : public utility::pointer::ReferenceCount
{
public: // typedefs


	typedef protocols::forge::components::VarLengthBuildCAP VarLengthBuildCAP;


public: //constructors

	RemodelConstraintGenerator();

	//RemodelConstraintGenerator( RemodelConstraintGenerator const & rval );


public:

	//virtual
	//RemodelConstraintGeneratorOP
	//clone() const = 0;

	void
	add_remodel_constraints_to_pose(
		core::pose::Pose & pose );

	void
	remove_remodel_constraints_from_pose(
		core::pose::Pose & pose ) const;

	virtual
	void
	generate_remodel_constraints(
		core::pose::Pose const & pose ) = 0;

	VarLengthBuildCAP
	vlb() const;

	void
	set_vlb( VarLengthBuildCAP vlb );

	void
	set_seqmap( core::id::SequenceMappingCOP seqmap );

	core::id::SequenceMappingCOP
	seqmap() const;

protected:

	void
	add_constraint(	core::scoring::constraints::ConstraintCOP cst );

	void
	add_constraints(	core::scoring::constraints::ConstraintCOPs csts );

	void
	clear_constraints();

private:

	utility::vector1< core::scoring::constraints::ConstraintCOP > remodel_csts_;

	core::id::SequenceMappingCOP seqmap_;

	VarLengthBuildCAP vlb_;


}; //class RemodelConstraintGenerator


} //namespace remodel
} //namespace forge
} //namespace protocols




#endif
