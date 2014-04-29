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
/// @modified Tom Linsky, tlinsky@uw.edu



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
#include <protocols/moves/Mover.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace remodel {


/// @brief pure virtual base class
class RemodelConstraintGenerator : public protocols::moves::Mover
{
public: // typedefs
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~RemodelConstraintGenerator();

	///@brief generates constraints and adds them to the pose
	virtual
	void apply( core::pose::Pose & pose );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	typedef protocols::forge::components::VarLengthBuildAP VarLengthBuildAP;


public: //constructors

	RemodelConstraintGenerator();

	RemodelConstraintGenerator( RemodelConstraintGenerator const & rval );


public:

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

	/// @brief Pose-specific setup routines go here
	virtual void
	init( core::pose::Pose const & ) {};

	VarLengthBuildAP
	vlb() const;

	void
	set_vlb( VarLengthBuildAP vlb );

	std::string
	id() const;

	void
	set_id( std::string const & id );

	void
	set_seqmap( core::id::SequenceMappingCOP seqmap );

	core::id::SequenceMappingCOP
	seqmap() const;

	// Undefined, commenting out to fix PyRosetta build  core::scoring::constraints::ConstraintCOPs constraints() const;

	static core::scoring::constraints::ConstraintCOPs const
	lookup_stored_constraints( std::string const & id );

protected:

	void
	add_constraint(	core::scoring::constraints::ConstraintCOP cst );

	void
	add_constraints(	core::scoring::constraints::ConstraintCOPs csts );

	void
	clear_constraints();

	void
	clear_stored_constraints();

	void
	store_constraints();

private:

	std::string id_;

	core::scoring::constraints::ConstraintCOPs csts_;

	core::id::SequenceMappingCOP seqmap_;

	VarLengthBuildAP vlb_;

	static std::map<std::string,core::scoring::constraints::ConstraintCOPs> cst_map_;

}; //class RemodelConstraintGenerator


} //namespace remodel
} //namespace forge
} //namespace protocols




#endif
