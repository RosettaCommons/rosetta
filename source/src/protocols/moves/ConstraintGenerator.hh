// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/ConstraintGenerator.hh
///
/// @brief  abstract base class that generates constraints
/// @author Florian Richter, floric@u.washington.edu, april 2009
/// @modified Tom Linsky, tlinsky@uw.edu


#ifndef INCLUDED_protocols_moves_ConstraintGenerator_hh
#define INCLUDED_protocols_moves_ConstraintGenerator_hh

//unit headers
#include <protocols/moves/ConstraintGenerator.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh> // WIN32 INCLUDE
#endif
#include <core/id/SequenceMapping.fwd.hh>
#include <protocols/moves/Mover.hh>

//utility headers
#include <utility/excn/EXCN_Base.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {

/// @brief pure virtual base class
class ConstraintGenerator : public protocols::moves::Mover
{
public: // constructors / destructors

	ConstraintGenerator();

	virtual ~ConstraintGenerator();

public: //overridden virtuals

	/// @brief generates constraints and adds them to the pose
	virtual
	void apply( core::pose::Pose & pose );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

public: // virtuals to be implemented by subclasses
	virtual core::scoring::constraints::ConstraintCOPs
	generate_constraints( core::pose::Pose const & pose ) = 0;

public:

	void
	add_constraints_to_pose( core::pose::Pose & pose );

	void
	remove_constraints_from_pose( core::pose::Pose & pose ) const;

	/// @brief Pose-specific setup routines go here
	virtual void
	init( core::pose::Pose const & ) {};

	std::string
	id() const;

	void
	set_id( std::string const & id );

	core::scoring::constraints::ConstraintCOPs const &
	constraints() const;

protected:
	void
	clear_stored_constraints();

	void
	store_constraints( core::scoring::constraints::ConstraintCOPs const & csts );

private:
	std::string id_;
	core::scoring::constraints::ConstraintCOPs csts_;
}; //class ConstraintGenerator

/// @brief class for looking up stored constraint sets
class ConstraintSetManager : public utility::SingletonBase< ConstraintSetManager > {
public:
	typedef std::map< std::string, core::scoring::constraints::ConstraintCOPs > ConstraintsMap;

public:
	ConstraintSetManager();

	static ConstraintSetManager *
	create_singleton_instance();

	bool
	constraints_exist( std::string const & name ) const;

	core::scoring::constraints::ConstraintCOPs const &
	retreive_constraints( std::string const & name ) const;

	void
	store_constraints( std::string const & name, core::scoring::constraints::ConstraintCOPs const & constraints );

	void
	remove_constraints( std::string const & name );

private:
	void
	print_valid_names( std::ostream & os ) const;

private:
	ConstraintsMap cst_map_;
};

class EXCN_RemoveCstsFailed : public utility::excn::EXCN_Base {
public:
	EXCN_RemoveCstsFailed():
		utility::excn::EXCN_Base()
	{}
	virtual void show( std::ostream & os ) const { os << "Remodel constraints somehow got lost along the way" << std::endl; }
};

} //namespace moves
} //namespace protocols


#endif
