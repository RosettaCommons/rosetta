// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/forge/constraints/NtoC_RCG.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky (tlinsky@uw.edu), Nov 2012

// Unit header
#include <protocols/forge/constraints/NtoC_RCG.hh>
#include <protocols/forge/constraints/NtoCCstGeneratorCreator.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/id/SequenceMapping.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

// Project headers
#include <core/chemical/ResidueType.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.forge.constraints.NtoC_RCG" );

namespace protocols {
namespace forge {
namespace constraints {

std::string
NtoCCstGeneratorCreator::keyname() const
{
	return NtoCCstGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
NtoCCstGeneratorCreator::create_mover() const {
	return protocols::moves::MoverOP( new NtoC_RCG() );
}

std::string
NtoCCstGeneratorCreator::mover_name()
{
	return "NtoCCstGenerator";
}

/// @brief
NtoC_RCG::NtoC_RCG():
	RemodelConstraintGenerator(),
	dist_( 11.0 ),
	coef_( 1.0 )
{}

NtoC_RCG::NtoC_RCG( NtoC_RCG const & rval )
: RemodelConstraintGenerator( rval ),
	dist_( rval.dist_ ),
	coef_( rval.coef_ )
{}

/// @brief
NtoC_RCG::NtoC_RCG( Real const dist, Real const coef ):
	RemodelConstraintGenerator(),
	dist_( dist ),
	coef_( coef )
{}

/// @brief
NtoC_RCG::~NtoC_RCG() {}

void
NtoC_RCG::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	RemodelConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	set_weight( tag->getOption< core::Real >( "weight", coef_ ) );
	set_distance( tag->getOption< core::Real >( "dist", dist_ ) );
}

std::string
NtoC_RCG::get_name() const
{
	return NtoCCstGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
NtoC_RCG::fresh_instance() const
{
	return protocols::moves::MoverOP( new NtoC_RCG() );
}

protocols::moves::MoverOP
NtoC_RCG::clone() const
{
	return protocols::moves::MoverOP( new NtoC_RCG( *this ) );
}

/// @brief set weight
void
NtoC_RCG::set_weight( Real const coef )
{
	coef_ = coef;
}

/// @brief set distance of constraint
void
NtoC_RCG::set_distance( Real const dist )
{
	dist_ = dist;
}


/// @brief
void
NtoC_RCG::generate_remodel_constraints( Pose const & pose )
{
	using namespace core::scoring::constraints;

	std::string tag( "constraint_between_N_&_C_terminal_Calpha" );
	Real lb( 0.0 );
	Real ub( dist_ );
	Real sd( 1.0 );
	core::scoring::func::ScalarWeightedFuncOP cstfunc( new core::scoring::func::ScalarWeightedFunc( coef_, core::scoring::func::FuncOP( new BoundFunc( lb, ub, sd, tag ) ) ) );

	Size last_residue = protocols::toolbox::match_enzdes_util::get_last_protein_residue( pose );
	Size first_residue = protocols::toolbox::match_enzdes_util::get_first_protein_residue( pose );
	TR << "first residue in NtoC generation is:" << first_residue << " and last is:" << last_residue << " out of total=" << pose.total_residue() << std::endl;
	core::id::AtomID atom1( pose.residue_type( first_residue ).atom_index( "CA" ), first_residue );
	core::id::AtomID atom2( pose.residue_type( last_residue ).atom_index( "CA" ), last_residue );
	ConstraintOP const cst( new AtomPairConstraint( atom1, atom2, cstfunc ) );

	TR << "Constraints between N- and C- terminal: " << first_residue << "-" << last_residue << ", dist=" << dist_ << std::endl;

	this->add_constraint( cst );
} //generate constraints


} //namespace constraints
} //namespace forge
} //namespace protocols
