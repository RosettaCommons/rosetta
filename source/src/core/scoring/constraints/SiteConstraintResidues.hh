// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/constraints/SiteConstraintResidues.hh
/// @brief This class is an AmbiguousConstraint in which the set is comprised of AtomPairConstraints
/// @brief of an atom of interest in one chain versus the CA of of another sets of residues
/// @author Lei Shi (shilei@uw.edu)

#ifndef INCLUDED_core_scoring_constraints_SiteConstraintResidues_hh
#define INCLUDED_core_scoring_constraints_SiteConstraintResidues_hh

// Unit header

#include <core/scoring/constraints/SiteConstraintResidues.fwd.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>


#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/SequenceMapping.fwd.hh>


#include <core/scoring/EnergyMap.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>


//Utility Headers
#include <numeric/xyzVector.fwd.hh>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


class SiteConstraintResidues : public AmbiguousConstraint {
public:

	/// @brief Constructor
	SiteConstraintResidues();

	/// @brief Constructor
	SiteConstraintResidues( ConstraintCOPs const & cst_in ) ;

	///
	virtual ConstraintOP clone() const;
	virtual bool operator == ( Constraint const & ) const;
	virtual bool same_type_as_me( Constraint const & ) const;
	virtual std::string type() const;

	/// @brief read in constraint defiinition
	void
	read_def( std::istream& data, pose::Pose const& pose,func::FuncFactory const& func_factory );

	void show( std::ostream& out) const;
	void setup_csts( Size res1, std::string name, Size res2, Size res3, core::pose::Pose const & pose, func::FuncOP const & func );

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //SiteConstraintResidues

} //constraints
} //scoring
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_SiteConstraintResidues )
#endif // SERIALIZATION


#endif
