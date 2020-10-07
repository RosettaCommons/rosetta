// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/pack/guidance_scoreterms/sap/SapMathConstraint.hh
/// @brief A constraint that allows you to subtract and add other SapConstraints
/// @details
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_SapMathConstraint_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_SapMathConstraint_hh

#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>
#include <core/pack/guidance_scoreterms/sap/SapMathConstraint.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintHelper.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

#include <string>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

class SapMathConstraint : public core::scoring::aa_composition_energy::SequenceConstraint {

public: //Constructor, destructor, copy, clone:

	/// @brief Constructor
	SapMathConstraint();

	/// @brief Copy constructor
	SapMathConstraint( SapMathConstraint const & src );

	/// @brief Destructor
	~SapMathConstraint() override;

	/// @brief Clone operator
	scoring::constraints::ConstraintOP clone() const override;

	SapMathConstraint& operator=( SapMathConstraint const & other );

	bool operator == ( scoring::constraints::Constraint const & /*other*/ ) const override;

	bool
	same_type_as_me( scoring::constraints::Constraint const & other ) const override;


public: //Functions that actually do stuff:


	void add_constraint( Real weight, std::string const & name );

	void upper_bound( Real upper );

	void lower_bound( Real lower );

	void penalty_per_unit( Real penalty );



	utility::vector1< std::pair< Real, SapConstraintHelperCOP > >
	parse_helpers( utility::vector1< SapConstraintHelperCOP > const & helpers) const;

	Real
	get_math_result( utility::vector1< std::pair< Real, SapConstraintHelperCOP > > const & parsed_helpers ) const;

	Real
	get_score( utility::vector1< std::pair< Real, SapConstraintHelperCOP > > const & parsed_helpers ) const;


	scoring::constraints::ConstraintOP remap_resid( core::id::SequenceMapping const & ) const override { return nullptr; }

private:
	// Member variables

	utility::vector1< std::pair< Real, std::string > > helper_weights_;
	Real upper_bound_;
	Real lower_bound_;
	Real penalty_per_unit_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //sap
} //guidance_scoreterms
} //pack
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_SapMathConstraint )
#endif // SERIALIZATION

#endif
