// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/pack/guidance_scoreterms/sap/SapConstraint.hh
/// @brief Apply the sap_score as a sequence constraint to the pose
/// @details
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraint_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraint_hh

#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraint.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>




#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

class SapConstraint : public core::scoring::aa_composition_energy::SequenceConstraint {

public: //Constructor, destructor, copy, clone:

	/// @brief Constructor
	SapConstraint( SapConstraintOptionsCOP const & options );

	/// @brief Copy constructor
	SapConstraint( SapConstraint const & src );

	/// @brief Destructor
	~SapConstraint() override;

	/// @brief Clone operator
	scoring::constraints::ConstraintOP clone() const override;

	SapConstraint& operator=( SapConstraint const & other );

	bool operator == ( scoring::constraints::Constraint const & /*other*/ ) const override;

	bool
	same_type_as_me( scoring::constraints::Constraint const & other ) const override;


public: //Functions that actually do stuff:

	SapConstraintOptionsOP
	get_options();

	SapConstraintOptionsCOP
	get_const_options() const;


	scoring::constraints::ConstraintOP remap_resid( core::id::SequenceMapping const & ) const override { return nullptr; }

private:
	// Member variables

	SapConstraintOptionsOP options_;

#ifdef    SERIALIZATION
public:
	SapConstraint() = default;
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
CEREAL_FORCE_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_SapConstraint )
#endif // SERIALIZATION

#endif
