// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh
/// @brief Options for SapConstraint
/// @details Contains all the options and a method to transform sap_score into sap_constraint
/// @author Brian Coventry (bcov@uw.edu)



#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintOptions_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintOptions_hh

// Unit headers
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapParameterOptions.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Project headers
#include <core/types.hh>
#include <string>

#include <utility/VirtualBase.hh> // AUTO IWYU For VirtualBase

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {


class SapConstraintOptions : public utility::VirtualBase {
public:

	SapConstraintOptions();

	SapConstraintOptionsOP clone() const;

	void sap_goal( Real goal );
	void sap_lb_goal( Real lb_goal );
	void packing_correction( Real adjustment );
	void penalty_per_sap( Real penalty );

	void full_accuracy_when_scoring( bool full_accuracy );

	void name( std::string const & name );

	SapParameterOptions & sap_parameter_options();
	SapParameterOptions const & sap_parameter_options() const;

	void score_selector( select::residue_selector::ResidueSelectorCOP const & sel );
	void sap_calculate_selector( select::residue_selector::ResidueSelectorCOP const & sel );
	void sasa_selector( select::residue_selector::ResidueSelectorCOP const & sel );

	void fast( bool fast );
	void lightning( bool lightning );

	Real sap_goal() const;
	Real sap_lb_goal() const;
	Real packing_correction() const;
	Real penalty_per_sap() const;

	bool full_accuracy_when_scoring() const;

	std::string name() const;

	select::residue_selector::ResidueSelectorCOP score_selector() const;
	select::residue_selector::ResidueSelectorCOP sap_calculate_selector() const;
	select::residue_selector::ResidueSelectorCOP sasa_selector() const;

	bool fast() const;

	bool lightning() const;

	Real transform_sap_to_score( Real score, bool approximate ) const;

	void sanity_check( ) const;

private:

	Real sap_goal_;
	Real sap_lb_goal_;
	Real packing_correction_;
	Real penalty_per_sap_;

	core::select::residue_selector::ResidueSelectorCOP score_selector_;
	core::select::residue_selector::ResidueSelectorCOP sap_calculate_selector_;
	core::select::residue_selector::ResidueSelectorCOP sasa_selector_;

	bool fast_;
	bool lightning_;

	bool full_accuracy_when_scoring_;
	std::string name_;

	SapParameterOptions sap_parameter_options_;


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
CEREAL_FORCE_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_SapConstraintOptions )
#endif // SERIALIZATION

#endif // INCLUDED_core_scoring_EtableEnergy_HH
