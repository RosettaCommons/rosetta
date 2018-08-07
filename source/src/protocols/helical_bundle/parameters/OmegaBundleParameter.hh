// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle/parameters/OmegaBundleParameter.hh
/// @brief Headers for a class for the omega paremeter, derived from the generic RealValuedParameter class.
/// @details Omega has a few additional options that can be configured, for pitch-copying vs. value-copying.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_helical_bundle_parameters_OmegaBundleParameter_hh
#define INCLUDED_protocols_helical_bundle_parameters_OmegaBundleParameter_hh

#include <protocols/helical_bundle/parameters/OmegaBundleParameter.fwd.hh>

// Core headers
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/Parameter.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace helical_bundle {
namespace parameters {

/// @brief A class for the omega paremeter, derived from the generic RealValuedParameter class.
/// @details Omega has a few additional options that can be configured, for pitch-copying vs. value-copying.
class OmegaBundleParameter : public core::conformation::parametric::RealValuedParameter {

public:

	/// @brief Default constructor.
	OmegaBundleParameter();

	/// @brief Copy constructor.
	OmegaBundleParameter(OmegaBundleParameter const & src);

	/// @details Destructor
	virtual ~OmegaBundleParameter();

	/// @brief Clone operator.
	/// @details Make a copy of this object and return a smart pointer to the copy.
	core::conformation::parametric::ParameterOP clone() const override;

	/// @brief Given another parameter of the same type, copy its value.  This does *not* set value_set_ to true.
	/// @details Performs type checking in debug mode.  Note that this version has options for copying pitch or
	/// copying the value directly.
	/// @returns Returns TRUE for failure, FALSE for success.
	bool copy_value_from_parameter( core::conformation::parametric::ParameterCOP other_parameter, core::conformation::parametric::ParametersCOP other_parameter_collection, core::conformation::parametric::ParametersCOP this_parameter_collection ) override;

	/// @brief Set whether we're copying pitch (true) or the twist per residue (false).
	inline void set_copies_pitch( bool const setting ) { pitch_copying_mode_ = setting; }

	/// @brief Get whether we're copying pitch (true) or the twist per residue (false).
	inline bool copies_pitch() const { return pitch_copying_mode_; }

public: //Parse functions

	/// @brief Return the XSD information for this parameter.
	/// @details Calls the equivalent function for RealValuedParameter, then adds additional settings for pitch copying.
	/// @note Must be implemented by derived classes.
	void provide_xsd_information( utility::tag::AttributeList & xsd, bool const provide_setting, bool const provide_copying, bool const provide_grid_sampling, bool const provide_perturbation ) const override;

	/// @brief Reset the copying options.
	/// @details This version resets pitch_copying_mode_, and calls the parent reset_copying_settings() function.
	void reset_copying_settings() override;

protected: //Functions that can be overridden by derived classes.

	/// @brief Parse a setting for this parameter (e.g. r0="12.5").
	/// @details Returns false if this parameter isn't provided.  Also parses copying.
	/// @note May be overridden by derived classes.
	bool parse_setting_only( utility::tag::TagCOP tag, bool const parse_copying_too ) override;

private: //Private functions

	/// @brief Given another parameter of the same type, copy its pitch.  This does *not* set value_set_ to true.
	/// @details Performs type checking in debug mode.
	/// @returns Returns TRUE for failure, FALSE for success.
	bool copy_pitch_from_parameter( core::conformation::parametric::ParameterCOP other_parameter, core::conformation::parametric::ParametersCOP other_parameter_collection, core::conformation::parametric::ParametersCOP this_parameter_collection );

	/// @brief Given a tag, parse out "pitch_from_helix" lines and configure this
	/// parameter accordingly.
	void set_pitch_copying_information( utility::tag::TagCOP tag );

	/// @brief Return the XSD information for copying the pitch value for this parameter from another.
	void provide_xsd_pitch_copying_information( utility::tag::AttributeList & xsd ) const;

private: //Private variables

	/// @brief When copying, do we copy pitch (true) or value (false)?
	bool pitch_copying_mode_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //protocols
} //helical_bundle
} //parameters

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_helical_bundle_parameters_OmegaBundleParameter )
#endif // SERIALIZATION

#endif //INCLUDED_protocols_helical_bundle_parameters_OmegaBundleParameter_hh





