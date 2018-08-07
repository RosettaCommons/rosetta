// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/MakeBundleHelix.hh
/// @brief  Headers for MakeBundleHelix.cc.  Builds a single helix as part of a helical bundle.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_MakeBundleHelix_hh
#define INCLUDED_protocols_helical_bundle_MakeBundleHelix_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/MakeBundleHelix.fwd.hh>
#include <protocols/helical_bundle/util.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.fwd.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/ContingentFilter.fwd.hh>
#include <protocols/filters/ContingentFilter.hh>

// Project Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace helical_bundle {

class MakeBundleHelix : public protocols::moves::Mover
{
public:
	//Typedefs:
	typedef core::conformation::parametric::Parameters Parameters;
	typedef core::conformation::parametric::ParametersOP ParametersOP;
	typedef core::conformation::parametric::ParametersSet ParametersSet;
	typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

	typedef protocols::helical_bundle::parameters::BundleParameters BundleParameters;
	typedef protocols::helical_bundle::parameters::BundleParametersOP BundleParametersOP;
	typedef protocols::helical_bundle::parameters::BundleParametersCOP BundleParametersCOP;
	typedef protocols::helical_bundle::parameters::BundleParametersSet BundleParametersSet;
	typedef protocols::helical_bundle::parameters::BundleParametersSetOP BundleParametersSetOP;
	typedef protocols::helical_bundle::parameters::BundleParametersSetCOP BundleParametersSetCOP;

public:

	/// @brief Default constructor.
	MakeBundleHelix();

	/// @brief Copy constructor.
	MakeBundleHelix( MakeBundleHelix const &src );

	/// @brief Initialization constructor: initializes this MakeBundleHelix mover with a BundleParametrizationCalculator.
	/// @details Input calculator is cloned.
	MakeBundleHelix( BundleParametrizationCalculatorCOP input_calculator );

	///@brief Destructor.
	~MakeBundleHelix() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Copy the parameter values for parameters that have not been set from the global parameters.
	/// @details This function should be called before apply().
	void copy_unset_params_from_globals( BundleParametrizationCalculatorCOP global_calculator );

	/// @brief Copy the parameter values for parameters that copy values from previous helices, from the previous helices.
	/// @details This function should be called before apply().
	/// @returns Returns true for failure, false for success.
	bool copy_params_from_previous_helices( core::pose::Pose const & prev_helices_pose );

	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;

	// Note that this mover is not intended to be configurable from RosettaScripts.  It is only
	// meant to be invoked from the MakeBundle and BundleGridSampler movers.  As such, it has no
	// parse_my_tag function.
	/*virtual void parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const &
	);*/

	/// @brief Set whether the input pose should be reset prior to building a helix.
	///
	void set_reset_pose( bool const reset_in) { reset_pose_ = reset_in; }

	/// @brief Returns a bool indicating whether the input pose will be reset prior to
	/// building a helix.
	inline bool reset_pose() const { return reset_pose_; }

	/// @brief Set the length of the helix, in residues.
	///
	void set_helix_length( core::Size const helix_length_in ) {
		runtime_assert_string_msg(helix_length_in >=2, "In protocols::helical_bundle::MakeBundleHelix::set_helix_length(): A helix must contain at least two residues.");
		helix_length_ = helix_length_in;
	}

	/// @brief Returns a bool indicating whether the input pose will be reset prior to
	/// building a helix.
	inline core::Size helix_length() const { return helix_length_; }

	/// @brief Set the residue type(s) (full name, not 3-letter code) that will make up the helix.
	/// @details If there is more than one residue per repeating unit in the minor helix, one residue
	/// name must be provided for each residue in the repeating unit.  Note that there is no check that
	/// the size of the residue_name_ vector matches the number of residues in the repeating unit until
	/// apply time.
	void set_residue_name(utility::vector1< std::string > const &names) { residue_name_=names; return; }

	/// @brief Get the name (full name, not 3-letter code) of one of the residue types in the repeating
	/// unit that will make up the helix.
	std::string const & residue_name( core::Size const repeat_index ) const {
		runtime_assert_string_msg(
			repeat_index > 0 && repeat_index <= residue_name_.size(),
			"Error in protocols::helical_bundle::MakeBundleHelix::residue_name(): The index is out of range."
		);
		return residue_name_[repeat_index];
	}

	/// @brief Set the minor helix parameters by reading them in from a file.
	///
	void set_minor_helix_params_from_file ( std::string const &filename );

	/// @brief Returns "true" if the last call to the apply function failed, false otherwise.
	///
	inline bool last_apply_failed() const { return last_apply_failed_; }

	/// @brief Non-const access to the calculator.
	BundleParametrizationCalculatorOP calculator_op() { return calculator_; }

	/// @brief Const access to the calculator.
	BundleParametrizationCalculatorCOP calculator_cop() const { return calculator_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE DATA                                                      //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Should the input pose be reset?
	/// @details Default true.
	bool reset_pose_;

	/// @brief Length of the helix, in residues.  Defaults to 10.
	///
	core::Size helix_length_;

	/// @brief Name (full-length, not 3-letter code) of the residue type(s) that the helix will be made of.
	/// @details Defaults to "ALA".  One residue must be specified for each in the repeating unit.
	utility::vector1< std::string > residue_name_;

	/// @brief Did the last apply fail?
	/// @details Initialized to "false"; "true" if the last apply failed.
	bool last_apply_failed_;

	/// @brief The BundleParametrizationCalculator object, which keeps track of parameter values and
	/// does the actual parametric math.
	BundleParametrizationCalculatorOP calculator_;

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Set whether the last call to the apply() function failed or not.
	/// @details Should only be called by the apply() function.
	void set_last_apply_failed( bool const val ) { last_apply_failed_=val; return; }

	/// @brief Add Crick parameter information to the Conformation object within the pose.
	/// @details This function updates the bundle_parameters_ object's links to residues within the pose,
	/// and then copies the owning pointer into the pose's Conformation object.
	void add_parameter_info_to_pose( core::pose::Pose &pose );

};

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_MakeBundleHelix_hh
