// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/MakeBarrelStrand.hh
/// @brief  Headers for MakeBarrelStrand.cc.  Builds a single strand as part of a beta-barrel.
/// @details The barrel is centred on the origin, with the barrel axis pointing along the
/// global z-axis.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_MakeBarrelStrand_hh
#define INCLUDED_protocols_beta_barrel_MakeBarrelStrand_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/beta_barrel/MakeBarrelStrand.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.hh>
#include <protocols/beta_barrel/parameters/BarrelParametersSet.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParametersSet.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <protocols/beta_barrel/BarrelParametrizationCalculator.fwd.hh>

// Scripter Headers
#include <protocols/moves/Mover.fwd.hh>

// Project Headers
#include <utility/vector1.hh>



///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace beta_barrel {

class MakeBarrelStrand : public protocols::moves::Mover
{
public:
	//Typedefs:
	typedef core::conformation::parametric::Parameters Parameters;
	typedef core::conformation::parametric::ParametersOP ParametersOP;
	typedef core::conformation::parametric::ParametersSet ParametersSet;
	typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

	typedef protocols::beta_barrel::parameters::BarrelParameters BarrelParameters;
	typedef protocols::beta_barrel::parameters::BarrelParametersOP BarrelParametersOP;
	typedef protocols::beta_barrel::parameters::BarrelParametersCOP BarrelParametersCOP;
	typedef protocols::beta_barrel::parameters::BarrelParametersSet BarrelParametersSet;
	typedef protocols::beta_barrel::parameters::BarrelParametersSetOP BarrelParametersSetOP;
	typedef protocols::beta_barrel::parameters::BarrelParametersSetCOP BarrelParametersSetCOP;

public:

	/// @brief Default constructor.
	MakeBarrelStrand();

	/// @brief Copy constructor.
	MakeBarrelStrand( MakeBarrelStrand const &src );

	/// @brief Initialization constructor: initializes this MakeBarrelStrand mover with a BarrelParametrizationCalculator.
	/// @details Input calculator is cloned.
	MakeBarrelStrand( BarrelParametrizationCalculatorCOP input_calculator );

	///@brief Destructor.
	~MakeBarrelStrand() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Copy the parameter values for parameters that have not been set from the global parameters.
	/// @details This function should be called before apply().
	void copy_unset_params_from_globals( BarrelParametrizationCalculatorCOP global_calculator );

	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;

	// Note that this mover is not intended to be configurable from RosettaScripts.  It is only
	// meant to be invoked from the MakeBarrel mover.  As such, it has no
	// parse_my_tag function.

	/// @brief Set whether the input pose should be reset prior to building a strand.
	///
	void set_reset_pose( bool const reset_in) { reset_pose_ = reset_in; }

	/// @brief Returns a bool indicating whether the input pose will be reset prior to
	/// building a strand.
	inline bool reset_pose() const { return reset_pose_; }

	/// @brief Set the length of the strand, in residues.
	///
	void set_strand_length( core::Size const strand_length_in ) {
		runtime_assert_string_msg(strand_length_in >=2, "In protocols::beta_barrel::MakeBarrelStrand::set_strand_length(): A strand must contain at least two residues.");
		strand_length_ = strand_length_in;
	}

	/// @brief Returns the length of the strand, in residues.
	///
	inline core::Size strand_length() const { return strand_length_; }

	/// @brief Set the residue type(s) (full name, not 3-letter code) that will make up the strand.
	/// @details If there is more than one residue per repeating unit in the minor helix, one residue
	/// name must be provided for each residue in the repeating unit.
	void set_residue_name(utility::vector1< std::string > const &names) { residue_name_=names; return; }

	/// @brief Get the name (full name, not 3-letter code) of one of the residue types in the repeating
	/// unit that will make up the strand.
	std::string const & residue_name( core::Size const repeat_index ) const {
		runtime_assert_string_msg(
			repeat_index > 0 && repeat_index <= residue_name_.size(),
			"Error in protocols::beta_barrel::MakeBarrelStrand::residue_name(): The index is out of range."
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
	BarrelParametrizationCalculatorOP calculator_op() { return calculator_; }

	/// @brief Const access to the calculator.
	BarrelParametrizationCalculatorCOP calculator_cop() const { return calculator_; }

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

	/// @brief Length of the strand, in residues.  Defaults to 7.
	///
	core::Size strand_length_;

	/// @brief Name (full-length, not 3-letter code) of the residue type(s) that the strand will be made of.
	/// @details Defaults to "ALA".  One residue must be specified for each in the repeating unit.
	utility::vector1< std::string > residue_name_;

	/// @brief Did the last apply fail?
	/// @details Initialized to "false"; "true" if the last apply failed.
	bool last_apply_failed_;

	/// @brief The BarrelParametrizationCalculator object, which keeps track of parameter values and
	/// does the actual parametric math.
	BarrelParametrizationCalculatorOP calculator_;

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Set whether the last call to the apply() function failed or not.
	/// @details Should only be called by the apply() function.
	void set_last_apply_failed( bool const val ) { last_apply_failed_=val; return; }

	/// @brief Add Crick parameter information to the Conformation object within the pose.
	/// @details This function updates the barrel_parameters object's links to residues within the pose,
	/// and then copies the owning pointer into the pose's Conformation object.
	void add_parameter_info_to_pose( core::pose::Pose &pose );

};

} //namespace beta_barrel
} //namespace protocols

#endif //INCLUDED_protocols_beta_barrel_MakeBarrelStrand_hh
