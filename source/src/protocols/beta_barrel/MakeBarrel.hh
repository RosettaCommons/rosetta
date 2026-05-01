// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/MakeBarrel.hh
/// @brief Headers for MakeBarrel.cc.  Builds a beta-barrel using the Crick parameters.
/// @details The barrel is centred on the origin, with the barrel axis pointing along the
/// global z-axis.  This mover calls the MakeBarrelStrand mover.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_MakeBarrel_hh
#define INCLUDED_protocols_beta_barrel_MakeBarrel_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/beta_barrel/MakeBarrel.fwd.hh>
#include <protocols/beta_barrel/MakeBarrelStrand.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.hh>
#include <protocols/beta_barrel/parameters/BarrelParametersSet.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParametersSet.hh>
#include <protocols/beta_barrel/BarrelParametrizationCalculator.fwd.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Project Headers
#include <utility/vector1.hh>



///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace beta_barrel {

class MakeBarrel : public protocols::moves::Mover
{
public: //Typedefs

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
	MakeBarrel();
	MakeBarrel( MakeBarrel const &src );
	~MakeBarrel() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	/// @brief Access the default calculator directly.
	/// @details Nonconst access is potentially dangerous!  Use with caution!
	inline BarrelParametrizationCalculatorOP default_calculator_nonconst() { return default_calculator_; }

	/// @brief Access the default calculator directly.
	/// @details Const access.
	inline BarrelParametrizationCalculatorCOP default_calculator() const { return utility::pointer::static_pointer_cast< BarrelParametrizationCalculator const >( default_calculator_ ); }

	/// @brief Set whether the input pose should be reset prior to building a barrel.
	///
	void set_reset_pose( bool const reset_in) { reset_pose_ = reset_in; }

	/// @brief Returns a bool indicating whether the input pose will be reset prior to
	/// building a barrel.
	bool reset_pose() const { return reset_pose_; }

	/// @brief Function to add a strand.
	/// @details This creates a MakeBarrelStrand mover that will be called at apply time.
	/// Note that this function assumes that default parameter values have been set already.
	void add_strand();

	/// @brief Non-const access to strand in the barrel.
	///
	protocols::beta_barrel::MakeBarrelStrandOP strand( core::Size const &index ) {
		runtime_assert_string_msg( index>0 && index<=make_barrel_strand_movers_.size(), "In MakeBarrel::strand(): index is out of range (i.e. it doesn't refer to an already-defined strand)." );
		return make_barrel_strand_movers_[index];
	}

	/// @brief Const access to strand in the barrel.
	///
	protocols::beta_barrel::MakeBarrelStrandCOP strand_cop( core::Size const &index ) const {
		runtime_assert_string_msg( index>0 && index<=make_barrel_strand_movers_.size(), "In MakeBarrel::strand_cop(): index is out of range (i.e. it doesn't refer to an already-defined strand)." );
		return make_barrel_strand_movers_[index];
	}

	/// @brief Returns the number of strands.
	///
	core::Size n_strands() const { return n_strands_; }

	/// @brief Set the number of strands.
	///
	void set_n_strands( core::Size const val ) { n_strands_ = val; }

	/// @brief Returns the shear number.
	///
	core::Size shear_number() const { return shear_number_; }

	/// @brief Set the shear number.
	///
	void set_shear_number( core::Size const val ) { shear_number_ = val; }

	/// @brief Returns whether the barrel is antiparallel.
	///
	bool antiparallel() const { return antiparallel_; }

	/// @brief Set antiparallel flag.
	///
	void set_antiparallel( bool const val ) { antiparallel_ = val; }

	/// @brief Set the default Crick params file name.
	/// @details Triggers a read from disk!
	void set_default_crick_params_file(std::string const &input_string);

	/// @brief Initialize the default calculator from the default Crick params file.
	/// @details Triggers a read from disk!
	void initialize_default_calculator_from_default_crick_params_file();

	/// @brief Set the Crick params file for a particular strand.
	/// @details Triggers a read from disk!
	void set_crick_params_file_for_strand( std::string const &filename, core::Size const strand_index );

	/// @brief Returns the default Crick params file name.
	///
	std::string const & default_crick_params_file() const;

	/// @brief Set the default residue name
	///
	void set_default_residue_name(utility::vector1< std::string > const &names);

	/// @brief Returns the default residue name for a particular position in the repeating unit.
	///
	std::string const & default_residue_name( core::Size const index_in_repeating_unit ) const;

	/// @brief Returns the default residue name vector.
	///
	utility::vector1 < std::string > const & default_residue_name() const;

	/// @brief Set the residue names for a given strand.
	void set_residue_name_for_strand( utility::vector1< std::string > const & names, core::Size const strand_index );

	/// @brief Set the default number of residues per strand
	///
	void set_default_strand_length(core::Size const &val);

	/// @brief Returns the default number of residues per strand.
	///
	core::Size default_strand_length() const;

	/// @brief Set the strand length, in residues, for the Nth strand.
	void set_strand_length_for_strand( core::Size const strand_length, core::Size const strand_index );

	/// @brief Set whether we're using degrees (true) or radians (false).
	///
	void set_use_degrees( bool const val=true );

	/// @brief Get whether we're using degrees (true) or radians (false).
	///
	bool use_degrees() const { return use_degrees_; }

	/// @brief Returns true or false based on whether the last call to the apply() function
	/// failed (true) or succeeded (false).
	bool last_apply_failed() const { return last_apply_failed_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //Function overrides needed for citation manager

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

private:

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE DATA                                                      //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Should the input pose be reset?
	/// @details Default true.
	bool reset_pose_;

	/// @brief The BarrelParameterizationCalculator object that will be used for setting up default parameters.
	/// @details This object gets cloned by each of the make_barrel_strand_movers_, and they use it as the initial (default) state before
	/// setting up their own parameters.
	BarrelParametrizationCalculatorOP default_calculator_;

	/// @brief A vector of owning pointers to the individual MakeBarrelStrand movers that will make each of
	/// the strands in the barrel.
	utility::vector1 < protocols::beta_barrel::MakeBarrelStrandOP > make_barrel_strand_movers_;

	/// @brief The number of strands in the barrel.
	///
	core::Size n_strands_;

	/// @brief The shear number of the barrel.
	///
	core::Size shear_number_;

	/// @brief Whether the barrel has antiparallel strand topology.
	///
	bool antiparallel_;

	/// DEFAULTS: Default values for strand parameters if no other values are set:

	/// @brief Default Crick params file name
	///
	std::string default_crick_params_file_;

	/// @brief Default residue name
	///
	utility::vector1< std::string > default_residue_name_;

	/// @brief Default number of residues in the strand.
	///
	core::Size default_strand_length_;

	/// @brief Are we using degrees (true) or radians (false)?  Default radians (false).
	///
	bool use_degrees_;

	/// @brief Did the last apply fail?
	///
	bool last_apply_failed_;

	/// @brief Have the default parameter values been set?
	/// @details If true, they cannot be set again.
	bool defaults_set_;

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Set whether the last apply failed.
	/// @details Called by the apply() function.
	void set_last_apply_failed (bool const val) { last_apply_failed_=val; return; }

};

} //namespace beta_barrel
} //namespace protocols

#endif //INCLUDED_protocols_beta_barrel_MakeBarrel_hh
