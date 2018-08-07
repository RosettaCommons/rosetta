// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/MakeBundle.hh
/// @brief Headers for MakeBundle.cc.  Builds a helical bundle using the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.  This mover calls the MakeBundleHelix mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_MakeBundle_hh
#define INCLUDED_protocols_helical_bundle_MakeBundle_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/MakeBundle.fwd.hh>
#include <protocols/helical_bundle/MakeBundleHelix.fwd.hh>
#include <protocols/helical_bundle/MakeBundleHelix.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.fwd.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>

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
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace helical_bundle {

class MakeBundle : public protocols::moves::Mover
{
public: //Typedefs

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
	MakeBundle();
	MakeBundle( MakeBundle const &src );
	~MakeBundle() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;


	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;


	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const &
	) override;

	/// @brief Access the default calculator directly.
	/// @details Nonconst access is potentially dangerous!  Use with caution!
	/// @note Used by BundleGridSampler to access the contained MakeBundle mover's calculator directly.
	inline BundleParametrizationCalculatorOP default_calculator_nonconst() { return default_calculator_; }

	/// @brief Access the default calculator directly.
	/// @details Const access.
	inline BundleParametrizationCalculatorCOP default_calculator() const { return utility::pointer::static_pointer_cast< BundleParametrizationCalculator const >( default_calculator_ ); }

	/// @brief Set symmetry and symmetry_copies options based on an input tag.
	///
	void set_symmetry_options_from_tag( utility::tag::TagCOP tag );

	/// @brief Set whether the input pose should be reset prior to building a helix.
	///
	void set_reset_pose( bool const reset_in) { reset_pose_ = reset_in; }

	/// @brief Returns a bool indicating whether the input pose will be reset prior to
	/// building a helix.
	bool reset_pose() const { return reset_pose_; }

	/// @brief Function to add a helix.
	/// @details This creates a MakeBundleHelix mover that will be called at apply time.  Note that this function assumes that default parameter values have been set already.  They
	/// cannot be set after calling this function.
	void add_helix();

	/// @brief Non-const access to helix in the bundle.
	///
	protocols::helical_bundle::MakeBundleHelixOP helix( core::Size const &index ) {
		runtime_assert_string_msg( index>0 && index<=make_bundle_helix_movers_.size(), "In MakeBundle::helix(): index is out of range (i.e. it doesn't refer to an already-defined helix)." );
		return make_bundle_helix_movers_[index];
	}

	/// @brief Const access to helix in the bundle.
	///
	protocols::helical_bundle::MakeBundleHelixCOP helix_cop( core::Size const &index ) const {
		runtime_assert_string_msg( index>0 && index<=make_bundle_helix_movers_.size(), "In MakeBundle::helix_cop(): index is out of range (i.e. it doesn't refer to an already-defined helix)." );
		return make_bundle_helix_movers_[index];
	}

	/// @brief Returns the number of helices.
	/// @details Note that this is not multiplied by the number of symmetry repeats.  That is,
	/// if the bundle has 3-fold symmetry and 2 helices are defined, there will be 6 helices in
	/// the final structure, but this function would still return 2.
	core::Size n_helices() const { return make_bundle_helix_movers_.size(); }

	/// @brief Set the symmetry of the bundle.
	/// @details See bundle_symmetry_ member variable for details.
	void set_symmetry( core::Size const symmetry_in ) { bundle_symmetry_ = symmetry_in; return; }

	/// @brief Get the symmetry of the bundle.
	/// @details See bundle_symmetry_ member variable for details.
	core::Size symmetry() const { return bundle_symmetry_; }

	/// @brief Set how many of the symmetry copies actually get generated.
	/// @details See bundle_symmetry_copies_ member variable for details.
	void set_symmetry_copies( core::Size const symmetry_copies_in ) { bundle_symmetry_copies_ = symmetry_copies_in; return; }

	/// @brief Get how many of the symmetry copies actually get generated.
	/// @details See bundle_symmetry_copies_ member variable for details.
	core::Size symmetry_copies() const { return bundle_symmetry_copies_; }

	/// @brief Set the default Crick params file name.
	/// @details Triggers a read from disk!
	void set_default_crick_params_file(std::string const &input_string);

	/// @brief Initialize the default calculator from the deafult Crick params file.
	/// @details Triggers a read from disk!
	void initialize_default_calculator_from_default_crick_params_file();

	/// @brief Set the Crick params file for a particular helix.
	/// @details Triggers a read from disk!
	void set_crick_params_file_for_helix( std::string const &filename, core::Size const helix_index );

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

	/// @brief Set the residue names for a given helix.
	void set_residue_name_for_helix( utility::vector1< std::string > const & names, core::Size const helix_index );

	/// @brief Set the default number of residues per helix
	///
	void set_default_helix_length(core::Size const &val);

	/// @brief Returns the default number of residues per helix.
	///
	core::Size default_helix_length() const;

	/// @brief Set the helix length, in residues, for the Nth helix.
	void set_helix_length_for_helix( core::Size const helix_length, core::Size const helix_index );

	/// @brief Set whether we're using degrees (true) or radians (false).
	///
	void set_use_degrees( bool const val=true );

	/// @brief Get whether we're using degrees (true) or radians (false).
	///
	bool use_degrees() const { return use_degrees_; }

	/// @brief Returns true or false based on whether the last call to the apply() function
	/// failed (true) or succeeded (false).
	/// @details The apply() function calls the private function set_last_apply_failed() to
	/// set this.
	bool last_apply_failed() const { return last_apply_failed_; }

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

	/// @brief The BundleParameterizationCalculator object that will be used for setting up default parameters.
	/// @details This object gets cloned by each of the make_bundle_helix_movers_, and they use it as the intial (default) state before
	/// setting up their own parameters.
	BundleParametrizationCalculatorOP default_calculator_;

	/// @brief A vector of owning pointers to the individual MakeBundleHelix movers that will make each of
	/// the helices in the bundle
	utility::vector1 < protocols::helical_bundle::MakeBundleHelixOP > make_bundle_helix_movers_;

	/// @brief The symmetry of the bundle.
	/// @details "0" or "1" are asymmetric.  "2" indicates 2-fold rotational symmetry, "3" indicates
	/// 3-fold rotational symmetry, and so forth.  Note that this need not correspond to the number of
	/// helices defined.  If 2 helices are defined with 3-fold symmetry, you end up with 6 helices.
	core::Size bundle_symmetry_;

	/// @brief The symmetry copies to generate.
	/// @details A value of 0 means to generate all copies.  Higher values mean to generate only the first N
	/// copies.  For example, if the symmetry were 16 but bundle_symmetry_copies_ were set to 4, only the
	/// first 4 symmetry repeats would be generated.
	core::Size bundle_symmetry_copies_;

	/// DEFAULTS: Default values for helix parameters if no other values are set:

	/// @brief Default Crick params file name
	///
	std::string default_crick_params_file_;

	/// @brief Default residue name
	///
	utility::vector1< std::string > default_residue_name_;

	/// @brief Default number of residues in the helix.
	///
	core::Size default_helix_length_;

	/// @brief Has the default number of residues in the helix been specified?
	///
	bool default_helix_length_set_;

	/// @brief Are we using degrees (true) or radians (false)?  Default radians (false).
	///
	bool use_degrees_;

	/// @brief Did the last apply fail?
	///
	bool last_apply_failed_;

	/// @brief Have the default parameter values been set?
	/// @details If true, they cannot be set again.  Note that once a helix is added, defaults_set_ becomes true.  This ensures a sensible order of
	/// operations when calling the MakeBundle mover from code, where configuration is not necessarily happening in the order specified in the
	/// parse_my_tag() function.
	bool defaults_set_;

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Set whether the last apply failed.
	/// @details Called by the apply() function.
	void set_last_apply_failed (bool const val) { last_apply_failed_=val; return; }

};

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_MakeBundle_hh
