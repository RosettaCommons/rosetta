// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/MakeBundle.hh
/// @brief  Headers for MakeBundle.cc.  Builds a helical bundle using the Crick parameters.
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
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/grid/CartGrid.hh>
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
		typedef protocols::helical_bundle::parameters::BundleParametersSet BundleParametersSet;
		typedef protocols::helical_bundle::parameters::BundleParametersSetOP BundleParametersSetOP;

public:
	MakeBundle();
	virtual ~MakeBundle();

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	///
	/// @brief Actually apply the mover to the pose.
	virtual void apply(core::pose::Pose & pose);

	virtual std::string get_name() const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const &
	);

	/// @brief Set whether the input pose should be reset prior to building a helix.
	///
	void set_reset_pose( bool const reset_in) { reset_pose_ = reset_in; }

	/// @brief Returns a bool indicating whether the input pose will be reset prior to
	/// building a helix.
	bool reset_pose() const { return reset_pose_; }

	/// @brief Function to add a helix.
	/// @details This creates a MakeBundleHelix mover that will be called at apply time.
	/// The new mover is only initialized if default values are provided by THIS mover.
	void add_helix();

	/// @brief Non-const access to helix in the bundle.
	///
	protocols::helical_bundle::MakeBundleHelixOP helix( core::Size const &index ) {
		runtime_assert_string_msg( index>0 && index<=make_bundle_helix_movers_.size(), "In MakeBundle::helix(): index is out of range (i.e. it doesn't refer to an already-defined helix)." );
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

	/// @brief Set the default r0 value (major helix radius)
	///
	void set_default_r0(core::Real const &val) { default_r0_=val; default_r0_set_=true; }

	/// @brief Returns true if and only if the default r0 value (major helix radius) has been set
	///
	bool default_r0_set() const { return default_r0_set_; }

	/// @brief Returns the default r0 value (major helix radius).
	///
	core::Real default_r0() const;

	/// @brief Set the default omega0 value (major helix rise per residue)
	///
	void set_default_omega0(core::Real const &val) { default_omega0_=val; default_omega0_set_=true; }

	/// @brief Returns true if and only if the default omega0 value (major helix rise per residue) has been set
	///
	bool default_omega0_set() const { return default_omega0_set_; }

	/// @brief Returns the default omega0 value (major helix rise per residue).
	///
	core::Real default_omega0() const;

	/// @brief Set the default delta_omega0 value (major helix turn)
	///
	void set_default_delta_omega0(core::Real const &val) { default_delta_omega0_=val; default_delta_omega0_set_=true; }

	/// @brief Returns true if and only if the default delta_omega0 value (major helix turn) has been set
	///
	bool default_delta_omega0_set() const { return default_delta_omega0_set_; }

	/// @brief Returns the default delta_omega0 value (major helix turn).
	///
	core::Real default_delta_omega0() const;


	/// @brief Set the default Crick params file name.
	///
	void set_default_crick_params_file(std::string const &input_string) { default_crick_params_file_=input_string; default_crick_params_file_set_=true; }

	/// @brief Returns true if and only if the default Crick params file name has been set.
	///
	bool default_crick_params_file_set() const { return default_crick_params_file_set_; }

	/// @brief Returns the default Crick params file name.
	///
	std::string default_crick_params_file() const;

	/// @brief Set the default omega1 value (minor helix turn per residue)
	///
	void set_default_omega1(core::Real const &val) { default_omega1_=val; default_omega1_set_=true; }

	/// @brief Returns true if and only if the default omega1 value (minor helix turn per residue) has been set
	///
	bool default_omega1_set() const { return default_omega1_set_; }

	/// @brief Returns the default omega1 value (minor helix turn per residue).
	///
	core::Real default_omega1() const;

	/// @brief Set the default z1 value (minor helix rise per residue)
	///
	void set_default_z1(core::Real const &val) { default_z1_=val; default_z1_set_=true; }

	/// @brief Returns true if and only if the default z1 value (minor helix rise per residue) has been set
	///
	bool default_z1_set() const { return default_z1_set_; }

	/// @brief Returns the default z1 value (minor helix rise per residue).
	///
	core::Real default_z1() const;

	/// @brief Set the default delta_omega1_all value (minor helix rotation)
	///
	void set_default_delta_omega1_all(core::Real const &val) { default_delta_omega1_all_=val; default_delta_omega1_all_set_=true; }

	/// @brief Returns true if and only if the default delta_omega1_all value (minor helix rotation) has been set
	///
	bool default_delta_omega1_all_set() const { return default_delta_omega1_all_set_; }

	/// @brief Returns the default delta_omega1_all value (minor helix rotation).
	///
	core::Real default_delta_omega1_all() const;

	/// @brief Set the default residue name
	///
	void set_default_residue_name(std::string const &name) { default_residue_name_=name; default_residue_name_set_=true; }

	/// @brief Returns true if and only if the default residue name has been set
	///
	bool default_residue_name_set() const { return default_residue_name_set_; }

	/// @brief Returns the default residue name.
	///
	std::string default_residue_name() const;

	/// @brief Set the default delta_t value (residue offset)
	///
	void set_default_delta_t(core::Real const &val) { default_delta_t_=val; default_delta_t_set_=true; }

	/// @brief Returns true if and only if the default delta_t value (residue offset) has been set
	///
	bool default_delta_t_set() const { return default_delta_t_set_; }

	/// @brief Returns the default delta_t value (residue offset).
	///
	core::Real default_delta_t() const;

	/// @brief Set the default invert value (should the helix be flipped?)
	///
	void set_default_invert(bool const &val) { default_invert_=val; default_invert_set_=true; }

	/// @brief Returns true if and only if the default invert value (should the helix be flipped?) has been set
	///
	bool default_invert_set() const { return default_invert_set_; }

	/// @brief Returns the default invert value (should the helix be flipped?)
	///
	bool default_invert() const;

	/// @brief Set the default number of residues per helix
	///
	void set_default_helix_length(core::Size const &val) { default_helix_length_=val; default_helix_length_set_=true; }

	/// @brief Returns true if and only if the default number of residues per helix has been set
	///
	bool default_helix_length_set() const { return default_helix_length_set_; }

	/// @brief Returns the default number of residues per helix.
	///
	core::Size default_helix_length() const;

	/// @brief Set the default for whether bond lengths should be set by the mover
	///
	void set_default_allow_bondlengths(bool const &val) { default_allow_bondlengths_=val; default_allow_bondlengths_set_=true; }

	/// @brief Returns true if and only if the default for whether bond lengths should be set by the mover has been set
	///
	bool default_allow_bondlengths_set() const { return default_allow_bondlengths_set_; }

	/// @brief Returns the default for whether bond lengths should be set by the mover
	///
	bool default_allow_bondlengths() const;

	/// @brief Set the default for whether bond lengths should be set by the mover
	///
	void set_default_allow_bondangles(bool const &val) { default_allow_bondangles_=val; default_allow_bondangles_set_=true; }

	/// @brief Returns true if and only if the default for whether bond lengths should be set by the mover has been set
	///
	bool default_allow_bondangles_set() const { return default_allow_bondangles_set_; }

	/// @brief Returns the default for whether bond lengths should be set by the mover
	///
	bool default_allow_bondangles() const;

	/// @brief Set the default for whether bond lengths should be set by the mover
	///
	void set_default_allow_dihedrals(bool const &val) { default_allow_dihedrals_=val; default_allow_dihedrals_set_=true; }

	/// @brief Returns true if and only if the default for whether bond lengths should be set by the mover has been set
	///
	bool default_allow_dihedrals_set() const { return default_allow_dihedrals_set_; }

	/// @brief Returns the default for whether bond lengths should be set by the mover
	///
	bool default_allow_dihedrals() const;

private:
////////////////////////////////////////////////////////////////////////////////
//          PRIVATE DATA                                                      //
////////////////////////////////////////////////////////////////////////////////

	/// @brief Should the input pose be reset?
	/// @details Default true.
	bool reset_pose_;

	/// @brief A vector of owning pointers to the individual MakeBundleHelix movers that will make each of
	/// the helices in the bundle
	utility::vector1 < protocols::helical_bundle::MakeBundleHelixOP > make_bundle_helix_movers_;

	/// @brief The symmetry of the bundle.
	/// @details "0" or "1" are asymmetric.  "2" indicates 2-fold rotational symmetry, "3" indicates
	/// 3-fold rotational symmetry, and so forth.  Note that this need not correspond to the number of
	/// helices defined.  If 2 helices are defined with 3-fold symmetry, you end up with 6 helices.
	core::Size bundle_symmetry_;


	/// DEFAULTS: Default values for helix parameters if no other values are set:

	/// @brief Default r0 value (major helix radius)
	///
	core::Real default_r0_;

	/// @brief Has the default r0 value (major helix radius) been specified?
	///
	bool default_r0_set_;

	/// @brief Default omega0 value (turn per residue of major helix)
	///
	core::Real default_omega0_;

	/// @brief Has the default omega0 value (turn per residue of major helix) been specified?
	///
	bool default_omega0_set_;

	/// @brief Default delta_omega0 value (rotation of major helix)
	///
	core::Real default_delta_omega0_;

	/// @brief Has the default delta_omega0 value (rotation of major helix) been specified?
	///
	bool default_delta_omega0_set_;

	/// @brief Default Crick params file name
	///
	std::string default_crick_params_file_;

	/// @brief Has the default Crick params file name been set?
	///
	bool default_crick_params_file_set_;

	/// @brief Default omega1 value (minor helix turn per residue)
	///
	core::Real default_omega1_;

	/// @brief Has the default omega1 value (minor helix turn per residue) been specified?
	///
	bool default_omega1_set_;

	/// @brief Default z1 value (minor helix rise per residue)
	///
	core::Real default_z1_;

	/// @brief Has the default z1 value (minor helix rise per residue) been specified?
	///
	bool default_z1_set_;

	/// @brief Default delta_omega1_all value (minor helix rotation)
	///
	core::Real default_delta_omega1_all_;

	/// @brief Has the default z1 value (minor helix rotation) been specified?
	///
	bool default_delta_omega1_all_set_;

	/// @brief Default residue name
	///
	std::string default_residue_name_;

	/// @brief Has the default residue name been specified?
	///
	bool default_residue_name_set_;

	/// @brief Default delta_t value (residue offset)
	///
	core::Real default_delta_t_;

	/// @brief Has the default delta_t value (residue offset) been specified?
	///
	bool default_delta_t_set_;

	/// @brief Default invert value (should the helix be flipped?)
	///
	bool default_invert_;

	/// @brief Has the default invert value been specified?
	///
	bool default_invert_set_;

	/// @brief Default number of residues in the helix.
	///
	core::Size default_helix_length_;

	/// @brief Has the default number of residues in the helix been specified?
	///
	bool default_helix_length_set_;

	/// @brief Default allow_bondlengths value (should the mover be allowed to set bond lengths?)
	///
	bool default_allow_bondlengths_;

	/// @brief Has the allow_bondlengths value been specified?
	///
	bool default_allow_bondlengths_set_;

	/// @brief Default allow_bondangles value (should the mover be allowed to set bond lengths?)
	///
	bool default_allow_bondangles_;

	/// @brief Has the allow_bondangles value been specified?
	///
	bool default_allow_bondangles_set_;

	/// @brief Default allow_dihedrals value (should the mover be allowed to set bond lengths?)
	///
	bool default_allow_dihedrals_;

	/// @brief Has the allow_dihedrals value been specified?
	///
	bool default_allow_dihedrals_set_;

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

};

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_MakeBundle_hh
