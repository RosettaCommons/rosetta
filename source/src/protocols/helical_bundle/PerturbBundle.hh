// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/PerturbBundle.hh
/// @brief  Headers for PerturbBundle.cc.  Perturbs a helical bundle by altering the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.  This mover calls the PerturbBundleHelix mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_PerturbBundle_hh
#define INCLUDED_protocols_helical_bundle_PerturbBundle_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/PerturbBundle.fwd.hh>
#include <protocols/helical_bundle/PerturbBundleHelix.fwd.hh>
#include <protocols/helical_bundle/PerturbBundleHelix.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <protocols/helical_bundle/PerturbBundleOptions.fwd.hh>

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

class PerturbBundle : public protocols::moves::Mover
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
	PerturbBundle();
	PerturbBundle( PerturbBundle const &src );
	virtual ~PerturbBundle();

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;


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

public:
////////////////////////////////////////////////////////////////////////////////
//          PUBLIC FUNCTIONS                                                  //
////////////////////////////////////////////////////////////////////////////////

	/// @brief Set which bundle parameters set will be used, if more than one is defined in the pose's Conformation object.
	/// @details  A value of n indicates that the nth bundle paramets set encountered will be perturbed.
	void set_bundleparametersset_index(core::Size const val) { bundleparametersset_index_=val; return; }

	/// @brief Returns a value indicating which bundle parameters set will be used, if more than one is defined in the pose's Conformation object.
	/// @details  A value of n indicates that the nth bundle paramets set encountered will be perturbed.
	core::Size bundleparametersset_index() const { return bundleparametersset_index_; }

	/// @brief Access the r0_ BundleOptions object, by index.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsOP r0( core::Size const index ) {
		runtime_assert_string_msg(index>0 && index<=r0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::r0() out of range.");
		return r0_[index];
	}

	/// @brief Access the r0_ BundleOptions object, by index.  This provides const access.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsCOP r0( core::Size const index ) const {
		runtime_assert_string_msg(index>0 && index<=r0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::r0() out of range.");
		return r0_[index];
	}

	/// @brief Access the omega0_ BundleOptions object, by index.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsOP omega0( core::Size const index ) {
		runtime_assert_string_msg(index>0 && index<=omega0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::omega0() out of range.");
		return omega0_[index];
	}

	/// @brief Access the omega0_ BundleOptions object, by index.  This provides const access.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsCOP omega0( core::Size const index ) const {
		runtime_assert_string_msg(index>0 && index<=omega0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::omega0() out of range.");
		return omega0_[index];
	}

	/// @brief Access the delta_omega0_ BundleOptions object, by index.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsOP delta_omega0( core::Size const index ) {
		runtime_assert_string_msg(index>0 && index<=delta_omega0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_omega0() out of range.");
		return delta_omega0_[index];
	}

	/// @brief Access the delta_omega0_ BundleOptions object, by index.  This provides const access.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsCOP delta_omega0( core::Size const index ) const {
		runtime_assert_string_msg(index>0 && index<=delta_omega0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_omega0() out of range.");
		return delta_omega0_[index];
	}

	/// @brief Access the delta_omega1_ BundleOptions object, by index.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsOP delta_omega1( core::Size const index ) {
		runtime_assert_string_msg(index>0 && index<=delta_omega1_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_omega1() out of range.");
		return delta_omega1_[index];
	}

	/// @brief Access the delta_omega1_ BundleOptions object, by index.  This provides const access.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsCOP delta_omega1( core::Size const index ) const {
		runtime_assert_string_msg(index>0 && index<=delta_omega1_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_omega1() out of range.");
		return delta_omega1_[index];
	}

	/// @brief Access the delta_t_ BundleOptions object, by index.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsOP delta_t( core::Size const index ) {
		runtime_assert_string_msg(index>0 && index<=delta_t_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_t() out of range.");
		return delta_t_[index];
	}

	/// @brief Access the delta_t_ BundleOptions object, by index.  This provides const access.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsCOP delta_t( core::Size const index ) const {
		runtime_assert_string_msg(index>0 && index<=delta_t_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_t() out of range.");
		return delta_t_[index];
	}
	
	/// @brief Access the z1_offset_ BundleOptions object, by index.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsOP z1_offset( core::Size const index ) {
		runtime_assert_string_msg(index>0 && index<=z1_offset_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::z1_offset() out of range.");
		return z1_offset_[index];
	}

	/// @brief Access the z1_offset_ BundleOptions object, by index.  This provides const access.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsCOP z1_offset( core::Size const index ) const {
		runtime_assert_string_msg(index>0 && index<=z1_offset_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::z1_offset() out of range.");
		return z1_offset_[index];
	}
	
	/// @brief Access the z0_offset_ BundleOptions object, by index.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsOP z0_offset( core::Size const index ) {
		runtime_assert_string_msg(index>0 && index<=z0_offset_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::z0_offset() out of range.");
		return z0_offset_[index];
	}

	/// @brief Access the z0_offset_ BundleOptions object, by index.  This provides const access.
	/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
	PerturbBundleOptionsCOP z0_offset( core::Size const index ) const {
		runtime_assert_string_msg(index>0 && index<=z0_offset_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::z0_offset() out of range.");
		return z0_offset_[index];
	}

	/// @brief Access the default_r0_ BundleOptions object.
	///
	PerturbBundleOptionsOP default_r0() { return default_r0_; }

	/// @brief Access the default_r0_ BundleOptions object (const-access).
	///
	PerturbBundleOptionsCOP default_r0() const { return default_r0_; }

	/// @brief Access the default_omega0_ BundleOptions object.
	///
	PerturbBundleOptionsOP default_omega0() { return default_omega0_; }

	/// @brief Access the default_omega0_ BundleOptions object (const-access).
	///
	PerturbBundleOptionsCOP default_omega0() const { return default_omega0_; }

	/// @brief Access the default_delta_omega0_ BundleOptions object.
	///
	PerturbBundleOptionsOP default_delta_omega0() { return default_delta_omega0_; }

	/// @brief Access the default_delta_omega0_ BundleOptions object (const-access).
	///
	PerturbBundleOptionsCOP default_delta_omega0() const { return default_delta_omega0_; }

	/// @brief Access the default_delta_omega1_ BundleOptions object.
	///
	PerturbBundleOptionsOP default_delta_omega1() { return default_delta_omega1_; }

	/// @brief Access the default_delta_omega1_ BundleOptions object (const-access).
	///
	PerturbBundleOptionsCOP default_delta_omega1() const { return default_delta_omega1_; }

	/// @brief Access the default_delta_t_ BundleOptions object.
	///
	PerturbBundleOptionsOP default_delta_t() { return default_delta_t_; }

	/// @brief Access the default_delta_t_ BundleOptions object (const-access).
	///
	PerturbBundleOptionsCOP default_delta_t() const { return default_delta_t_; }

	/// @brief Access the default_z1_offset_ BundleOptions object.
	///
	PerturbBundleOptionsOP default_z1_offset() { return default_z1_offset_; }

	/// @brief Access the default_z1_offset_ BundleOptions object (const-access).
	///
	PerturbBundleOptionsCOP default_z1_offset() const { return default_z1_offset_; }

	/// @brief Access the default_z0_offset_ BundleOptions object.
	///
	PerturbBundleOptionsOP default_z0_offset() { return default_z0_offset_; }

	/// @brief Access the default_z0_offset_ BundleOptions object (const-access).
	///
	PerturbBundleOptionsCOP default_z0_offset() const { return default_z0_offset_; }

	/// @brief Add options for a new helix
	/// @details Return value is the current total number of helices after the addition.
	core::Size add_helix( core::Size const helix_index );


private:
////////////////////////////////////////////////////////////////////////////////
//          PRIVATE DATA                                                      //
////////////////////////////////////////////////////////////////////////////////

	/// @brief Default options for perturbing r0.
	/// @details May be overridden on a helix-by-helix basis.
	PerturbBundleOptionsOP default_r0_;

	/// @brief Helix-by-helix options for perturbing r0.
	///
	PerturbBundleOptionsOPs r0_;

	/// @brief Default options for perturbing omega0.
	/// @details May be overridden on a helix-by-helix basis.
	PerturbBundleOptionsOP default_omega0_;

	/// @brief Helix-by-helix options for perturbing omega0.
	///
	PerturbBundleOptionsOPs omega0_;

	/// @brief Default options for perturbing delta_omega0.
	/// @details May be overridden on a helix-by-helix basis.
	PerturbBundleOptionsOP default_delta_omega0_;

	/// @brief Helix-by-helix options for perturbing delta_omega0.
	///
	PerturbBundleOptionsOPs delta_omega0_;

	/// @brief Default options for perturbing delta_omega1.
	/// @details May be overridden on a helix-by-helix basis.
	PerturbBundleOptionsOP default_delta_omega1_;

	/// @brief Helix-by-helix options for perturbing delta_omega1.
	///
	PerturbBundleOptionsOPs delta_omega1_;

	/// @brief Default options for perturbing delta_t.
	/// @details May be overridden on a helix-by-helix basis.
	PerturbBundleOptionsOP default_delta_t_;

	/// @brief Helix-by-helix options for perturbing delta_t.
	///
	PerturbBundleOptionsOPs delta_t_;
	
	/// @brief Default options for perturbing z1_offset.
	/// @details May be overridden on a helix-by-helix basis.
	PerturbBundleOptionsOP default_z1_offset_;
	
	/// @brief Helix-by-helix options for perturbing z1_offset.
	///
	PerturbBundleOptionsOPs z1_offset_;

	/// @brief Default options for perturbing z0_offset.
	/// @details May be overridden on a helix-by-helix basis.
	PerturbBundleOptionsOP default_z0_offset_;

	/// @brief Helix-by-helix options for perturbing z0_offset.
	///
	PerturbBundleOptionsOPs z0_offset_;

	/// @brief Which set of bundle parameters (if there exists more than one) should the mover alter?
	/// @details Defaults to 1.  Higher values indicate the nth set encountered in the ParametersSet list.
	core::Size bundleparametersset_index_;
	
private:
////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

	/// @brief Is a value in a list?
	///
	bool is_in_list( core::Size const val, utility::vector1 < core::Size> const &list ) const;

	/// @brief Confirms that a helix has not yet been defined.  Returns "true" if the helix
	/// has NOT been defined, false otherwise.
	bool helix_not_defined( core::Size const helix_index) const;

	/// @brief Perturb the helical bundle parameter values in the pose, subject to the options already set.
	/// @details Called by the apply() function.  Returns true for success, false for failure.
	bool perturb_values( BundleParametersSetOP params_set) const;

	/// @brief Rebuild the helical bundle conformation using the bundle parameter values in the pose.
	///
	void rebuild_conformation(
		core::pose::Pose &pose,
		BundleParametersSetOP params_set,
		core::Size const params_set_index,
		bool &failed
	) const;

	/// @brief Write out the perturbed Crick parameters.
	/// @details The "before" parameter determines whether this is a pre-perturbation
	/// report or a post-perturbation report.
	void write_report(
		BundleParametersSetOP params_set,
		bool const before
	) const;


}; //PerturbBundle class

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_PerturbBundle_hh
