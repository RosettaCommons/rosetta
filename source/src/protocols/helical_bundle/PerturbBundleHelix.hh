// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/PerturbBundleHelix.hh
/// @brief  Headers for PerturbBundleHelix.cc.  Reads the Crick parameters for a piece of a pose from
/// the input pose and sets the mainchain torsions accordingly.  The parameters are presumed to
/// have been perturbed by another mover.  This mover is intended to be called by the PerturbBundle
/// mover, which handles the perturbation of the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_PerturbBundleHelix_hh
#define INCLUDED_protocols_helical_bundle_PerturbBundleHelix_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/PerturbBundleHelix.fwd.hh>
#include <protocols/helical_bundle/util.hh>
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

class PerturbBundleHelix : public protocols::moves::Mover
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
	PerturbBundleHelix();
	~PerturbBundleHelix() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;


	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;


	/*virtual void parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const &
	);*/

public: //Setters:

	/// @brief Sets the index of the ParametersSet object in the pose that's being perturbed.
	/// @details This ParametersSet contains the Parameters object that defines the helix that's being perturbed.
	void set_parameters_set_index( core::Size const index ) { parameters_set_index_=index; return; }

	/// @brief Sets the index of the Parameters object in the ParametersSet object in the pose that's being perturbed.
	/// @details This Parameters object contains the Crick parameters that define the helix that's being perturbed.
	void set_parameters_index( core::Size const index ) { parameters_index_=index; return; }


public: //Getters:

	/// @brief Returns the index of the ParametersSet object in the pose that's being perturbed.
	/// @details This ParametersSet contains the Parameters object that defines the helix that's being perturbed.
	core::Size parameters_set_index() const { return parameters_set_index_; }

	/// @brief Returns the index of the Parameters object in the ParametersSet object in the pose that's being perturbed.
	/// @details This Parameters object contains the Crick parameters that define the helix that's being perturbed.
	core::Size parameters_index() const { return parameters_index_; }

	bool last_apply_failed() const { return last_apply_failed_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE DATA                                                      //
	////////////////////////////////////////////////////////////////////////////////

private:

	/// @brief The index of the ParametersSet object in the pose that's being perturbed.
	/// @details This ParametersSet contains the Parameters object that defines the helix that's being perturbed.
	core::Size parameters_set_index_;

	/// @brief The index of the Parameters object in the ParametersSet object in the pose that's being perturbed.
	/// @details This Parameters object contains the Crick parameters that define the helix that's being perturbed.
	core::Size parameters_index_;

	/// @brief Stores whether the last apply attempt failed.
	///
	bool last_apply_failed_;

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

private:

	void set_last_apply_failed(bool const status) { last_apply_failed_=status; return;}

};

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_PerturbBundleHelix_hh
