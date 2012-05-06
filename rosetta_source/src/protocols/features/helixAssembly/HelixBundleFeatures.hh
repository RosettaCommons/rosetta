// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file HelixBundleFeatures.hh
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_helixAssembly_HelixBundleFeatures_hh
#define INCLUDED_protocols_features_helixAssembly_HelixBundleFeatures_hh

//Unit
#include <protocols/features/helixAssembly/HelixBundleFeatures.fwd.hh>

//External
#include <boost/uuid/uuid.hpp>

//Protocols
#include <protocols/features/FeaturesReporter.hh>

//Devel
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility
#include <utility/vector1.hh>

namespace protocols {
namespace features {
namespace helixAssembly {

class HelixBundleFeatures : public protocols::features::FeaturesReporter {

public:

	HelixBundleFeatures();

	void init_from_options();

	virtual
	std::string
	type_name() const  {
		return "HelixBundleFeatures";
	}

	///@brief return sql statements that sets up the appropriate tables
	///to contain the features.
	virtual
	std::string
	schema() const;

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	///@brief collect all the feature data for the pose
	virtual
	core::Size
	report_features(
					core::pose::Pose const & pose,
					utility::vector1<bool> const & relevant_residues,
					boost::uuids::uuid struct_id,
					utility::sql_database::sessionOP db_session
					);

	utility::vector1<HelicalFragment> get_full_helices(boost::uuids::uuid struct_id, utility::sql_database::sessionOP db_session);

	bool checkHelixContacts(core::pose::Pose const & pose, HelicalFragment helix_1, HelicalFragment helix_2, HelicalFragment helix_3);

private:

	core::Real helix_cap_dist_cutoff_;
	core::Real helix_contact_dist_cutoff_;
	core::Size min_helix_size_;

};

} //namespace helixAssembly
} //namespace features
} //namespace protocols

#endif
