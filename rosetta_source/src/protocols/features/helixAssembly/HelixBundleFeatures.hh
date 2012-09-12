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
#include <protocols/simple_filters/InterfaceSasaFilter.hh>

//Devel
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility
#include <utility/vector1.hh>

namespace protocols {
namespace features {
namespace helixAssembly {

class HelixBundleFeatures : public protocols::features::FeaturesReporter
{

public:

	HelixBundleFeatures();

	void init_from_options();

	virtual
	std::string
	type_name() const  {
		return "HelixBundleFeatures";
	};

	///@brief generate the table schemas and write them to the database
	virtual
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session
	) const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose
	);

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;
	
	void
	generate_comparison_list(
		utility::vector1<HelicalFragmentOP> all_helix_fragments,
		core::Size prev_index,
		utility::vector1<HelicalFragmentOP> fragment_list
	);

	bool
	check_cap_distances(
		core::pose::Pose const & pose,
		utility::vector1<HelicalFragmentOP> frag_set
	);

	///@brief collect all the feature data for the pose
	virtual
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		boost::uuids::uuid struct_id,
		utility::sql_database::sessionOP db_session
	);

	utility::vector1<HelicalFragmentOP>
	get_helix_fragments(
		boost::uuids::uuid struct_id,
		utility::sql_database::sessionOP db_session
	);

	///@brief create a bundle-pose from the combination of fragments
	/// and record the "interface" SASA for each helix against the
	/// rest of the bundle
	void
	record_helix_sasas(
		core::pose::Pose const & pose,
		utility::vector1<HelicalFragmentOP> & bundle_fragments
	);

private:

	//Maximum allowable distance, in angstroms, between any two helix ends
	core::Real helix_cap_dist_cutoff_;
	
	//number of helices in the bundle
	core::Size bundle_size_;
	
	//number of residues in each helix
	core::Size helix_size_;
	
	//List of all sublists of size bundle_size from the total number of helix fragments
	utility::vector1< utility::vector1<HelicalFragmentOP> > comparison_list_;
	
	protocols::simple_filters::InterfaceSasaFilter sasa_filter_;

};

} //namespace helixAssembly
} //namespace features
} //namespace protocols

#endif
