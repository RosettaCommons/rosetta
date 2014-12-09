// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelixBundleFeatures.hh
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_helixAssembly_HelixBundleFeatures_hh
#define INCLUDED_protocols_features_helixAssembly_HelixBundleFeatures_hh

//Unit
#include <protocols/features/helixAssembly/HelixBundleFeatures.fwd.hh>

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

struct FragmentPair
{
	FragmentPair(HelicalFragmentOP frag1, HelicalFragmentOP frag2):
	fragment_1(frag1),
	fragment_2(frag2)
	{}

	HelicalFragmentOP fragment_1;
	HelicalFragmentOP fragment_2;

	core::Real end_1_distance;
	core::Real end_2_distance;
	core::Real fa_attr;
	core::Real fa_fraction;
	core::Real crossing_angle;
};

typedef std::map< std::pair<core::Size, core::Size>, FragmentPair> PairMap;

class HelixBundleFeatures : public protocols::features::FeaturesReporter
{

public:

	HelixBundleFeatures();

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
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose
	);

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	bool
	overlapping(
		HelicalFragmentOP const & fragment_1,
		HelicalFragmentOP const & fragment_2
	);

	///@brief collect all the feature data for the pose
	virtual
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

	utility::vector1<HelicalFragmentOP>
	get_helix_fragments(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

	void
	calc_pc_and_com(
		core::pose::Pose const & pose,
		HelicalFragmentOP fragment
	);

	///@brief create a bundle-pose from the combination of fragments
	/// and record the "interface" SASA for each helix against the
	/// rest of the bundle
	void
	record_helix_sasas(
		core::pose::Pose const & pose,
		std::set<HelicalFragmentOP> const & frag_set
	);

	PairMap
	get_helix_pairs(
		core::pose::Pose const & pose,
		utility::vector1<HelicalFragmentOP> helix_fragments
	);

	///@brief calculate the shared fa_attr for each pair of helices
	/// in the bundle
	void
	calc_fa_energy(
		core::pose::Pose const & pose,
		FragmentPair & fragment_pair
	);

	///@brief calculate the crossing angles of the helix fragment in the bundle set
	void
	calc_crossing_angles(
		core::pose::Pose const & pose,
		FragmentPair & fragment_pair
	);

private:

	//Maximum allowable distance, in angstroms, between any two helix ends
	core::Real helix_cap_dist_cutoff_;

	//number of helices in the bundle
	core::Size bundle_size_;

	//number of residues in each helix
	core::Size helix_size_;

	//minimum residue_normalized fa_attr between two helix fragments in order to be considered
	//interacting
	core::Real min_per_residue_fa_attr_;

	//minimum fraction of
	core::Real min_interacting_set_fraction_;

	//maximum degrees off parallel for crossing angle
	core::Real max_degrees_off_parallel_;

	protocols::simple_filters::InterfaceSasaFilter sasa_filter_;

	core::scoring::ScoreFunctionOP scorefxn_;
};

} //namespace helixAssembly
} //namespace features
} //namespace protocols

#endif
