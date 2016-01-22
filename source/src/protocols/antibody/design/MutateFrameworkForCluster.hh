// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/MutateFrameworkForCluster.hh
/// @brief Mutates Framework regions after insertion of a particular cluster
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_MutateFrameworkForCluster_hh
#define INCLUDED_protocols_antibody_design_MutateFrameworkForCluster_hh

#include <protocols/antibody/design/MutateFrameworkForCluster.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>

#include <map>

// Forward

namespace protocols {
namespace antibody {
namespace design {

struct MutantPosition{
	AntibodyNumberingSchemeEnum
		numbering_scheme_;

	std::string
		pdb_position_;

	utility::vector1<bool>
		mutants_allowed_;

};

/// @brief Mutates framework positions to what is needed for a particular cluster
/// if that cluster has needed mutations.
///
/// @details Will use clusters from pose data cache if there, otherwise,
/// will use clusters set from AntibodyInfo.
///
/// Use set_cdr or set_cdrs to limit to a particular CDR or set of CDRs.
///
///
class MutateFrameworkForCluster : public protocols::moves::Mover {
public:

	/// @brief Constructor for RosettaScripts only.
	MutateFrameworkForCluster();

	MutateFrameworkForCluster(AntibodyInfoCOP ab_info);

	MutateFrameworkForCluster(MutateFrameworkForCluster const & src);

	virtual ~MutateFrameworkForCluster();

	std::string
	get_name() const;

	//virtual moves::MoverOP
	//clone() const;

	virtual moves::MoverOP
	fresh_instance() const;

	virtual void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const &,
		moves::Movers_map const &,
		Pose const & pose
	);

	virtual void
	apply(core::pose::Pose & pose);


public:

	void
	set_cdr_only(CDRNameEnum const & cdr);

	void
	set_cdrs(utility::vector1<bool> const & cdrs);

	/// @brief Set the distance for the packing shell.
	void
	set_pack_shell(core::Real const pack_shell);

	void
	set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn);

	/// @brief Set custom data instead of loading it from the database.
	void
	set_custom_data(std::map<clusters::CDRClusterEnum, utility::vector1<MutantPosition> > const & mutant_info);

	/// @brief Get data used to do the design.
	std::map<clusters::CDRClusterEnum, utility::vector1<MutantPosition> >
	get_data();

	/// @brief Set to keep the current amino acid(s) for the position in the mutation set.
	/// Each position is not 100 percent a particular mutant for a cluster, so we allow this by default and hope Rosetta can choose properly.
	void
	keep_current( bool keep_current);

public:

	/////////// Helper functions //////////////

	bool
	has_framework_dependant_clusters(core::pose::Pose const & pose);

	bool
	has_framework_dependant_cluster(core::pose::Pose const & pose, CDRNameEnum const cdr);

	/// @brief Get all framework dependant clusters in list
	utility::vector1< clusters::CDRClusterEnum >
	framework_dependant_clusters();

	/// @brief Get all framework positions that would be required by clusters regardless of whether that cluster is in the pose.
	utility::vector1< bool >
	framework_dependant_positions(core::pose::Pose const & pose);

	/// @brief Get all framework positions that would be required by a cluster regardless of whether that cluster is in the pose.
	utility::vector1< bool >
	framework_dependant_positions(core::pose::Pose const & pose, clusters::CDRClusterEnum const cluster);

	/// @brief Get all framework mutations for a particular framework dependent position, regardless of whether that cluster is in the pose.
	utility::vector1<bool>
	framework_dependant_mutations(core::pose::Pose const & pose, clusters::CDRClusterEnum const cluster, core::Size const resnum);

private:

	void
	set_defaults();

	void
	load_data();

private:

	AntibodyInfoCOP ab_info_;
	core::scoring::ScoreFunctionCOP scorefxn_;

	std::map<clusters::CDRClusterEnum, utility::vector1<MutantPosition> > mutant_info_;
	utility::vector1<bool> cdrs_;
	core::Real pack_shell_;
	bool keep_current_;



};

} //design
} //antibody
} //protocols



#endif //INCLUDED_protocols_antibody_design_MutateFrameworkForCluster_hh







