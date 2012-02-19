// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   HelixAssemblyMover.hh
///
/// @brief
/// @author Tim jacobs

#ifndef INCLUDED_apps_pilot_tjacobs_HelixAssemblyMover_hh
#define INCLUDED_apps_pilot_tjacobs_HelixAssemblyMover_hh

//Package
#include <devel/helixAssembly/HelixAssemblyJob.hh>

//Core
#include <core/pose/Pose.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

//Protocols
#include <protocols/moves/Mover.hh>
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//mover definition
//class HelixAssemblyMover : public protocols::moves::Mover {
class HelixAssemblyMover {
public:

        ///@brief
        HelixAssemblyMover(protocols::features::helixAssembly::HelicalFragment query_frag_1, protocols::features::helixAssembly::HelicalFragment query_frag_2);
        ~HelixAssemblyMover();

        virtual std::string get_name() const {
          return "HelixAssemblyMover";
        }

        core::scoring::ScoreFunctionOP get_scorefxn() const;
        protocols::features::helixAssembly::HelicalFragment get_query_frag_1() const;
        protocols::features::helixAssembly::HelicalFragment get_query_frag_2() const;
        core::Real get_helix_cap_dist_cutoff() const;
        core::Real get_helix_contact_dist_cutoff() const;
        core::Real get_helix_pair_rmsd_cutoff() const;
        core::Size get_minimum_helix_contacts() const;
        std::string get_query_structure_path() const;
        std::string get_query_structure_string() const;
        core::Real get_single_helix_rmsd_cutoff() const;
        void set_scorefxn(const core::scoring::ScoreFunctionOP & scorefxn_);
        void set_query_frag_1(const protocols::features::helixAssembly::HelicalFragment & frag_1_);
        void set_query_frag_2(const protocols::features::helixAssembly::HelicalFragment & frag_2_);
        void set_helix_cap_dist_cutoff(core::Real helix_cap_dist_cutoff_);
        void set_helix_contact_dist_cutoff(core::Real helix_contact_dist_cutoff_);
        void set_helix_pair_rmsd_cutoff(core::Real helix_pair_rmsd_cutoff_);
        void set_minimum_helix_contacts(core::Size minimum_helix_contacts_);
        void set_query_structure_path(std::string query_structure_path_);
        void set_query_structure_string(std::string query_structure_string_);
        void set_single_helix_rmsd_cutoff(core::Real single_helix_rmsd_cutoff_);

        void init_from_options();

        void combinePoses(core::pose::Pose & pose1, const core::pose::Pose & pose2);

        utility::vector1<protocols::features::helixAssembly::HelicalFragment> findHelices(const core::pose::Pose & pose);

        utility::vector1<protocols::features::helixAssembly::HelicalFragment> findFragmentMatches(core::pose::Pose const & search_structure,
            core::pose::Pose const & query_structure, protocols::features::helixAssembly::HelicalFragment query_fragment,
            utility::vector1<protocols::features::helixAssembly::HelicalFragment> all_helices);

        bool checkHelixContacts(const core::pose::Pose & query_structure, std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> helix_pair,
            protocols::features::helixAssembly::HelicalFragment helix_to_check);

        bool closenessCheck(const core::Distance maxRange, const core::Distance end1Dist, const core::Distance end2Dist,
            const core::pose::Pose & search_structure, protocols::features::helixAssembly::HelicalFragment search_frag_1, protocols::features::helixAssembly::HelicalFragment search_frag_2);

        utility::vector1<protocols::features::helixAssembly::HelicalFragment> findPartnerHelices(core::pose::Pose const & search_structure,
            std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> query_match, utility::vector1<protocols::features::helixAssembly::HelicalFragment> all_helices,
            bool first_round, bool direction_needed);

        std::map<core::id::AtomID, core::id::AtomID> getFragmentMap(const core::pose::Pose & pose_1,
            const core::pose::Pose & pose_2, protocols::features::helixAssembly::HelicalFragment pose_1_fragment, protocols::features::helixAssembly::HelicalFragment pose_2_fragment);

        std::map<core::id::AtomID, core::id::AtomID> getFragmentPairMap(const core::pose::Pose & pose_1,
            const core::pose::Pose & pose_2, const std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> & pose_1_fragments,
            const std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> & pose_2_fragments);

        void superimposeBundles(const core::pose::Pose & query_structure, core::pose::Pose & results_structure, std::pair< protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> matching_pair);

        core::Real bb_score(core::pose::Pose & pose, core::Size unique_chain_num, core::scoring::ScoreFunctionOP & scorefxn);

        void removeDuplicateFragmentPairs(const core::pose::Pose & pose,
        		utility::vector1< std::pair<protocols::features::helixAssembly::HelicalFragment,protocols::features::helixAssembly::HelicalFragment> > & helix_pairs);

        void removeDuplicateFragments(const core::pose::Pose & pose, const utility::vector1<protocols::features::helixAssembly::HelicalFragment> & all_helix_fragments,
        		utility::vector1<protocols::features::helixAssembly::HelicalFragment> & helix_fragments);

        std::vector<HelixAssemblyJob> apply(HelixAssemblyJob & job);


private:

        core::scoring::ScoreFunctionOP scorefxn_;
        protocols::features::helixAssembly::HelicalFragment query_frag_1_;
        protocols::features::helixAssembly::HelicalFragment query_frag_2_;
        std::string query_structure_path_;
        std::string query_structure_string_;
        core::Real single_helix_rmsd_cutoff_;
        core::Real helix_pair_rmsd_cutoff_;
        core::Real helix_cap_dist_cutoff_;
        core::Real helix_contact_dist_cutoff_;
        core::Size minimum_helix_contacts_;

}; //end HelixAssemblyMover

#endif

