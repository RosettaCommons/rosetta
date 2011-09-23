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

//mover definition
//class HelixAssemblyMover : public protocols::moves::Mover {
class HelixAssemblyMover {
public:

        ///@brief
        HelixAssemblyMover();

        virtual ~HelixAssemblyMover();

        virtual std::string get_name() const {
          return "HelixAssemblyMover";
        }

        core::scoring::ScoreFunctionOP get_scorefxn() const;
        HelicalFragment get_query_frag_1() const;
        HelicalFragment get_query_frag_2() const;
        core::Real get_helix_cap_distance_cutoff() const;
        core::Real get_helix_contact_distance_cutoff() const;
        core::Real get_helix_pair_rmsd_cutoff() const;
        core::Size get_minimum_helix_contacts() const;
        std::string get_query_structure_path() const;
        std::string get_query_structure_string() const;
        core::Real get_single_helix_rmsd_cutoff() const;
        void set_scorefxn(core::scoring::ScoreFunctionOP scorefxn_);
        void set_query_frag_1(HelicalFragment frag_1_);
        void set_query_frag_2(HelicalFragment frag_2_);
        void set_helix_cap_distance_cutoff(core::Real helix_cap_distance_cutoff_);
        void set_helix_contact_distance_cutoff(core::Real helix_contact_distance_cutoff_);
        void set_helix_pair_rmsd_cutoff(core::Real helix_pair_rmsd_cutoff_);
        void set_minimum_helix_contacts(core::Size minimum_helix_contacts_);
        void set_query_structure_path(std::string query_structure_path_);
        void set_query_structure_string(std::string query_structure_string_);
        void set_single_helix_rmsd_cutoff(core::Real single_helix_rmsd_cutoff_);

        void init_from_options();

        core::pose::Pose combinePoses(const core::pose::Pose & pose1, const core::pose::Pose & pose2);

        utility::vector1<std::pair<core::Size,core::Size> > findHelices(const core::pose::Pose & pose);

        utility::vector1<HelicalFragment> findFragmentMatches(core::pose::Pose const & search_structure,
            core::pose::Pose const & query_structure, HelicalFragment query_fragment,
            utility::vector1< std::pair< core::Size,core::Size > > helix_endpts);

        bool checkHelixContacts(const core::pose::Pose & query_structure, const core::pose::Pose & fragment1,
            const core::pose::Pose & fragment2, HelicalFragment helix_to_check);

        bool closenessCheck(const core::Distance maxRange, const core::Distance end1Dist, const core::Distance end2Dist,
            const core::pose::Pose & search_structure, HelicalFragment search_frag_1, HelicalFragment search_frag_2);

        utility::vector1<HelicalFragment> findPartnerHelices(core::pose::Pose const & search_structure,
            core::pose::Pose const & fragment1, core::pose::Pose const & fragment2,
            std::pair<HelicalFragment, HelicalFragment> helix_pair, utility::vector1< std::pair< core::Size,core::Size > > helix_endpts,
            bool first_round, bool direction_needed);

        std::map<core::id::AtomID, core::id::AtomID> getFragmentMap(const core::pose::Pose & pose_1,
            const core::pose::Pose & pose_2, HelicalFragment pose_1_fragment, HelicalFragment pose_2_fragment);

        std::map<core::id::AtomID, core::id::AtomID> getFragmentPairMap(const core::pose::Pose & pose_1,
            const core::pose::Pose & pose_2, const std::pair<HelicalFragment, HelicalFragment> & pose_1_fragments,
            const std::pair<HelicalFragment, HelicalFragment> & pose_2_fragments);

        void superimposeBundles(core::pose::Pose & query_structure, const core::pose::Pose & results_structure);

        core::Real bb_score(core::pose::Pose & pose, core::Size unique_chain_num, core::scoring::ScoreFunctionOP & scorefxn);

        std::vector<HelixAssemblyJob> apply(HelixAssemblyJob & job);


private:

        core::scoring::ScoreFunctionOP scorefxn_;
        HelicalFragment query_frag_1_;
        HelicalFragment query_frag_2_;
        std::string query_structure_path_;
        std::string query_structure_string_;
        core::Real single_helix_rmsd_cutoff_;
        core::Real helix_pair_rmsd_cutoff_;
        core::Real helix_cap_distance_cutoff_;
        core::Real helix_contact_distance_cutoff_;
        core::Size minimum_helix_contacts_;

}; //end HelixAssemblyMover

#endif

