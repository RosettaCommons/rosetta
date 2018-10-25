// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file
/// @brief Parses the VALL by secondary structure after it is read in hdf5 format
/// @author TJ Brunette <tjbrunette@uw.edu>

#ifndef INCLUDED_protocols_indexed_structure_store_SSHashedFragmentStore_hh
#define INCLUDED_protocols_indexed_structure_store_SSHashedFragmentStore_hh

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/indexed_structure_store/FragmentStore.fwd.hh>
#include <protocols/indexed_structure_store/SSHashedFragmentStore.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/types.hh>
#include <map>
#include <set>

// Utility Headers
#include <utility/SingletonBase.hh>


namespace protocols {
namespace indexed_structure_store {

class SSHashedFragmentStore : public utility::SingletonBase< SSHashedFragmentStore >
{
public:
	friend class utility::SingletonBase< SSHashedFragmentStore >;
	void set_threshold_distance(numeric::Real threshold_distance);
	void init_SS_stub_HashedFragmentStore();
	numeric::Size get_valid_resid(core::pose::Pose const & pose,int resid);
	std::set<std::string> potential_valid_ss_strings(std::string frag_ss);
	numeric::Real max_rmsd_in_region(core::pose::Pose const & pose, utility::vector1<numeric::Size> resids);
	numeric::Real lookback_account_for_dssp_inaccuracy(core::pose::Pose const & pose, numeric::Size resid,bool find_closest,numeric::Real rms_threshold);
	numeric::Real lookback_account_for_dssp_inaccuracy(core::pose::Pose const & pose, numeric::Size resid,std::string frags_ss, bool find_closest,numeric::Real rms_threshold);
	numeric::Real lookback_account_for_dssp_inaccuracy(core::pose::Pose const & pose, numeric::Size resid,std::string frag_ss, numeric::Real & match_rmsd, numeric::Size & match_index, numeric::Size & match_ss_index);
	numeric::Real lookback(core::pose::Pose const & pose, numeric::Size resid);
	numeric::Real lookback(core::pose::Pose const & pose, numeric::Size resid,std::string frag_ss,bool find_closest);
	void lookback_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, char resTypeBeforeLoop,char resTypeAfterLoop ,numeric::Size loop_length, numeric::Real & match_rmsd, numeric::Size & match_index, numeric::Size & match_ss_index);
	void lookback_uncached_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, numeric::Size stub_match_ss_index, numeric::Size loop_length, numeric::Real & match_rmsd, numeric::Size & match_index);
	std::vector< numeric::xyzVector<numeric::Real> > get_fragment_coordinates(numeric::Size db_index,numeric::Size match_index);
	// vector<FragmentLookupResult> get_N_fragments(std::string abego_string,numeric::Size topNFrags);
	// vector<FragmentLookupResult> get_topN_fragments(std::string selectionType,numeric::Size topNFrags, pose::Pose const & pose, numeric::Size resid,numeric::Real rms_threshold,std::string fragAbegoStr);
	// vector<FragmentLookupResult> get_fragments_below_rms(pose::Pose const & pose, numeric::Size resid,numeric::Real rms_threshold);
	// vector<FragmentLookupResult> get_fragments_below_rms(pose::Pose const & pose, numeric::Size resid,numeric::Real rms_threshold,std::string fragAbegoStr);
	void get_hits_below_rms(core::pose::Pose const & pose, numeric::Size resid, numeric::Real rms_threshold, utility::vector1< std::vector<numeric::Real> > & hits_cen, utility::vector1<numeric::Real> & hits_rms, utility::vector1<std::string> & hits_aa);
	protocols::indexed_structure_store::FragmentStoreOP get_fragment_store(numeric::Size db_index);
	protocols::indexed_structure_store::FragmentStoreOP get_fragment_store();
	std::map<numeric::Size, protocols::indexed_structure_store::FragmentStoreOP> get_hashed_fragment_store();
	numeric::Size get_fragment_length();
private:
	SSHashedFragmentStore();
	std::map<numeric::Size, protocols::indexed_structure_store::FragmentStoreOP> SSHashedFragmentStore_;
	std::map<std::string, utility::vector1<numeric::Size> > SS_stub_HashedFragmentStoreIndex_;
	std::vector<bool> generate_subset_residues_to_compare(numeric::Size loop_length,numeric::Size fragment_length,bool match_tail);
	static SSHashedFragmentStore * create_singleton_instance();
};

}
}
#endif
