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

#ifndef INCLUDED_core_indexed_structure_store_SSHashedFragmentStore_hh
#define INCLUDED_core_indexed_structure_store_SSHashedFragmentStore_hh

#include <utility/pointer/ReferenceCount.hh>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>
#include <core/indexed_structure_store/SSHashedFragmentStore.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <map>
#include <set>

// Utility Headers
#include <utility/SingletonBase.hh>


namespace core {
namespace indexed_structure_store {

class SSHashedFragmentStore : public utility::SingletonBase< SSHashedFragmentStore >
{
public:
	friend class utility::SingletonBase< SSHashedFragmentStore >;
	void set_threshold_distance(Real threshold_distance);
	void init_SS_stub_HashedFragmentStore();
	Size get_valid_resid(core::pose::Pose const pose,int resid);
	std::set<std::string> potential_valid_ss_strings(std::string frag_ss);
	Real max_rmsd_in_region(pose::Pose const pose, utility::vector1<Size> resids);
	Real lookback_account_for_dssp_inaccuracy(pose::Pose const pose, Size resid,bool find_closest,Real rms_threshold);
	Real lookback_account_for_dssp_inaccuracy(pose::Pose const pose, Size resid,std::string frags_ss, bool find_closest,Real rms_threshold);
	Real lookback_account_for_dssp_inaccuracy(pose::Pose const pose, Size resid,std::string frag_ss, Real & match_rmsd, Size & match_index, Size & match_ss_index);
	Real lookback(pose::Pose const pose, Size resid);
	Real lookback(pose::Pose const pose, Size resid,std::string frag_ss,bool find_closest);
	void lookback_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, char resTypeBeforeLoop,char resTypeAfterLoop ,Size loop_length, Real & match_rmsd, Size & match_index, Size & match_ss_index);
	void lookback_uncached_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, Size stub_match_ss_index, Size loop_length, Real & match_rmsd, Size & match_index);
	std::vector< numeric::xyzVector<numeric::Real> > get_fragment_coordinates(Size db_index,Size match_index);
	// vector<FragmentLookupResult> get_N_fragments(std::string abego_string,Size topNFrags);
	// vector<FragmentLookupResult> get_topN_fragments(std::string selectionType,Size topNFrags, pose::Pose const pose, Size resid,Real rms_threshold,std::string fragAbegoStr);
	// vector<FragmentLookupResult> get_fragments_below_rms(pose::Pose const pose, Size resid,Real rms_threshold);
	// vector<FragmentLookupResult> get_fragments_below_rms(pose::Pose const pose, Size resid,Real rms_threshold,std::string fragAbegoStr);
	void get_hits_below_rms(pose::Pose const pose, Size resid, Real rms_threshold, utility::vector1< std::vector<Real> > & hits_cen, utility::vector1<Real> & hits_rms, utility::vector1<std::string> & hits_aa);
	core::indexed_structure_store::FragmentStoreOP get_fragment_store(Size db_index);
	core::indexed_structure_store::FragmentStoreOP get_fragment_store();
	std::map<Size, core::indexed_structure_store::FragmentStoreOP> get_hashed_fragment_store();
	Size get_fragment_length();
private:
	SSHashedFragmentStore();
	std::map<Size, core::indexed_structure_store::FragmentStoreOP> SSHashedFragmentStore_;
	std::map<std::string, utility::vector1<Size> > SS_stub_HashedFragmentStoreIndex_;
	std::vector<bool> generate_subset_residues_to_compare(Size loop_length,Size fragment_length,bool match_tail);
	static SSHashedFragmentStore * create_singleton_instance();
};

}
}
#endif
