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
/// @brief A wrapper to allow me to load the hashed fragment store on the datacache.
/// @author TJ Brunette <tjbrunette@uw.edu>

#ifndef INCLUDED_core_indexed_structure_store_ABEGOHashedFragmentStore_hh
#define INCLUDED_core_indexed_structure_store_ABEGOHashedFragmentStore_hh

#include <utility/pointer/ReferenceCount.hh>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>
#include <core/indexed_structure_store/ABEGOHashedFragmentStore.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <map>
#include <set>

// Utility Headers
#include <utility/SingletonBase.hh>


namespace core
{
namespace indexed_structure_store
{
using utility::vector1;

class ABEGOHashedFragmentStore : public utility::SingletonBase< ABEGOHashedFragmentStore >
{
public:
	friend class utility::SingletonBase< ABEGOHashedFragmentStore >;
	void set_threshold_distance(Real threshold_distance);
	void generate_ss_stub_to_abego();
	void lookback_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, std::string stub_ss, Real & match_rmsd, Size & match_index, Size & match_abego);
	void lookback_uncached_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, std::string short_stub_ss,std::string full_stub_ss, Real & match_rmsd, Size & match_index, Size & match_abego);
	Real lookback(pose::Pose const pose, Size resid);
	Real lookback(pose::Pose const pose, Size resid,std::string fragAbegoStr);
	std::vector< numeric::xyzVector<numeric::Real> > lookback_xyz(pose::Pose const pose, Size resid);
	std::vector< numeric::xyzVector<numeric::Real> > get_fragment_coordinates(Size match_abego,Size match_index);
	std::vector<FragmentLookupResult> get_N_fragments(std::string abego_string,Size topNFrags);
	std::vector<FragmentLookupResult> get_topN_fragments(std::string selectionType,Size topNFrags, pose::Pose const pose, Size resid,Real rms_threshold,std::string fragAbegoStr);
	std::vector<FragmentLookupResult> get_fragments_below_rms(pose::Pose const pose, Size resid,Real rms_threshold);
	std::vector<FragmentLookupResult> get_fragments_below_rms(pose::Pose const pose, Size resid,Real rms_threshold,std::string fragAbegoStr);
	core::indexed_structure_store::FragmentStoreOP get_fragment_store(std::string fragment_abego);
	core::indexed_structure_store::FragmentStoreOP get_fragment_store(Size base5index);
	core::indexed_structure_store::FragmentStoreOP get_fragment_store();
	std::string get_abego_string(pose::Pose const pose, Size resid);
	Size get_fragment_length();
private:
	ABEGOHashedFragmentStore();
	std::map<Size, core::indexed_structure_store::FragmentStoreOP> ABEGOHashedFragmentStore_;
	std::map<std::string, vector1<Size> > ss_stub_to_abego_;
	void init_ss_stub_to_abego();
	std::vector<bool>  generate_ss_subset_match(FragmentStoreOP fragStoreOP, std::string stub_ss);
	bool valid_ss_stub_abego_match(std::string ss_stub,std::string abego_string);
	std::set<std::string> get_ss_stubs_per_fragmentStoreOP(Size base5ABEGOindex, core::indexed_structure_store::FragmentStoreOP fragment_storeOP);
	static ABEGOHashedFragmentStore * create_singleton_instance();
};

}
}
#endif
