// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file       protocols/frag_picker/fragment_movers/FragSetFromH5Mover.hh
///
/// @brief      Pick fragments from the H5 database.
/// @details    Allows based on ABEGO def,SSdef,PredictedSS,9mers,3mers pushes to stack
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note

#ifndef INCLUDED_protocols_frag_picker_frag_movers_FragSetFromH5Mover_hh
#define INCLUDED_protocols_frag_picker_frag_movers_FragSetFromH5Mover_hh

// Project Headers
#include <basic/datacache/DataMap.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>

#include <core/indexed_structure_store/ABEGOHashedFragmentStore.hh>
#include <core/pose/Pose.hh>


#include <protocols/moves/Mover.hh>
#include <core/io/external/PsiPredInterface.fwd.hh>
#include <protocols/ss_prediction/SS_predictor.fwd.hh>

// C++ Headers
#include <string>
#include <map>

namespace protocols {
namespace frag_picker {
namespace frag_movers {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace fragment;
using utility::vector1;

class FragSetFromH5Mover : public protocols::moves::Mover {
public:
	FragSetFromH5Mover();
	void convert_ss_to_abego_helper(char tmpChar, vector1<std::string> & abego_strings);
	vector1<std::string> convert_ss_to_abegos(std::string ss);
	vector1<vector1<Real> > get_ss_prediction(const core::pose::Pose & pose);
	bool fragSet_needs_update();
	core::pose::Pose get_selected_pose(const core::pose::Pose pose);
	vector1 <std::string> get_ss_strings_for_residue_from_ssPred(vector1<vector1<Real> > ss_prediction,Size position);
	vector1 <std::string> get_abego_strings_for_residue_from_ssPred(vector1<vector1<Real> > ss_prediction,Size res,Size fragLength);
	vector1 <std::string> get_abego_strings_for_residue(const core::pose::Pose pose,Size res,Size fragLength);
	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	moves::MoverOP clone() const { return moves::MoverOP( new FragSetFromH5Mover( *this ) ); }
	virtual void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::indexed_structure_store::ABEGOHashedFragmentStore * ABEGOHashedFragmentStore_;
	basic::datacache::DataMap datamap_;
	std::string abego_;
	std::string ss_;
	std::string fragSetName_;
	bool use_svm_;
	bool use_psipred_;
	protocols::ss_prediction::SS_predictorOP ss_predictor_;
	core::io::external::PsiPredInterfaceOP psipred_interface_;
	bool use_pose_;
	bool use_ssPred_;
	Real ssPred_cutoff_;
	bool fragSetUnchanged_;
	Size nFrag_;
	std::string chain_;
};

}//frag_movers
}//frag_picker
}//protocols

#endif
