// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/pose_length_moves/NearNativeLoopCloser.hh
/// @brief

#ifndef INCLUDED_protocols_pose_length_moves_NearNativeLoopCloser_hh
#define INCLUDED_protocols_pose_length_moves_NearNativeLoopCloser_hh


#include <protocols/moves/Mover.hh>

#include <protocols/pose_length_moves/NearNativeLoopCloser.fwd.hh>

#include <protocols/indexed_structure_store/SSHashedFragmentStore.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// C++ Headers
#include <string>
// Utility Headers
#include <core/types.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/xyzVector.fwd.hh>


namespace protocols {
namespace pose_length_moves {


class PossibleLoop : public utility::VirtualBase {
public:
	PossibleLoop(int resAdjustmentBeforeLoop, int resAdjustmentAfterLoop,core::Size loopLength,core::Size resBeforeLoop, core::Size resAfterLoop, char resTypeBeforeLoop, char resTypeAfterLoop, core::Size insertedBeforeLoopRes, core::Size insertedAfterLoopRes, core::pose::PoseOP fullLengthPoseOP, core::pose::PoseOP  orig_atom_type_fullLengthPoseOP, std::string fragment_store_path, std::string fragment_store_format, std::string fragment_store_compression,core::Size numb_stubs_to_consider);
	~PossibleLoop() override;
	void evaluate_distance_closure();
	void generate_stub_rmsd(core::Real stubRmsdThreshold);
	//void generate_uncached_stub_rmsd();
	void generate_output_pose(bool output_closed, bool ideal_loop, core::Real rms_threshold,std::string allowed_loop_abegos,std::string closure_type);
	void setup_finalPose_copy_labels();
	void setup_finalPose_copy_rotamers();
	core::pose::PoseOP get_finalPoseOP();
	core::Real get_stubRMSD(){return stub_rmsd_top_match_;}
	//core::Real get_uncached_stubRMSD(){return uncached_stub_rmsd_;}
	core::Real get_final_RMSD(){return final_rmsd_;}
	bool outputed(){return outputed_;}
	void outputed(bool outputed){outputed_=outputed;}
	bool get_below_distance_threshold(){return below_distance_threshold_;}
	void label_loop(std::string label);
	std::string get_description();

private:
	void trimRegion(core::pose::PoseOP & poseOP, core::Size resStart, core::Size resStop);
	void extendRegion(bool towardCTerm, core::Size resStart, core::Size numberAddRes,core::pose::PoseOP & poseOP);
	void generate_overlap_range(core::Size & front_overlap, core::Size & back_overlap);
	bool check_loop_abego(core::pose::PoseOP & poseOP, core::Size loopStart, core::Size loopStop, std::string allowed_loop_abegos,core::Real current_rmsd);
	void assign_phi_psi_omega_from_lookback(core::Size db_index, core::Size fragment_index, core::pose::PoseOP & poseOP,core::Size front_overlap_res_length);
	std::vector<core::Real> get_center_of_mass( const core::Real* coordinates, int number_of_atoms);
	void output_fragment_debug(std::vector< numeric::xyzVector<numeric::Real> > coordinates, std::string filename);
	void add_coordinate_csts_from_lookback(core::Size stub_ss_index_match, core::Size fragment_index, core::Size pose_residue, bool match_stub_alone, core::Size front_overlap_res_length, core::pose::PoseOP & poseOP);
	void add_dihedral_csts_from_lookback(core::Size stub_ss_index_match,core::Size fragment_index,core::Size pose_residue,core::pose::PoseOP & poseOP);
	core::Size get_valid_resid(int resid,core::pose::Pose const pose);
	std::vector< numeric::xyzVector<numeric::Real> > get_coordinates_from_pose(core::pose::PoseOP const poseOP,core::Size resid,core::Size length);
	core::Real rmsd_between_coordinates(std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates,std::vector< numeric::xyzVector<numeric::Real> > coordinates);
	bool kic_closure(core::scoring::ScoreFunctionOP scorefxn_, core::pose::PoseOP & poseOP,core::Size firstLoopRes, core::Size lastLoopRes,core::Size numb_kic_cycles);
	void minimize_loop(core::scoring::ScoreFunctionOP scorefxn,bool ideal_loop,core::pose::PoseOP & poseOP);
	core::Size insertedBeforeLoopRes_;
	core::Size insertedAfterLoopRes_;
	int resAdjustmentBeforeLoop_;
	int resAdjustmentAfterLoop_;
	char resTypeBeforeLoop_;
	char resTypeAfterLoop_;
	core::Size loopLength_;
	core::Size resBeforeLoop_;
	core::Size resAfterLoop_;
	core::Size fullLength_resBeforeLoop_;
	bool below_distance_threshold_;
	core::Size numb_stubs_to_consider_;
	utility::vector1<indexed_structure_store::BackboneStub> stubVector_;
	core::Real stub_rmsd_top_match_;
	bool outputed_;
	core::pose::PoseOP original_atom_type_fullLengthPoseOP_;
	core::pose::PoseOP fullLengthPoseOP_;
	core::pose::PoseOP finalPoseOP_;
	protocols::indexed_structure_store::SSHashedFragmentStoreOP SSHashedFragmentStoreOP_;
	core::Real final_rmsd_;
	std::string fragment_store_path_;
	std::string fragment_store_format_;
	std::string fragment_store_compression_;

};

struct BackboneStubRMSDComparator {
	bool operator()(const PossibleLoopOP left, const PossibleLoopOP right) const {
		return(left->get_stubRMSD() < right->get_stubRMSD());
	}
};

struct FinalRMSDComparator {
	bool operator()(const PossibleLoopOP left, const PossibleLoopOP right) const {
		return(left->get_final_RMSD() < right->get_final_RMSD());
	}
};

class NearNativeLoopCloser : public protocols::moves::Mover {
public:
	NearNativeLoopCloser();
	NearNativeLoopCloser(int resAdjustmentStartLow,int resAdjustmentStartHigh,int resAdjustmentStopLow,int resAdjustmentStopHigh,int resAdjustmentStartLow_sheet,int resAdjustmentStartHigh_sheet,int resAdjustmentStopLow_sheet,int resAdjustmentStopHigh_sheet,core::Size loopLengthRangeLow, core::Size loopLengthRangeHigh,core::Size resBeforeLoop,core::Size resAfterLoop,std::string const & chainBeforeLoop, std::string const & chainAfterLoop,core::Real rmsThreshold, core::Real max_vdw_change, bool idealExtension,bool ideal, bool output_closed, std::string closure_type="lookback",std::string allowed_loop_abegos="",std::string label_loop="", std::string fragment_store_path="",std::string fragment_store_format="",std::string fragment_store_compression="",core::Size numb_stubs_to_consider=1);
	moves::MoverOP clone() const override { return utility::pointer::make_shared< NearNativeLoopCloser >( *this ); }
	core::Real close_loop(Pose & pose);
	void apply( Pose & pose ) override;
	void combine_chains(Pose & pose);
	void switch_chain_order(Pose & pose, utility::vector1<core::Size> new_chain_order);
	void extendRegion(bool towardCTerm, core::Size resStart, char neighborResType, core::Size numberAddRes,core::pose::PoseOP & poseOP);
	core::pose::PoseOP create_maximum_length_pose(char resTypeBeforeLoop, char resTypeAfterLoop, core::pose::Pose pose);
	utility::vector1<PossibleLoopOP> create_potential_loops(core::pose::Pose pose);
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap ) override;
	core::pose::PoseOP get_additional_output() override;
	core::pose::PoseOP get_additional_output_with_rmsd(core::Real & return_rmsd);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	int resAdjustmentStartLow_;
	int resAdjustmentStartHigh_;
	int resAdjustmentStopLow_;
	int resAdjustmentStopHigh_;
	int resAdjustmentStartLow_sheet_;
	int resAdjustmentStartHigh_sheet_;
	int resAdjustmentStopLow_sheet_;
	int resAdjustmentStopHigh_sheet_;
	core::Size loopLengthRangeLow_;
	core::Size loopLengthRangeHigh_;
	core::Size resBeforeLoop_;
	core::Size resAfterLoop_;
	std::string chainBeforeLoop_;
	std::string chainAfterLoop_;
	core::Real rmsThreshold_;
	bool idealExtension_;
	bool output_closed_;
	bool output_all_;
	core::Size numb_outputed_;
	core::Size max_number_of_results_;
	bool top_outputed_;
	core::Real max_vdw_change_;
	bool ideal_;
	std::string closure_type_;
	std::string pose_name_;
	std::string label_loop_;
	protocols::indexed_structure_store::SSHashedFragmentStoreOP SSHashedFragmentStoreOP_;
	utility::vector1<PossibleLoopOP> possibleLoops_;
	std::string allowed_loop_abegos_;
	std::string fragment_store_path_;
	std::string fragment_store_format_;
	std::string fragment_store_compression_;
	core::Size numb_stubs_to_consider_;
};


} // pose_length_moves
} // protocols

#endif
