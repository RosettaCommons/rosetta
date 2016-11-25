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
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>
#include <protocols/generalized_kinematic_closure/perturber/GeneralizedKICperturber.hh>

#include <protocols/pose_length_moves/NearNativeLoopCloser.fwd.hh>

#include <core/indexed_structure_store/SSHashedFragmentStore.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

// C++ Headers
#include <string>
#include <map>
// Utility Headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <numeric/xyzVector.hh>


namespace protocols {
namespace pose_length_moves {

class PossibleLoop : public utility::pointer::ReferenceCount {
public:
	PossibleLoop(int resAdjustmentBeforeLoop, int resAdjustmentAfterLoop,Size loopLength,Size resBeforeLoop, Size resAfterLoop, char resTypeBeforeLoop, char resTypeAfterLoop, Size insertedBeforeLoopRes, Size insertedAfterLoopRes, core::pose::PoseOP fullLengthPoseOP, core::pose::PoseOP  orig_atom_type_fullLengthPoseOP);
	~PossibleLoop();
	void evaluate_distance_closure();
	void generate_stub_rmsd();
	void generate_uncached_stub_rmsd();
	void generate_output_pose(bool output_closed, bool ideal_loop, core::Real rms_threshold,std::string closure_type);
	core::pose::PoseOP get_finalPoseOP();
	core::Real get_stubRMSD(){return stub_rmsd_match_;}
	core::Real get_uncached_stubRMSD(){return uncached_stub_rmsd_;}
	core::Real get_final_RMSD(){return final_rmsd_;}
	bool outputed(){return outputed_;}
	void outputed(bool outputed){outputed_=outputed;}
	bool get_below_distance_threshold(){return below_distance_threshold_;}
	std::string get_description();

private:
	void trimRegion(core::pose::PoseOP & poseOP, Size resStart, Size resStop);
	void extendRegion(bool towardCTerm, Size resStart, Size numberAddRes,core::pose::PoseOP & poseOP);
	void assign_phi_psi_omega_from_lookback(Size db_index, Size fragment_index, core::pose::PoseOP & poseOP);
	std::vector<core::Real> get_center_of_mass(core::Real* coordinates, int number_of_atoms);
	void output_fragment_debug(std::vector< numeric::xyzVector<numeric::Real> > coordinates, std::string filename);
	void add_coordinate_csts_from_lookback(Size stub_ss_index_match, Size fragment_index, Size pose_residue, bool match_stub_alone, core::pose::PoseOP & poseOP);
	void add_dihedral_csts_from_lookback(Size stub_ss_index_match,Size fragment_index,Size pose_residue,core::pose::PoseOP & poseOP);
	Size get_valid_resid(int resid,core::pose::Pose const pose);
	std::vector< numeric::xyzVector<numeric::Real> > get_coordinates_from_pose(core::pose::PoseOP const poseOP,Size resid,Size length);
	core::Real rmsd_between_coordinates(std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates,std::vector< numeric::xyzVector<numeric::Real> > coordinates);
	bool kic_closure(core::scoring::ScoreFunctionOP scorefxn_, core::pose::PoseOP & poseOP,Size firstLoopRes, Size lastLoopRes,Size numb_kic_cycles);
	void minimize_loop(core::scoring::ScoreFunctionOP scorefxn,bool ideal_loop,core::pose::PoseOP & poseOP);
	Size insertedBeforeLoopRes_;
	Size insertedAfterLoopRes_;
	int resAdjustmentBeforeLoop_;
	int resAdjustmentAfterLoop_;
	char resTypeBeforeLoop_;
	char resTypeAfterLoop_;
	Size loopLength_;
	Size resBeforeLoop_;
	Size resAfterLoop_;
	Size fullLength_resBeforeLoop_;
	bool below_distance_threshold_;
	core::Real stub_rmsd_match_;
	Size stub_index_match_;
	Size stub_ss_index_match_;
	core::Real uncached_stub_rmsd_;
	Size uncached_stub_index_;
	bool outputed_;
	core::pose::PoseOP original_atom_type_fullLengthPoseOP_;
	core::pose::PoseOP fullLengthPoseOP_;
	core::pose::PoseOP finalPoseOP_;
	core::indexed_structure_store::SSHashedFragmentStore * SSHashedFragmentStore_;
	core::Real final_rmsd_;
};

struct StubRMSDComparator {
	bool operator()(const PossibleLoopOP left, const PossibleLoopOP right) const {
		return(left->get_uncached_stubRMSD() < right->get_uncached_stubRMSD());
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
	NearNativeLoopCloser(int resAdjustmentStartLow,int resAdjustmentStartHigh,int resAdjustmentStopLow,int resAdjustmentStopHigh,int resAdjustmentStartLow_sheet,int resAdjustmentStartHigh_sheet,int resAdjustmentStopLow_sheet,int resAdjustmentStopHigh_sheet,Size loopLengthRangeLow, Size loopLengthRangeHigh,Size resBeforeLoop,Size resAfterLoop,
		char chainBeforeLoop, char chainAfterLoop,core::Real rmsThreshold, core::Real max_vdw_change, bool idealExtension,bool ideal, bool output_closed, std::string closure_type="lookback");
	moves::MoverOP clone() const override { return moves::MoverOP( new NearNativeLoopCloser( *this ) ); }
	core::Real close_loop(Pose & pose);
	void apply( Pose & pose ) override;
	void combine_chains(Pose & pose);
	void extendRegion(bool towardCTerm, Size resStart, char neighborResType, Size numberAddRes,core::pose::PoseOP & poseOP);
	core::pose::PoseOP create_maximum_length_pose(char resTypeBeforeLoop, char resTypeAfterLoop, core::pose::Pose pose);
	utility::vector1<PossibleLoopOP> create_potential_loops(core::pose::Pose pose);
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
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
	Size loopLengthRangeLow_;
	Size loopLengthRangeHigh_;
	Size resBeforeLoop_;
	Size resAfterLoop_;
	char chainBeforeLoop_;
	char chainAfterLoop_;
	core::Real rmsThreshold_;
	bool idealExtension_;
	bool output_closed_;
	bool output_all_;
	bool top_outputed_;
	core::Real max_vdw_change_;
	bool ideal_;
	std::string closure_type_;
	std::string pose_name_;
	core::indexed_structure_store::SSHashedFragmentStore * SSHashedFragmentStore_;
	utility::vector1<PossibleLoopOP> possibleLoops_;
};


} // pose_length_moves
} // protocols

#endif
