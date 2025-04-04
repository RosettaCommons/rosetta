// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_cryst_cryst_movers_hh
#define INCLUDED_protocols_cryst_cryst_movers_hh


#include <utility>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/movemap/MoveMapFactory.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/id/AtomID.fwd.hh> // AUTO IWYU For AtomID
#include <utility/vector1.hh> // AUTO IWYU For vector1


namespace protocols {
namespace cryst {


class ReportGradientsMover : public moves::Mover {
public:
	ReportGradientsMover() :
		Mover(), verbose_(false) {}

	ReportGradientsMover(core::scoring::ScoreFunctionOP sfin) :
		Mover(), verbose_(false), score_function_(std::move(sfin)) {}

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< ReportGradientsMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;

	// compute gradients
	core::Real compute(core::pose::Pose & pose );

	// helper function normalizes gradients for a single atom
	core::Real normalization(core::pose::Pose & pose, core::id::AtomID atmid, core::scoring::ScoreFunctionOP sfxn );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	bool verbose_;
	core::scoring::ScoreFunctionOP score_function_;
	std::string outfile_;
};

class SetCrystWeightMover : public moves::Mover {
public:
	SetCrystWeightMover() :
		Mover(), autoset_wt_(true), cartesian_(false), weight_(0.0), weight_scale_(1.0), weight_min_(1.0) {}

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< SetCrystWeightMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	bool autoset_wt_;
	bool cartesian_;
	core::Real weight_, weight_scale_, weight_min_;

	core::scoring::ScoreFunctionOP score_function_;
	core::scoring::ScoreFunctionOP score_function_ref_;
	core::select::movemap::MoveMapFactoryOP mmf_;
};


class RecomputeDensityMapMover : public moves::Mover {
public:
	RecomputeDensityMapMover() : Mover(), keep_sidechains_(true) {}

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< RecomputeDensityMapMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	bool keep_sidechains_;
};

///@brief This class loads a density map into global data.
/// Each apply will regenerate this data and force a reload of the data by default.
class LoadDensityMapMover : public moves::Mover {
public:

	///@brief Default Constructor
	LoadDensityMapMover() : Mover(), mapfile_("") {}

	///@brief Constructor with mapfile
	LoadDensityMapMover( std::string const & mapfile);

	////////////////////////////////////////
	///////////// Setters //////////////////
	////////////////////////////////////////

	///@brief Set a mapfile to be loaded
	void
	set_mapfile( std::string const & mapfile);

	///@brief Set sidechain scaling in density map to this factor (default 1.0)
	void
	set_sc_scale( core::Real const sc_scale);

	///@brief Set window in the density map to this value (default 3).  Used if using window scoring.
	void
	set_window( core::Size const window);

	////////////////////////////////////////
	///////////// Getters //////////////////
	////////////////////////////////////////

	///@brief Get the mapfile path set for this mover.
	std::string
	get_mapfile();

	///@brief Get the sc_scale set for this mover.
	core::Real
	get_sc_scale();

	///@brief Get the window set for this mover.
	core::Size
	get_window();

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< LoadDensityMapMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;


	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	std::string mapfile_;
	core::Real sc_scale_;
	core::Size window_;
};


class FitBfactorsMover : public moves::Mover {
public:
	FitBfactorsMover() : Mover(), adp_strategy_("individual"), b_min_(5.0), b_max_(5.0) {}

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< FitBfactorsMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	void randomize_bs( core::pose::Pose & pose );

	std::string adp_strategy_;
	core::Real b_min_, b_max_;
};

class UpdateSolventMover : public moves::Mover {
public:
	UpdateSolventMover() : Mover(), update_mask_(true), update_fcalc_(true), optimize_mask_(false), optimize_params_(false) {}

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< UpdateSolventMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	bool update_mask_ , update_fcalc_, optimize_mask_, optimize_params_;
};


class TagPoseWithRefinementStatsMover : public moves::Mover {
public:
	TagPoseWithRefinementStatsMover() : Mover(), tag_(""), dump_pose_(false),report_grads_(false) {}

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< TagPoseWithRefinementStatsMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string tag_;
	bool dump_pose_,report_grads_;
};

class SetRefinementOptionsMover : public moves::Mover {
public:
	SetRefinementOptionsMover() :
		Mover(), res_high_(0.0), res_low_(0.0), sharpen_b_(0.0), twin_law_(""), algo_(""), target_(""), map_type_(""), setmap_type_(false)
	{}

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< SetRefinementOptionsMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Real res_high_, res_low_;
	core::Real sharpen_b_;
	std::string twin_law_, algo_, target_, map_type_;
	utility::vector1<std::string> cif_files_;
	bool setmap_type_;
};


}
}

#endif
