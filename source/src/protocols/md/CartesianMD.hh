// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/md/CartesianMD.hh
/// @brief   Cartesian MD
/// @details
/// @author  Hahnbeom Park
#ifndef INCLUDED_protocols_md_CartesianMD_hh
#define INCLUDED_protocols_md_CartesianMD_hh

#include <protocols/md/CartesianMD.fwd.hh>
#include <protocols/md/MDBase.hh>
#include <protocols/md/Rattle.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/types.hh>

#ifdef WIN32
#include <time.h>
#include <WinSock2.h>
#else
#include <sys/time.h>
#endif

namespace protocols {
namespace md {

using namespace core::optimization;

class CartesianMD : public protocols::md::MDBase {

public:
	//constructor
	CartesianMD( );

	//constructor
	CartesianMD( core::pose::Pose const &pose,
		core::scoring::ScoreFunctionCOP sfxn,
		core::kinematics::MoveMapCOP movemap = nullptr );

	CartesianMD( core::pose::Pose const & pose,
		core::scoring::ScoreFunction const &sfxn );

	CartesianMD( core::pose::Pose const & pose,
		core::scoring::ScoreFunction const &sfxn,
		core::kinematics::MoveMap const &movemap );

	//destructor
	~CartesianMD() override;

	// From Mover
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new CartesianMD() ); };
	protocols::moves::MoverOP clone() const override;
	void apply( core::pose::Pose & pose ) override;
	// XRW TEMP  std::string get_name() const override;

	// From MDbase
	//void set_movemap(
	//core::pose::Pose const &,
	//	core::kinematics::MoveMapCOP movemap) override;
	//core::kinematics::MoveMapOP movemap() const override { return movemap_; }

	void use_rattle( bool const value );

	Multivec get_current_eqxyz() const;
	void update_restraint( core::pose::Pose & pose,
		core::optimization::CartesianMinimizerMap const &min_map );

	void cst_on_pose_simple( core::pose::Pose &pose ) const;

	void cst_on_pose_dynamic( core::pose::Pose &pose,
		core::optimization::Multivec const &ref_xyz,
		core::optimization::Multivec const &curr_eqxyz,
		core::optimization::Multivec &prv_eqxyz,
		core::optimization::CartesianMinimizerMap const &min_map ) const;


	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & ) override;

	void
	parse_opts(
		TagCOP tag,
		basic::datacache::DataMap & data,
		//Filters_map const &,
		//protocols::moves::Movers_map const &,
		Pose const & pose );

	void
	parse_movemap(
		TagCOP tag,
		basic::datacache::DataMap & data,
		//Filters_map const &,
		//protocols::moves::Movers_map const &,
		Pose const & pose );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	utility::vector1< core::pose::Pose >
	dump_poses( core::pose::Pose const &pose_ref ) const;

private:
	void get_native_info( core::pose::Pose const &pose );

	void do_initialize( core::pose::Pose &pose );

	// deprecated
	void Berendsen_Integrator( core::pose::Pose & pose,
		core::optimization::CartesianMinimizerMap &min_map );

	void VelocityVerlet_Integrator( core::pose::Pose & pose,
		core::optimization::CartesianMinimizerMap &min_map,
		md::Rattle & rattle,
		bool const update_score = false );

	void do_minimize( core::pose::Pose &pose,
		core::optimization::MinimizerOptions const &options,
		bool const &show_energy );

	void do_MD( core::pose::Pose & pose,
		core::Size const &nstep,
		core::Real const &temp0 = 300,
		bool const &initialize = false );

	void initialize_velocity( core::Real const &temperature );

	void report_MD( core::pose::Pose &pose,
		core::optimization::CartesianMinimizerMap const &min_map,
		bool const report_trj );

private:

	core::optimization::CartesianMinimizerMap min_map_;
	timeval inittime_;
	bool use_rattle_;

	core::pose::Pose native_;
	bool native_given_;
	std::map< Size, Size > native_resmap_;

}; //class

} //namespace md
} //namespace protocols

#endif
