// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/md/CartesianMD.hh
/// @brief   Cartesian MD
/// @detailed
/// @author  Hahnbeom Park
#ifndef INCLUDED_protocols_md_CartesianMD_hh
#define INCLUDED_protocols_md_CartesianMD_hh

#include <protocols/md/CartesianMD.fwd.hh>
#include <protocols/md/MDBase.hh>
#include <protocols/md/Rattle.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizerMap.hh>

#ifdef WIN32
#include <time.h>
#include <WinSock2.h>
#else
#include <sys/time.h>
#endif

namespace protocols {
namespace md {

class CartesianMD : public protocols::md::MDBase {

public:
	//constructor
	CartesianMD( );

	//constructor
	CartesianMD( core::pose::Pose const &pose, 
	 core::scoring::ScoreFunctionCOP sfxn,
	 core::kinematics::MoveMapCOP movemap );

	void
	init( );
	
	//destructor
	virtual ~CartesianMD();

	// From Mover
	virtual	protocols::moves::MoverOP	fresh_instance() const { return new CartesianMD(); };
	virtual protocols::moves::MoverOP clone() const;
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	// From MDbase
	virtual void set_movemap(
		core::pose::Pose const &,
		core::kinematics::MoveMapCOP movemap);

	virtual core::kinematics::MoveMapOP movemap() const { return movemap_; }

	void use_rattle( bool const value ){
		use_rattle_ = value;
		if( use_rattle_ ){
			dt_ = 0.002;
		} else {
			dt_ = 0.001;
		}
	}

	virtual
	void parse_my_tag(
		TagPtr const,
		protocols::moves::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	void 
	parse_opts(
  	TagPtr const tag,
	  protocols::moves::DataMap & data,
  	//Filters_map const &,
  	//protocols::moves::Movers_map const &,
  	Pose const & pose );

	void 
	parse_movemap(
  	TagPtr const tag,
  	protocols::moves::DataMap & data,
  	//Filters_map const &,
  	//protocols::moves::Movers_map const &,
  	Pose const & pose );

	utility::vector1< pose::Pose > 
	dump_poses( pose::Pose const &pose_ref ) const;

private:
	void get_native_info( pose::Pose const &pose );
	
	void do_initialize(
		core::pose::Pose &pose, 
		core::Real const &temp0 );

	// deprecated
	void Berendsen_Integrator( core::pose::Pose & pose, 
														 core::optimization::CartesianMinimizerMap &min_map );
	
	void VelocityVerlet_Integrator( core::pose::Pose & pose,
																	core::optimization::CartesianMinimizerMap &min_map,
																	md::Rattle & rattle );

	void do_minimize( pose::Pose &pose,
										core::optimization::MinimizerOptions const &options,
										bool const &show_energy );

	void do_MD( core::pose::Pose & pose,
							core::Size const &nstep,
							core::Real const &temp0 = 300, 
							bool const &initialize = false );

	void initialize_velocity( core::Real const &temperature );

	void report_MD( core::pose::Pose &pose );

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
