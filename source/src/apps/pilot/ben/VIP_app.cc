// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/public/RosettaVIP/VIP_app.cc
/// @brief RosettaVIP protocol. Locates buried voids
/// using RosettaHoles and attempts point mutants (using GOE) on a fixed backbone to reduce void
/// volumes and improve packing.
/// @author ben (bborgo@genetics.wustl.edu)


#include <basic/options/option.hh>

//protocols library (Movers)
#include <protocols/vip/VIP_Mover.hh>
#include <basic/options/keys/cp.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility/options/keys/OptionKey.hh>
#include <core/io/pdb/file_data.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_moves/ScoreMover.hh>

//utilities
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>

//local options
namespace core{ namespace options{ namespace OptionKeys{
basic::options::BooleanOptionKey const minimize_sidechains("minimize_sidechains");
}}}//basic::options::OptionKeys

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	devel::init(argc, argv);

	core::pose::Pose in_pose;
	core::io::pdb::build_pose_from_pdb_as_is(
		in_pose,
		option[ OptionKeys::in::file::s ]().vector().front());

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	protocols::simple_moves::ScoreMover scoreme = protocols::simple_moves::ScoreMover( scorefxn );

	bool iterate = true;
	core::Size it = 1;

if ( option[ cp::ncycles ] ) {

	core::Size ncycles = option[ cp::ncycles ];
	core::pose::Pose out_pose;

	while( it <= ncycles ){
        	scoreme.apply(in_pose);
		core::Real old_energy = in_pose.energies().total_energy();

        	protocols::vip::VIP_Mover();
        	protocols::vip::VIP_Mover vip_mover;
        		vip_mover.set_initial_pose( in_pose );
        		vip_mover.apply();

		out_pose = vip_mover.get_final_pose();
        	core::Real new_energy = vip_mover.get_final_energy();

        if( new_energy < old_energy ){
                old_energy = new_energy;
                in_pose = out_pose;
                it++;}}
                out_pose.dump_pdb( option[ cp::output ] );}

else{
	while( iterate == true ){
        	scoreme.apply(in_pose);
		core::Real old_energy = in_pose.energies().total_energy();

		protocols::vip::VIP_Mover();
		protocols::vip::VIP_Mover vip_mover;
			vip_mover.set_initial_pose( in_pose );
			vip_mover.apply();

		core::pose::Pose out_pose = vip_mover.get_final_pose();
		core::Real new_energy = vip_mover.get_final_energy();

	if( new_energy < old_energy ){
		old_energy = new_energy;
		in_pose = out_pose;
		iterate = true;}
	else{
		out_pose.dump_pdb( option[cp::output] );
		iterate = false;}}}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
