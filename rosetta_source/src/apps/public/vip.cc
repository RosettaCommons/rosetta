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
#include <basic/options/keys/cp.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//protocols library (Movers)
#include <protocols/vip/VIP_Mover.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/string_util.hh>
#include <protocols/simple_moves/ScoreMover.hh>

#include <basic/Tracer.hh>
static basic::Tracer TR("VIP");

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
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	devel::init(argc, argv);

	core::pose::Pose in_pose;
	core::io::pdb::build_pose_from_pdb_as_is(
		in_pose,
		option[ OptionKeys::in::file::s ]().vector().front());

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score12_full" );
	protocols::simple_moves::ScoreMover scoreme = protocols::simple_moves::ScoreMover( scorefxn );

//	bool iterate = true;
	core::Size it = 1;
	core::Size const ncycles = option[ cp::ncycles ];

	core::pose::Pose out_pose;
	bool not_finished( true );

	while( not_finished ){
		scoreme.apply(in_pose);
		core::Real old_energy = in_pose.energies().total_energy();

		protocols::vip::VIP_Mover();
		protocols::vip::VIP_Mover vip_mover;
		vip_mover.set_initial_pose( in_pose );
		vip_mover.apply();

		out_pose = vip_mover.get_final_pose();
		core::Real new_energy = vip_mover.get_final_energy();

		TR << "Comparing new energy " << new_energy << " with old energy " << old_energy << std::endl;
		bool improved( new_energy < old_energy ? true : false );

		if( improved ){
			// Print out the accepted mutation

			for( core::Size j = 1 ; j <= in_pose.total_residue() ; ++j ) {
				if( out_pose.residue(j).name() != in_pose.residue(j).name() ) {
					core::Size pdb_position( out_pose.pdb_info()->number(j) );
					char pdb_chain( out_pose.pdb_info()->chain(j) );
        	TR << "Accepting mutation at position " << pdb_position << " chain " << pdb_chain <<
												" from " << in_pose.residue(j).name() << " to " <<
												out_pose.residue(j).name() << std::endl;

					if( option[ cp::print_intermediate_pdbs ] ) {
						std::string pdb_file = "vip_iter_" + utility::to_string( it ) + ".pdb";
						TR << "Dumping intermediate pdb file " << pdb_file << std::endl;
						out_pose.dump_pdb( pdb_file );
					}

					break;
				}
			}

			old_energy = new_energy;
			in_pose = out_pose;
			it++;
		} else {
			TR << "Rejecting attempted mutation - finished!" << std::endl;
		}

		not_finished = ( ncycles == 0 ? improved : it <= ncycles );
	}
	out_pose.dump_pdb( option[ cp::output ] );
}
