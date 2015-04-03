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
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//protocols library (Movers)
#include <protocols/vip/VIP_Mover.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/string_util.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/MiniRelax.hh>


#include <basic/Tracer.hh>
static thread_local basic::Tracer TR( "VIP" );

//utilities
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

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

	core::Real initial_E = 0.0;
	core::Real old_energy = 0.0;

	core::pose::Pose in_pose;
	core::import_pose::pose_from_pdb(
		in_pose,
		option[ OptionKeys::in::file::s ]().vector().front());

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::relax_sfxn] );
	protocols::simple_moves::ScoreMover scoreme = protocols::simple_moves::ScoreMover( scorefxn );

	// Perform a pre-relaxation of the starting point

	{
		std::string rmover = option[ cp::relax_mover ];
		if( rmover == "relax" ){
			protocols::relax::RelaxProtocolBaseOP relaxmover( new protocols::relax::FastRelax( scorefxn , 15 ) );
			relaxmover->apply(in_pose);
		} else if( rmover == "classic_relax" ){
			protocols::relax::RelaxProtocolBaseOP relaxmover( new protocols::relax::ClassicRelax( scorefxn ) );
			relaxmover->apply(in_pose);
		} else if( rmover == "cst_relax" ){
			protocols::relax::RelaxProtocolBaseOP cstrelaxmover( new protocols::relax::MiniRelax( scorefxn ) );
			cstrelaxmover->apply(in_pose);
		}


	}

	//	bool iterate = true;
	core::Size it = 1;
	core::Size const ncycles = option[ cp::ncycles ];
	core::Size const max_failures = option[ cp::max_failures ];

	bool use_unrelaxed_mutants( option[ cp::use_unrelaxed_starting_points ] );
	bool easy_acceptance( option[ cp::easy_vip_acceptance ] );

	// How many failures at a given level to tolerate before failing
	core::Size current_failures( 0 );

	core::pose::Pose stored_unrelaxed_pose;
	core::pose::Pose out_pose;
	bool not_finished( true );

	while( not_finished ){
		TR << "Entering VIP loop with it # " << it << std::endl;
		scoreme.apply(in_pose);
		if( it == 1 ) {
			old_energy = in_pose.energies().total_energy();
			initial_E = old_energy;
		}

		protocols::vip::VIP_Mover();
		protocols::vip::VIP_Mover vip_mover;
		if( ( it > 1 ) && use_unrelaxed_mutants && ( stored_unrelaxed_pose.total_residue() > 0 ) ) {
			vip_mover.set_initial_pose( stored_unrelaxed_pose );
			//TR << "Recovering stored pose with total residues = " << stored_unrelaxed_pose.total_residue() << std::endl;
		} else {
			vip_mover.set_initial_pose( in_pose );
		}
		vip_mover.set_iteration( it );
		vip_mover.apply();

		out_pose = vip_mover.get_final_pose();
		core::Real new_energy = vip_mover.get_final_energy();

		bool improved( new_energy < old_energy ? true : false );

		// Check for bogus final_pose
		if( out_pose.total_residue() == 0 ) {
			improved = false;
		} else {
			if( use_unrelaxed_mutants ) {
				stored_unrelaxed_pose = vip_mover.get_unrelaxed_pose();
			}
			TR << "Comparing new energy " << new_energy << " with old energy " << old_energy << std::endl;
		}

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

				vip_mover.set_energy_to_beat( easy_acceptance ? initial_E : new_energy );
				vip_mover.set_use_stored_energy( true );

					break;
				}
			}

			old_energy = ( easy_acceptance ? initial_E : new_energy );
			in_pose = out_pose;
			it++;
			current_failures = 0;
		} else {
			current_failures++;
			TR << "Rejecting attempted mutation!" << std::endl;
		}


		bool done_improving( (current_failures >= max_failures) && !improved );

		// if the requested ncycles is 0, quit as soon as improvement stops.  If ncycles
		// is set at a non-zero value, quit either after improvement stops or the number of
		// cycles has been hit.
		not_finished = ( ncycles == 0 ? !done_improving : (it <= ncycles && !done_improving) );
	}
	out_pose.dump_pdb( option[ cp::output ] );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
