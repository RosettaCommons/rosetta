// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_AddDeleteMonteCarlo
/// @brief SWA Monte Carlo -- generalization of Monte Carlo to permit addition/deletion of RNA residues
/// @detailed
/// @author Rhiju Das

#include <protocols/stepwise/monte_carlo/rna/RNA_AddDeleteMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/mover/AddOrDeleteMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_O2PrimeMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_TorsionMover.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/conformation/Residue.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>


using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
// RNA_AddDeleteMonteCarlo (a prototype for SWA Monte Carlo) -- framework for
//  sampling RNA through torsion moves,
//  and moves that delete or add residues at chain termini.
//  This probably should be deprecated soon -- could be captured by
//  StepWiseMonteCarlo
//////////////////////////////////////////////////////////////////////////

static numeric::random::RandomGenerator RG(29111);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.stepwise.monte_carlo.rna.RNA_AddDeleteMonteCarlo" ) ;

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	RNA_AddDeleteMonteCarlo::RNA_AddDeleteMonteCarlo(  mover::AddOrDeleteMoverOP rna_add_or_delete_mover,
																										 RNA_TorsionMoverOP     rna_torsion_mover,
																										 RNA_O2PrimeMoverOP      rna_o2prime_mover,
																										 core::scoring::ScoreFunctionOP scorefxn ):
		rna_add_or_delete_mover_( rna_add_or_delete_mover ),
		rna_torsion_mover_( rna_torsion_mover ),
		rna_o2prime_mover_( rna_o2prime_mover ),
		scorefxn_( scorefxn ),
		num_cycles_( 50000 ),
		output_period_( 5000 ),
		sample_range_small_( 5.0 ),
		sample_range_large_( 50.0 ),
		kT_( 0.5 ),
		do_add_delete_( true ),
		silent_file_( "" ),
		silent_file_data_( new core::io::silent::SilentFileData )
	{
		initialize_next_suite_atoms(); // used in RMSD calculations.
	}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  RNA_AddDeleteMonteCarlo::~RNA_AddDeleteMonteCarlo()
  {}

  //////////////////////////////////////////////////////////////////////////
  void
  RNA_AddDeleteMonteCarlo::apply( core::pose::Pose & pose )
	{

		using namespace protocols::moves;
		using namespace core::pose::full_model_info;

		MonteCarloOP monte_carlo = new MonteCarlo( pose, *scorefxn_, kT_ );

		//bool accepted( true );
		std::string move_type;

		for (Size count = 1; count <= num_cycles_; count++) {

			Real const random_number = RG.uniform();

			move_type = "";

			utility::vector1< Size > moving_res_list = get_moving_res_from_full_model_info( pose );

			if ( (random_number < 0.01 && do_add_delete_) || moving_res_list.size() == 0 /*got to add something!*/ ) {
				rna_add_or_delete_mover_->apply( pose, move_type );
			}

			if ( move_type.size() == 0 /*no move yet!*/ && random_number  < 0.8 ){

				Real const random_number2 = RG.uniform();
				if ( random_number2  < 0.5 ){
					move_type = "sml";
					rna_torsion_mover_->apply( pose, move_type, sample_range_small_ );
				} else {
					move_type = "lrg";
					rna_torsion_mover_->apply( pose, move_type, sample_range_large_ );
				}
			} else{
				rna_o2prime_mover_->apply( pose, move_type );
			}

			/*accepted =*/ monte_carlo->boltzmann( pose, move_type );

			//Real const current_score = (*scorefxn_)( pose );

			if ( count % output_period_ == 0  || count == num_cycles_ ) {
				std::cout << "On " << count << " of " << num_cycles_ << " trials." << std::endl;
				output_silent_file( pose, count );
			}

		}

		monte_carlo->show_counters();
		monte_carlo->recover_low( pose );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_AddDeleteMonteCarlo::initialize_next_suite_atoms(){
		// could initialize this once somewhere else.
		next_suite_atoms_.clear();
		next_suite_atoms_.push_back(" P  ");
		next_suite_atoms_.push_back(" OP2");
		next_suite_atoms_.push_back(" OP1");
		next_suite_atoms_.push_back(" O5'");
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Following assumed that pose is already properly aligned to native pose!
	// stuff into RNA_AddDeleteMonteCarlo.
	void
	RNA_AddDeleteMonteCarlo::output_silent_file( pose::Pose & pose, Size const count ){

		using namespace core::io::silent;
		using namespace core::conformation;
		using namespace ObjexxFCL;
		using namespace core::pose::full_model_info;

		utility::vector1< Size > const & working_res_list = get_moving_res_from_full_model_info( pose );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );

		// useful to keep track of what's a working residue and what's not.
		utility::vector1< bool > is_working_res;
		for ( Size i = 1; i <= pose.total_residue(); i++ ) is_working_res.push_back( false );
		for ( Size n = 1; n <= working_res_list.size(); n++ ) is_working_res[ working_res_list[n] ] = true;

		std::string const tag = "S_"+lead_zero_string_of(count,6);
		BinarySilentStruct s( pose, tag );

		Pose const & native_pose = *get_native_pose();

		Real dev( 0.0 );
		Real rmsd( 0.0 );
		Size natoms( 0 );

		for ( Size n = 1; n <= working_res_list.size(); n++ ){
			Size const i      = working_res_list[ n ];
			Size const i_full = res_list[ i ];

			//		std::cout << "NATIVE CHECK: " <<  i << ' ' << i_full << std::endl;

			Residue const & rsd        = pose.residue( i );
			Residue const & rsd_native = native_pose.residue( i_full );

			if ( rsd.aa() != rsd_native.aa() )			std::cout << "mismatch:   pose " << i << ' ' << rsd.aa() << "   native " << i_full << rsd_native.aa() << std::endl;
			runtime_assert( rsd.aa() == rsd_native.aa() );

			for ( Size j = 1; j <= rsd.nheavyatoms(); j++ ){
				if ( rsd.is_virtual( j ) ) continue;  // note virtual phosphates on 5'-ends of loops.
				std::string atom_name = rsd.atom_name( j );
				//			std::cout << "RMSD: " << i << atom_name << std::endl;
				if ( !rsd_native.has( atom_name ) ) continue;
				Size const j_full = rsd_native.atom_index( atom_name );
				dev += ( rsd_native.xyz( j_full ) - rsd.xyz( j ) ).length_squared();
				natoms++;
			}

			// also add in atoms in next suite, if relevant (and won't be covered later in rmsd calc.)
			if ( i < pose.total_residue() && !is_working_res[ i+1 ] &&  !pose.fold_tree().is_cutpoint(i) ){
				Size const i_next      = i+1;
				Size const i_next_full = res_list[ i+1 ];
				runtime_assert( i_next_full == i_full + 1 ); //better be a connection in both the pose & native pose!

				Residue const & rsd_next        = pose.residue( i_next );
				Residue const & rsd_next_native = native_pose.residue( i_next_full );
				runtime_assert( rsd_next.aa() == rsd_next_native.aa() );

				//std::cout << "RMSD:  num_suite_atoms " << next_suite_atoms_.size() << std::endl;
				for (Size k = 1; k <= next_suite_atoms_.size(); k++ ){
					std::string atom_name = next_suite_atoms_[ k ];
					//std::cout << "RMSD: " << i+1 << atom_name << std::endl;
					runtime_assert( rsd_next.has( atom_name ) );
					runtime_assert( rsd_next_native.has( atom_name ) );
					Size const j = rsd_next.atom_index( atom_name );
					Size const j_full = rsd_next_native.atom_index( atom_name );
					dev += ( rsd_next_native.xyz( j_full ) - rsd_next.xyz( j ) ).length_squared();
					natoms++;
				}

			}

		}
		if ( natoms > 0 ) rmsd = std::sqrt( dev / static_cast<Real>( natoms ) );

		s.add_energy( "rms",rmsd );
		s.add_string_value( "count", ObjexxFCL::format::I(9,count) );

		std::string built_res = "";
		if ( working_res_list.size() == 0 ) {
			built_res = "-";
		} else {
			built_res += string_of( working_res_list[1] );
			for ( Size n = 2; n <= working_res_list.size(); n++ ) built_res += "-"+string_of( res_list[ working_res_list[n] ] );
		}
		s.add_string_value( "built_res", built_res);

		silent_file_data_->write_silent_struct( s, silent_file_, false /*score_only*/ );
	}



	std::string
	RNA_AddDeleteMonteCarlo::get_name() const {
		return "RNA_AddDeleteMonteCarlo";
	}

} //rna
} //monte_carlo
} //stepwise
} //protocols
