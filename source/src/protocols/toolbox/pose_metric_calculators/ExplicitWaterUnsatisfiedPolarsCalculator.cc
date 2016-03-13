// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin ExplicitWaterUnsatisfiedPolarsCalculator
///
/// @brief This Calculator tries to solvate all polar groups by
/// 				docking explicit TP5 water molecules, then counts
///					unsatisfied hydrogen bonds using the same criteria in 
///					BuriedUnsatisfiedHydrogenBondCalculator
///				
///
///
///
///
///
///
///
///
///
///
///
///
///
/// @author Chris King - dr.chris.king@gmail.com
/// 
///
/// @last_modified 4.8.2011
/////////////////////////////////////////////////////////////////////////
#include <protocols/toolbox/pose_metric_calculators/ExplicitWaterUnsatisfiedPolarsCalculator.hh>

#include <numeric/random/random.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

//protocol headers
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RB_geometry.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <basic/MetricValue.hh>

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

	////////////////////////////////////////////////
	// danger USING ////////////////////////////////
	using namespace core;
	using namespace basic;
	using namespace id;
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace kinematics;
	using namespace scoring;
	using namespace options;
	using namespace basic::options::OptionKeys;
	using namespace optimization;
	using utility::vector1;

static thread_local basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.ExplicitWaterUnsatisfiedPolarsCalculator" );

	///@brief default constructor sets shell_cutoff to 4.0. 
	ExplicitWaterUnsatisfiedPolarsCalculator::ExplicitWaterUnsatisfiedPolarsCalculator( ScoreFunctionOP scorefxn ):
		scorefxn_( scorefxn ),
		shell_cutoff_( 4.0 ),
		all_unsat_polars_( 0 )
	{

	}


	///@brief constructor user defined water shell max distance
	ExplicitWaterUnsatisfiedPolarsCalculator::ExplicitWaterUnsatisfiedPolarsCalculator( ScoreFunctionOP scorefxn, Real shell_cutoff ):
		scorefxn_( scorefxn ),
		shell_cutoff_( shell_cutoff ),
		all_unsat_polars_( 0 )
	{

	}





void ExplicitWaterUnsatisfiedPolarsCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const{
	if ( key == "all_unsat_polars" ) { 
			basic::check_cast( valptr, &all_unsat_polars_, "all_unsat_polars expects to return a Size" );
			(static_cast<basic::MetricValue<Size> *>(valptr))->set( all_unsat_polars_ );
	}else {
			basic::Error() << "ExplicitWaterUnsatisfiedPolarsCalculator cannot compute the requested metric " << key << std::endl;
			utility_exit();
	}
}


std::string ExplicitWaterUnsatisfiedPolarsCalculator::print( std::string const & key ) const{
  if ( key == "all_unsat_polars" ) {
    return utility::to_string( all_unsat_polars_ );
	}
	basic::Error() << "ExplicitWaterUnsatisfiedPolarsCalculator cannot compute metric " << key << std::endl;
	  utility_exit();
	  return "";

}

void
append_rsd_by_jump_near_atom(
  pose::Pose & pose,
  Size seqpos, 
  Size atomno,
  conformation::Residue new_rsd,
  Size new_atomno,
  Real dist_min,
  Real dist_max
) 
{
  typedef  numeric::xyzMatrix< Real > Matrix;
  
  Residue rsd( pose.residue( seqpos ) );
  //append by jump from seqpos atomno to new_rsd atom 1, maybe make random downstream atom?
  pose.append_residue_by_jump( new_rsd, seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), true );
  Size new_seqpos( pose.total_residue() );
  Size jump_number( pose.fold_tree().num_jump() );
  Jump jump( pose.jump( jump_number ) );

  //set jump distance as random val from dist_min to dist_max
  Real jump_dist( dist_min + numeric::random::rg().uniform() * ( dist_max - dist_min ) );
  jump.random_trans( jump_dist );
  //set jump rotation as random matrix
  jump.set_rotation( protocols::geometry::random_reorientation_matrix( 360, 360 ) );

  //and set jump in pose at our new jump
  pose.set_jump( jump_number, jump );

}

void
dock_waters_to_atom(
  pose::Pose & pose,
  ScoreFunctionOP scorefxn,
  Size seqpos,
  Size atomno,
  conformation::Residue wat_rsd,
  Size new_atomno,
  Real dist_min,
  Real dist_max
)
{
  typedef  numeric::xyzMatrix< Real > Matrix;

  //attempt appending in 20 random orientations
	Size const max_attempt_per_atom( 20 );
	Size const max_wat_per_atom( 5 ); 
	Size n_wat( 0 ); 
  pose::Pose start_pose( pose );
  protocols::moves::MonteCarloOP mc_create( new protocols::moves::MonteCarlo( pose, *scorefxn, 0.8 ) );
  for( Size i = 1; i <= max_attempt_per_atom; ++i ){
		if( n_wat >= max_wat_per_atom ) break;

    pose = start_pose;
		append_rsd_by_jump_near_atom( pose, seqpos, atomno, wat_rsd, new_atomno, dist_min, dist_max );

    Size jump_number( pose.fold_tree().num_jump() );
    //gaussian perturbations to RB dofs
    protocols::moves::MonteCarloOP mc_dock( new protocols::moves::MonteCarlo( pose, *scorefxn, 0.8 ) );
    for( Size i = 1; i <= 10; ++i ){
      Jump jump( pose.jump( jump_number ) );
      jump.gaussian_move( 1, 0.05, 90.0 );
      pose.set_jump( jump_number, jump );
      mc_dock->boltzmann( pose );
    }
    mc_dock->recover_low( pose );
		//if accept, minimize
    if( mc_create->boltzmann( pose ) ){
				MoveMapOP mm = new MoveMap;
				mm->set_jump( pose.fold_tree().num_jump(), true );
				protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, false );
				min_mover->apply( pose );
				++n_wat;
		}
  }
  mc_create->recover_low( pose );
}

//count residue unsat polar atoms, pass back bool vector of unsat atoms
void
find_res_unsat_polars(
    Pose const pose,
    Size const seqpos,
		vector1< bool > & atm_is_unsat
    )
{
  using namespace core::scoring::hbonds;
  HBondSet hbset;
  hbset.setup_for_residue_pair_energies( pose, false, false );
  conformation::Residue rsd( pose.residue( seqpos ) );
  //count n hbonds for donors and acceptors
	vector1< Size > n_atom_hbonds( rsd.natoms(), 0 );
	for( Size i = 1; i <= hbset.nhbonds(); ++i ){
		HBond hb( hbset.hbond( i ) );
		if( hb.don_res() == seqpos ){
			Size don_atm( rsd.atom_base( hb.don_hatm() ) );
			n_atom_hbonds[ don_atm ] += 1;
		}
		else if( hb.acc_res() == seqpos ){
			n_atom_hbonds[ hb.acc_atm() ] += 1;
		}
	}
  //calc unsat hbonds a la BuriedUnsatisfied calculator
  for( Size atm = 1; atm <= rsd.nheavyatoms(); ++atm ){
    if( !( rsd.atom_type( atm ).is_acceptor() || rsd.atom_type( atm ).is_donor() ) ) continue;
    //we need at least one hbond
    if( n_atom_hbonds[ atm ] == 0 ){
      atm_is_unsat[ atm ] = true; continue;
    }
    //may need > 1 hbond, see BuriedUnsatCalc
    //behavior copied from Florian!
    std::string atom_type( rsd.type().atom_type( atm ).name() );
    Size satisfac_cut = 3;
    if( atom_type == "OH" || atom_type == "OCbb" || atom_type == "S" ){
      satisfac_cut = 2;
    }
    Size bonded_heavyatoms = rsd.n_bonded_neighbor_all_res( atm )
      - rsd.type().number_bonded_hydrogens( atm );
    if( bonded_heavyatoms + n_atom_hbonds[ atm ] < satisfac_cut ){
      atm_is_unsat[ atm ] = true; 
    }
  }
}

///@brief this should just be called "compute"
void ExplicitWaterUnsatisfiedPolarsCalculator::recompute( core::pose::Pose const & pose_in ){
	ScoreFunctionOP scorefxn( scorefxn_ );
	all_unsat_polars_ = 0;
  ResidueTypeSet const & rsd_set( pose_in.residue( 1 ).residue_type_set() );

	//for each residue
	for( core::Size seqpos = 1; seqpos <= pose_in.total_residue(); ++seqpos ){
			//check each residue seperately
			Pose pose( pose_in );
			Residue rsd( pose.residue( seqpos ) ); 

			//first check for unsat polars in residue
			//we only need to try solvating unsat atoms
			vector1< bool > atm_is_unsat( rsd.natoms(), false );
			find_res_unsat_polars( pose, seqpos, atm_is_unsat );

			Real min_dist( 1.5 );
			ResidueTypeSet const & rsd_set( rsd.residue_type_set() );
			ResidueOP wat_rsd( ResidueFactory::create_residue( rsd_set.name_map( "TP5" ) ) ); 
			//turn off solvation score
			Real fa_sol_wt( scorefxn->get_weight( fa_sol ) ); 
			scorefxn->set_weight( fa_sol, 0.0 );

			//do polar hydrogens
			for ( chemical::AtomIndices::const_iterator hnum  = rsd.Hpos_polar().begin(),
				hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {  
				Size const iatom( *hnum );
				if( !atm_is_unsat[ iatom ] ) continue;
				//append water molecule from iatom to water oxygen (atom 1)
				dock_waters_to_atom( pose, scorefxn, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff_ );
			}
			//then do acceptors
			for ( chemical::AtomIndices::const_iterator anum  = rsd.accpt_pos().begin(),
				anume = rsd.accpt_pos().end(); anum != anume; ++anum ){   
				Size const iatom( *anum );
				if( !atm_is_unsat[ iatom ] ) continue;
				//append water molecule from iatom to water 
				dock_waters_to_atom( pose, scorefxn, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff_ );
			}

			//now check again for unsat hbonds in residue
			find_res_unsat_polars( pose, seqpos, atm_is_unsat );
			for( Size iatom = 1; iatom <= atm_is_unsat.size(); ++iatom ){
				if( atm_is_unsat[ iatom ] ) ++all_unsat_polars_;
			}
	}


}


}
}
}
