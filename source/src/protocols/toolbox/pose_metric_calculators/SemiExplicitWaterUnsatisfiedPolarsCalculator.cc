// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin SemiExplicitWaterUnsatisfiedPolarsCalculator
///
/// @brief
/// How many unsatisfied polars are there?
///
/// @detailed
/// Buried unsatisfied polar hbonds are destabilizing for proteins. It is good to have less.
/// In a study of 2299 high resolution crystal structures of 1.5A or better, there was an average
/// 71 unsatisfied buried polar hbonds. The normalized average (normalized against aa #) was 0.30 (unpublished).
/// To get this piece of code to work, you must first load in your pdb. Then, you need the following lines:
///
///	core::pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator = new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
///	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );
///
///	core::pose::metrics::PoseMetricCalculatorOP unsat_calculator = new protocols::toolbox::pose_metric_calculators::SemiExplicitWaterUnsatisfiedPolarsCalculator("num_hbonds");
///	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "unsat", unsat_calculator );
///
/// This segment of code sets everything up to be used in the calculator. To use this on your protein, you simply need to
/// write the following: pose.print_metric("unsat", "all_unsat_polars");
///
/// @author
/// Chris King
/// Steven Combs - comments
///
/// @last_modified 9/2/2011
/////////////////////////////////////////////////////////////////////////
/// @file   core/pose/metrics/SemiExplicitWaterUnsatisfiedPolarsCalculator.cc
/// @brief  number of hbonds calculator class
/// @author Chris King, templated off of BuriedUnsatisfiedPolarsCalculator by Florian Richter

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/SemiExplicitWaterUnsatisfiedPolarsCalculator.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>

// AUTO-REMOVED #include <numeric/random/random.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

//protocol headers
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>
// AUTO-REMOVED #include <protocols/rigid/RB_geometry.hh>


// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <basic/MetricValue.hh>


#include <cassert>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <numeric/constants.hh>



using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;
using utility::vector1;
using std::string;
using utility::to_string;
using numeric::constants::f::pi;
using numeric::constants::f::pi_2;

static thread_local basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.SemiExplicitWaterUnsatisfiedPolarsCalculator" );

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

SemiExplicitWaterUnsatisfiedPolarsCalculator::SemiExplicitWaterUnsatisfiedPolarsCalculator(
  std::string hbond_calc,
	scoring::ScoreFunctionOP scorefxn,
	core::Real semiexpl_water_cutoff
) :
	hb_database_( core::scoring::hbonds::HBondDatabase::get_database( choose_hbond_parameter_set() ) ),
	all_unsat_polars_( 0 ),
	special_region_unsat_polars_(0),
	semiexpl_water_cutoff_( semiexpl_water_cutoff ),
	name_of_hbond_calc_( hbond_calc ),
	scorefxn_( scorefxn )
{
  atom_unsat_.clear();
  residue_unsat_polars_.clear();
  atom_semiexpl_score_.clear();
  residue_semiexpl_score_.clear();
  special_region_.clear();
  assert_calculators();

}

SemiExplicitWaterUnsatisfiedPolarsCalculator::SemiExplicitWaterUnsatisfiedPolarsCalculator(
  std::string hbond_calc,
	scoring::ScoreFunctionOP scorefxn,
  std::set< core::Size > const & special_region,
	core::Real semiexpl_water_cutoff
) :
	hb_database_( core::scoring::hbonds::HBondDatabase::get_database( choose_hbond_parameter_set() ) ),
	all_unsat_polars_(0),
	special_region_unsat_polars_(0),
	semiexpl_water_cutoff_( semiexpl_water_cutoff ),
	name_of_hbond_calc_( hbond_calc ),
	scorefxn_( scorefxn ),
	special_region_( special_region )
{
  atom_unsat_.clear();
  residue_unsat_polars_.clear();
  atom_semiexpl_score_.clear();
  residue_semiexpl_score_.clear();
  assert_calculators();
}


void
SemiExplicitWaterUnsatisfiedPolarsCalculator::assert_calculators()
{
	if( !CalculatorFactory::Instance().check_calculator_exists( name_of_hbond_calc_ ) ){
		if( name_of_hbond_calc_ != "default" ) TR << "Attention: couldn't find the specified hbond calculator ( " << name_of_hbond_calc_ << " ), instantiating default one." << std::endl;
		name_of_hbond_calc_ = "unsat_calc_default_hbond_calc";
		if( !CalculatorFactory::Instance().check_calculator_exists( name_of_hbond_calc_ ) ){
			CalculatorFactory::Instance().register_calculator( name_of_hbond_calc_, new NumberHBondsCalculator() );
		}
	}
}


void
SemiExplicitWaterUnsatisfiedPolarsCalculator::lookup(
  std::string const & key,
  basic::MetricValueBase * valptr
) const
{

   if ( key == "all_unsat_polars" ) {
     basic::check_cast( valptr, &all_unsat_polars_, "all_unsat_polars expects to return a Size" );
     (static_cast<basic::MetricValue<Size> *>(valptr))->set( all_unsat_polars_ );

   } else if ( key == "special_region_unsat_polars" ) {
     basic::check_cast( valptr, &special_region_unsat_polars_, "special_region_unsat_polars expects to return a Size" );
     (static_cast<basic::MetricValue<Size> *>(valptr))->set( special_region_unsat_polars_ );

   } else if ( key == "atom_unsat" ) {
     basic::check_cast( valptr, &atom_unsat_, "atom_unsat expects to return a id::AtomID_Map< bool >" );
     (static_cast<basic::MetricValue<id::AtomID_Map< bool > > *>(valptr))->set( atom_unsat_ );

   } else if ( key == "atom_semiexpl_score" ) {
     basic::check_cast( valptr, &atom_semiexpl_score_, "atom_semiexpl_score expects to return a id::AtomID_Map< Real >" );
     (static_cast<basic::MetricValue<id::AtomID_Map< Real > > *>(valptr))->set( atom_semiexpl_score_ );

   } else if ( key == "residue_unsat_polars" ) {
     basic::check_cast( valptr, &residue_unsat_polars_, "residue_unsat_polars expects to return a utility::vector1< Size >" );
     (static_cast<basic::MetricValue<utility::vector1< Size > > *>(valptr))->set( residue_unsat_polars_ );

   } else if ( key == "residue_semiexpl_score" ) {
     basic::check_cast( valptr, &residue_semiexpl_score_, "residue_semiexpl_score expects to return a utility::vector1< Real >" );
     (static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set( residue_semiexpl_score_ );

   } else {
     basic::Error() << "SemiExplicitWaterUnsatisfiedPolarsCalculator cannot compute the requested metric " << key << std::endl;
     utility_exit();
   }

} //lookup



std::string
SemiExplicitWaterUnsatisfiedPolarsCalculator::print( std::string const & key ) const
{

  if ( key == "all_unsat_polars" ) {
    return utility::to_string( all_unsat_polars_ );
  } else if ( key == "special_region_unsat_polars" ) {
    return utility::to_string( special_region_unsat_polars_ );
  } else if ( key == "atom_Hbonds" ) {
    basic::Error() << "id::AtomID_Map< bool > has no output operator, for metric " << key << std::endl;
    utility_exit();
  } else if ( key == "residue_unsat_polars" ) {
    return utility::to_string( residue_unsat_polars_ );
  }

  basic::Error() << "SemiExplicitWaterUnsatisfiedPolarsCalculator cannot compute metric " << key << std::endl;
  utility_exit();
  return "";

} //print

//check for clashes based on atom distances
//only check atoms in AtomID vector
//should be way faster than calculating entire score
bool
fast_clash_check(
  Pose const & pose,
  vector1< id::AtomID > const check_atids,
  Real const clash_dist_cut
)
{
  Real const clash_dist2_cut( clash_dist_cut * clash_dist_cut );
  for( Size iatid = 1; iatid <= check_atids.size(); ++iatid ){
    Vector const at1_xyz( pose.xyz( check_atids[ iatid ] ) );
    for( Size res2 = 1; res2 <= pose.total_residue(); ++res2 ){
      for( Size at2 = 1; at2 <= pose.residue( res2 ).natoms(); ++at2 ){
        //skip virtual atoms!
        if( pose.residue( res2 ).atom_type( at2 ).lj_wdepth() == 0.0 ) continue;
        id::AtomID atid2( at2, res2 );
        //skip if atid2 is in check_atids
        bool skip_at2( false );
        for( Size jatid = 1; jatid <= check_atids.size(); ++jatid ){
          if( atid2 == check_atids[ jatid ] ){ skip_at2 = true; break; }
        }
        if( skip_at2 ) continue;
        Real const dist2( at1_xyz.distance_squared( pose.xyz( atid2 ) ) );
        if( dist2 < clash_dist2_cut ){
          //TR_unsat << "CLASH!: " << check_atids[ iatid ] << " - " << atid2 <<
          //   " = " << dist2 << std::endl;
          return true;
        }
      }
    }
  }
  return false;
}



Real
SemiExplicitWaterUnsatisfiedPolarsCalculator::semiexpl_water_hbgeom_score(
  pose::Pose pose,
  scoring::ScoreFunctionOP scorefxn,
  Size seqpos,
  Size atomno,
  conformation::Residue new_rsd,
  Size new_atomno
)
{
  using namespace id;
  using namespace conformation;
  using namespace scoring;
  using namespace scoring::hbonds;
  using namespace kinematics;

  Residue rsd( pose.residue( seqpos ) ); 
  Pose ref_pose( pose );
  //store hbond e before adding water
  /*Real rsd_hbond_e_before(
			pose.energies().residue_total_energies( seqpos )[ hbond_bb_sc ] +
			pose.energies().residue_total_energies( seqpos )[ hbond_sc ]
  );*/

  //append by jump from seqpos atomno to new_rsd atom 1, maybe make random downstream atom?
  pose.append_residue_by_jump( new_rsd, seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), true );

  Size new_seqpos( pose.total_residue() );
  Size jump_number( pose.fold_tree().num_jump() );
  Jump jump( pose.jump( jump_number ) ); 
  //store water atom ids for clash check
  vector1< id::AtomID > clash_check_atids;
  for( Size iat = 1; iat <= new_rsd.natoms(); ++iat ){
    clash_check_atids.push_back( id::AtomID( iat, new_seqpos ) ); 
  }

  //which is acceptor, donor?
  Size aatm( 0 ), acc_pos( 0 ), hatm( 0 ), don_pos( 0 ); 
  bool wat_is_acc( false );
  //water is acceptor
  if( rsd.atom_type( atomno ).is_polar_hydrogen() &&
    new_rsd.atom_type( new_atomno ).is_acceptor() ){
    don_pos = seqpos;
    hatm = atomno;
    acc_pos = new_seqpos;
    aatm = new_atomno;
    wat_is_acc = true;
  }
  //or water is donor
  else if( new_rsd.atom_type( new_atomno ).is_polar_hydrogen() &&
    rsd.atom_type( atomno ).is_acceptor() ){
    don_pos = new_seqpos;
    hatm = new_atomno;
    acc_pos = seqpos;
    aatm = atomno;
  }
  else{ utility_exit_with_message( "ERROR: res " + to_string( seqpos ) + " atom " + to_string( atomno ) +
    " res " + to_string( new_seqpos ) + " atom " + to_string( new_atomno ) + " is not HB don/acc pair!!\n" );
  }

  //now get their base atoms to get datm and batm
  Size datm( pose.residue( don_pos ).atom_base( hatm ) );
  Size batm( pose.residue( acc_pos ).atom_base( aatm ) );
  Size b2atm( pose.residue( acc_pos ).abase2( aatm ) );
  Size dbatm( pose.residue( don_pos ).atom_base( datm ) ); //hpol base2

  //add vrt res so final torsion exists
  chemical::ResidueTypeSet const & rsd_set( rsd.residue_type_set() );
  conformation::ResidueOP vrt_rsd( conformation::ResidueFactory::create_residue( rsd_set.name_map( "VRT" ) ) );
  pose.append_residue_by_jump( *vrt_rsd, pose.total_residue() );
  FoldTree f_jump( pose.fold_tree() );
  //just min the new jump
  //MoveMapOP mm = new MoveMap;
  //mm->set_jump( jump_number, true );
  //protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm, scorefxn, "dfpmin", 0.01, true );

  //new naive fold tree
  FoldTree f_rot( pose.total_residue() );
  //switch to chem bond so can use bond angle defs directly
  f_rot.new_chemical_bond( seqpos, new_seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), pose.total_residue() - 2 );
  pose.fold_tree( f_rot );

	Size water_hb_states_tot( 0 );
	Size water_hb_states_good( 0 );
	hbonds::HBEvalTuple hbe_type( datm, pose.residue( don_pos ), aatm, pose.residue( acc_pos ) );
	//granularity of enumeration
	//Size steps( 7 );
	//must import from hbonds/constants.hh
	// in case you're wondering, we're sampling in angle space instead of cos(angle)
	// space so we sample more equally in polar coord space (not more densely when cos(pi-angle)->1)
	Real AHdist_min(MIN_R), AHdist_max(MAX_R); Size AHdist_steps( 5 );
//  Real cosBAH_min( MIN_xH ), cosBAH_max( MAX_xH ); Size cosBAH_steps( 5 );
//  Real cosAHD_min( MIN_xD ), cosAHD_max( MAX_xD ); Size cosAHD_steps( 5 );
  Real BAHang_min( pi - std::acos( MIN_xH ) ), BAHang_max( pi - std::acos( MAX_xH ) ); Size BAHang_steps( 5 );
  Real AHDang_min( pi - std::acos( MIN_xD ) ), AHDang_max( pi - std::acos( MAX_xD ) ); Size AHDang_steps( 5 );
//  Real B2BAHchi_min( 0 ), B2BAHchi_max( pi_2 ); Size B2BAHchi_steps( 1 );
  Real B2BAHchi_min( 0 ), B2BAHchi_max( pi_2 ); Size B2BAHchi_steps( 10 );
  Real BAHDchi_min( 0 ), BAHDchi_max( pi_2 ); Size BAHDchi_steps( 3 );
	//gimbal lock… <sigh>
  for( Real AHdist = AHdist_min; AHdist <= AHdist_max;
      AHdist += (AHdist_max - AHdist_min)/static_cast<Real>(AHdist_steps-1)){
//		for( Real cosBAH = cosBAH_min - 0.0001; cosBAH <= cosBAH_max;
//				cosBAH += (cosBAH_max - cosBAH_min)/static_cast<Real>(cosBAH_steps-1)){
//			for( Real cosAHD = cosAHD_min - .0001; cosAHD <= cosAHD_max;
//					cosAHD += (cosAHD_max - cosAHD_min)/static_cast<Real>(cosAHD_steps-1)){
		for( Real BAHang = BAHang_min - 0.0001; BAHang <= BAHang_max;
				BAHang += (BAHang_max - BAHang_min)/static_cast<Real>(BAHang_steps-1)){
			for( Real AHDang = AHDang_min - .0001; AHDang <= AHDang_max;
					AHDang += (AHDang_max - AHDang_min)/static_cast<Real>(AHDang_steps-1)){
				//this loop should be longer than the water internal orientation loops
				for( Real B2BAHchi = B2BAHchi_min; B2BAHchi <= B2BAHchi_max;
						B2BAHchi += (B2BAHchi_max - B2BAHchi_min)/static_cast<Real>(B2BAHchi_steps-1)){
					for( Real BAHDchi = BAHDchi_min; BAHDchi <= BAHDchi_max;
							BAHDchi += (BAHDchi_max - BAHDchi_min)/static_cast<Real>(BAHDchi_steps-1)){

//						Real BAHang( pi - std::acos( cosBAH ) );
//						Real AHDang( pi - std::acos( cosAHD ) );
						Real cosBAH( std::cos( pi - BAHang ) );
						Real cosAHD( std::cos( pi - AHDang ) );
            //call hbonds::hbond_compute_energy( ... ) with these angles
            // and skip if is zero hb energy
						// this allows computing exactly what frac of actual HB phase space
						// can be filled w/ a favorable water molecule
						Real hb_energy( 0.0 );
						scoring::hbonds::hbond_compute_energy( *( hb_database_ ),
							scorefxn->energy_method_options().hbond_options(),
							hbe_type, AHdist, cosAHD, cosBAH, B2BAHchi, hb_energy );
						//TR << hb_energy << std::endl;
						if( hb_energy >= 0.0 ) continue; //was not actually an hbond

						++water_hb_states_tot;

            //reset chem  bond ftree
            pose.fold_tree( f_rot );

						//Real water_ang( 0.0 );
						//store water internal bond angle
						//WARNING: hard-coded for TP3 water (H1-O-H2)
//            if( wat_is_acc ){
//							water_ang = pose.conformation().bond_angle( AtomID( 2, acc_pos ),
//									AtomID( 1, acc_pos  ),
//									AtomID( 3, acc_pos  ) );
//						}
            pose.conformation().set_bond_angle( AtomID( batm, acc_pos ),
                AtomID( aatm, acc_pos  ),
                AtomID( hatm, don_pos  ),
                BAHang );
//            if( wat_is_acc ){
//							//need to fix water internal bond angle
//							pose.conformation().set_bond_angle( AtomID( 2, acc_pos ),
//									AtomID( 1, acc_pos  ),
//									AtomID( 3, don_pos  ),
//									water_ang );
//						}

            pose.conformation().set_bond_angle( AtomID( aatm, acc_pos ),
                AtomID( hatm, don_pos  ),
                AtomID( datm, don_pos  ),
                AHDang );

            pose.conformation().set_torsion_angle( AtomID( batm, acc_pos ),
                AtomID( aatm, acc_pos  ),
                AtomID( hatm, don_pos  ),
                AtomID( datm, don_pos  ),
                BAHDchi );
            if( wat_is_acc ){
              //need to redefine hbond chi torsion if water acceptor (hbond chi undefined)
              pose.conformation().set_torsion_angle( AtomID( dbatm, don_pos ),
                  AtomID( datm, don_pos  ),
                  AtomID( hatm, don_pos  ),
                  AtomID( aatm, acc_pos  ),
                  B2BAHchi );
            }
            else{
              pose.conformation().set_torsion_angle( AtomID( b2atm, acc_pos ),
                  AtomID( batm, acc_pos  ),
                  AtomID( aatm, acc_pos  ),
                  AtomID( hatm, don_pos  ),
                  B2BAHchi );
            }

            pose.conformation().set_bond_length( AtomID( aatm, acc_pos ),
                AtomID( hatm, don_pos  ),
                AHdist );

            //do fast clash check, OH hbonds are only 0.8A!
						if( fast_clash_check( pose, clash_check_atids, 0.8 ) ) continue;

            pose.fold_tree( f_jump );
            //minimize the new jump, too slow!
            //min_mover->apply( pose );
            scorefxn->score( pose );

            //get score
            Real wat_score( pose.energies().residue_total_energies( new_seqpos ).dot( scorefxn->weights() ) );

						if( wat_score <= 0.0 ){
							++water_hb_states_good;
							//TR_unsat << "AHdist: " << AHdist << " ";
							//TR_unsat << "\tBAHang: " << BAHang << " ";
							//TR_unsat << "\tBAHDchi: " << BAHDchi << " ";
							//TR_unsat << "\twat_score: " << wat_score << std::endl;
							//TR_unsat << pose.energies().total_energies().weighted_to_string( scorefxn->weights() )
							//	+ " total_score: " + to_string( pose.energies().total_energies()[ total_score ] );
							//pose.dump_pdb( "watest." + to_string( seqpos ) + "." + to_string( atomno ) + "." + to_string( Size( numeric::conversions::degrees( AHDang ) ) ) + "." + to_string( Size( numeric::conversions::degrees( BAHang ) ) ) + "." + to_string( Size( numeric::conversions::degrees( BAHDchi ) ) ) + "." + to_string( Size( numeric::conversions::degrees( B2BAHchi ) ) ) + ".pdb" );
						}
          }
        }
      }
    }
  }
	return ( static_cast< Real >( water_hb_states_good ) /
		static_cast< Real >( water_hb_states_tot ) );
}

/// @brief this should just be caled "compute"
/// @brief for non-hbonded polar atoms, attempt docking single water residues, then count number still unsatisfied
void
SemiExplicitWaterUnsatisfiedPolarsCalculator::recompute( Pose const & in_pose )
{

	all_unsat_polars_ = 0;
	special_region_unsat_polars_ = 0;
	//Real min_dist( core::scoring::hbonds::MIN_R );
	//Real shell_cutoff( core::scoring::hbonds::MAX_R );
	chemical::ResidueTypeSet const & rsd_set( in_pose.residue( 1 ).residue_type_set() );
	conformation::ResidueOP wat_rsd( conformation::ResidueFactory::create_residue( rsd_set.name_map( "TP3" ) ) );
	Size wat_O_at( 1 ); //warning! hardcoded for TP3 water!
	Size wat_H1_at( 2 ); //warning! hardcoded for TP3 water!
	//turn off solvation score
	scoring::ScoreFunctionOP scorefxn( scorefxn_ );
	scorefxn->set_weight( scoring::fa_sol, 0.0 );


	if( in_pose.total_residue() != residue_unsat_polars_.size() ){
			residue_unsat_polars_.resize( in_pose.total_residue() );
			atom_unsat_.resize( in_pose.total_residue() );
	}

  basic::MetricValue< id::AtomID_Map< Size > > atom_hbonds_dry;
  in_pose.metric( name_of_hbond_calc_, "atom_Hbonds", atom_hbonds_dry );

  for( Size i = 1; i <= in_pose.total_residue(); ++i ){

		residue_unsat_polars_[i] = 0;
    conformation::Residue const & rsd = in_pose.residue( i );

		for( Size at = 1; at <= rsd.natoms(); ++at){
			//reset pose after each atom
			pose::Pose pose( in_pose );

			core::id::AtomID atid( at, i );
			bool this_atom_unsat(false);

			//counting acceptors and donors
			if( !( rsd.atom_type( at ).is_acceptor() || rsd.atom_type( at ).is_donor() ) ){
				atom_unsat_.set( atid, this_atom_unsat );
				continue;
			}

			//how about instead we iter over heavy atoms, then simply check all H's is is donor
			//can just avg al child H's water hbond fraction
			//just add this vale into the bonded+hbonds value to gauge 
			//this is not the best way, but is the way that is comparable w/ BuriedUnsatCalc...
			Size satisfac_cut = satisfaction_cutoff( rsd.type().atom_type( at ).name() );
			Size bonded_heavyatoms = rsd.n_bonded_neighbor_all_res( at ) - rsd.type().number_bonded_hydrogens( at );
			Size n_hbonds_needed( satisfac_cut - bonded_heavyatoms );

			//check for unsat hbonds first
			Size n_atom_hbonds( atom_hbonds_dry.value()[ atid ] );
			if( n_atom_hbonds >= n_hbonds_needed ){
				atom_unsat_.set( atid, this_atom_unsat );
				continue;
			}

			//could atom make enough water hbonds?
			Real semiexpl_water_score( 0.0 );
			if( rsd.atom_type( at ).is_acceptor() ){
				//hmmmm, acc score should count for > max 1
				//n_hbonds_needed is a good proxy
				semiexpl_water_score += semiexpl_water_hbgeom_score(
					pose, scorefxn, i, at, *wat_rsd, wat_H1_at );
				TR << rsd.name3() << " " << to_string( i ) << " " << rsd.atom_name( at ) <<
					" semiexpl_wat_score: " << semiexpl_water_score << std::endl;
				if( semiexpl_water_score >= semiexpl_water_cutoff_ ){
					n_atom_hbonds += n_hbonds_needed;
				}
			}
			if( rsd.atom_type( at ).is_donor() ){
				//H atom score should only count for 1 hbond
				//but to mimic bur unsat hbonds logic, make it so that
				//we just avg hbgeom score and see if above cutoff
				Size n_bonded_H = rsd.type().attached_H_end( at ) - rsd.type().attached_H_begin( at );
				for( Size h_at = rsd.type().attached_H_begin( at );
						h_at<= rsd.type().attached_H_end( at ); h_at++){
					//get avg wat score of all attached H's
					semiexpl_water_score += ( semiexpl_water_hbgeom_score(
								pose, scorefxn, i, h_at, *wat_rsd, wat_O_at ) / n_bonded_H );
				}
				TR << rsd.name3() << " " << to_string( i ) << " " << rsd.atom_name( at ) <<
					" semiexpl_wat_score: " << semiexpl_water_score << std::endl;
				if( semiexpl_water_score >= semiexpl_water_cutoff_ ){
					n_atom_hbonds += n_hbonds_needed;
				}
			}

			//debug				
			//TR << "HBONDS_WET: res\t" << i << "\t" << rsd.atom_name( at ) << "\t" << atom_hbonds_wet.value()[ atid ] << std::endl;
			//store semiexpl score
			atom_semiexpl_score_.set( atid, semiexpl_water_score );
			residue_semiexpl_score_[ i ] += semiexpl_water_score;

			//now check if still unsat
			if( n_atom_hbonds < n_hbonds_needed ){

				all_unsat_polars_++;
				residue_unsat_polars_[i]++;
				this_atom_unsat = true;

				if( special_region_.find( i ) != special_region_.end() ) special_region_unsat_polars_++;

				//debug				
				//TR << "UNSAT_POLARS: res\t" << i << "\t" << rsd.atom_name( at ) << "\t";
				//if( this_atom_unsat ) TR << 1 << std::endl;
				//else TR << 0 << std::endl;
			}
			atom_unsat_.set( atid, this_atom_unsat );

		} //atom
	} //residue
} //recompute

//sum of n bonded heavyatoms and n hbonds must reach this val
core::Size
SemiExplicitWaterUnsatisfiedPolarsCalculator::satisfaction_cutoff( std::string atom_type )
{

	//according to jk, buried hydroxyls are often seen making only one hydrogen bond. also, ether oxygens often are bad h-bond acceptors
	if( atom_type == "OH" ) return 2;

	//backbone oxygens also only have one h-bbond in most secondary structure elements
	else if (atom_type == "OCbb") return 2;

	else if( atom_type ==  "S")	return 2;

	//everything else we expect to have 3 bonded/h-bonded neighbours to count as satisfied
	else return 3;


}

} //namespace pose_metric_calculators
} //namespace toolbox
} //namespace protocols
