#include <protocols/mpi_refinement/StructAvrgMover.hh>
#include <protocols/mpi_refinement/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/relax/AtomCoordinateCstMover.hh>

#include <core/scoring/rms_util.hh>
#include <numeric/model_quality/rms.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

using namespace core;

namespace protocols{
namespace mpi_refinement{

static basic::Tracer TR("MPI.LHR.StructAvrger");

// Constructor
StructAvrgMover::StructAvrgMover( pose::Pose const &pose,
				  protocols::wum::SilentStructStore const &decoys,
				  bool const minimize )
{
  set_default();
  minimize_ = minimize;

  TR << "Received " << decoys.size() << " decoys as input." << std::endl;

  // convert silentstore to pose vectors
  pose::Pose ref_pose( pose ); // copy
  poses_.resize( 0 );
  decoys.get_pose( 0, ref_pose );

  protocols::wum::SilentStructStore::const_iterator iter;
  for( iter = decoys.begin(); iter != decoys.end(); ++iter ){
    pose::Pose pose_tmp( ref_pose );
    (*iter)->fill_pose( pose_tmp );
    poses_.push_back( pose_tmp );
  }
}

StructAvrgMover::~StructAvrgMover(){}

void
StructAvrgMover::set_default()
{
  kT_ = 1.0;
  shave_frac_ = 0.0;
  niter_ = 1;
  mulfactor_ = 1.8;
  sfxn_ = scoring::ScoreFunctionFactory::create_score_function( "talaris2013" );
  sfxn_->set_weight( scoring::cart_bonded, 0.5 );
  // Very strong weight to fix the structure!
  sfxn_->set_weight( scoring::coordinate_constraint, 10.0 );
}

void
StructAvrgMover::apply( pose::Pose &pose )
{

  // Get coordinate fluctuations
  pose::Pose ref_pose( pose );

  utility::vector1< utility::vector1< Real > > deviation;
  deviation.resize( ref_pose.total_residue() );
  for( Size i_pose = 1; i_pose <= poses_.size(); ++i_pose ){
    add_deviations( ref_pose, poses_[i_pose], deviation );
  }

  CAvar_ = calculate_variations( deviation );

  // Iterate by shaving outliers; by default niter = 1, no shaving
  for( core::Size iter = 1; iter <= niter_; ++iter ){
    TR.Debug << "Avrg, iter " << iter << std::endl;
    pose = weighted_average( poses_, sfxn_, ref_pose, 
			     CAvar_, true );

    if( shave_frac_ > 0.0 ) shave_poses( poses_, pose, shave_frac_ );
  }
}

utility::vector1< std::pair< Size, Size > >
StructAvrgMover::predict_region( utility::vector1< core::Real > const CAvar, 
				 utility::vector1< bool > &is_region,
				 core::Real const min_fluc,
				 core::Real const frac_base ) const
{

  // First export to std::vector to measurely sort
  std::vector< core::Real > CAvar_sortable;
  for( Size i = 1; i <= CAvar.size(); ++i ){
    CAvar_sortable.push_back( CAvar[i] );
  }

  // Sort and get baseline for overall fluctuation
  std::sort( CAvar_sortable.begin(), CAvar_sortable.end() );
  Size imax_base = (Size)(frac_base*CAvar.size()) - 1;

  core::Real fluc = 0.0;
  for( Size i = 0; i < imax_base; ++i ) fluc += CAvar_sortable[i];

  fluc /= (core::Real)(imax_base);
  core::Real const fluc_cut = mulfactor_*fluc;

  TR << "RMSF: Fluc cut is set as " << fluc_cut << std::endl;

  // Scan over residue, and set as variable region 
  // if satisfy both condition CArmsd > fluc*mulfactor and CArmsd > min_fluc
  is_region.resize( CAvar.size(), false );
  for( Size ires = 1; ires <= CAvar.size(); ++ires ){
    if( CAvar[ires] > min_fluc && CAvar[ires] > fluc_cut ) is_region[ires] = true;
  }

  utility::vector1< std::pair< Size, Size > > var_region;

  // 1. Sandwich gap and then Remove outlier peak
  {
    utility::vector1< bool > is_region_org( is_region );
    for( Size ires = 2; ires <= CAvar.size()-1; ++ires ){
      if( !is_region_org[ires] && is_region_org[ires-1] && is_region_org[ires+1] )
	is_region[ires] = true;
    }
  }

  {
    utility::vector1< bool > is_region_org( is_region );
    for( Size ires = 2; ires <= CAvar.size()-1; ++ires ){
      if( is_region_org[ires] && !is_region_org[ires-1] && !is_region_org[ires] )
	is_region[ires] = false;
    }
  }
    
  //2. Extend by 2res
  utility::vector1< bool > is_region_org( is_region );
  for( Size ires = 3; ires <= CAvar.size()-2; ++ires ){
    if( is_region_org[ires] ){
      is_region[ires-2] = true; is_region[ires-1] = true;
      is_region[ires+1] = true; is_region[ires+2] = true;
    }
  }

  //3.Split into region
  Size var_last = -1;
  for( Size ires = 1; ires <= CAvar.size(); ++ires ){
    if( !is_region[ires] ) continue;
    
    Size i_seg = var_region.size();
    if( ires - var_last > 1 ){ // Start over with a new piece
      var_region.push_back( std::pair< Size, Size > (ires, ires) );
    } else { // Extension
      var_region[i_seg].second = ires;
    }
    var_last = ires;
  }

  // 4. update 
  for( Size ires = 1; ires <= CAvar.size(); ++ires ) is_region[ires] = false;
  for( Size ireg = 1; ireg <= var_region.size(); ++ireg ){
    for( Size ires = var_region[ireg].first; ires <= var_region[ireg].second; ++ires )
      is_region[ires] = true;
  }

  return var_region;
}

void
StructAvrgMover::report_dev( core::pose::Pose const & ref_pose ) const 
{

  using namespace ObjexxFCL::format;

  //bool report_native( false );
  //if( native_given_ ) report_native = true;
  utility::vector1< bool > is_region;
  utility::vector1< std::pair< Size, Size > > var_region = predict_region( CAvar_, is_region );

  // report region
  core::Size nres_region( 0 );
  TR << "RMSF: Variable region: ";
  for( core::Size ireg = 1; ireg <= var_region.size(); ++ireg ){
    Size startres = ref_pose.pdb_info()->number( var_region[ireg].first );
    Size endres = ref_pose.pdb_info()->number( var_region[ireg].second );
    Size dres = var_region[ireg].second - var_region[ireg].first;
    // prevent from weird Cterms...
    if( var_region[ireg].second == ref_pose.total_residue() && endres == 1 )
      endres = ref_pose.pdb_info()->number( var_region[ireg].second - 1 );

    if( startres > endres ) endres = startres + dres - 1;
    TR << " " << startres << "-" << endres;
    nres_region += endres - startres + 1;
  }
  TR << std::endl;
  TR << "RMSF: Frac var region: " << I(4,nres_region) << " " << I(4,CAvar_.size());
  TR << " " << F(8,4,(core::Real)(nres_region)/(core::Real)(CAvar_.size()) ) << std::endl;

  // report full log
  TR << "RMSF: " << std::setw(4) << "Res" << std::setw(4) << "PDB";
  TR << std::setw(4) << "Var" << std::setw(9) << "Fluc";
  TR << std::endl;
  for( Size ires = 1; ires <= CAvar_.size(); ++ires ){
    Size resno = ref_pose.pdb_info()->number( ires );
    // sometimes weird thing can happend at Cterm...
    if( ires == ref_pose.total_residue() && resno == 1 ) continue;
    std::string is_reg; is_region[ires] ? is_reg = "1" : is_reg = "0";
    TR << "RMSF:" << I(4,ires) << " ";
    TR << I(4,resno) << "   " << is_reg;
    TR << " " << F(8,3,CAvar_[ires]);
    TR << std::endl;
  }
}

utility::vector1< core::Real >
StructAvrgMover::calculate_variations( utility::vector1< utility::vector1< core::Real > > const deviation ){
  utility::vector1< core::Real > CAvars( deviation.size(), 0.0 );

  for( Size ires = 1; ires <= deviation.size(); ++ires ){
    utility::vector1< core::Real > const &dev_res = deviation[ires];

    if( dev_res.size() == 0 ) continue;

    core::Real rmsd( 0.0 );
    for( Size i = 1; i <= dev_res.size(); ++i ){
      rmsd += dev_res[i]*dev_res[i];
    }

    rmsd /= (core::Real)(dev_res.size());
    rmsd = std::sqrt( rmsd );

    CAvars[ires] = rmsd;
  }

  // Smoothen
  core::Size const winsize( 9 );
  Size const pos_shift( (winsize - 1)/2 );
  utility::vector1< core::Real > const CAvars0( CAvars );
  
  utility::vector1< core::Real > w( winsize );
  for( Size i_pos = 1; i_pos <= pos_shift+1; ++i_pos ){
    w[winsize-i_pos] = i_pos*0.1;
    w[i_pos] = i_pos*0.1;
  }

  //Clean
  CAvars.resize( CAvars0.size(), 0.0 );

  for( Size i_res = 1; i_res <= CAvars.size(); ++i_res ){
    core::Real valsum( 0.0 );
    core::Real wsum( 0.0 );
    for( Size i_w = 1; i_w <= winsize; ++i_w ){
      Size const resno( i_res+i_w-pos_shift-1);
      if( resno < 1 || resno > CAvars.size() ) continue;

      core::Real wval = w[i_w]*CAvars0[resno];
      valsum += wval;
      wsum += w[i_w];
    }
    CAvars[i_res] = valsum/wsum;
  } 

  return CAvars;
}

void
StructAvrgMover::add_deviations( pose::Pose ref_pose, 
				 pose::Pose pose,
				 utility::vector1< utility::vector1< core::Real > > &deviation
				 )
{
  // First superimpose
  scoring::calpha_superimpose_pose( pose, ref_pose );

  for( Size ires = 1; ires <= ref_pose.total_residue(); ++ ires ){
    if( !ref_pose.residue( ires ).has(" CA ") ) continue;

    Vector const &xyz1 = ref_pose.residue(ires).xyz( "CA" );
    Vector const &xyz2 = pose.residue(ires).xyz( "CA" );
    Real const dis = xyz1.distance( xyz2 );
    deviation[ires].push_back( dis );
  }

}

pose::Pose
StructAvrgMover::weighted_average( utility::vector1< pose::Pose > &poses,
				   scoring::ScoreFunctionCOP sfxn,
				   pose::Pose const &pose_ref,
				   utility::vector1< Real > const , //CAvar,
				   bool const weighted
				   )
{

  pose::Pose pose_avrg( pose_ref );
  std::map< id::AtomID, Vector > avrg_crd; 

  utility::vector1< std::string > atoms_copy;
  atoms_copy.push_back( " N  " );
  atoms_copy.push_back( " CA " );
  atoms_copy.push_back( " C  " );
  atoms_copy.push_back( " O  " );
  atoms_copy.push_back( " H  " );
  atoms_copy.push_back( " HA " );
  atoms_copy.push_back( "1HA " );
  atoms_copy.push_back( "2HA " );
  atoms_copy.push_back( " OXT" );
  //atoms_copy.push_back( "CB" );
  
  //First get atomIDs available
  utility::vector1< id::AtomID > ids;
  for( Size ires = 1; ires <= pose_ref.total_residue(); ++ires ){
    conformation::Residue const &rsd = pose_avrg.residue(ires);
    for( Size iatm = 1; iatm <= atoms_copy.size(); ++iatm ){
      if( rsd.has( atoms_copy[iatm] ) ){
	id::AtomID id( rsd.atom_index(atoms_copy[iatm]), ires );
	ids.push_back( id );
	avrg_crd[id] = Vector( 0.0 );
      }
    }
  }

  utility::vector1< Real > weights( poses.size(), (Real)(1.0/poses.size()) );
  if( weighted ){
    // get min/max E
    std::vector< Real > scores;
    utility::vector1< Real > scores_cp;
    for( Size ipose = 1; ipose <= poses.size(); ++ipose ){
      Real score = sfxn->score( poses[ipose] );
      scores.push_back( score );
      scores_cp.push_back( score );
    }

    std::sort( scores.begin(), scores.end() );
    
    Real const Emin( scores[0] ); 
    Real const denominator( scores[ scores.size() ] - Emin );

    // Then get normalized weight
    Real prob_sum( 0.0 );
    for( Size ipose = 1; ipose <= scores_cp.size(); ++ipose ){
      Real const Ess = scores_cp[ipose];
      
      Real const normalizedE( -(Ess - Emin)/denominator ); // Emin = 0, Emax = -1
      Real const prob = exp( normalizedE/kT_ );
      prob_sum += prob;
      weights[ipose] = prob;
    }

    for( Size i = 1; i <= weights.size(); ++ i ){
      weights[i] /= prob_sum;
      TR.Debug << "Set pose " << std::setw(4) << i << " with score ";
      TR.Debug << std::setw(8) << scores_cp[i] << " weight as: ";
      TR.Debug << weights[i] << std::endl;
    }
  }

  // Finally 
  for( Size ipose = 1; ipose <= poses.size(); ++ipose ){
    pose::Pose &pose_tmp = poses[ipose]; 

    // Superimpose to ref
    scoring::calpha_superimpose_pose( pose_tmp, pose_ref );

    for( Size i_id = 1; i_id <= ids.size(); ++i_id ){
      id::AtomID const &id = ids[i_id];
      // H
      std::string atmname = pose_tmp.residue( id.rsd() ).atom_name( id.atomno() );
      Vector const &xyz = pose_tmp.xyz( id );

      avrg_crd[id] += weights[ipose] * xyz;
    }
  }

  // Export into pose_avrg
  for( Size i_id = 1; i_id <= ids.size(); ++ i_id ){
    id::AtomID const &id = ids[i_id];
    pose_avrg.set_xyz( id, avrg_crd[id] );
  }

  // Get closest pose to borrow sidechain 
  //std::vector< Real > rmsds;
  utility::vector0< Real > rmsds;
  utility::vector1< Real > rmsds_cp;

  for( Size i_pose = 1; i_pose <= poses.size(); ++i_pose ){
    Real const rmsd = scoring::CA_rmsd( poses[i_pose], pose_avrg );
    rmsds.push_back( rmsd );
    rmsds_cp.push_back( rmsd );
  }

  // Sort to get rmsd cut
  std::sort( rmsds.begin(), rmsds.end() );
  Size closest( 1 );

  // just doing explicitly because below breaks due to crazy NULL definition... 
  //Size closest = (core::Size)(rmsds_cp.index_of( rmsds[0] ));
  for( core::Size i = 1; i <= rmsds_cp.size(); ++i ){
    if( rmsds_cp[i]-rmsds[0] < 1.0e6 ){
      closest = i;
      break;
    }
  }
  pose::Pose &pose_close = poses[closest];

  // Copy
  for( Size ires = 1; ires <= pose_avrg.total_residue(); ++ires ){
    Size start = pose_avrg.residue(ires).type().first_sidechain_atom();
    for( Size iatm = start; iatm <= pose_avrg.residue(ires).natoms(); ++iatm ){
      id::AtomID id( iatm, ires );
      pose_avrg.set_xyz( id, pose_close.residue( ires ).xyz( iatm ) );
    }
  }

  // minimize
  if( minimize_ ){
    protocols::relax::AtomCoordinateCstMover coord_cst_mover;
    core::scoring::ScoreFunctionOP scorefxn_loc = sfxn->clone();
    if( (*scorefxn_loc)[ core::scoring::cart_bonded ] == 0.0 )
      scorefxn_loc->set_weight( core::scoring::cart_bonded, 0.5 );
    if( (*scorefxn_loc)[ core::scoring::pro_close ] > 0.0 )
      scorefxn_loc->set_weight( core::scoring::pro_close, 0.0 );

    scorefxn_loc->set_weight( core::scoring::coordinate_constraint, 10.0 ); //10.0 will a
    coord_cst_mover.cst_sd( 1.0 );
    coord_cst_mover.apply( pose_avrg );
    ramp_minpack_pose( pose_avrg, scorefxn_loc ); // this is cartmin!
  }
  
  return pose_avrg;
}

void
StructAvrgMover::shave_poses( utility::vector1< pose::Pose > &poses,
			      pose::Pose const &avrg_pose,
			      Real const frac
			      ){

  std::vector< Real > rmsds;
  utility::vector1< Real > rmsds_cp;

  for( Size i_pose = 1; i_pose <= poses.size(); ++i_pose ){
    Real const rmsd = scoring::CA_rmsd( poses[i_pose], avrg_pose );
    rmsds.push_back( rmsd );
    rmsds_cp.push_back( rmsd );
  }

  // Sort to get rmsd cut
  std::sort( rmsds.begin(), rmsds.end() );
  Size const icut = (Size)( (1.0 - frac)*rmsds.size() );
  Real const rmsdcut = rmsds[icut];

  utility::vector1< pose::Pose > poses_shaved;
  for( Size i_pose = 1; i_pose <= poses.size(); ++i_pose ){
    if( rmsds_cp[i_pose] < rmsdcut )
      poses_shaved.push_back( poses[i_pose] );
  }

  TR << "N? " << poses.size() << " " << poses_shaved.size() << std::endl;
  poses = poses_shaved;
}


} 
}
