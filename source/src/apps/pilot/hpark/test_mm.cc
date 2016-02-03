#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>

//#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <fstream>
#include <cmath>

OPT_1GRP_KEY(Integer, fpd, resno)
OPT_1GRP_KEY(String, fpd, prefix)

using namespace core;

utility::vector1< Real > set_chival( Size nchi, Size &nbin )
{
  utility::vector1< Real > chival;

  if( nchi == 1 ){
    Real chibin( 5.0 );
    Size binsize = (int) (360.0/chibin);

    Real chi1 = 0.0 + chibin*nbin;

    chival.push_back( chi1 );
    chival.push_back( 180.0 );
    chival.push_back( 180.0 );
    chival.push_back( 180.0 );

  } else if(nchi == 2 ){
    Real chibin( 5.0 );
    Size binsize = (int) (360.0/chibin);

    Size ichi1 = (nbin%binsize);
    Size ichi2 = (nbin/binsize);
    
    Real chi1 = 0.0 + chibin*ichi1;
    Real chi2 = 0.0 + chibin*ichi2;

    chival.push_back( chi1 );
    chival.push_back( chi2 );
    chival.push_back( 180.0 );
    chival.push_back( 180.0 );
    
  } else if( nchi == 3 ){
    Real chibin( 10.0 );
    Size binsize = (int) (360.0/chibin);

    Size ichi1 = (nbin%binsize);
    Size ichi2 = (nbin/binsize)%binsize;
    Size ichi3 = (nbin/(binsize*binsize));
    
    Real chi1 = 0.0 + chibin*ichi1;
    Real chi2 = 0.0 + chibin*ichi2;
    Real chi3 = 0.0 + chibin*ichi3;

    chival.push_back( chi1 );
    chival.push_back( chi2 );
    chival.push_back( chi3 );
    chival.push_back( 180.0 );

  } else if( nchi == 4 ){
    Real chibin( 15.0 );
    Size binsize = (int) (360.0/chibin);

    Size ichi1 = (nbin%binsize);
    Size ichi2 = (nbin/binsize)%binsize;
    Size ichi3 = (nbin/(binsize*binsize))%binsize;
    Size ichi4 = (nbin/(binsize*binsize*binsize));
    
    Real chi1 = 0.0 + chibin*ichi1;
    Real chi2 = 0.0 + chibin*ichi2;
    Real chi3 = 0.0 + chibin*ichi3;
    Real chi4 = 0.0 + chibin*ichi4;

    chival.push_back( chi1 );
    chival.push_back( chi2 );
    chival.push_back( chi3 );
    chival.push_back( chi4 );
  }

  nbin++;
  return chival;
}

void scan_energy( pose::Pose &pose, 
		  int const &resno,
		  Size const &nchi,
		  scoring::ScoreFunctionOP sfxn,
		  std::string const &prefix )
{

  pose::PoseOP poseOP = new pose::Pose;
  *poseOP = pose;
  conformation::Residue const & ires( pose.residue( resno ) ); 

  Size nbin( 0 );
  Size max_nbin;

  if( nchi == 0 ){
    return;
  } else if( nchi == 1 ){
    max_nbin = (int) (360.0/5.0);
  } else if( nchi == 2 ){
    max_nbin = (int) ((360.0/5.0)*(360.0/5.0));
  } else if( nchi == 3 ){
    max_nbin = (int) ((360.0/10.0)*(360.0/10.0)*(360.0/10.0));
  } else if( nchi == 4 ){
    max_nbin = (int) ((360.0/15.0)*(360.0/15.0)*(360.0/15.0)*(360.0/15.0));
  }
  
  utility::vector1< bool > is_scoringres( pose.total_residue(), false );
  is_scoringres[ resno ] = true;

  // Production

  std::ofstream myfile;
  myfile.open( prefix.c_str() );

  while( nbin <= max_nbin ){
    // Perturb chi angle
    utility::vector1< Real > chival = set_chival( nchi, nbin );
    //for( Size ichi = 1; ichi <= nchi; ++ichi ){
    for( Size ichi = 1; ichi <= nchi; ++ichi ){
      poseOP->set_chi( ichi, resno, chival[ ichi ] );
    }

    // Score
    Real score_dumm = sfxn->score( *poseOP );
    Real score = sfxn->get_sub_score( *poseOP, is_scoringres );

    /*
    core::scoring::EnergyMap emap( sfxn->weights() );
    sfxn->get_sub_score( pose_tmp, is_scoringres, emap );
    */
    
    if( chival[1] <= 0.01 ){ myfile << std::endl; }

    // Report
    //std::cout << prefix;
    
    myfile << std::setw(6) << chival[1] << " ";
    myfile << std::setw(6) << chival[2] << " ";
    myfile << std::setw(6) << chival[3] << " ";
    myfile << std::setw(6) << chival[4] << " ";
    myfile << std::setw(10) << score << std::endl;

    //std::cout << std::setw(10) << emap[ scoring::mm_twist ] << std::endl;
  }
  myfile.close();
}

void
scan_rotamer( pose::Pose const &pose, 
	      scoring::ScoreFunctionOP const &sfxn,
	      Size const &ir ){

  pack::dunbrack::RotamerLibrary const & rotamer_library_ 
    = * pack::dunbrack::RotamerLibrary::get_instance();

  pack::rotamers::SingleResidueRotamerLibraryCAP 
    residue_rotamer_library( rotamer_library_.get_rsd_library(pose.residue( ir ).type()) );

  pack::dunbrack::SingleResidueDunbrackLibraryCAP 
    dun_rotlib( dynamic_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const * >
	       ( residue_rotamer_library.get() ));

  utility::vector1< pack::dunbrack::DunbrackRotamerSampleData > sample_data 
    = dun_rotlib->get_all_rotamer_samples( pose.phi( ir ), pose.psi( ir ) );

  utility::vector1< bool > is_scoringres( pose.total_residue(), false );
  is_scoringres[ ir ] = true;

  for( Size j = 1; j <= sample_data.size(); j++ ){
    
    pose::Pose pose_tmp( pose );

    Size nchis = sample_data[j].nchi();
    pack::dunbrack::Real4 chis = sample_data[j].chi_mean();

    for( Size ichi = 1; ichi <= nchis; ++ichi ){
      pose_tmp.set_chi( ichi, ir, chis[ ichi ] );
    }
    
    sfxn->score( pose_tmp );
    Real score = sfxn->get_sub_score( pose_tmp, is_scoringres );
    
    std::cout << std::setw(6) << j << std::setw(6) << sample_data.size();
    std::cout << " " << std::setw(12) << sample_data[j].probability();
    
    for ( Size jj = 1; jj <= nchis; ++jj ){
      if ( chis[ jj ] < 0.0 ) { chis[ jj ] += 360.0;}
      std::cout << " " << std::setw(10) << chis[ jj ];
    }

    std::cout << " " << std::setw(15) << score;
    std::cout << std::endl;
  }
}

void
scan_neutral( pose::Pose const &pose, 
	      scoring::ScoreFunctionOP sfxn,
	      std::string const &prefix ){
  
  // I, L, F, W, Y, N, D, H
  utility::vector1<std::string> aastr( 8 );
  utility::vector1<std::string> aastr1( 8 );
  aastr[1] = "ILE"; aastr[2] = "LEU";  aastr[3] = "PHE";  aastr[4] = "TRP";  aastr[5] = "TYR";  aastr[6] = "ASN";  aastr[7] = "ASP";  aastr[8] = "HIS";
  aastr1[1] = "I"; aastr1[2] = "L";  aastr1[3] = "F";  aastr1[4] = "W";  aastr1[5] = "Y";  aastr1[6] = "N";  aastr1[7] = "D";  aastr1[8] = "H";

  std::stringstream ss;

  chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

  for( Size ii = 1; ii <= 8; ++ii ){
    conformation::ResidueOP rsd( conformation::ResidueFactory::create_residue( rsd_set->name_map( aastr[ ii ] ) ) );

    pose::Pose pose_tmp( pose );

    pose_tmp.replace_residue( 2, *rsd, false );

    pose_tmp.set_phi( 2, -130.0 );
    pose_tmp.set_psi( 2, 90.0 );

    Size nchi( rsd->nchi() );

    ss.str( std::string() );
    ss.clear();
    ss << prefix << "." << aastr1[ii] << ".N.dat";
    std::string filename = ss.str();
    
    scan_energy( pose_tmp, 2, nchi, sfxn, filename );
  }

}

int main( int argc, char * argv [] )
{
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  //NEW_OPT( fpd::resno, "resno", 0 );
  NEW_OPT( fpd::prefix, "prefix", "" );

  devel::init(argc, argv);

  //Size resno( option[ fpd::resno ]() );
  std::string prefix( option[ fpd::prefix ]() );

  // Score setup

  scoring::ScoreFunctionOP score_in 
    = core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

  chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

  pose::Pose pose;
  import_pose::pose_from_file( pose, *rsd_set, option[ in::file::s ](1) , core::import_pose::PDB_file);

  std::stringstream ss;

  scan_neutral( pose, score_in, prefix );
  return 0;

  // for residue-wise scanning
  for (int resno = 2; resno < pose.total_residue(); ++resno ){
    Size nchi( pose.residue(resno).nchi() );
    if (nchi == 0) continue;

    std::cout << "Chi angle energy scanning on Residue " << resno;
    std::cout << " " << pose.residue(resno).name();
    std::cout << ", nchi " << nchi;
    std::cout << ", phi/psi " << pose.phi(resno) << " " << pose.psi(resno) << std::endl;
    
    std::cout << "Probability and chi angles: " << std::endl;

    scan_rotamer( pose, score_in, resno );
    scan_energy( pose, resno, nchi, score_in, prefix );
  }

  return 0;
}
