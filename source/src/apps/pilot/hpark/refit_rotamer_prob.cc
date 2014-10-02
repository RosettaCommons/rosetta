#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

//#include <iomanip>
#include <iostream>
#include <fstream>

OPT_1GRP_KEY(Boolean, fpd, verb)
OPT_1GRP_KEY(Boolean, fpd, rotscan)

namespace myspace{

using namespace core;

struct DataStruct
{
  core::Real phi;
  core::Real psi;
  std::string aa;
  bool sized;
  bool scored;
  utility::vector1< utility::vector1< core::Real > > scores;
  utility::vector1< pack::dunbrack::Real4 > chis;
  utility::vector1< Real > prob;
  utility::vector1< Size > nclash;
};

Real 
periodic( Real ang ){
  if( ang > 180.0 ){
    return ang - 360.0;
  } else {
    return ang;
  }
}

Real 
nonperiodic( Real ang ){
  if( ang < 0.0 ){
    return ang + 360.0;
  } else {
    return ang;
  }
}

Size 
find_nearest_rot( utility::vector1< Real > chi,  
		  utility::vector1< pack::dunbrack::Real4 > rotchis ){
  
  Size minid( 0 );
  Size nchi( chi.size() );
  Real mindev( 1.0e8 );
  for( Size i = 1; i <= rotchis.size(); ++i ){
    Real dev( 0.0 );
    for( Size j = 1; j <= nchi; ++j ){
      Real chidiff = chi[j] - rotchis[i][j];
      if( chidiff > 360.0 ){
	chidiff -= 360.0;
      } else if (chidiff < -360.0){
	chidiff += 360.0;
      }
      dev += chidiff*chidiff;
    }
    if( dev < mindev ){ 
      minid = i;
      mindev = dev;
    }
  }

  return minid;
}

Real 
stdev( utility::vector1< Real > const & vals,
       Real const & mean ){

  Real tmp = 0.0;
  for( Size i = 1; i <= vals.size(); ++i ){
    tmp += (vals[i] - mean)*(vals[i] - mean);
  }
  tmp /= vals.size();
  return std::sqrt( tmp );
}

void
report_data( utility::vector1< DataStruct > const &data,
	     bool report_verb ){
  
  Size n( data.size() );
  std::string aa_prv = "";
  // Iter over phi/psi bin
  for( Size i = 1; i <= n; ++i ){
    DataStruct const &datum ( data[i] );
    if( datum.aa.length() < 3 || (!report_verb && !datum.scored)) continue;

    Size nrot( datum.scores.size() );
    
    if( aa_prv.compare( datum.aa ) != 0 ){
      std::cout << "RESIDUE " << datum.aa << std::endl;
    }

    std::cout << std::endl;

    for( Size j = 1; j <= nrot; ++j ){
      Size const ndata( datum.scores[j].size() );

      Real score_mean( 0.0 );
      Real minscore( 999.0 );
      Real std( 0.0 );

      for( Size k = 1; k <= ndata; ++k ){
	score_mean += datum.scores[j][k];
	minscore = std::min( datum.scores[j][k], minscore );
      }
      if( ndata > 0 ) {
	score_mean /= ndata;
	std = stdev( datum.scores[j], score_mean );
      }

      //Real zcut = score_mean + 2.0*stdev( datum.scores[j], score_mean );

      std::cout << " " << std::setw(6) << i;
      std::cout << " " << std::setw(3) << datum.aa;
      std::cout << " " << std::setw(6) << std::setprecision(4) << datum.phi;
      std::cout << " " << std::setw(6) << std::setprecision(4) << datum.psi;
      std::cout << " " << std::setw(6) << ndata;
      std::cout << " " << std::setw(4) << j;
      std::cout << " " << std::setw(12) << std::setprecision(5) << score_mean;
      std::cout << " " << std::setw(12) << std::setprecision(5) << minscore;
      std::cout << " " << std::setw(12) << std::setprecision(3) << std;
      std::cout << " " << std::setw(8) << std::setprecision(3) << datum.prob[j];
      std::cout << " " << std::setw(6) << std::setprecision(3) << datum.chis[j][1];
      std::cout << " " << std::setw(6) << std::setprecision(3) << datum.chis[j][2];
      std::cout << " " << std::setw(6) << std::setprecision(3) << datum.chis[j][3];
      std::cout << " " << std::setw(6) << std::setprecision(3) << datum.chis[j][4];
      std::cout << std::endl;
    }
    aa_prv = datum.aa;
  }
}

Size aa_to_id( std::string const &aa ){
  if      (aa.compare("CYS") == 0){ return 1; }
  else if (aa.compare("ASP") == 0){ return 2; }
  else if (aa.compare("GLU") == 0){ return 3; }
  else if (aa.compare("PHE") == 0){ return 4; }
  else if (aa.compare("HIS") == 0){ return 5; }
  else if (aa.compare("ILE") == 0){ return 6; }
  else if (aa.compare("LYS") == 0){ return 7; }
  else if (aa.compare("LEU") == 0){ return 8; }
  else if (aa.compare("MET") == 0){ return 9; }
  else if (aa.compare("ASN") == 0){ return 10; }
  else if (aa.compare("PRO") == 0){ return 11; }
  else if (aa.compare("GLN") == 0){ return 12; }
  else if (aa.compare("ARG") == 0){ return 13; }
  else if (aa.compare("SER") == 0){ return 14; }
  else if (aa.compare("THR") == 0){ return 15; }
  else if (aa.compare("VAL") == 0){ return 16; }
  else if (aa.compare("TRP") == 0){ return 17; }
  else if (aa.compare("TYR") == 0){ return 18; }
  else { return 0; }
}

std::string 
id_to_aa( Size const &id ){
  if      (id == 1 ){ return "CYS"; }
  else if (id == 2 ){ return "ASP"; }
  else if (id == 3 ){ return "GLU"; }
  else if (id == 4 ){ return "PHE"; }
  else if (id == 5 ){ return "HIS"; }
  else if (id == 6 ){ return "ILE"; }
  else if (id == 7 ){ return "LYS"; }
  else if (id == 8 ){ return "LEU"; }
  else if (id == 9 ){ return "MET"; }
  else if (id == 10){ return "ASN"; }
  else if (id == 11){ return "PRO"; }
  else if (id == 12){ return "GLN"; }
  else if (id == 13){ return "ARG"; }
  else if (id == 14){ return "SER"; }
  else if (id == 15){ return "THR"; }
  else if (id == 16){ return "VAL"; }
  else if (id == 17){ return "TRP"; }
  else if (id == 18){ return "TYR"; }
  else { return "NAN"; }
}

void
scan_Evalue( pose::Pose &pose, 
	     scoring::ScoreFunctionOP sfxn,
	     //scoring::ScoreFunctionOP sfxn_vdw,
	     utility::vector1< DataStruct > &data,
	     bool const rotscan
	     ){


  // Set index
  Size phibin = (Size)(nonperiodic(pose.phi( 3 ) + 5.0)/10.0) + 1;
  Size psibin = (Size)(nonperiodic(pose.psi( 3 ) + 5.0)/10.0) + 1;
  Size bbid = phibin + 36*(psibin-1);
  Size aaid = aa_to_id( pose.residue( 3 ).type().name() );

  // Set scoring mask
  DataStruct &datum = data[ (aaid-1)*36*36 + bbid ];
  datum.scored = true;
  
  if( rotscan ){

    for( Size j = 1; j <= datum.chis.size(); ++j ){
      pose::Pose pose_tmp( pose );

      for( Size k = 1; k <= 4; ++k )
	pose_tmp.set_chi( k, 3, datum.chis[j][k] );

      sfxn->score( pose_tmp );
      Real score = pose_tmp.energies().residue_total_energy( 3 );
      
      datum.scores[j].push_back( score );
      if( score > 0.0 ){
	datum.nclash[j] += 1;
      }
    }

  } else {
    // Find nearest rotamer index
    Size j = find_nearest_rot( pose.residue( 3 ).chi(), datum.chis );
    if( j == 0 ){
      //std::cout << ;
      return;
    }
    
    // score
    sfxn->score( pose );
    Real score = pose.energies().residue_total_energy( 3 );
    
    // add only if there is no clash
    if( score < 10.0 ){
      //datum.scores[j].push_back( std::min( 10.0, score ) );
      datum.scores[j].push_back( score );
    }
  }

}

utility::vector1< DataStruct >
initialize_data(){

  utility::vector1< DataStruct > data;
  Size k; 
  Real phi;
  Real psi;
  data.resize(18*36*36);

  pack::dunbrack::RotamerLibrary const & rotlib = * pack::dunbrack::RotamerLibrary::get_instance();

  chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  
  for( Size iaa = 1; iaa <= 18; ++ iaa ){
    // Setup rotamer lib
    chemical::ResidueType const rsdtype( rsd_set->name_map( id_to_aa(iaa) ) );

    pack::dunbrack::SingleResidueRotamerLibraryCAP 
      residue_rotamer_library( rotlib.get_rsd_library( rsdtype ) );
    pack::dunbrack::SingleResidueDunbrackLibraryCAP 
      dun_rotlib( dynamic_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const * >
		  ( residue_rotamer_library.get() ));

    // Iter over bins
    for( Size iphi = 1; iphi <= 36; ++ iphi ){
      for( Size ipsi = 1; ipsi <= 36; ++ ipsi ){
	k = (iaa-1)*36*36 + (ipsi-1)*36 + iphi;

	DataStruct &datum( data[k] );

	datum.aa = id_to_aa( iaa );
	datum.phi = periodic(iphi*10.0);
	datum.psi = periodic(ipsi*10.0);

	utility::vector1< pack::dunbrack::DunbrackRotamerSampleData > sample_data 
	  = dun_rotlib->get_all_rotamer_samples( periodic((iphi-1)*10.0), periodic((ipsi-1)*10.0) );
	
	datum.scores.resize( sample_data.size() );
	for ( Size j = 1; j <= sample_data.size(); ++j ){
	  datum.chis.push_back( sample_data[j].chi_mean() );
	  datum.prob.push_back( sample_data[j].probability() );
	  datum.nclash.push_back( 0 );
	}
	
      }
    }
  }
  return data;
}

}//myspace

int main( int argc, char * argv [] )
{
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace myspace;

  devel::init(argc, argv);

  NEW_OPT( fpd::verb, "verb", false );
  NEW_OPT( fpd::rotscan, "rotscan", false );

  bool verb( false );
  if( option[ fpd::verb] ){
    verb = option[ fpd::verb ]();
  }

  chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  // Score setup
  scoring::ScoreFunctionOP score_in 
    = core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );
  scoring::ScoreFunctionOP score_vdw;
    //    = core::scoring::ScoreFunctionFactory::create_score_function( option[score::we );

  utility::vector1< DataStruct > data = initialize_data();

  // silent
  io::silent::SilentFileData sfd;
  sfd.read_file( *(option[ in::file::silent ]().begin()) );
  
  pose::Pose pose;
  Size i( 0 );
  for( io::silent::SilentFileData::iterator iter = sfd.begin(); iter != sfd.end(); ++iter ){
    i++;
    if( i%1000 == 0 ){
      std::cerr << "Scanning " << i << std::endl;
    }
    iter->fill_pose( pose, *rsd_set );
    scan_Evalue( pose, score_in, data, option[ fpd::rotscan ]() );
  }

  report_data( data, verb );

  return 0;
}
