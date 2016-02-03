#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <iostream>
#include <fstream>

OPT_1GRP_KEY(String, fpd, input)

namespace myspace{

using namespace core;

struct InputData
{
  std::string pdbname;
  Size resno;
  
  std::string aa;
  utility::vector1< Real > phi;
  utility::vector1< Real > psi;
  utility::vector1< Real > omg;
  Size fraglen;
  Size index;
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

utility::vector1< InputData > 
read_input( std::string const filename ){

  utility::io::izstream instream;
  instream.open( filename );
  
  std::string fileline;
  std::string line;
  utility::vector1< InputData > inputs;

  while( instream ) {
    InputData input;

    getline( instream, line );
    std::istringstream linestream( line );
    if( line.size() < 2 ) break;
    if( line[0] == '#' ) continue;

    linestream >> input.aa >> input.fraglen;

    input.phi.resize( input.fraglen );
    input.psi.resize( input.fraglen );
    input.omg.resize( input.fraglen, 180.0 );

    for( Size ipos = 1; ipos <= input.fraglen; ++ipos )
      linestream >> input.phi[ipos] >> input.psi[ipos];

    input.index = inputs.size() + 1;

    inputs.push_back( input );
  }

  return inputs;
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

pose::Pose
pose_from_scratch( InputData input ){
  pose::Pose pose;
  std::string sequence;
  Size poseoff = (input.fraglen-1)/2;

  for(Size i = 1; i <= poseoff; ++i) sequence += "A";
  sequence += input.aa;
  for(Size i = 1; i <= poseoff; ++i) sequence += "A";

  pose::make_pose_from_sequence( pose, sequence,
	*( chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" ))  );
  
  for( Size pos = 1; pos <= input.fraglen; pos++ ) {
    pose.set_phi( pos, input.phi[pos] );
    pose.set_psi( pos, input.psi[pos] );
    pose.set_omega( pos, input.omg[pos]  );
  }

  return pose;
}

void
scan_Evalue( pose::Pose &pose, 
	     scoring::ScoreFunctionOP sfxn,
	     chemical::ResidueTypeSetCAP rsd_set,
	     pack::dunbrack::RotamerLibrary const &rotlib,
	     InputData input
	     ){

  Size cenpos( input.fraglen/2 + 1 );
  scoring::EnergyMap const wts = sfxn->weights();
  scoring::ScoreTypes scores = sfxn->get_nonzero_weighted_scoretypes();

  std::cout << pose.residue(cenpos).name();
  printf(" %6.1f %6.1f\n",
	 pose.phi(cenpos), pose.psi(cenpos));

  printf(" j n %6s %6s %6s %6s %8s", "chi1", "chi2", "chi3", "chi4", "score");
  for( Size i = 1; i <= scores.size(); ++i ){
    if( wts[scores[i]] > 0.001 ) std::cout << " " << std::setw(10) << scoring::name_from_score_type( scores[i] );
  }
  printf("\n");

  // Rotlib setup
  chemical::ResidueType const rsdtype( rsd_set->name_map( pose.residue(cenpos).name() ) );

  pack::rotamers::SingleResidueRotamerLibraryCAP 
    residue_rotamer_library( rotlib.get_rsd_library( rsdtype ) );
  pack::dunbrack::SingleResidueDunbrackLibraryCAP 
    dun_rotlib( dynamic_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const * >
		( residue_rotamer_library.get() ));

  utility::vector1< pack::dunbrack::DunbrackRotamerSampleData > sample_data 
    = dun_rotlib->get_all_rotamer_samples( pose.phi(cenpos), pose.psi(cenpos) );

  //pose.dump_pdb("test.pdb");

  // Iter over rotamer
  for( Size j = 1; j <= sample_data.size(); ++j ){
    pose::Pose pose_tmp( pose );
    
    pack::dunbrack::Real4 chis = sample_data[j].chi_mean();
    for( Size k = 1; k <= pose_tmp.residue( cenpos ).chi().size(); ++k )
      pose_tmp.set_chi( k, cenpos, chis[k] );
    
    sfxn->score( pose_tmp );

    Real const score = pose_tmp.energies().residue_total_energy( cenpos );
    scoring::EnergyMap const emap = pose_tmp.energies().residue_total_energies( cenpos );
    
    printf("%2d %2d %6.1f %6.1f %6.1f %6.1f %8.3f", 
	   j, sample_data.size(),
	   chis[1], chis[2], chis[3], chis[4], score );

    for( Size i = 1; i <= scores.size(); ++i ){
      //if( scores[i] == scoring::rama || scores[i] == scoring::p_aa_pp_offset || scores[i] == scoring::ref ||
      //	  scores[i] == scoring::cart_bonded || scores[i] == dslf_fa13 || scores[i] == scoring::omega ||
      //	  scores[i] == scoring::fa_dun_dev ) continue;
    if( wts[scores[i]] > 0.001 )
      printf(" %10.3f",wts[scores[i]]*emap[scores[i]] );
    }
    printf("\n");

    std::stringstream pdbname("");
    pdbname << "input" << input.index << "." << j << ".pdb";
    pose_tmp.dump_pdb( pdbname.str() );
  }

}
} // end myspace

int main( int argc, char * argv [] )
{
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace myspace;

  NEW_OPT( fpd::input, "input", "" );

  devel::init(argc, argv);

  chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  pack::dunbrack::RotamerLibrary const &rotlib = * pack::dunbrack::RotamerLibrary::get_instance();

  // Score setup
  scoring::ScoreFunctionOP score_in 
    = core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

  utility::vector1<InputData> inputs = read_input( option[ fpd::input ]() );
  
  for( Size i = 1; i <= inputs.size(); ++i ){
    pose::Pose pose;
    //import_pose::pose_from_file( pose, *rsd_set, inputs[i].pdbname , core::import_pose::PDB_file); 
    std::cout << "Reading input " << i << std::endl;
    pose = pose_from_scratch( inputs[i] );
    scan_Evalue( pose, score_in, rsd_set, rotlib, inputs[i] );
  }

  return 0;
}
