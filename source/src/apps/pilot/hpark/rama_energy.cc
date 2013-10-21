#include <devel/init.hh>

#include <basic/options/option.hh>
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
#include <sstream>

using namespace core;

pose::Pose
pose_from_scratch( chemical::AA aa ){
  pose::Pose pose;
  std::stringstream s( "" );
  std::string sequence;

  s << "A";
  s << chemical::oneletter_code_from_aa( aa );
  s << "A";
  s >> sequence;

  pose::make_pose_from_sequence( pose, sequence,
	*( chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" ))  );

  return pose;
}

void
scan_Evalue( chemical::AA aa,
	     scoring::ScoreFunctionOP sfxn,
	     chemical::ResidueTypeSetCAP rsd_set
	     ){

  pose::Pose pose( pose_from_scratch( aa ) );

  scoring::EnergyMap const wts = sfxn->weights();
  scoring::ScoreTypes scores = sfxn->get_nonzero_weighted_scoretypes();

  printf("AA  %6s %6s %8s %8s %8s\n", "phi", "psi", "Erama", "Epaapp", "Mindun" );

  pack::dunbrack::RotamerLibrary const &rotlib = pack::dunbrack::RotamerLibrary::get_instance();

  pack::dunbrack::SingleResidueRotamerLibraryCAP 
    residue_rotamer_library( rotlib.get_rsd_library( pose.residue( 2 ).type() ) );
  pack::dunbrack::SingleResidueDunbrackLibraryCAP 
    dun_rotlib( dynamic_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const * >
		( residue_rotamer_library.get() ));

  for( int ii = -18; ii != 18; ++ii ){
    for( int jj = -18; jj != 18; ++jj ){
      Real const phi( ii*10.0 );
      Real const psi( jj*10.0 );

      pose.set_phi( 2, phi );
      pose.set_psi( 2, psi );
      pose.set_omega( 2, 180.0  );

      if( pose.residue( 2 ).chi().size() > 0 ){
	utility::vector1< pack::dunbrack::DunbrackRotamerSampleData > sample_data 
	  = dun_rotlib->get_all_rotamer_samples( phi, psi );
	pack::dunbrack::Real4 chis = sample_data[1].chi_mean();
	for( Size k = 1; k <= pose.residue( 2 ).chi().size(); ++k )
	  pose.set_chi( k, 2, chis[k] );
      }

      sfxn->score( pose );

      Real const score = pose.energies().residue_total_energy( 2 );
      scoring::EnergyMap const emap = pose.energies().residue_total_energies( 2 );
    
      //if( emap[ scoring::rama ] > 3.0 ) continue;

      std::cout << std::setw(3) << aa;
      printf(" %6.1f %6.1f %8.3f %8.3f %8.3f\n", 
	     phi, psi, 
	     wts[ scoring::rama ]*emap[ scoring::rama ], 
	     wts[ scoring::p_aa_pp ]*emap[ scoring::p_aa_pp ],
	     wts[ scoring::fa_dun ]*emap[ scoring::fa_dun ]
	     );
    }
  }

}

int main( int argc, char * argv [] )
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  devel::init(argc, argv);

  chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

  // Score setup
  scoring::ScoreFunctionOP score_in 
    = core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

  std::list< chemical::AA >::const_iterator it;
  for( it = rsd_set->aas_defined_begin(); it != rsd_set->aas_defined_end(); ++it ){
    if( *it > chemical::num_canonical_aas ) break;

    scan_Evalue( *it, score_in, rsd_set );
  }

  return 0;
}
