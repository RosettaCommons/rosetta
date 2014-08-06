#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/pose/PDBInfo.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>

#include <utility/io/izstream.hh>
#include <iostream>
#include <fstream>

int main( int argc, char * argv [] )
{
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  devel::init(argc, argv);

  chemical::ResidueTypeSetCAP rsd_set_fa
    = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  chemical::ResidueTypeSetCAP rsd_set_cen
    = chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

  // Score setup
  scoring::ScoreFunctionOP score_in 
    = core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

  // silent
  io::silent::SilentFileData sfd;
  sfd.read_file( *(option[ in::file::silent ]().begin()) );

  protocols::moves::MoverOP tofa 
    = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );
  
  pose::Pose pose;
  Size i( 0 );

  std::string outprefix = option[ out::prefix ]();

  for( io::silent::SilentFileData::iterator iter = sfd.begin(); iter != sfd.end(); ++iter ){
    iter->fill_pose( pose, *rsd_set_cen );
    tofa->apply( pose );

    score_in->score( pose );

    io::silent::SilentStructOP ss = 
      io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");

    ss->fill_struct( pose );
    ss->set_decoy_tag( pose.pdb_info()->name() );
    ss->energies_from_pose( pose );
    sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
  }
  return 0;
}
