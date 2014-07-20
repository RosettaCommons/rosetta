//id
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Pose
#include <core/pose/PDB_Info.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Scoring
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <devel/init.hh>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

void
get_resmap( pose::Pose const &pose, 
	    pose::Pose const &native,
	    std::map< Size, Size > &resmap
	    //id::AtomID_Map< id::AtomID > atommap
	    )
{

  for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
    Size ii_pdb( pose.pdb_info()->number( ii ) );
    //Size i_ca1 = pose.residue(ii).atom_index(" CA ");
    //id::AtomID id1( i_ca1, ii );

    for ( Size jj = 1; jj <= native.total_residue(); ++jj ) {
      Size jj_pdb( native.pdb_info()->number( jj ) );
      if( ii_pdb == jj_pdb ){
	resmap[ii] = jj;
	//Size i_ca2 = native.residue(jj).atom_index(" CA ");
	//id::AtomID id2( i_ca2, jj );
	//atommap[ id1 ] =  id2;
	break;
      }
    }
  }
}

void
decompose( pose::Pose const &pose0, 
	   pose::Pose const &pose,
	   pose::Pose const &native,
	   std::string const &pdbname,
	   std::map< Size, Size > &resmap,
	   bool const verbose )
{
  
  Size nloop( 0 );
  Size ncore( 0 );
  Real core_GDT( 0.0 );
  Real loop_GDT( 0.0 );

  for( Size ires = 1; ires <= pose.total_residue(); ++ires ){

    if( resmap.count( ires ) == 0 ) continue;

    Size jres = resmap[ires];
    Vector dcrd = pose.residue( ires ).xyz(" CA ") - native.residue( jres ).xyz(" CA ");
    Real dist = std::sqrt( dcrd.dot( dcrd ) );

    Real val( 0 );
    if( dist < 1 ){
      val = 4;
    } else if ( dist < 2 ){
      val = 3;
    } else if ( dist < 4 ){
      val = 2;
    } else if ( dist < 8 ){
      val = 1;
    }

    if( pose0.secstruct( ires ) == 'L' ){
      nloop++;
      loop_GDT += val;
    } else {
      ncore++;
      core_GDT += val;
    }

    if( verbose )
      printf("%3d %3d %1a %8.3f %1d\n",
	     ires, jres, pose.secstruct( ires ), dist, val );
  }

  Real tot_GDT = loop_GDT + core_GDT;
  Real nres = (Real)(nloop+ncore);
  tot_GDT *= 25.0/nres;

  core_GDT *= 25.0/nres;
  loop_GDT *= 25.0/nres;

  printf("%30s %s %6.2f %6.2f %6.2f %2d %2d\n",
	 pdbname.c_str(),
	 "Decomposed TotGDT/coreGDT/loopGDT/ncore/nloop: ",
	 tot_GDT, core_GDT, loop_GDT, 
	 ncore, nloop,
	 ncore*100.0/nres, nloop*100.0/nres
	 );
}

int main( int argc, char *argv [] ){

  devel::init(argc, argv);

  core::pose::Pose pose0, pose_work, native;
  core::chemical::ResidueTypeSetCAP rsd_set
    = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  
  Size npdb = option[ in::file::s ]().size();
  std::string const pose_init = option[ in::file::s ](1);
  
  core::import_pose::pose_from_pdb( native, *rsd_set, option[ in::file::native ]() ); 
  core::import_pose::pose_from_pdb( pose0, *rsd_set, pose_init ); 

  std::map< Size, Size > resmap;
  get_resmap( pose0, native, resmap );

  // Run dssp first
  core::scoring::dssp::Dssp dssp( pose0 );
  dssp.insert_ss_into_pose( pose0 );
  
  //decompose
  for( Size i = 1; i <= npdb; ++i ){
    std::string const pdbname = option[ in::file::s ](i);
    core::import_pose::pose_from_pdb( pose_work, *rsd_set, pdbname ); 
    decompose( pose0, pose_work, native, pdbname, resmap, false );
  }

  return 0;
}

