#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/pose/PDBInfo.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <iostream>
#include <fstream>

OPT_1GRP_KEY(String, fpd, countfile)
OPT_1GRP_KEY(String, fpd, prefix)

namespace myspace{

using namespace core;

struct DataStruct
{
  core::Real phi;
  core::Real psi;
  std::string aa;
  utility::vector1< pack::dunbrack::Real4 > chis;
  utility::vector1< Real > prob;
  //utility::vector1< bool > contact_exist;
  bool contact_exist;
  utility::vector1< utility::vector1< io::silent::SilentStructOP > > silents;
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

bool
is_contacting( conformation::Residue const &rsd1,
	       conformation::Residue const &rsd2 )
{
  Real D2CUT = 9.0;
  
  for( Size iatm = 1; iatm <= rsd1.natoms(); ++iatm ){
    if( (rsd1.atom_is_backbone( iatm ) || rsd1.atom_name( iatm ).compare(" CB ") == 0 ||
	 rsd1.atom_name( iatm ).compare(" HB ") == 0 ||
	 rsd1.atom_name( iatm ).compare("1HB ") == 0 || rsd1.atom_name( iatm ).compare("2HB ") == 0) ) continue;

    for( Size jatm = 1; jatm <= rsd2.natoms(); ++jatm ){
      if( !(rsd2.atom_is_backbone( jatm ) || rsd2.atom_name( jatm ).compare(" CB ") == 0) ) continue;
      
      Vector dxyz( rsd1.xyz(iatm) - rsd2.xyz(jatm) );
      Real d2 = dxyz.dot_product(dxyz);
      if( d2 < D2CUT ){
	return true;
      }
    }
  }

  return false;
}

Size 
find_nearest_rot( utility::vector1< Real > chi,  
		  utility::vector1< pack::dunbrack::Real4 > rotchis )
{
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
    if( dev < mindev ) { minid = i; }
  }

  return minid;
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
  else { return "NULL"; }
}

void
get_data_index( pose::Pose const & pose,
		Size const ires,
		utility::vector1< DataStruct > const &data,
		Size &data_index,
		Size &rot_index ){

  rot_index = 0;
  data_index = 0;
  conformation::Residue const &rsd( pose.residue( ires ) );

  Size aaid = aa_to_id( rsd.name() );
  if( aaid == 0){
    std::cerr << "Passing " << rsd.name() << std::endl;
    return;
  }

  // Set index
  Size phibin = (Size)(nonperiodic(pose.phi( ires ) + 5.0)/10.0) + 1;
  Size psibin = (Size)(nonperiodic(pose.psi( ires ) + 5.0)/10.0) + 1;
  Size bbid = phibin + 36*(psibin-1);
  
  // Return
  data_index = (aaid-1)*36*36 + bbid;
  rot_index = find_nearest_rot( pose.residue( ires ).chi(), data[data_index].chis );
  return;
}
 
void
scan_contact( pose::Pose const &pose, 
	     Size const &winsize,
	     utility::vector1< DataStruct > &data
	     ){

  Size data_index;
  Size rot_index;

  for( Size ires = 1+winsize; ires <= pose.total_residue()-winsize; ++ires ){
    std::string resname( pose.residue( ires ).name() );
    if( resname.compare("GLY") == 0 || resname.compare("ALA") == 0) continue;

    get_data_index( pose, ires, data, data_index, rot_index );
    if( data_index == 0 ) continue;

    DataStruct &datum = data[ data_index ];

    if( datum.contact_exist ) continue;

    for( Size k = 0; k <= winsize*2; ++k ){
      if( k == winsize ) continue;
      Size kres = ires-winsize+k;
      
      if( is_contacting( pose.residue( ires ), pose.residue( kres ) )){
	datum.contact_exist = true;
	std::cout << "Contact Found at " << ires << " " << data_index << " " << pose.residue( ires ).name() << 
	  " " << datum.phi << " " << datum.psi << " " << rot_index << std::endl;
	break;
      }
    }
  }
}

void
collect_silent( pose::Pose &pose, 
		Size const &winsize,
		utility::vector1< DataStruct > &data
		){

  Size data_index;
  Size rot_index;
  Size ncollect( 0 );

  core::pose::Pose pose5;

  // Is it better to use ALA? or native sequence?
  std::string sequence( 5, 'A' );
  core::pose::make_pose_from_sequence( pose5, sequence,
	*( chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" ))  );

  for( Size ires = 1+winsize; ires <= pose.total_residue()-winsize; ++ires ){
    std::string resname( pose.residue( ires ).name() );
    if( resname.compare("GLY") == 0 || resname.compare("ALA") == 0) continue;

    get_data_index( pose, ires, data, data_index, rot_index );

    DataStruct &datum = data[ data_index ];

    //std::cout << "data?" << data_index << " " << datum.contact_exist << " " << false << std::endl;

    if( !datum.contact_exist ) continue;

    ncollect ++;
    // mutate pose5
    pose5.replace_residue( 3, pose.residue( ires ), true );
    for( Size pos = 1; pos <= 5; pos++ ) {
      pose5.set_phi( pos, pose.phi( ires-3+pos ) );
      pose5.set_psi( pos, pose.psi( ires-3+pos ) );
      pose5.set_omega( pos, pose.omega( ires-3+pos ) );
    }
    
    std::stringstream tag;
    tag << pose.pdb_info()->name() << "." << ires;
    
    io::silent::SilentStructOP ss = 
      io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");

    ss->fill_struct( pose5 );
    ss->set_decoy_tag( tag.str() );
    datum.silents[rot_index].push_back( ss );

  }
  std::cerr << "Collected "<< ncollect << " structures from " << pose.pdb_info()->name() << std::endl;
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

    for( Size iphi = 1; iphi <= 36; ++ iphi ){
      for( Size ipsi = 1; ipsi <= 36; ++ ipsi ){
	k = (iaa-1)*36*36 + (ipsi-1)*36 + iphi;

	DataStruct &datum( data[k] );

	datum.aa = id_to_aa( iaa );
	datum.phi = periodic(iphi*10.0);
	datum.psi = periodic(ipsi*10.0);
	datum.contact_exist = false;

	utility::vector1< pack::dunbrack::DunbrackRotamerSampleData > sample_data 
	  = dun_rotlib->get_all_rotamer_samples( periodic((iphi-1)*10.0), periodic((ipsi-1)*10.0) );

	datum.silents.resize( sample_data.size() );
	for ( Size j = 1; j <= sample_data.size(); ++j ){
	  datum.chis.push_back( sample_data[j].chi_mean() );
	  datum.prob.push_back( sample_data[j].probability() );
	}
	
      }
    }
  }

  return data;
}

void
write_silent( utility::vector1< DataStruct > &data,
	      std::string const prefix ){
  
  Size MIN_DATA_SIZE( 0 );

  for( Size data_index = 1; data_index <= data.size(); ++data_index ){
    DataStruct &datum = data[ data_index ];
    if( !datum.contact_exist ) continue;

    for( Size rotid = 1; rotid <= datum.silents.size(); ++ rotid ){
      //if ( datum.silents[rotid].size() < MIN_DATA_SIZE ) continue;

      std::stringstream ssname;
      ssname << datum.aa << "/" << prefix << "." << datum.aa << "." << (Size)(nonperiodic(datum.phi)) << 
	"." << (Size)(nonperiodic(datum.psi)) << "." << rotid << ".out";

      io::silent::SilentFileData sfd;
      for( Size i_ss = 1; i_ss <= datum.silents[rotid].size(); ++i_ss ){
	io::silent::SilentStructOP ss( datum.silents[rotid][i_ss] );
	sfd.write_silent_struct( *ss, ssname.str() );
      }
    }
  }
}

// This is just for reading already prepared db
void
read_contact( std::string const infile, 
	      utility::vector1< DataStruct > &data ){
  
  utility::io::izstream instream( infile );
  //instream.open( infile );

  std::string aaname;
  Size data_index, aaid, phibin, psibin, ndat;
  Size nread = 0;
  
  while (instream) {
    //instream >> aaname >> data_index >> phibin >> psibin >> ndat;
    instream >> data_index;
    aaid = aa_to_id( aaname );

    data[data_index].contact_exist = true;
    nread ++;

    // Set index
    //Size phibin = (Size)(nonperiodic(phi( ires ) + 5.0)/10.0) + 1;
    //Size psibin = (Size)(nonperiodic(psi( ires ) + 5.0)/10.0) + 1;
    //Size bbid = phibin + 36*(psibin-1);
    //data_index = (aaid-1)*36*36 + bbid;
  }

  std::cerr << "Read " << nread << " contacts from " << infile << std::endl;
}

}//myspace

int main( int argc, char * argv [] )
{
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace myspace;

  NEW_OPT( fpd::countfile, "countfile", "" );
  NEW_OPT( fpd::prefix, "prefix", "" );

  devel::init(argc, argv);

  Size window_size( 2 );

  chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

  utility::vector1< DataStruct > data = initialize_data();

  utility::vector1< std::string > pdbs;
  utility::io::izstream list( option[ in::file::l ](1).name().c_str() );
  while ( list ) {
    std::string pdbname;
    list >> pdbname;
    if( pdbname.length() > 4){
      pdbs.push_back( pdbname );
    }
  }

  if( option[fpd::countfile].user() ){
    // Read from datafile
    read_contact( option[fpd::countfile](), data );

  } else {
    //Scan contact
    for( Size i = 1; i <= pdbs.size(); ++i ){
      std::cerr << "Scanning " << i << " th pdb: " << pdbs[i] << std::endl; 
      pose::Pose pose;
      import_pose::pose_from_pdb( pose, *rsd_set, pdbs[i] );
      scan_contact( pose, window_size, data );
    }
  }

  for( Size i = 1; i <= pdbs.size(); ++i ){
    std::cerr << "collecting " << i << " th pdb: " << pdbs[i] << std::endl; 
    pose::Pose pose;
    import_pose::pose_from_pdb( pose, *rsd_set, pdbs[i] );
    collect_silent( pose, window_size, data );
  }

  std::string outprefix = option[ fpd::prefix ]();
  write_silent( data, outprefix );

  return 0;
}
