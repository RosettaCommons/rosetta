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
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

//#include <iomanip>
#include <iostream>
#include <fstream>

OPT_1GRP_KEY(Boolean, fpd, eval_full)
OPT_1GRP_KEY(Boolean, fpd, dump_pdb)
OPT_1GRP_KEY(String,  fpd, aa_dump)

using namespace core;
using namespace core::scoring;

// Predeclaration
namespace myspace{

struct DatumStruct
{
  utility::vector1< Real > Emap;
  Real Esum;
  std::string tag;
};

typedef std::map< std::string, utility::vector1< DatumStruct > > DatastoreType;

class DataClass
{

public:
  DataClass( bool const by_aa, bool const score_full )
  {
    std::map< std::string, utility::vector1< Real > > Evector;
    compname_.resize( 0 );
    compname_.push_back( fa_atr );
    compname_.push_back( fa_rep );
    compname_.push_back( facts_elec );
    compname_.push_back( facts_solv );
    compname_.push_back( facts_sasa );
    compname_.push_back( fa_dun );
    compname_.push_back( rama );
    compname_.push_back( p_aa_pp );

    Ncomp_ = compname_.size();
    if( by_aa ) by_aa_ = true;
    if( score_full ) score_full_ = true;
  }

  ~DataClass(){}
  
  void scan_Evalue( pose::Pose &pose,
		    scoring::ScoreFunctionOP sfxn,
		    Size const resno,
		    std::string const tag,
		    std::string const mode );

  void scan_Evalue_and_write( pose::Pose &pose,
			      scoring::ScoreFunctionOP sfxn,
			      Size const resno,
			      std::string const trg,
			      std::string const aa1,
			      std::string const aa2 );

  void report_by_aa( std::string const mode );

  void open_trgfile( std::string const filename ){ 

    trgfile.open( filename.c_str() ); 
    trgfile << " Trg res aa1 aa2";
    trgfile << std::setw(10) << "Esum";
    for( Size icomp = 1; icomp <= Ncomp_; ++icomp )
      trgfile << " " << std::setw(10) << compname_[icomp];
    trgfile << std::endl;
  }

  void set_dumppdb( bool dump_pdb, std::string aa_dump );
  
  void close_trgfile(){ trgfile.close(); }
  bool by_aa(){ return by_aa_; }

private:
  inline 
  utility::vector1< Real >
  Emap_to_Evector( scoring::EnergyMap const emap ){
    utility::vector1< Real > Evector( Ncomp_ );
    for( Size i = 1; i <= Ncomp_; ++i )
      Evector[i] = emap[ compname_[i] ];

    return Evector;
  }

private:
  Size Ncomp_;
  bool by_aa_;
  bool score_full_;
  bool dump_pdb_;
  std::string aa2_dump_;

  utility::vector1< core::scoring::ScoreType > compname_;
  std::ofstream trgfile;

  DatastoreType data_self_, data_far_;
};

  //}
  //namespace myspace{

void
DataClass::scan_Evalue( pose::Pose &pose, 
			scoring::ScoreFunctionOP sfxn,
			Size const resno,
			std::string const tag,
			std::string const mode
			){

  DatumStruct datum;
  scoring::EnergyMap emap;

  // score
  sfxn->score( pose );

  if( score_full_ ){
    emap = pose.energies().total_energies();
    datum.Esum = pose.energies().total_energy();
  } else {
    emap  = pose.energies().residue_total_energies( resno );
    datum.Esum = pose.energies().residue_total_energy( resno );
  }

  // weighted Emap
  for( Size i_comp = 1; i_comp <= Ncomp_; ++i_comp ){
    scoring::ScoreType scorename( compname_[i_comp] );
    emap[ scorename ] *= sfxn->get_weight( scorename );
  }

  datum.Emap = Emap_to_Evector( emap );
  datum.tag = tag;

  std::string const aatype( pose.residue( resno ).name3() );

  if( mode == "self" ){
    data_self_[aatype].push_back( datum );
  } else if( mode == "far" ){
    data_far_[aatype].push_back( datum );
  }
  
}

void
DataClass::scan_Evalue_and_write( pose::Pose &pose, 
				  scoring::ScoreFunctionOP sfxn,
				  Size const resno,
				  std::string const trg,
				  std::string const aa1,
				  std::string const aa2
				  ){

  // score
  sfxn->score( pose );

  DatumStruct datum;
  scoring::EnergyMap emap;

  if( score_full_ ){
    emap = pose.energies().total_energies();
    datum.Esum = pose.energies().total_energy();
  } else {
    emap  = pose.energies().residue_total_energies( resno );
    datum.Esum = pose.energies().residue_total_energy( resno );
  }

  // weighted Emap
  for( Size i_comp = 1; i_comp <= Ncomp_; ++i_comp ){
    scoring::ScoreType scorename( compname_[i_comp] );
    emap[ scorename ] *= sfxn->get_weight( scorename );
  }
  datum.Emap = Emap_to_Evector( emap );
  
  trgfile << trg << " " << std::setw(3) << resno << " " << aa1 << " " << aa2;
  trgfile << " " << std::setw(10) << datum.Esum;
  for( Size icomp = 1; icomp <= Ncomp_; ++icomp )
    trgfile << " " << std::setw(10) << datum.Emap[icomp];
  trgfile << std::endl;


  if( report_pdb_failed_ && aa2.compare( aa2_dump_ ) ){
    pose.dump_pdb( tag );
  }
}

void 
DataClass::set_dumppdb( bool dump_pdb, std::string aa_dump ){
  dump_pdb_ = dump_pdb;
  aa2_dump_ = aa_dump;
}

void 
DataClass::report_by_aa( std::string const mode ){
  
  if( !by_aa_ ) return;
  //std::map< std::string, utility::vector1< utility::vector1< Real > > >::iterator it;
  
  DatastoreType data;
  if( mode.compare( "self" ) == 0 ){
    data = data_self_;
  } else if( mode.compare( "far" ) == 0 ){
    data = data_far_;
  }

  DatastoreType::iterator it;
  for( it = data.begin(); it != data.end(); ++it ){
    std::string aa( it->first );
    utility::vector1< DatumStruct > const &aa_data = it->second;

    std::string filename = aa+"."+mode+".score";
    std::cout << "filename? " << filename << std::endl;
    std::ofstream aafile( filename.c_str() );
    
    //std::cout << "Size for " << aa << " " << aa_emap.size() << std::endl;

    // Header
    aafile << std::setw(10) << "Esum";
    for( Size icomp = 1; icomp <= Ncomp_; ++icomp )
      aafile << " " << std::setw(10) << compname_[icomp];
    aafile << " Tag" << std::endl;
    
    // Body
    for( Size i = 1; i <= aa_data.size(); ++i ){
      aafile << std::setw(10) << aa_data[i].Esum;
      for( Size icomp = 1; icomp <= Ncomp_; ++icomp )
	aafile << " " << std::setw(10) << aa_data[i].Emap[icomp];
      
      aafile << aa_data[i].tag << std::endl;
    }

    aafile.close();
  }

}

}//myspace

int main( int argc, char * argv [] )
{
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace myspace;

  devel::init(argc, argv);

  //NEW_OPT( fpd::rotscan, "rotscan", false );
  NEW_OPT( fpd::eval_full, "eval_full", true );
  NEW_OPT( fpd::aa_dump,   "aa_dump",   "" );
  NEW_OPT( fpd::dump_pdb,  "eval_full", true );

  chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  // Score setup
  scoring::ScoreFunctionOP score_in 
    = core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

  //DataClass data( false, false );
  bool by_aa( false );
  bool eval_full( option[ fpd::eval_full ]() );
  std::string aa_dump( option[ fpd::aa_dump ]() );
  bool dump_pdb( option[ fpd::dump_pdb ]() );

  // Full structure evaluation, dump pdb
  DataClass data( by_aa, eval_full ); 
  data.set_dumppdb( dump_pdb, aa_dump );

  Size const nsilents( option[ in::file::silent ]().size() );

  // Scan score
  pose::Pose pose;
  for( Size i_file = 1; i_file <= nsilents; ++i_file ){
    io::silent::SilentFileData sfd;
    std::string filename( option[ in::file::silent ](i_file) );

    sfd.read_file( filename );
    utility::vector1< std::string > tags( sfd.tags() );

    std::string const trg( tags[1], 0, 4 );
    data.open_trgfile( trg+".mut.sc" );

    std::cout << "Read " << tags.size() << " from " << filename << std::endl;

    Size i( 0 );
    for( io::silent::SilentFileData::iterator iter = sfd.begin(); iter != sfd.end(); ++iter ){

      // Get resno/mutation state of interest
      i++;
      std::string resinfo( utility::string_split( tags[i], '.' )[2], 3 );

      std::string resnostr( resinfo, 3, resinfo.size()-3 );

      Size const resno( atoi(resnostr.c_str()) );
      std::string aa1( resinfo, 0, 3 );
      std::string aa2( resinfo, resinfo.size()-3, resinfo.size() );

      std::string mode;
      if( aa1.compare( aa2 ) == 0 ){
	mode = "self";
      } else {
	mode = "far";
      }

      // Run scoring
      iter->fill_pose( pose, *rsd_set );

      if( data.by_aa() ){
	data.scan_Evalue( pose, score_in, resno, tags[i], mode );
      } else {
	data.scan_Evalue_and_write( pose, score_in, resno, trg, aa1, aa2 );
	
      }

    }
    
    data.close_trgfile();
  }

  data.report_by_aa( "self" );
  data.report_by_aa( "far" );

  return 0;
}
