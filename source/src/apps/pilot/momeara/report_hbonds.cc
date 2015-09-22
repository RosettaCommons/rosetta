// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/momeara/report_hbonds.cc
///
/// @brief  report all data about hydrogen bonds
/// @author Matthew O'Meara


#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.fwd.hh>


#include <boost/algorithm/string.hpp>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>


#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/types.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/Tracer.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/ShakeStructureMover.hh>
#include <protocols/relax_protocols.hh>
// utility headers
#include <utility/io/izstream.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

// c++ headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

// namespaces
using namespace std;
using namespace boost;
using namespace core;
  using namespace conformation;
  using namespace kinematics;
  using namespace chemical;
  using namespace pose;
  using namespace protocols::jd2;
  using namespace scoring;
    using namespace hbonds;


using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;
using utility::vector1;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.momeara.HBondReporter" );


OPT_1GRP_KEY( String, HBondReporter, output )
OPT_1GRP_KEY( Boolean, HBondReporter, relax )

OPT_1GRP_KEY( String, HBondReporter, relevant_chains )
//OPT_1GRP_KEY( String, HBondReporter, job_data_fname )


class HBondReporter : public protocols::moves::Mover {

public:
  HBondReporter():
    allowNonProtein_( true ),
    scfxn( get_score_function() ),
    relevant_chains_( "*" ),
    hb_database_( HBondDatabase::get_database() ){
}

  virtual ~HBondReporter(){};


  void load_job_data( Pose & pose ){

    vector< string > lines, tokens;
    stringstream chains;
    if (!pose.pdb_info()){
      // this is true if, say, the structure came from a silent file
      return;
    }

    if ( basic::options::option[ basic::options::OptionKeys::HBondReporter::relevant_chains ].active() ){

      utility::io::izstream relevant_chains_file( basic::options::option[ basic::options::OptionKeys::HBondReporter::relevant_chains ]() );


      if ( !relevant_chains_file ){
	TR.Error << " Cannot open job data file "<< basic::options::option[ basic::options::OptionKeys::HBondReporter::relevant_chains ]() << std::endl;
	return;
      }
      string line;
      FileName pose_filename  = pose.pdb_info()->name();
      pose_filename.to_local_name();

      string chains;
      while( getline( relevant_chains_file, line ) ){
	//TR << "For line '"<< line <<"'" << std::endl;
	split(tokens, line, is_any_of("\t") );
	//TR << "\t tokens[0]="<< tokens[0] << " tokens[1]=" <<tokens[1] <<std::endl;
	//TR << "\t local pdb_filename " << pose_filename << std::endl;
	if( tokens[0]  == pose_filename ){
	  relevant_chains_ = tokens[1];
	  TR << "Restricting pdb '"<< pose.pdb_info()->name() << "' to chains " << relevant_chains_ << std::endl;

	  //	  frag3_file = tokens[2];

	  return;
	}
      }
    }

    TR << "Job data for PDB '"<< pose.pdb_info()->name() << "' not found!" << std::endl;
  }

  string pose_name(Pose& pose){
    //silent files and pdbs set the name of the pose differently
    string name = "No_Name_Found";
    if (pose.pdb_info()){
      name = pose.pdb_info()->name();
    } else if ( pose.data().has( datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
      name = static_cast< basic::datacache::CacheableString const & >
	( pose.data().get( datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
    } else {
      name = JobDistributor::get_instance()->current_job()->input_tag();
    }
    return name;
  }


  void get_aa_counts( Pose & pose ){
    vector1<Size> aa_counts(num_canonical_aas, 0);
    for( Size i=1; i <= pose.total_residue(); ++i ){
      if( pose.residue( i ).is_protein() ) {
	aa_counts[ pose.residue( i ).aa() ] += 1;
      }
    }

    TR << "Amino acid counts:";
    for( vector1<Size>::iterator i = aa_counts.begin(), i_end = aa_counts.end(); i != i_end; ++i ){
      TR << " " << *i;
    }
    TR << std::endl;
  }


  Real bfactor( Pose & pose, Size resNum, Size atm, bool is_backbone) {
    if (!pose.pdb_info()){
      return 9999;
    }

    if ( is_backbone )
      return pose.pdb_info()->temperature( resNum, atm);
    else {
      Real max_temp = -1000;
      for ( Size i=1; i <= pose.residue( resNum ).natoms(); i++)
	max_temp = max( max_temp, pose.pdb_info()->temperature( resNum, i ) );
      return max_temp;
    }
  }


  void decoyify( Pose& pose ) {

    MoveMapOP mm( new MoveMap );
    mm->clear();


    if( relevant_chains_ == "*" || relevant_chains_.empty() ){
      TR << "Information about relevant chains not found for " << pose_name(pose) << "!" << std::endl;
      mm->set_bb( true );
      mm->set_chi( true );
    } else if ( !relevant_chains_.empty() ) {
      for( Size j=1; j<=pose.total_residue(); ++j ){
	if( relevant_chains_.find( pose.pdb_info()->chain( j ), 0 ) != std::string::npos ){
	  mm->set_chi( j, true );
	  mm->set_bb( j, true );
	}
      }
    }

    protocols::relax::FastRelax relax(scfxn);
    //protocols::relax::ClassicRelax relax( scfxn, mm );

    // apply relax mover to pose
    TR << "Initial score: " << (*scfxn)(pose)<< std::endl;
    TR << "Relaxing...";
    relax.apply(pose);
    TR << "Done!!!" << std::endl;
    TR << "Decoy score: " << (*scfxn)(pose) << std::endl;
  }

virtual
void
apply( Pose& pose ){


    // proof of concept

//  TR << "Before updating the options system" << std::endl;
//    TR << basic::options::option[ basic::options::OptionKeys::in::file::frag3 ]() << std::endl;
//    TR << basic::options::option[ basic::options::OptionKeys::in::file::frag9 ]() << std::endl;
//    TR << basic::options::option[ basic::options::OptionKeys::in::file::native ]() << std::endl;
//
//    int argc = 3;
//    char * argv [] = {"refresh_options", "-frag3", "new_frag3" };
//    basic::options::option.load( argc, argv, false );
//    TR << "After updating the options system" << std::endl;
//    TR << basic::options::option[ basic::options::OptionKeys::in::file::frag3 ]() << std::endl;
//    TR << basic::options::option[ basic::options::OptionKeys::in::file::frag9 ]() << std::endl;
//    TR << basic::options::option[ basic::options::OptionKeys::in::file::native ]()<< std::endl;


  clock_t start_time = clock();

  TR<< "getting hbonds for " << pose_name(pose) << std::endl;

  load_job_data( pose );

  if ( basic::options::option[ basic::options::OptionKeys::HBondReporter::relax ].value() ) {
    decoyify( pose );
  }

  ofstream fout;
  if ( basic::options::option[ basic::options::OptionKeys::HBondReporter::output ].user() ) {
    fout.open( basic::options::option[ basic::options::OptionKeys::HBondReporter::output ]().c_str(), std::ios::app );
    TR << "Writing to file " << basic::options::option[ basic::options::OptionKeys::HBondReporter::output ]() << std::endl;
  } else {
    FileName pose_filename  = pose_name(pose);
    pose_filename.to_local_name();
    stringstream output_filename;
    output_filename << pose_filename << "_hbonds.txt";
    fout.open( output_filename.str().c_str(), std::ios::out );
    TR << "Writing to file " << output_filename.str() << std::endl;
  }


  HBondSet set1;
  (*scfxn)(pose);
  pose.update_residue_neighbors();
  set1.setup_for_residue_pair_energies( pose, false, false );

  TR << "Number of hydrogen bonds found: " << set1.nhbonds() << std::endl;

  // compute sasa pack
  Real const probe_radius(1.4);
  id::AtomID_Map< core::Real > atom_sasa;
  utility::vector1< core::Real > residue_sasa;
  scoring::calc_per_atom_sasa( pose, atom_sasa, residue_sasa, probe_radius);


  // compute dssp
  core::scoring::dssp::Dssp dssp( pose );
  dssp.insert_ss_into_pose( pose );

  core::Real hbond_energies = 0; // to get average energy for structure;
  core::Size nhbonds_used = 0;
  for (Size i = 1; i<= set1.nhbonds(); i++) {
    HBond bond = set1.hbond( i );

    //need to access donor and acc Residues as well as ints
    Size hatm = bond.don_hatm();
    HBEvalType type = bond.eval_type();
    Size accResNum = bond.acc_res();
    Size donResNum = bond.don_res();
    //get acc and donor residues from sequence numbers
    Residue accRes = pose.residue( accResNum );
    Residue donRes = pose.residue( donResNum );
    Size const datm( donRes.atom_base( hatm ) );
    Vector const & hatm_xyz( donRes.atom( hatm ).xyz() );
    //donRes.atom() is of type Atom
    Vector const & datm_xyz( donRes.atom( datm ).xyz() );
    Size aatm = bond.acc_atm();
    Size const base( accRes.atom_base( aatm ) );
    Size const base_of_base( accRes.atom_base( base ) );

    Size const base2( accRes.abase2( aatm ) );
    Size const dbase( donRes.atom_base( datm ) );
    Vector const & aatm_xyz = accRes.atom( aatm ).xyz();
    Vector const & base_xyz = accRes.atom( base ).xyz();
    Vector const & base_of_base_xyz = accRes.atom( base_of_base ).xyz();

    Vector const & base2_xyz = accRes.atom( base2 ).xyz();
    Vector const & dbase_xyz = donRes.atom( dbase ).xyz();

    Size donType = donRes.atom( datm ).type();
    Size accType = accRes.atom( aatm ).type();
    std::string donResName = donRes.name3();
    std::string accResName = accRes.name3();
    bool isProtein = true;
    if ( ! donRes.is_protein() || ! accRes.is_protein() ) {
      isProtein = false;
    }
    if ( !( (set1.allow_hbond(i) && isProtein) || (set1.allow_hbond(i) && allowNonProtein_) ) ) continue;


    const AtomType & donAtomType = donRes.atom_type( datm );
    const AtomType & accAtomType = accRes.atom_type( aatm );
    const std::string donElement = donAtomType.element();
    const std::string accElement = accAtomType.element();
    const std::string & donAtomName = donRes.atom_name( datm );
    const std::string & accAtomName = accRes.atom_name( aatm );
    const std::string & hAtomName = donRes.atom_name( hatm );
    int donElemNum, accElemNum;	//usual atomic numbers from periodic table

    char donCh, accCh;
    int donPdbResNum, accPdbResNum;
    char donPdbICode, accPdbICode;

    if (pose.pdb_info()){

      donCh = pose.pdb_info()->chain( donResNum );
      accCh = pose.pdb_info()->chain( accResNum );

      donPdbResNum = pose.pdb_info()->number( donResNum );
      accPdbResNum = pose.pdb_info()->number( accResNum );

      donPdbICode = pose.pdb_info()->icode( donResNum );
      accPdbICode = pose.pdb_info()->icode( accResNum );
    } else {
      donCh = 'A';
      accCh = 'A';

      donPdbResNum = donResNum;
      accPdbResNum = accResNum;

      donPdbICode = ' ';
      accPdbICode = ' ';
    }

    if( relevant_chains_ != "*" ){
      if( ( relevant_chains_.find( donCh, 0 ) == std::string::npos ) ||
	  (relevant_chains_.find( accCh, 0 ) == std::string::npos ) ) continue;
    }

    Size const donNbrs( set1.nbrs( donResNum ) );
    Size const accNbrs( set1.nbrs( accResNum ) );

    Real const donSasa( atom_sasa( donResNum, datm ) );
    Real const accSasa( atom_sasa( accResNum, aatm ) );

    Real donBfactor, accBfactor;
    if (!pose.pdb_info()){
      donBfactor = 9999;
      accBfactor = 9999;
    } else {
      donBfactor = pose.pdb_info()->temperature( donResNum, datm );
      accBfactor = pose.pdb_info()->temperature( accResNum, aatm );
    }

    char const & donSS = pose.secstruct( donResNum );
    char const & accSS = pose.secstruct( accResNum );


    HBond::Deriv deriv = bond.deriv();
    Real energy = bond.energy();
    hb_energy_deriv( database, type, datm_xyz, hatm_xyz, aatm_xyz, base_xyz, base2_xyz, energy, true, deriv );
    Real weight = bond.weight();


    hbond_energies +=  energy;


    if ( ! donElement.compare( "O" ) ) {
      donElemNum = 8;
    } else if ( ! donElement.compare("N") ) {
      donElemNum = 7;
    } else if ( ! donElement.compare("S") ) {
      donElemNum = 16;
    }	else {
      donElemNum = -1;
    }

    if ( ! accElement.compare( "O" ) ) {
      accElemNum = 8;
    } else if ( ! accElement.compare("N") ) {
      accElemNum = 7;
    }	else if ( ! donElement.compare("S") ) {
      accElemNum = 16;
    } else {
      accElemNum = -1;
    }

    /*
      std::string donor_back;
      if ( bond.don_hatm_is_protein_backbone() ){
      donor_back = "donBK";
      }	else {
      donor_back = "donSC";
      }
      std::string acc_back;
      if ( bond.acc_atm_is_protein_backbone() ){
      acc_back = "accBK";
      }	else {
      acc_back = "accSC";
      }
    */
    //The type of the hbond
    fout << type << "\t";

    // coordinates of atoms
    fout << datm_xyz.x()         << "\t" << datm_xyz.y()         << "\t" << datm_xyz.z()         << "\t";
    fout << hatm_xyz.x()         << "\t" << hatm_xyz.y()         << "\t" << hatm_xyz.z()         << "\t";
    fout << aatm_xyz.x()         << "\t" << aatm_xyz.y()         << "\t" << aatm_xyz.z()         << "\t";
    fout << base_xyz.x()         << "\t" << base_xyz.y()         << "\t" << base_xyz.z()         << "\t";
    fout << base_of_base_xyz.x() << "\t" << base_of_base_xyz.y() << "\t" << base_of_base_xyz.z() << "\t";
    fout << base2_xyz.x()        << "\t" << base2_xyz.y()        << "\t" << base2_xyz.z()        << "\t";
    fout << dbase_xyz.x()        << "\t" << dbase_xyz.y()        << "\t" << dbase_xyz.z()        << "\t";
    // chemical information about donor and acceptor
    fout << donCh << "\t" << donResName << "\t" << donResNum << "\t" << donType << "\t" << donAtomName << "\t" << donElemNum << "\t";
    fout << accCh << "\t" << accResName << "\t" << accResNum << "\t" << accType << "\t" << accAtomName << "\t" << accElemNum << "\t";


    // misc parameters
    fout << donNbrs << "\t" << donSasa << "\t" << donBfactor << "\t" << donSS << "\t";
    fout << accNbrs << "\t" << accSasa << "\t" << accBfactor << "\t" << accSS << "\t";

    // Energy and environmental weight
    fout << energy << "\t" << weight << "\t";

    fout << hAtomName << "\t";
    fout << donPdbResNum << "\t" << donPdbICode << "\t";
    fout << accPdbResNum << "\t" << accPdbICode << "\t";


    // coordinates of closest reduce placed hydrogen
    fout << "0\t0\t0\t";

    // First and second derivative
    fout << deriv.first[  0 ] << "\t" << deriv.first[  1 ] << "\t" << deriv.first[  2 ] << "\t";
    fout << deriv.second[ 0 ] << "\t" << deriv.second[ 1 ] << "\t" << deriv.second[ 2 ] << "\t";

    // PDB name
    fout << pose_name(pose);
    fout << "\n";

    nhbonds_used++;

  }	//close for
  fout.close();

  get_aa_counts( pose );


  TR << "Average hbond energy = " << hbond_energies << "/" << nhbonds_used <<"="<< hbond_energies / nhbonds_used << std::endl;
  clock_t stoptime = clock();
  TR << "Timing: " << ((double) stoptime - start_time )/CLOCKS_PER_SEC << "s" << std::endl;
  TR << "DONE getting hbonds for " << pose_name(pose) << std::endl;

}


private:
  bool allowNonProtein_;
  core::scoring::ScoreFunctionOP scfxn;
  std::string relevant_chains_;
  HBondDatabaseCOP hb_database_;
};

typedef utility::pointer::owning_ptr< HBondReporter > HBondReporterOP;

int main( int argc, char* argv[] )
{
	try{
	  NEW_OPT( HBondReporter::output, "Location of dumped hbond data", "hbonds.txt");
	  NEW_OPT( HBondReporter::relax, "Perform relaxation before dumping hbonds", false);
	  NEW_OPT( HBondReporter::relevant_chains, "relevant_chains", "*");
	//NEW_OPT( HBondReporter::job_data_fname, "pdb indexed table with job related information", "");

	  devel::init(argc, argv);
	  JobDistributor::get_instance()->go(new HBondReporter);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
