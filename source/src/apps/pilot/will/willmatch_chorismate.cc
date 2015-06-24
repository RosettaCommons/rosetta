// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/willmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/will_util.ihh>


using core::Real;
using core::Size;
using core::pose::Pose;
using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using std::string;
using utility::vector1;
using numeric::min;
using core::import_pose::pose_from_pdb;
using basic::options::option;
using numeric::min;
using numeric::max;

typedef utility::vector1<core::Real> Reals;
typedef utility::vector1<core::Size> Sizes;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

static thread_local basic::Tracer TR( "willmatch_chorismate" );


void myoptH(Pose & pose, ScoreFunctionOP sf) {
  add_lower_terminus_type_to_pose_residue(pose,1);
  add_upper_terminus_type_to_pose_residue(pose,pose.n_residue());
  core::pack::optimizeH(pose,*sf);
  remove_lower_terminus_type_from_pose_residue(pose,1);
  remove_upper_terminus_type_from_pose_residue(pose,pose.n_residue());
}


void run() {
  using namespace basic::options::OptionKeys;
  using namespace core::id;

  core::chemical::ResidueTypeSetCOP cen_residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ) );
  core::chemical::ResidueTypeSetCOP  fa_residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
  //Real tmpdis = option[willmatch::max_dis_metal]();  // unused ~Labonte
  // Real const MXDSMTL = tmpdis*tmpdis;
  // Real const MXAGMTL = option[willmatch::max_ang_metal]();
  // Real const MATCH_OVERLAP_DOT = cos(numeric::conversions::radians(option[willmatch::match_overlap_ang]()));


  core::io::silent::SilentFileData sfd;

  vector1<string> infiles;
  if( option[in::file::l].user() ) {
    utility::io::izstream in(option[in::file::l]()[1]);
    string tmp;
    while(in >> tmp) infiles.push_back(tmp);
  } else if(option[in::file::s].user()) {
    infiles = option[in::file::s]();
  } else {
    utility_exit_with_message("no input!");
  }

  std::string startfile = "";
  Size startrsd = 0;
  {
    utility::io::izstream inz(option[out::file::o]()+"/willmatch_chorismate.progress");
    while( inz >> startfile >> startrsd ) ;
    inz.close();
    if(startfile != "") TR<<"continuing from "<<startfile<<" "<<startrsd<<std::endl;
  }
  utility::io::ozstream oprogress(option[out::file::o]()+"/willmatch_chorismate.progress");

  for(Size ifile = 1; ifile <= infiles.size(); ifile++) {
    string infile = infiles[ifile];
    if( startfile != "" && startfile != infile ) continue; // CHECKPOINT
    Pose in_fa;
    pose_from_pdb(in_fa, *fa_residue_set,infile);
    for(Size ir = 1; ir <= in_fa.n_residue(); ++ir) {
      if(in_fa.residue(ir).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(in_fa,ir);
      if(in_fa.residue(ir).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(in_fa,ir);
    }
    Pose native = in_fa;
    Size nres = in_fa.n_residue();
    core::chemical::ResidueType const & rtala( in_fa.residue(1).residue_type_set().name_map("ALA") );
    // core::chemical::ResidueType const & rtasp( in_fa.residue(1).residue_type_set().name_map("ASP") );
    // core::chemical::ResidueType const & rtglu( in_fa.residue(1).residue_type_set().name_map("GLU") );
    core::chemical::ResidueType const & rtarg( in_fa.residue(1).residue_type_set().name_map("ARG") );
    for(Size i = 1; i <= nres; ++i) core::pose::replace_pose_residue_copying_existing_coordinates(in_fa,i,rtala);
    Pose const fa_pose = in_fa;
    ImplicitFastClashCheck clashcheck(fa_pose,basic::options::option[basic::options::OptionKeys::willmatch::clash_dis]());

    ScoreFunctionOP sf = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
    sf->set_weight(core::scoring::fa_dun,1.0);

    Real chi1incr = option[willmatch::chi1_increment]();
    Real chi2incr = option[willmatch::chi2_increment]();
    vector1<Real> CHI1,CHI2;
    Real ofst1 = floor(chi1incr/2.0)+0.5, ofst2 = floor(chi2incr/2.0)+0.5;

    for(Real i = ofst1; i < 360.0; i+= chi1incr) CHI1.push_back( (i>180.0 ? i-360.0 : i));
    for(Real i = ofst2; i < 360.0; i+= chi2incr) CHI2.push_back( (i>180.0 ? i-360.0 : i));

    //utility::io::ozstream bhhout(option[out::file::o]()+"/"+utility::file_basename(infile)+".bhh_match");
    // setup HIS residues for checking
    Pose hse,hsd,glu,arg,asp;
    core::pose::make_pose_from_sequence(glu,"E"       ,*fa_residue_set,false);
    core::pose::make_pose_from_sequence(arg,"R"       ,*fa_residue_set,false);
    core::pose::make_pose_from_sequence(asp,"D"       ,*fa_residue_set,false);
    core::pose::remove_lower_terminus_type_from_pose_residue(glu,1);
    core::pose::remove_lower_terminus_type_from_pose_residue(arg,1);
    core::pose::remove_lower_terminus_type_from_pose_residue(asp,1);
    core::pose::remove_upper_terminus_type_from_pose_residue(glu,1);
    core::pose::remove_upper_terminus_type_from_pose_residue(arg,1);
    core::pose::remove_upper_terminus_type_from_pose_residue(asp,1);

    vector1<Size> scanres;
    if(option[willmatch::residues].user()) {
      TR<<"input scanres!!!!!!"<<std::endl;
      scanres = option[willmatch::residues]();
    } else {
      for(Size i = 1; i <= in_fa.n_residue(); ++i) {
        if(!in_fa.residue(i).has("N" )) { continue; }
        if(!in_fa.residue(i).has("CA")) { continue; }
        if(!in_fa.residue(i).has("C" )) { continue; }
        if(!in_fa.residue(i).has("O" )) { continue; }
        if(!in_fa.residue(i).has("CB")) { continue; }
        if(in_fa.residue(i).name3()=="PRO") { continue; }
        if(in_fa.residue(i).name3()=="GLY") { continue; }
        scanres.push_back(i);
      }
    }

    core::pose::Pose pose(native);

    vector1<Real> natsasa; {core::id::AtomID_Map<Real> atom_sasa;core::scoring::calc_per_atom_sasa( native, atom_sasa, natsasa, 2.0, false );}
    vector1<Size> cbnbrs(native.n_residue()); {
      for(vector1<Size>::const_iterator rit = scanres.begin(); rit != scanres.end(); ++rit) {
        cbnbrs[*rit] = 0;
				Vec P1 = pose.residue(*rit).xyz("CA");
				if(pose.residue(*rit).has("CB")) P1 = pose.residue(*rit).xyz("CB");
        for(vector1<Size>::const_iterator cit = scanres.begin(); cit != scanres.end(); ++cit) {
					Vec P2 = pose.residue(*cit).xyz("CA");
					if(pose.residue(*cit).has("CB")) P2 = pose.residue(*cit).xyz("CB");
          if( P1.distance_squared(P2) < 81.0 ) cbnbrs[*rit]++;
        }
        cbnbrs[*rit]--; // self nbr
      }
    }


    Real cbcb_thresh2 = 11.5*11.5;
    Real cbcg_thresh2 = 10.8*10.8;
    // Real cbcd_thresh2 = 9.3*9.3;
    Real cgcg_thresh2 = 9.9*9.9;
    Real cgcd_thresh2 = 8.5*8.5;
    Real cdcd_thresh2 = 7.0*7.0;
    core::pack::rotamers::SingleResidueRotamerLibraryCOP rlib( core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rtarg ) );
    core::pack::dunbrack::RotamerLibraryScratchSpace scratch;
    Size dhit=0,ehit=0;
    for(vector1<Size>::const_iterator rit = scanres.begin(); rit != scanres.end(); ++rit) {
      //if(natsasa[*rit] > 0.1) continue;
      if(cbnbrs[*rit]<15.0) continue;
      pose.replace_residue(*rit,arg.residue(1),true);
      for(vector1<Size>::const_iterator cit = scanres.begin(); cit != scanres.end(); ++cit) {
        if(*cit== *rit) continue;
        //if(natsasa[*cit] > 0.1) continue;
        if(cbnbrs[*cit]<15.0) continue;
        if( pose.residue(*rit).xyz(5).distance_squared( pose.residue(*cit).xyz(5) ) > cbcb_thresh2 ) continue;
        std::cerr<<*rit<<" "<<*cit<<" "<<dhit<<" "<<ehit<<std::endl;
        pose.replace_residue(*cit,glu.residue(1),true);
        for(Size ich1 = 1; ich1 <= CHI1.size(); ich1++) {
          pose.set_chi(1,*rit,CHI1[ich1]);
          if( pose.residue(*rit).xyz(6).distance_squared( pose.residue(*cit).xyz(5) ) > cbcg_thresh2 ) continue;
          if( !clashcheck.clash_check( pose.residue(*rit).xyz(6),*rit) ) continue;
          for(Size jch1 = 1; jch1 <= CHI1.size(); jch1++) {
            pose.set_chi(1,*cit,CHI1[jch1]);
            Vec const CG = pose.residue(*cit).xyz(6);
            if( pose.residue(*rit).xyz(6).distance_squared( CG ) > cgcg_thresh2 ) continue;
            if( !clashcheck.clash_check( CG,*cit) ) continue;
            for(Size ich2 = 1; ich2 <= CHI2.size(); ich2++) {
              pose.set_chi(2,*rit,CHI2[ich2]);
              if( pose.residue(*rit).xyz(7).distance_squared( CG ) > cgcd_thresh2 ) continue;
              if( !clashcheck.clash_check( pose.residue(*rit).xyz(7),*rit) ) continue;
              for(Size jch2 = 1; jch2 <= CHI2.size(); jch2++) {
                pose.set_chi(2,*cit,CHI2[jch2]);
                Vec const CD = pose.residue(*cit).xyz(7);
                if( pose.residue(*rit).xyz(7).distance_squared( CD ) > cdcd_thresh2 ) continue;
                if( !clashcheck.clash_check( CD,*cit) ) continue;
                for(Size ich3 = 1; ich3 <= CHI2.size(); ich3++) {
                  pose.set_chi(3,*rit,CHI2[ich3]);
                  if( !clashcheck.clash_check( pose.residue(*rit).xyz(8),*rit) ) continue;
                  for(Size ich4 = 1; ich4 <= CHI2.size(); ich4++) {
                    pose.set_chi(4,*rit,CHI2[ich4]);
                    Vec const CZ = pose.residue(*rit).xyz(9);
                    Real const dg = CZ.distance_squared(CG);
                    Real const dd = CZ.distance_squared(CD);
                    if( 17.64 < dg && 17.64 < dd ) continue;
                    if( !clashcheck.clash_check( pose.residue(*rit).xyz(9),*rit) ) continue;
                    if( !clashcheck.clash_check( pose.residue(*rit).xyz(10),*rit) ) continue;
                    if( !clashcheck.clash_check( pose.residue(*rit).xyz(11),*rit) ) continue;
                    Vec CB = pose.residue(*cit).xyz(5);
                    Vec NZ = pose.residue(*rit).xyz(8);
										Vec NH1 = pose.residue(*rit).xyz(10);
                    Vec oz  = (CZ-NZ).normalized();
                    Vec oz2 = (CZ-NH1).normalized();
                    if( jch2==1 && 14.44 < dg && dg < 17.64 ) { // match CG -> ASP
                      Vec og = (CG-CB).normalized();
                      if( og.dot(oz ) > -0.95 && og.dot(oz2) > -0.95 ) continue;
                      if( og.dot(oz2) > -0.95 && projperp(oz ,CG).distance( projperp(oz ,CZ) ) > 0.3 ) continue;
                      if( og.dot(oz ) > -0.95 && projperp(oz2,CG).distance( projperp(oz2,CZ) ) > 0.3 ) continue;
                      if( rlib->rotamer_energy(pose.residue(*rit),scratch) > 8.0 ) continue;
                      dhit++;
                      TR<<"HIT D "<<*rit<<" "<<*cit<<" "<<CHI1[ich1]<<" "<<CHI2[ich2]<<" "<<CHI2[ich3]<<" "<<CHI2[ich4]<<" "<<CHI1[jch1]<<std::endl;
                      TR.flush();
											Pose tmp = pose;
											tmp.replace_residue(*cit,asp.residue(1),true);
											//											Mat R = rotation_matrix_degrees(og,ang);


											utility::io::ozstream o("D_"+str(*rit)+"-"+str(*cit)+"-"+str(ich1)+"-"+str(jch1)+"-"+str(ich2)+"-"+str(ich3)+"-"+str(ich4)+".pdb.gz");
											Size ano = 1;
											core::io::pdb::dump_pdb_residue(tmp.residue(*rit),ano,o);
											core::io::pdb::dump_pdb_residue(tmp.residue(*cit),ano,o);
											o.close();
                    }
                    if( 14.44 < dd && dd < 17.64 ) { // match CD -> GLU
                      Vec od = (CD-CG).normalized();
                      if( od.dot(oz ) > -0.95 && od.dot(oz2) > -0.95 ) continue;
                      if( od.dot(oz2) > -0.95 && projperp(oz ,CD).distance( projperp(oz ,CZ) ) > 0.3 ) continue;
                      if( od.dot(oz ) > -0.95 && projperp(oz2,CD).distance( projperp(oz2,CZ) ) > 0.3 ) continue;
                      if( rlib->rotamer_energy(pose.residue(*rit),scratch) > 8.0 ) continue;
                      ehit++;
                      TR<<"HIT E "<<*rit<<" "<<*cit<<" "<<CHI1[ich1]<<" "<<CHI2[ich2]<<" "<<CHI2[ich3]<<" "<<CHI2[ich4]<<" "<<CHI1[jch1]<<" "<<CHI2[jch2]<<std::endl;
                      TR.flush();
											utility::io::ozstream o("E_"+str(*rit)+"-"+str(*cit)+"-"+str(ich1)+"-"+str(jch1)+"-"+str(ich2)+"-"+str(jch2)+"-"+str(ich3)+"-"+str(ich4)+".pdb.gz");
											Size ano = 1;
											core::io::pdb::dump_pdb_residue(pose.residue(*rit),ano,o);
											core::io::pdb::dump_pdb_residue(pose.residue(*cit),ano,o);
											o.close();
                    }
                    //                Real rdun = rlib->rotamer_energy(R,scratch);
                    // if(rdun < 3.0) {
                    //  utility::io::ozstream out("test.pdb");
                    //  Size N = 1;
                    //  core::io::pdb::dump_pdb_residue(R,N,out);
                    //  out.close();
                    //
                  }
                }
              }
            }
          }
        }
      }
    }
    utility_exit_with_message("test backup");
  }

}


int main (int argc, char *argv[]) {

	try {

  devel::init(argc,argv);
  run();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
