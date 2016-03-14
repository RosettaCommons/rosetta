// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/*
WISH LIST
	buried unsats
		scan for rotamers that cna make H-bonds
		detect and penalize missing BB density
	strand pairs
	centriod score components ?
p	disulfide-compatible positions
	low-res sc? maybe patchdock-like complimentarity measure?
	statistically "good" BB xforms? all of by good contacts?
IGNORE
	contacts by SS ?? avg. degree probably sufficient
DONE
	termini distance
	contacts weight by avg. deg.
*/

#include <protocols/sic_dock/SICFast.hh>
#include <protocols/sic_dock/RigidScore.hh>
#include <protocols/sic_dock/util.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/sasa.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <devel/init.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>

#include <numeric/xyzTransform.hh>


// #include <time.h>
// #ifdef __MACH__
// #include <mach/mach_time.h>
// #endif
// Real time_highres() {
// #ifdef __MACH__
//   mach_timebase_info_data_t info;
//   mach_timebase_info(&info);
//   return mach_absolute_time() / 1000000000.0;
//   //uint64_t duration = mach_absolute_time();
//   //duration *= info.numer;
//   //duration /= info.denom;
// #else
//   timespec tp;
//   clock_gettime(CLOCK_REALTIME, &tp);
//   return tp.tv_sec + tp.tv_nsec/1000000000.0;
// #endif
//   return 0;
// }


#ifdef USE_OPENMP
#include <omp.h>
#endif
using core::Size;
using core::Real;
using core::pose::Pose;
using std::string;
using utility::vector1;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using ObjexxFCL::format::RJ;
using numeric::min;
using numeric::max;
using std::cout;
using std::cerr;
using std::endl;
typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;
typedef numeric::xyzVector<Real> Vecf;
typedef numeric::xyzMatrix<Real> Matf;
static THREAD_LOCAL basic::Tracer TR( "symdock_enum" );
OPT_1GRP_KEY( Real , tcdock, intra )
OPT_1GRP_KEY( Real , tcdock, intra1 )
OPT_1GRP_KEY( Real , tcdock, intra2 )
OPT_1GRP_KEY( Real , tcdock, termini_weight )
OPT_1GRP_KEY( Real , tcdock, termini_cutoff )
OPT_1GRP_KEY( Real , tcdock, termini_cutoff_short )
OPT_1GRP_KEY( Integer , tcdock, termini_trim )
OPT_1GRP_KEY( Integer , tcdock, nsamp1 )
OPT_1GRP_KEY( Integer , tcdock, topx )
OPT_1GRP_KEY( Integer , tcdock, peak_grid_size)
OPT_1GRP_KEY( Integer , tcdock, peak_grid_smooth)
OPT_1GRP_KEY( Boolean , tcdock, reverse )
OPT_1GRP_KEY( Boolean , tcdock, dump_pdb )
OPT_1GRP_KEY( Boolean , tcdock, dump_gz )
OPT_1GRP_KEY( Boolean , tcdock, dump_pdb_primary_sub )
OPT_1GRP_KEY( Boolean , tcdock, separate_components )
OPT_1GRP_KEY( IntegerVector , tcdock, dump_pdb_grid )
OPT_1GRP_KEY( Real , tcdock, grid_delta_rot )
OPT_1GRP_KEY( Real , tcdock, grid_delta_disp )
OPT_1GRP_KEY( Integer , tcdock, grid_size_rot )
OPT_1GRP_KEY( Integer , tcdock, grid_size_disp )
OPT_1GRP_KEY( IntegerVector , tcdock, dump_peak_grids )
OPT_1GRP_KEY( IntegerVector , tcdock, justone )
OPT_1GRP_KEY( FileVector, tcdock, I2 )
OPT_1GRP_KEY( FileVector, tcdock, I3 )
OPT_1GRP_KEY( FileVector, tcdock, I5 )
OPT_1GRP_KEY( FileVector, tcdock, O2 )
OPT_1GRP_KEY( FileVector, tcdock, O3 )
OPT_1GRP_KEY( FileVector, tcdock, O4 )
OPT_1GRP_KEY( FileVector, tcdock, T2 )
OPT_1GRP_KEY( FileVector, tcdock, T3 )
OPT_1GRP_KEY( Boolean , tcdock, usen1 )
OPT_1GRP_KEY( Boolean , tcdock, usec1 )
OPT_1GRP_KEY( Boolean , tcdock, require_exposed_termini )
OPT_1GRP_KEY( Boolean , tcdock, dry_run )
OPT_1GRP_KEY( Real , tcdock, term_min_expose )
OPT_1GRP_KEY( Real , tcdock, term_max_angle )
OPT_1GRP_KEY( String , tcdock, clash_atoms )
OPT_1GRP_KEY( Real , tcdock, redundant_angle )
// OPT_1GRP_KEY( Real , tcdock, redundant_dist )
OPT_1GRP_KEY( Boolean , tcdock, fast_stage_one )
OPT_1GRP_KEY( Integer , tcdock, max_linker_len )
OPT_1GRP_KEY( Integer , tcdock, linker_lookup_radius )
OPT_1GRP_KEY( Integer , tcdock, max_res )
OPT_1GRP_KEY( Boolean , tcdock, cb_weight_secstruct )
OPT_1GRP_KEY( Boolean , tcdock, cb_weight_average_degree )
OPT_1GRP_KEY( Real    , tcdock, trim_floppy_termini )
OPT_1GRP_KEY( Boolean , tcdock, debug )


void register_options() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		OPT( in::file::s );
		NEW_OPT( tcdock::intra            ,"weight for intra contacts"           ,  1.0   );
		NEW_OPT( tcdock::intra1           ,"weight for comp1 intra contacts"     ,  1.0   );
		NEW_OPT( tcdock::intra2           ,"weight for conp2 intra contacts"     ,  1.0   );
		NEW_OPT( tcdock::nsamp1           ,"check around X top lowres hits"      , 2000   );
		NEW_OPT( tcdock::topx             ,"output top X hits"                   , 10     );
		NEW_OPT( tcdock::peak_grid_size   ,"peak detect grid size (2*N+1)"       , 24     );
		NEW_OPT( tcdock::peak_grid_smooth ,"peak detect grid smooth (0+)"        ,  1     );
		NEW_OPT( tcdock::reverse          ,"reverse one component"               , false  );
		NEW_OPT( tcdock::dump_pdb         ,"dump pdbs"                           , false  );
		NEW_OPT( tcdock::dump_gz          ,"dump gz pdbs"                        , true   );
		NEW_OPT( tcdock::dump_pdb_primary_sub ,"dump pdb primary subunit only"   , false  );
		NEW_OPT( tcdock::separate_components ,"different chains per subsub"   , true  );
		NEW_OPT( tcdock::dump_pdb_grid    ,"dump pdb design grids"               , -1  );
		NEW_OPT( tcdock::grid_delta_rot   ,"grid samp size theta"                ,  0.5   );
		NEW_OPT( tcdock::grid_delta_disp  ,"grid samp size r"                    ,  0.5   );
		NEW_OPT( tcdock::grid_size_rot    ,"grid nsamp theta +- (2*N+1 samp)"    ,  1     );
		NEW_OPT( tcdock::grid_size_disp   ,"grid samp r"                         ,  5     );
		NEW_OPT( tcdock::dump_peak_grids  ,"dump specific grids grids"           , -1     );
		NEW_OPT( tcdock::justone          ,"dump info on one structure"          , -1     );
		NEW_OPT( tcdock::I2               ,"file(s) for icos 2fold dimer"        , ""     );
		NEW_OPT( tcdock::I3               ,"file(s) for icos 3fold trimer"       , ""     );
		NEW_OPT( tcdock::I5               ,"file(s) for icos 5fold pentamer"     , ""     );
		NEW_OPT( tcdock::O2               ,"file(s) for octa 2fold dimer"        , ""     );
		NEW_OPT( tcdock::O3               ,"file(s) for octa 3fold trimer"       , ""     );
		NEW_OPT( tcdock::O4               ,"file(s) for octa 4fold tetramer"     , ""     );
		NEW_OPT( tcdock::T2               ,"file(s) for tetr 2fold dimer"        , ""     );
		NEW_OPT( tcdock::T3               ,"file(s) for tetr 3fold trimer"       , ""     );
		NEW_OPT( tcdock::termini_cutoff   ,"tscore = w*max(0,cut-x)"             , 20.0   );
		NEW_OPT( tcdock::termini_cutoff_short   ,"tscore = w*max(0,cut-x)"       ,  0.0   );
		NEW_OPT( tcdock::termini_weight   ,"tscore = w*max(0,cut-x)"             ,  0.0   );
		NEW_OPT( tcdock::termini_trim     ,"trim termini up to for termini score",  0     );
		NEW_OPT( tcdock::usec1            ,"use comp1 cterm"                     , true  );
		NEW_OPT( tcdock::usen1            ,"use comp1 nterm"                     , true  );
		NEW_OPT( tcdock::require_exposed_termini,""         , false );
		NEW_OPT( tcdock::dry_run          ,"no calculations, just setup"         , false );
		NEW_OPT( tcdock::term_max_angle ,""             , 45.0  );
		NEW_OPT( tcdock::term_min_expose ,""            , 0.1  );
		NEW_OPT( tcdock::clash_atoms     ,""            , "BB"  );
		NEW_OPT( tcdock::redundant_angle ,""            , 0  );
		// NEW_OPT( tcdock::redundant_dist  ,""            , 0  );
		NEW_OPT( tcdock::fast_stage_one  ,"faster stage one, may miss some"   , false  );
		NEW_OPT( tcdock::max_linker_len  ,""   ,  10    );
		NEW_OPT( tcdock::linker_lookup_radius  ,""   , 1 );
		NEW_OPT( tcdock::max_res  ,""   , 99999 );
		NEW_OPT( tcdock::cb_weight_secstruct     ,""         , true );
		NEW_OPT( tcdock::cb_weight_average_degree,""         , true );
		NEW_OPT( tcdock::trim_floppy_termini ,"", 0.0 );
		NEW_OPT( tcdock::debug,""         , false );


}

inline
numeric::xyzTransform<core::Real> const
toXform(core::kinematics::Stub const & stub){
	return numeric::xyzTransform<core::Real>(stub.M,stub.v);
}

template<typename T> inline T sqr(T x) { return x*x; }
inline Real sigmoid( Real const & sqdist, Real const & start, Real const & stop ) {
	if( sqdist > stop*stop ) {
		return 0.0;
	} else if( sqdist < start*start ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return (stop-dist)/(stop-start);
		//return sqr(1.0	- sqr( (dist - start) / (stop - start) ) );
	}
}

void dump_points_pdb(utility::vector1<Vecf> const & p, std::string fn) {
	using namespace ObjexxFCL::format;
	std::ofstream o(fn.c_str());
	for(Size i = 1; i <= p.size(); ++i) {
		std::string rn = "VIZ";
		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x())<<F(8,3,p[i].y())<<F(8,3,p[i].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}
void dump_points_pdb(utility::vector1<Vecf> const & p, Vec t, std::string fn) {
	using namespace ObjexxFCL::format;
	std::ofstream o(fn.c_str());
	for(Size i = 1; i <= p.size(); ++i) {
		std::string rn = "VIZ";
		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x()+t.x())<<F(8,3,p[i].y()+t.y())<<F(8,3,p[i].z()+t.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}
// void dump_points_pdb(numeric::geometry::hashing::xyzStripeHashWithMeta<Real>::float4 const *p, int n, Vec t, std::string fn) {
// 	using namespace ObjexxFCL::format;
// 	std::ofstream o(fn.c_str());
// 	for(Size i = 0; i < n; ++i) {
// 		std::string rn = "VIZ";
// 		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x-t.x())<<F(8,3,p[i].y-t.y())<<F(8,3,p[i].z-t.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
// 	}
// 	o.close();
// }
// void dump_points_pdb(numeric::geometry::hashing::xyzStripeHashWithMeta<Real>::float4 const *p, int n, std::string fn) {
// 	using namespace ObjexxFCL::format;
// 	std::ofstream o(fn.c_str());
// 	for(Size i = 0; i < n; ++i) {
// 		std::string rn = "VIZ";
// 		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x)<<F(8,3,p[i].y)<<F(8,3,p[i].z)<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
// 	}
// 	o.close();
// }

void xform_pose( core::pose::Pose & pose, core::kinematics::Stub const & s, Size sres=1, Size eres=0 ) {
  if(eres==0) eres = pose.n_residue();
  for(Size ir = sres; ir <= eres; ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
    }
  }
}
void xform_pose_rev( core::pose::Pose & pose, core::kinematics::Stub const & s ) {
  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s.global2local(pose.xyz(aid)) );
    }
  }
}


void trans_pose( Pose & pose, Vecf const & trans, Size start=1, Size end=0 ) {
	if(0==end) end = pose.n_residue();
	for(Size ir = start; ir <= end; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, pose.xyz(aid) + (Vec)trans );
		}
	}
}
void rot_pose( Pose & pose, Mat const & rot, Size start=1, Size end=0 ) {
	if(0==end) end = pose.n_residue();
	for(Size ir = start; ir <= end; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}
void rot_pose( Pose & pose, Mat const & rot, Vecf const & cen, Size start=1, Size end=0 ) {
	trans_pose(pose,-cen,start,end);
	rot_pose(pose,rot,start,end);
	trans_pose(pose,cen,start,end);
}
void rot_pose( Pose & pose, Vecf const & axis, Real const & ang, Size start=1, Size end=0 ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),start,end);
}
void rot_pose( Pose & pose, Vecf const & axis, Real const & ang, Vecf const & cen, Size start=1, Size end=0 ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen,start,end);
}
void alignaxis(Pose & pose, Vecf newaxis, Vecf oldaxis, Vecf cen = Vecf(0,0,0) ) {
	newaxis.normalize();
	oldaxis.normalize();
	if(fabs(newaxis.dot(oldaxis)) < 0.9999) {
		Vecf axis = newaxis.cross(oldaxis).normalized();
		Real ang = (Real)-acos(numeric::max((Real)-1.0,numeric::min((Real)1.0,newaxis.dot(oldaxis))))*(Real)180.0/numeric::constants::f::pi;
		rot_pose(pose,axis,ang,cen);
	}
}
inline Vec projperp(Vec const & u, Vec const & v) {
  return v - projection_matrix(u)*v;
}
void prune_cb_pairs(vector1<Vecf> & cba, vector1<Vecf> & cbb, vector1<Real> & wa_in, vector1<Real> & wb_in, Real CTD2) {
	vector1<Vecf> a,b;
	vector1<Real> wa,wb;
	vector1<Real>::const_iterator iwa = wa_in.begin();
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia,++iwa) {
		vector1<Real>::const_iterator iwb = wb_in.begin();
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib,++iwb) {
			if( ib->distance_squared(*ia) < CTD2 ) {
				a.push_back(*ia);
				b.push_back(*ib);
				wa.push_back(*iwa);
				wb.push_back(*iwb);
			}
		}
	}
	cba = a;
	cbb = b;
	wa_in = wa;
	wb_in = wb;
}
int neighbor_count(Pose const &pose, int ires, Real distance_threshold=10.0) {
	core::conformation::Residue const resi( pose.residue( ires ) );
	Size resi_neighbors( 0 );
	for(Size jres = 1; jres <= pose.n_residue(); ++jres) {
		core::conformation::Residue const resj( pose.residue( jres ) );
		Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
		if( distance <= distance_threshold ){
			++resi_neighbors;
		}
	}
	return resi_neighbors;
}
Real brute_mindis(vector1<Vecf> const & pa, vector1<Vecf> const & pb, Vecf const ofst) {
	Real mindis = 9e9;
	for(vector1<Vecf>::const_iterator i = pa.begin(); i != pa.end(); ++i) {
		for(vector1<Vecf>::const_iterator j = pb.begin(); j != pb.end(); ++j) {
			mindis = min( mindis, i->distance_squared(*j+ofst) );
		}
	}
	return mindis;
}
struct LMAX {
	Real score,radius;
	int icmp2,icmp1,iori;
	LMAX() :	score(0),radius(0),icmp2(0),icmp1(0),iori(0) {}
	LMAX(Real _score, Real _radius, int _icmp2, int _icmp1, int _iori) :
		score(_score),radius(_radius),icmp2(_icmp2),icmp1(_icmp1),iori(_iori) {}
};
int compareLMAX(const LMAX a,const LMAX b) {
	return a.score > b.score;
}


struct TCDock {
	vector1<Real> cmp2mnpos_,cmp1mnpos_,cmp2mnneg_,cmp1mnneg_,cmp2dspos_,cmp1dspos_,cmp2dsneg_,cmp1dsneg_;
	ObjexxFCL::FArray2D<Real> cmp2cbpos_,cmp1cbpos_,cmp2cbneg_,cmp1cbneg_;
	ObjexxFCL::FArray3D<Real> gradii,gscore;
	Vecf cmp1axs_,cmp2axs_;
	Real alpha_,sin_alpha_,tan_alpha_;
	Real cmp1diapos_,cmp1dianeg_,cmp2diapos_,cmp2dianeg_;
	Real cmp1scale_,cmp2scale_;
	core::pose::Pose cmp1in_,cmp2in_;
	vector1<Vecf> cmp1pts_,cmp1cbs_,cmp2pts_,cmp2cbs_;
	vector1<Real> cmp1wts_,cmp2wts_;
	core::id::AtomID_Map<Real> cmp1wtmap_,cmp2wtmap_;
	string cmp1name_,cmp2name_,cmp1type_,cmp2type_,symtype_;
	int cmp1nangle_,cmp2nangle_,cmp1nsub_,cmp2nsub_;
	std::map<string,Vecf> axismap_;
	protocols::sic_dock::SICFast sic_;
	bool abort_;
	core::id::AtomID_Map<core::Real> clashmap1_,clashmap2_;
	protocols::sic_dock::RigidScoreCOP rigid_sfxn_;
	protocols::sic_dock::LinkerScoreCOP lnscore_;
	Size conf_count_;
	Real start_time_;

	TCDock(
		Pose const & pose1,
		Pose const & pose2,
		string cmp1pdb,
		string cmp2pdb,
		string cmp1type,
		string cmp2type
	):
		cmp1type_(cmp1type),
		cmp2type_(cmp2type),
		abort_(false),
		conf_count_(0)
	{
		#ifdef USE_OPENMP
		cout << "OMP info: " << num_threads() << " " << thread_num() << endl;
		#endif
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		cmp1in_ = pose1;
		cmp2in_ = pose2;
		// core::chemical::ResidueTypeSetCAP crs=core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD);
		// core::import_pose::pose_from_file(cmp1in_,*crs,cmp1pdb, core::import_pose::PDB_file);
		// core::import_pose::pose_from_file(cmp2in_,*crs,cmp2pdb, core::import_pose::PDB_file);

		cmp1diapos_=0.0,cmp1dianeg_=0.0,cmp2diapos_=0.0,cmp2dianeg_=0.0;
		for(Size i = 1; i <= cmp1in_.n_residue(); ++i) {
			int const natom = (cmp1in_.residue(i).name3()=="GLY") ? 4 : 5;
			for(int j = 1; j <= natom; ++j) {
				Vec const & v( cmp1in_.residue(i).xyz(j) );
				cmp1diapos_ = max(cmp1diapos_, v.z());
				cmp1dianeg_ = max(cmp1dianeg_,-v.z());
			}
		}
		for(Size i = 1; i <= cmp2in_.n_residue(); ++i) {
			int const natom = (cmp2in_.residue(i).name3()=="GLY") ? 4 : 5;
			for(int j = 1; j <= natom; ++j) {
				Vec const & v( cmp2in_.residue(i).xyz(j) );
				cmp2diapos_ = max(cmp2diapos_, v.z());
				cmp2dianeg_ = max(cmp2dianeg_,-v.z());
			}
		}

		if(cmp1type[1]=='2') { make_dimer   (cmp1in_); cmp1nsub_ = 2; }
		if(cmp1type[1]=='3') { make_trimer  (cmp1in_); cmp1nsub_ = 3; }
		if(cmp1type[1]=='4') { make_tetramer(cmp1in_); cmp1nsub_ = 4; }
		if(cmp1type[1]=='5') { make_pentamer(cmp1in_); cmp1nsub_ = 5; }
		if(cmp2type[1]=='2') { make_dimer   (cmp2in_); cmp2nsub_ = 2; }
		if(cmp2type[1]=='3') { make_trimer  (cmp2in_); cmp2nsub_ = 3; }
		if(cmp2type[1]=='4') { make_tetramer(cmp2in_); cmp2nsub_ = 4; }
		if(cmp2type[1]=='5') { make_pentamer(cmp2in_); cmp2nsub_ = 5; }
		if(option[tcdock::reverse]()) rot_pose(cmp1in_,Vecf(0,1,0),180.0);
		core::scoring::dssp::Dssp(cmp1in_).insert_ss_into_pose(cmp1in_);
		core::scoring::dssp::Dssp(cmp2in_).insert_ss_into_pose(cmp2in_);
		if(option[tcdock::debug]()) cmp1in_.dump_pdb("cmp1in.pdb");
		if(option[tcdock::debug]()) cmp2in_.dump_pdb("cmp2in.pdb");

		// for(int i = 1; i <= cmp1in_.n_residue(); ++i) cout << cmp1in_.secstruct(i); cout << endl;
		// for(int i = 1; i <= cmp2in_.n_residue(); ++i) cout << cmp2in_.secstruct(i); cout << endl;

		bool nt1good=1,nt2good=1,ct1good=1,ct2good=1;
		// protocols::sic_dock::termini_exposed(cmp1in_,nt1good,ct1good);
		// protocols::sic_dock::termini_exposed(cmp2in_,nt2good,ct2good);
		cout << cmp1pdb << " nterm 1: " << nt1good << endl;
		cout << cmp1pdb << " cterm 1: " << ct1good << endl;
		cout << cmp2pdb << " nterm 2: " << nt2good << endl;
		cout << cmp2pdb << " cterm 2: " << ct2good << endl;
		if(!option[tcdock::usec1].user()) {
			if( ct1good && nt2good ) option[tcdock::usec1](true);
			else                     option[tcdock::usec1](false);
		}
		if(!option[tcdock::usen1].user()) {
			if( ct2good && nt1good ) option[tcdock::usen1](true);
			else                     option[tcdock::usen1](false);
		}
		cout << "TERMINI " << cmp1pdb << " " << cmp2pdb << " " << option[tcdock::usec1]() << " " << option[tcdock::usen1]() << endl;
		if( !option[tcdock::usec1]() && !option[tcdock::usen1]() ) {
			if(option[tcdock::require_exposed_termini]()) {
				std::cout << "no compatible good termini" << std::endl;
				abort_ = true;
				return;
			}
		}


		cmp1name_ = utility::file_basename(cmp1pdb);
		cmp2name_ = utility::file_basename(cmp2pdb);
		if(cmp1name_.substr(cmp1name_.size()-3)==".gz" ) cmp1name_ = cmp1name_.substr(0,cmp1name_.size()-3);
		if(cmp2name_.substr(cmp2name_.size()-3)==".gz" ) cmp2name_ = cmp2name_.substr(0,cmp2name_.size()-3);
		if(cmp1name_.substr(cmp1name_.size()-4)==".pdb") cmp1name_ = cmp1name_.substr(0,cmp1name_.size()-4);
		if(cmp2name_.substr(cmp2name_.size()-4)==".pdb") cmp2name_ = cmp2name_.substr(0,cmp2name_.size()-4);
		symtype_ = cmp2type.substr(0,1);

		std::map<string,int> comp_nangle;
		comp_nangle["I2"] = 180; comp_nangle["O2"] = 180; comp_nangle["T2"] = 180;
		comp_nangle["I3"] = 120; comp_nangle["O3"] = 120; comp_nangle["T3"] = 120;
		comp_nangle["I5"] =  72;
		comp_nangle["O4"] =  90;
		axismap_["I2"] = Vecf( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized();
		axismap_["I3"] = Vecf( 0.93417235896272, 0.00000000000000, 0.35682208977309 ).normalized();
		axismap_["I5"] = Vecf( 0.85065080835204, 0.52573111211914, 0.00000000000000 ).normalized();

		axismap_["O2"] = Vecf( 1.00000000000000, 1.00000000000000, 0.00000000000000 ).normalized();
		axismap_["O3"] = Vecf( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized();
		axismap_["O4"] = Vecf( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized();

		axismap_["T2"] = Vecf( 0.81649657940871, 0.00000000000000, 0.57735027133783 ).normalized();
		axismap_["T3"] = Vecf( 0.00000000000000, 0.00000000000000, 1.00000000000000 ).normalized();

		cmp1axs_ = axismap_[cmp1type_];
		cmp2axs_ = axismap_[cmp2type_];
		cmp1nangle_ = comp_nangle[cmp1type_];
		cmp2nangle_ = comp_nangle[cmp2type_];
		cout << cmp1type_ << " " << cmp2type_ << " nang1 " << cmp1nangle_ << " nang2 " << cmp2nangle_ << endl;

		alpha_ = angle_degrees(cmp1axs_,Vecf(0,0,0),cmp2axs_);
		sin_alpha_ = sin(numeric::conversions::radians(alpha_));
		tan_alpha_ = tan(numeric::conversions::radians(alpha_));

		alignaxis(cmp1in_,cmp1axs_,Vec(0,0,1),Vec(0,0,0));
		alignaxis(cmp2in_,cmp2axs_,Vec(0,0,1),Vec(0,0,0));

		// COMPUTE WEIGHTS BB, AND CBS
		using core::id::AtomID;
		core::pose::initialize_atomid_map(cmp1wtmap_,cmp1in_, 0.0);
		core::pose::initialize_atomid_map(cmp2wtmap_,cmp2in_, 0.0);
		core::pose::initialize_atomid_map(clashmap1_,cmp1in_,-1.0);
		core::pose::initialize_atomid_map(clashmap2_,cmp2in_,-1.0);
		for(Size i12 = 0; i12 < 2; ++i12){
			Pose const & ptmp( i12?cmp1in_:cmp2in_ );
			for(Size i = 1; i <= ptmp.n_residue(); ++i) {
				if(ptmp.residue(i).has("CB")) {
					Real wt = 1.0;
					if( option[tcdock::cb_weight_average_degree]()                          ) wt *= min(1.0,(Real)neighbor_count(ptmp,i)/20.0);
					if( option[tcdock::cb_weight_secstruct     ]() && ptmp.secstruct(i)=='L') wt /= 3.0; //TODO make option somehow
					(i12?cmp1cbs_:cmp2cbs_).push_back(Vecf(ptmp.residue(i).xyz("CB")));
					(i12?cmp1wts_:cmp2wts_).push_back(wt);
					(i12?cmp1wtmap_:cmp2wtmap_)[AtomID(ptmp.residue(i).atom_index("CB"),i)] = wt;
				}
				Size natom = 0;
				if( "BB"    == option[tcdock::clash_atoms]() ) natom = (ptmp.residue(i).name3()=="GLY") ? 4 : 5;
				if( "HEAVY" == option[tcdock::clash_atoms]() ) natom = ptmp.residue(i).nheavyatoms();
				for(Size j = 1; j <= natom; ++j) {
					Vec const & v( ptmp.residue(i).xyz(j) );
					(i12?clashmap1_:clashmap2_)[AtomID(j,i)] = ptmp.residue(i).atom_type(j).lj_radius();
					(i12?cmp1pts_:cmp1pts_).push_back(v);
				}
			}
		}

		sic_.init(cmp1in_,cmp2in_,clashmap1_,clashmap2_);
		Real CTD(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]());
		Real CLD(basic::options::option[basic::options::OptionKeys::sicdock::clash_dis]());

		TR << "creating score functions" << std::endl;
		protocols::sic_dock::JointScoreOP  jtscore = new protocols::sic_dock::JointScore;
		protocols::sic_dock::CBScoreCOP cbscore = new protocols::sic_dock::CBScore(cmp1in_,cmp2in_,CLD,CTD,cmp1wtmap_,cmp2wtmap_);
		jtscore->add_score(cbscore ,1.0);
		if(option[tcdock::max_linker_len]() > 0){
			string lnscore_tag = cmp1name_+"_"+cmp2name_+"_"+symtype_+(option[tcdock::reverse]()?"R":"F")+"_";
			lnscore_ = new protocols::sic_dock::LinkerScore(cmp1in_,cmp2in_,option[tcdock::max_linker_len](),option[tcdock::linker_lookup_radius](),lnscore_tag);
			jtscore->add_score(lnscore_,1.0);
		}
		rigid_sfxn_ = jtscore;

		cmp1mnpos_.resize(cmp1nangle_,0.0);
		cmp2mnpos_.resize(cmp2nangle_,0.0);
		cmp1mnneg_.resize(cmp1nangle_,0.0);
		cmp2mnneg_.resize(cmp2nangle_,0.0);
		cmp1dspos_.resize(cmp1nangle_,0.0);
		cmp2dspos_.resize(cmp2nangle_,0.0);
		cmp1dsneg_.resize(cmp1nangle_,0.0);
		cmp2dsneg_.resize(cmp2nangle_,0.0);
		cmp1cbpos_.dimension(cmp1nangle_,200,0.0);
		cmp2cbpos_.dimension(cmp2nangle_,200,0.0);
		cmp1cbneg_.dimension(cmp1nangle_,200,0.0);
		cmp2cbneg_.dimension(cmp2nangle_,200,0.0);
		gradii.dimension(cmp2nangle_,cmp1nangle_,360,-9e9);
		gscore.dimension(cmp2nangle_,cmp1nangle_,360,-9e9);


		core::scoring::methods::RG_Energy_Fast rgcalc;
		Real tmp1 = rgcalc.calculate_rg_score(cmp1in_);
		Real tmp2 = rgcalc.calculate_rg_score(cmp2in_);
		cmp1scale_ = 2.0*tmp1/(tmp1+tmp2);
		cmp2scale_ = 2.8*tmp2/(tmp1+tmp2);
	}
	virtual ~TCDock() {
	}
	int num_threads() {
		#ifdef USE_OPENMP
			return omp_get_max_threads();
		#else
			return 1;
		#endif
	}
	int thread_num() {
		#ifdef USE_OPENMP
			return omp_get_thread_num();
		#else
			return 0;
		#endif
	}
	Vecf swap_axis(vector1<Vecf> & pts,string type) {
		Vecf other_axis;
		if( type[1] == '2' ) {
			Mat r = rotation_matrix_degrees( axismap_[type[0]+string("3")],120.0);
			other_axis = r*axismap_[type];
			for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) *i = r*(*i);
		} else {
			Mat r = rotation_matrix_degrees( axismap_[type[0]+string("2")],180.0);
			other_axis = r*axismap_[type];
			for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) *i = r*(*i);
		}
		return other_axis.normalized();
	}
	Matf swap_axis_rotation(string type) {
		if( type[1] == '2' ) {
			return rotation_matrix_degrees( axismap_[type[0]+string("3")],120.0);
		} else {
			return rotation_matrix_degrees( axismap_[type[0]+string("2")],180.0);
		}
		return Matf::identity();
	}
	void precompute_intra() {
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		Real const CONTACT_D	= basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]();
		Real const CLASH_D    = basic::options::option[basic::options::OptionKeys::sicdock::	clash_dis]();
		Real const CONTACT_D2 = sqr(CONTACT_D);
		// Real const CLASH_D2	= sqr(CLASH_D);
		// compute high/low min dis for pent and cmp1 here, input to sicfast and don't allow any below
		cout << "precomputing one-component interactions every 1°" << endl;
		for(int i12 = 0; i12 < 2; ++i12) {
			Vecf const axis = (i12?cmp1axs_:cmp2axs_);
			vector1<Real>             & cmpwts  ( i12?cmp1wts_  :cmp2wts_   );
			vector1<Real>             & cmpmnpos( i12?cmp1mnpos_:cmp2mnpos_ );
			vector1<Real>             & cmpmnneg( i12?cmp1mnneg_:cmp2mnneg_ );
			vector1<Real>             & cmpdspos( i12?cmp1dspos_:cmp2dspos_ );
			vector1<Real>             & cmpdsneg( i12?cmp1dsneg_:cmp2dsneg_ );
			ObjexxFCL::FArray2D<Real> & cmpcbpos( i12?cmp1cbpos_:cmp2cbpos_ );
			ObjexxFCL::FArray2D<Real> & cmpcbneg( i12?cmp1cbneg_:cmp2cbneg_ );

			protocols::sic_dock::SICFast sic;
			sic.init( i12? cmp1in_   :cmp2in_,
				      i12? cmp1in_   :cmp2in_   ,
				      i12? clashmap1_:clashmap2_,
				      i12? clashmap1_:clashmap2_);
			protocols::sic_dock::CBScore sfxn(i12?cmp1in_:cmp2in_,i12?cmp1in_:cmp2in_,CLASH_D,CONTACT_D);

			for(int ipn = 0; ipn < 2; ++ipn) {
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int icmp = 0; icmp < (i12?cmp1nangle_:cmp2nangle_); icmp+=1) {
					vector1<Vecf> ptsA,cbA;
					if(i12) get_cmp1(icmp,ptsA,cbA); else get_cmp2(icmp,ptsA,cbA);
					vector1<Vecf> ptsB(ptsA),cbB(cbA);
					assert( cbA.size() == cmpwts.size() && cbB.size() == cmpwts.size() );
					/*              */ swap_axis(ptsB,i12?cmp1type_:cmp2type_);
					Vecf const axis2 = swap_axis( cbB,i12?cmp1type_:cmp2type_);
					Vecf const sicaxis = ipn ? (axis2-axis).normalized() : (axis-axis2).normalized();
					Real score = 0;
					core::kinematics::Stub xa,xb;
					xa.M = rotation_matrix_degrees(axis2,(Real)icmp) * swap_axis_rotation(i12?cmp1type_:cmp2type_);
					xb.M = rotation_matrix_degrees(axis ,(Real)icmp);
					Real const d = -protocols::sic_dock::slide_into_contact_and_score_DEPRICATED(sic,sfxn,xa,xb,sicaxis,score);
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for cmppos! "+ObjexxFCL::string_of(icmp));
					if(fabs(d) > 9e8){
						utility_exit_with_message("precompute_intra error");
					}
					(ipn?cmpmnpos:cmpmnneg)[icmp+1] = (ipn?-1.0:1.0) * d/2.0/sin( angle_radians(axis2,Vecf(0,0,0),axis)/2.0 );
					(ipn?cmpdspos:cmpdsneg)[icmp+1] = (ipn?-1.0:1.0) * d;
					for(vector1<Vecf>::iterator iv = cbB.begin(); iv != cbB.end(); ++iv) *iv = (*iv) - d*sicaxis;
					vector1<Real> wA(cmpwts),wB(cmpwts);
					prune_cb_pairs(cbA,cbB,wA,wB,CONTACT_D2); // mutates inputs!!!!!!!!!!
					// Real lastcbc = 9e9;
					for(int i = 1; i <= 200; ++i) {
						Real cbc = 0.0;
						vector1<Real>::const_iterator iwa=wA.begin(), iwb=wB.begin();
						for( vector1<Vecf>::const_iterator ia=cbA.begin(), ib=cbB.begin(); ia != cbA.end(); ++ia,++ib,++iwa,++iwb) {
							cbc += sigmoid( ia->distance_squared(*ib) , CLASH_D, CONTACT_D ) * (*iwa) * (*iwa);
						}
						// assert(lastcbc >= cbc);
						// lastcbc = cbc;
						(ipn?cmpcbpos:cmpcbneg)(icmp+1,i) = cbc;
						if(cbc==0.0) break;
						for(vector1<Vecf>::iterator iv = cbA.begin(); iv != cbA.end(); ++iv) *iv = (*iv) + (ipn?0.1:-0.1)*axis;
						for(vector1<Vecf>::iterator iv = cbB.begin(); iv != cbB.end(); ++iv) *iv = (*iv) + (ipn?0.1:-0.1)*axis2;
					}
				}
			}
		}
	}
	void dump_onecomp() {
		using	namespace	basic::options;
		using	namespace	basic::options::OptionKeys;
		using utility::file_basename;
		cout << "dumping 1D stats: " << option[out::file::o]()+"/"+cmp2name_+"_POS_1D.dat" << endl;
		{ utility::io::ozstream out(option[out::file::o]()+"/"+cmp2name_+"_POS_1D.dat");
			for(int i = 1; i <= cmp2nangle_; ++i) out << i << " " << cmp2dspos_[i] << " " << cmp2cbpos_(i,1) << " " << cmp2cbpos_(i,2) << " " << cmp2cbpos_(i,3) << " " << cmp2cbpos_(i,4) << endl;
			out.close(); }
		{ utility::io::ozstream out(option[out::file::o]()+"/"+cmp2name_+"_NEG_1D.dat");
			for(int i = 1; i <= cmp2nangle_; ++i) out << i << " " << cmp2dsneg_[i] << " " << cmp2cbneg_(i,1) << " " << cmp2cbneg_(i,2) << " " << cmp2cbneg_(i,3) << " " << cmp2cbneg_(i,4) << endl;
			out.close(); }
		{ utility::io::ozstream out(option[out::file::o]()+"/"+cmp1name_+"_POS_1D.dat");
			for(int i = 1; i <= cmp1nangle_; ++i) out << i << " " << cmp1dspos_[i] << " " << cmp1cbpos_(i,1) << " " << cmp1cbpos_(i,2) << " " << cmp1cbpos_(i,3) << " " << cmp1cbpos_(i,4) << endl;
			out.close(); }
		{ utility::io::ozstream out(option[out::file::o]()+"/"+cmp1name_+"_NEG_1D.dat");
			for(int i = 1; i <= cmp1nangle_; ++i) out << i << " " << cmp1dsneg_[i] << " " << cmp1cbneg_(i,1) << " " << cmp1cbneg_(i,2) << " " << cmp1cbneg_(i,3) << " " << cmp1cbneg_(i,4) << endl;
			out.close(); }
		utility_exit_with_message("1COMP!");
	}
	void get_cmp1(Real acmp1, vector1<Vecf> & pts, vector1<Vecf> & cbs ) {
		Matf R = rotation_matrix_degrees( cmp1axs_, acmp1 );
		pts = cmp1pts_;
		cbs = cmp1cbs_;
		for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) (*i) = R*(*i);
		for(vector1<Vecf>::iterator i = cbs.begin(); i != cbs.end(); ++i) (*i) = R*(*i);
	}
	void get_cmp2(Real acmp2, vector1<Vecf> & pts, vector1<Vecf> & cbs ){
		Matf R = rotation_matrix_degrees( cmp2axs_, acmp2 );
		pts = cmp2pts_;
		cbs = cmp2cbs_;
		for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) (*i) = R*(*i);
		for(vector1<Vecf>::iterator i = cbs.begin(); i != cbs.end(); ++i) (*i) = R*(*i);
	}
	void make_dimer(core::pose::Pose & pose) {
		// cerr << "make_dimer" << endl;
		core::pose::Pose t2(pose);
		rot_pose(t2,Vecf(0,0,1),180.0);
		for(Size i = 1; i <= t2.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
		// pose.dump_pdb("dimer.pdb");
	}
	void make_trimer(core::pose::Pose & pose) {
		// cerr << "make_trimer" << endl;
		core::pose::Pose t2(pose),t3(pose);
		rot_pose(t2,Vecf(0,0,1),120.0);
		rot_pose(t3,Vecf(0,0,1),240.0);
		for(Size i = 1; i <= t2.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
		for(Size i = 1; i <= t3.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t3.residue(i),1); else pose.append_residue_by_bond(t3.residue(i));
		// pose.dump_pdb("trimer.pdb");
	}
	void make_tetramer(core::pose::Pose & pose) {
		// cerr << "make_tetramer" << endl;
		core::pose::Pose t2(pose),t3(pose),t4(pose);
		rot_pose(t2,Vecf(0,0,1), 90.0);
		rot_pose(t3,Vecf(0,0,1),180.0);
		rot_pose(t4,Vecf(0,0,1),270.0);
		for(Size i = 1; i <= t2.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
		for(Size i = 1; i <= t3.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t3.residue(i),1); else pose.append_residue_by_bond(t3.residue(i));
		for(Size i = 1; i <= t4.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t4.residue(i),1); else pose.append_residue_by_bond(t4.residue(i));
	}
	void make_pentamer(core::pose::Pose & pose) {
		// cerr << "make_pentamer" << endl;
		core::pose::Pose t2(pose),t3(pose),t4(pose),t5(pose);
		rot_pose(t2,Vecf(0,0,1), 72.0);
		rot_pose(t3,Vecf(0,0,1),144.0);
		rot_pose(t4,Vecf(0,0,1),216.0);
		rot_pose(t5,Vecf(0,0,1),288.0);
		for(Size i = 1; i <= t2.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
		for(Size i = 1; i <= t3.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t3.residue(i),1); else pose.append_residue_by_bond(t3.residue(i));
		for(Size i = 1; i <= t4.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t4.residue(i),1); else pose.append_residue_by_bond(t4.residue(i));
		for(Size i = 1; i <= t5.n_residue(); ++i) if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()) pose.append_residue_by_jump(t5.residue(i),1); else pose.append_residue_by_bond(t5.residue(i));
	}
	// TODO refactor this to NOT use poses
	void
	get_best_sub1_contact_delta_rotations(
		Real icmp1, Real dcmp1, Real icmp2, Real dcmp2,
		Real & out_delta_ang1,  Real & out_delta_ang2
	){
		Real const contact_dis = basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]();
		Real const contact_dis2 = contact_dis*contact_dis;
		// std::cerr << p1.residue(1).xyz("CB") << endl << cmp1cbs_[1] << endl;
		// std::cerr << p2.residue(1).xyz("CB") << endl << cmp2cbs_[1] << endl << endl;

		// pick rotations with best num contacts
		core::kinematics::Stub x1(rotation_matrix_degrees(cmp1axs_,icmp1),dcmp1*cmp1axs_);
		core::kinematics::Stub x2(rotation_matrix_degrees(cmp2axs_,icmp2),dcmp2*cmp2axs_);
		// xform_pose(p1,x1);
		// xform_pose(p2,x2);

				// std::cerr << p1.residue(1).xyz("CB") << endl << x1.local2global(cmp1cbs_[1]) << endl;
				// std::cerr << p2.residue(1).xyz("CB") << endl << x2.local2global(cmp2cbs_[1]) << endl << endl;

		// cerr << cmp1cbs_.size() << " " << p1.n_residue() << std::endl;
		// cerr << cmp2cbs_.size() << " " << p2.n_residue() << std::endl;

		int best1=0,best2=0,bestcount=0;
		for(int ir1 = 0; ir1 < cmp1nsub_; ++ir1){
			for(int ir2 = 0; ir2 < cmp2nsub_; ++ir2){
				int ccount = 0;
				for(Size ir=1; ir<=cmp1cbs_.size()/cmp1nsub_; ++ir){ Vec const cb1( x1.local2global(cmp1cbs_[ir]) );
				for(Size jr=1; jr<=cmp2cbs_.size()/cmp2nsub_; ++jr){ Vec const cb2( x2.local2global(cmp2cbs_[jr]) );
					if( cb1.distance_squared(cb2) < contact_dis2 ) ccount++;
					// if( p1.residue(ir).xyz(2).distance_squared(p2.residue(jr).xyz(2)) < dist2 ) ccount++;
				}}
				if(ccount > bestcount){
					bestcount = ccount;
					best1 = ir1;
					best2 = ir2;
				}

				// std::cerr << p1.residue(1).xyz("CB") << endl << x1.local2global(cmp1cbs_[1]) << endl;
				// std::cerr << p2.residue(1).xyz("CB") << endl << x2.local2global(cmp2cbs_[1]) << endl << endl;

				// rot_pose(p2,rotation_matrix_degrees(cmp2axs_,(Real)cmp2nangle_));
				x2.M = rotation_matrix_degrees(cmp2axs_,(Real)cmp2nangle_) * x2.M;
			}
			// rot_pose(p1,rotation_matrix_degrees(cmp1axs_,(Real)cmp1nangle_));
			x1.M = rotation_matrix_degrees(cmp1axs_,(Real)cmp1nangle_) * x1.M;

		}
		out_delta_ang1 = (Real)best1*(Real)cmp1nangle_;
		out_delta_ang2 = (Real)best2*(Real)cmp2nangle_;
	}

	void dump_pdb(int icmp2, int icmp1, int iori, string fname, int idx, bool dumpsym=true, bool sepcomp=true, bool termini=true) {
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using ObjexxFCL::string_of;
		Real d,dcmp2,dcmp1,icbc,cmp1cbc,cmp2cbc;
		dock_get_geom(icmp2,icmp1,iori,d,dcmp2,dcmp1,icbc,cmp2cbc,cmp1cbc);
		Pose p1,p2,symm;
		#ifdef USE_OPENMP
		#pragma omp critical
		#endif
		{
			p1.append_residue_by_jump(cmp1in_.residue(1),1,"","",false);
			for(Size i = 2; i <= cmp1in_.n_residue()/cmp1nsub_; ++i) {
				if(p1.residue(i-1).is_terminus()||p1.residue(i-1).is_ligand())
				     p1.append_residue_by_jump(cmp1in_.residue(i),1);
				else p1.append_residue_by_bond(cmp1in_.residue(i));
			}
			if(termini && !p1.residue(       1      ).is_lower_terminus()) add_lower_terminus_type_to_pose_residue(p1,      1       );
			if(termini && !p1.residue(p1.n_residue()).is_upper_terminus()) add_upper_terminus_type_to_pose_residue(p1,p1.n_residue());
			p2.append_residue_by_jump(cmp2in_.residue(1),1,"","", true ); // other chain iff not dumping sym complex (too many chains)
			for(Size i = 2; i <= cmp2in_.n_residue()/cmp2nsub_; ++i) {
				if(p2.residue(p2.n_residue()).is_terminus()||p2.residue(p2.n_residue()).is_ligand())
				     p2.append_residue_by_jump(cmp2in_.residue(i),1);
				else p2.append_residue_by_bond(cmp2in_.residue(i));
			}
			if(termini && !p2.residue(       1      ).is_lower_terminus()) add_lower_terminus_type_to_pose_residue(p2,      1       );
			if(termini && !p2.residue(p2.n_residue()).is_upper_terminus()) add_upper_terminus_type_to_pose_residue(p2,p2.n_residue());
		}

		core::kinematics::Stub x1(rotation_matrix_degrees(cmp1axs_,(Real)icmp1),dcmp1*cmp1axs_);
		core::kinematics::Stub x2(rotation_matrix_degrees(cmp2axs_,(Real)icmp2),dcmp2*cmp2axs_);
		xform_pose(p1,x1);
		xform_pose(p2,x2);

		// make joint structure
		{
			symm.append_residue_by_jump(p1.residue(1),1,"","",false);
			for(Size i = 2; i <= p1.n_residue(); ++i) {
				if(symm.residue(i-1).is_terminus()||symm.residue(i-1).is_ligand()) symm.append_residue_by_jump(p1.residue(i),1);
				else                                symm.append_residue_by_bond(p1.residue(i));
			}
			symm.append_residue_by_jump(p2.residue(1),1,"","", sepcomp ); // other chain iff not dumping sym complex (too many chains)
			for(Size i = 2; i <= p2.n_residue(); ++i) {
				if(symm.residue(symm.n_residue()).is_terminus()||symm.residue(symm.n_residue()).is_ligand()) symm.append_residue_by_jump(p2.residue(i),1);
				else                                             symm.append_residue_by_bond(p2.residue(i));
			}
		}
		// symm.dump_pdb("test_"+string(dumpsym?"presym":"nosym")+".pdb");
		if(dumpsym) {
			option[symmetry::symmetry_definition]("input/"+cmp1type_.substr(0,1)+"_tcdock.sym");
			core::pose::symmetry::make_symmetric_pose(symm);
			// symm.dump_pdb("test_sym.pdb");
		}
		core::io::pdb::dump_pdb(symm,option[out::file::o]()+"/"+fname);


		if(lnscore_){
			lnscore_->dump_linkers(toXform(x1),toXform(x2),option[out::file::o]()+"/"+fname);
		}

		// find contacts

		vector1<int> dumpg = option[tcdock::dump_pdb_grid]();
		if( std::find(dumpg.begin(),dumpg.end(),idx)!=dumpg.end() ){
			Pose asym;
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			{
				asym.append_residue_by_jump(p1.residue(1),1);
				for(Size i = 2; i <= p1.n_residue(); ++i) {
				  if(asym.residue(i-1).is_terminus()||asym.residue(i-1).is_ligand()) asym.append_residue_by_jump(p1.residue(i),1);
					else                                asym.append_residue_by_bond(p1.residue(i));
				}
				asym.append_residue_by_jump(p2.residue(1),1);
				for(Size i = 2; i <= p2.n_residue(); ++i) {
				  if(asym.residue(asym.n_residue()).is_terminus()||asym.residue(asym.n_residue()).is_ligand()) asym.append_residue_by_jump(p2.residue(i),1);
					else                                             asym.append_residue_by_bond(p2.residue(i));
				}
			}
			Real pr1=0.0,pr2=0.0;
			for(vector1<Vecf>::const_iterator i = cmp1pts_.begin(); i != cmp1pts_.end(); ++i) pr1 = max(pr1,projperp(cmp1axs_,*i).length());
			for(vector1<Vecf>::const_iterator i = cmp2pts_.begin(); i != cmp2pts_.end(); ++i) pr2 = max(pr2,projperp(cmp2axs_,*i).length());
			int const na = option[tcdock::grid_size_rot]();
			int const nr = option[tcdock::grid_size_disp]();
			Real const da1 = numeric::conversions::degrees(asin(option[tcdock::grid_delta_rot]()*1.2/pr1));
			Real const da2 = numeric::conversions::degrees(asin(option[tcdock::grid_delta_rot]()*1.2/pr1));
			Real const dr1 = -option[tcdock::grid_delta_disp]()*dcmp1/d;
			Real const dr2 = -option[tcdock::grid_delta_disp]()*dcmp2/d;
			rot_pose(asym,cmp1axs_,-(na+1)*da1,               1, p1.n_residue());
			rot_pose(asym,cmp2axs_, (na  )*da2,p1.n_residue()+1,asym.n_residue());
			for(int ia1 = 0; ia1 < 2*na+1; ++ia1) {
				rot_pose(asym,cmp1axs_,da1,1, p1.n_residue());
				rot_pose(asym,cmp2axs_,-(2*na+1)*da2,p1.n_residue()+1,asym.n_residue());
				for(int ia2 = 0; ia2 < 2*na+1; ++ia2) {
					rot_pose(asym,cmp2axs_,da2,p1.n_residue()+1,asym.n_residue());
					for(int ir = 0; ir < nr; ++ir) {
						#ifdef USE_OPENMP
						#pragma omp critical
						#endif
						core::io::pdb::dump_pdb(asym,option[out::file::o]()+"/"+fname+"_"+string_of(ia1)+"_"+string_of(ia2)+"_"+string_of(ir)+".pdb"+(option[tcdock::dump_gz]()?".gz":""));
						trans_pose(asym,cmp1axs_*dr1,1, p1.n_residue());
						trans_pose(asym,cmp2axs_*dr2,p1.n_residue()+1,asym.n_residue());
					}
					trans_pose(asym,cmp1axs_*-nr*dr1,1, p1.n_residue());
					trans_pose(asym,cmp2axs_*-nr*dr2,p1.n_residue()+1,asym.n_residue());
				}
			}
		}
	}
	Real __dock_base__(int icmp2,int icmp1,int iori,Real&dori,Real&dcmp2,Real&dcmp1,Real&icbc,Real&cmp2cbc,Real&cmp1cbc,bool cache=true){
		if(!cache || gradii(icmp2+1,icmp1+1,iori+1)==-9e9 ) {
			using basic::options::option;
			using namespace basic::options::OptionKeys;

			Vecf sicaxis = rotation_matrix_degrees( cmp1axs_.cross(cmp2axs_) ,(Real)iori) * (cmp1axs_+cmp2axs_).normalized();
			core::kinematics::Stub xa,xb;
			xa.M = rotation_matrix_degrees(cmp1axs_,(Real)icmp1);
			xb.M = rotation_matrix_degrees(cmp2axs_,(Real)icmp2);


			// Real const d = -slide_into_contact_and_score(sic_,*rigid_sfxn_,xa,xb,sicaxis,icbc);
			Real const d = -sic_.slide_into_contact_DEPRICATED(xa,xb,sicaxis);
			xa.v -= d*sicaxis;

			if(fabs(d) > 9e8){
				gscore(icmp2+1,icmp1+1,iori+1) = 0.0;
				gradii(icmp2+1,icmp1+1,iori+1) = 9e9;
				return 0.0;
			}
			if(d > 0){
				// TR << "WARNING slide in wrong direction!" << std::endl;
				gscore(icmp2+1,icmp1+1,iori+1) = 0.0;
				gradii(icmp2+1,icmp1+1,iori+1) = 9e9;
				return 0.0;
			}

			dori = d;
			Real const theta=(Real)iori;
			Real const gamma=numeric::conversions::radians(theta-alpha_/2.0);
			Real const sin_gamma=sin(gamma), cos_gamma=cos(gamma), x=d*sin_gamma, y=d*cos_gamma, w=x/sin_alpha_, z=x/tan_alpha_;
			dcmp2 = y+z;  dcmp1 = w;
						//
						// cerr << icmp1 << " " << icmp2 << " " << iori << " " << d << endl;
						// cerr << dcmp1 << " " << dcmp2 << endl;
						// utility_exit_with_message("FOO");
						//
			if( w > 0 ) {
				Real const cmp2mn = cmp2mnpos_[icmp2+1];
				Real const cmp1mn = cmp1mnpos_[icmp1+1];
				if( dcmp1 < cmp1mn || dcmp2 < cmp2mn ) {
					Real const dmncmp2 = cmp2mn/(cos_gamma+sin_gamma/tan_alpha_);
					Real const dmncmp1 = cmp1mn*sin_alpha_/sin_gamma;
					gradii(icmp2+1,icmp1+1,iori+1) = min(dmncmp2,dmncmp1);
					gscore(icmp2+1,icmp1+1,iori+1) = 0.0;
					return -12345;
				}
				int dp = (int)(dcmp2-cmp2mn)*10+1;
				int dt = (int)(dcmp1-cmp1mn)*10+1;
				if(dp < 1) utility_exit_with_message("bad");
				if(dt < 1) utility_exit_with_message("bad");
				if( 0 < dp && dp <= 200 ) cmp1cbc = cmp2cbpos_(icmp2+1,dp);
				if( 0 < dt && dt <= 200 ) cmp2cbc = cmp1cbpos_(icmp1+1,dt);
			} else {
				Real const cmp2mn = cmp2mnneg_[icmp2+1];
				Real const cmp1mn = cmp1mnneg_[icmp1+1];
				if( dcmp1 > cmp1mn || dcmp2 > cmp2mn ) {
					Real const dmncmp2 = cmp2mn/(cos_gamma+sin_gamma/tan_alpha_);
					Real const dmncmp1 = cmp1mn*sin_alpha_/sin_gamma;
					gradii(icmp2+1,icmp1+1,iori+1) = min(dmncmp2,dmncmp1);
					gscore(icmp2+1,icmp1+1,iori+1) = 0.0;
					return -12345;
				}
				int dp = (int)(-dcmp2+cmp2mn)*10+1;
				int dt = (int)(-dcmp1+cmp1mn)*10+1;
				if(dp < 1) utility_exit_with_message("bad");
				if(dt < 1) utility_exit_with_message("bad");
				if( 0 < dp && dp <= 200 ) cmp1cbc = cmp2cbneg_(icmp2+1,dp);
				if( 0 < dt && dt <= 200 ) cmp2cbc = cmp1cbneg_(icmp1+1,dt);
			}
			if(icbc!=-12345.0) {
				icbc = rigid_sfxn_->score(toXform(xa),toXform(xb));
				gscore(icmp2+1,icmp1+1,iori+1)  = icbc;
				gscore(icmp2+1,icmp1+1,iori+1) += option[tcdock::intra]()*option[tcdock::intra2]()*cmp2cbc;
				gscore(icmp2+1,icmp1+1,iori+1) += option[tcdock::intra]()*option[tcdock::intra1]()*cmp1cbc;
				gscore(icmp2+1,icmp1+1,iori+1) += termini_score(dcmp1,icmp1,dcmp2,icmp2);
			}
			gradii(icmp2+1,icmp1+1,iori+1) = d;
		}
		return gscore(icmp2+1,icmp1+1,iori+1);
	}
	void   dock_no_score(int icmp2,int icmp1,int iori){
		Real dori,dcmp2,dcmp1,icbc,cmp2cbc,cmp1cbc;
		icbc = -12345.0;
		__dock_base__(icmp2,icmp1,iori,dori,dcmp2,dcmp1,icbc,cmp2cbc,cmp1cbc,true);
	}
	Real dock_score(int icmp2,int icmp1,int iori,Real&icbc,Real&cmp2cbc,Real&cmp1cbc){
		Real dori,dcmp2,dcmp1;
		return __dock_base__(icmp2,icmp1,iori,dori,dcmp2,dcmp1,icbc,cmp2cbc,cmp1cbc,false);
	}
	Real dock_score(int icmp2,int icmp1,int iori){
		Real dori,dcmp2,dcmp1,icbc,cmp2cbc,cmp1cbc;
		return __dock_base__(icmp2,icmp1,iori,dori,dcmp2,dcmp1,icbc,cmp2cbc,cmp1cbc,false);
	}
	Real dock_get_geom(int icmp2,int icmp1,int iori,Real&dori,Real&dcmp2,Real&dcmp1,Real&icbc,Real&cmp2cbc,Real&cmp1cbc){
		return __dock_base__(icmp2,icmp1,iori,dori,dcmp2,dcmp1,icbc,cmp2cbc,cmp1cbc,false);
	}
	Real min_termini_dis(Real d1, Real a1, Real d2, Real a2){
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		bool usen1 = option[tcdock::usen1]();
		bool usec1 = option[tcdock::usec1]();
		Real td = 9e9;
		for(int t1 = 0; t1 <= option[tcdock::termini_trim](); ++t1) {
			Vecf n1_0 = cmp1in_.xyz(core::id::AtomID(1,1+t1));
			for(int t2 = 0; t2 <= option[tcdock::termini_trim](); ++t2) {
				Vecf n2_0 = cmp2in_.xyz(core::id::AtomID(1,1+t2));
				for(int t3 = 0; t3 <= option[tcdock::termini_trim](); ++t3) {
					Vecf c1_0 = cmp1in_.xyz(core::id::AtomID(3,cmp1in_.n_residue()-t3));
					for(int t4 = 0; t4 <= option[tcdock::termini_trim](); ++t4) {
						Vecf c2_0 = cmp2in_.xyz(core::id::AtomID(3,cmp2in_.n_residue()-t4));
						for(int i1 = 0; i1 < cmp1nsub_; ++i1){
							Vecf n1 = d1*cmp1axs_ + rotation_matrix_degrees(cmp1axs_,a1 + 360.0/(Real)cmp1nsub_*(Real)i1) * n1_0;
							Vecf c1 = d1*cmp1axs_ + rotation_matrix_degrees(cmp1axs_,a1 + 360.0/(Real)cmp1nsub_*(Real)i1) * c1_0;
							for(int i2 = 0; i2 < cmp2nsub_; ++i2){
								Vecf n2 = d2*cmp2axs_ + rotation_matrix_degrees(cmp2axs_,a2 + 360.0/(Real)cmp2nsub_*(Real)i2) * n2_0;
								Vecf c2 = d2*cmp2axs_ + rotation_matrix_degrees(cmp2axs_,a2 + 360.0/(Real)cmp2nsub_*(Real)i2) * c2_0;
								if(usen1) td = min( n1.distance(c2), td );
								if(usec1) td = min( c1.distance(n2), td );
							}
						}
					}
				}
			}
		}
		return td;
	}
	Real min_termini_proj(Real a1, Real a2, Vec sicaxis){
		using basic::options::option;
		using namespace basic::options::OptionKeys;
//		return 0.0;
		bool usen1 = option[tcdock::usen1]();
		bool usec1 = option[tcdock::usec1]();
		Real td = 9e9;
		int t1=1;
		Vecf n1_0 = cmp1in_.xyz(core::id::AtomID(1,1+t1));
		int t2=1;
		Vecf n2_0 = cmp2in_.xyz(core::id::AtomID(1,1+t2));
		int t3=0;
		Vecf c1_0 = cmp1in_.xyz(core::id::AtomID(3,cmp1in_.n_residue()-t3));
		int t4=0;
		Vecf c2_0 = cmp2in_.xyz(core::id::AtomID(3,cmp2in_.n_residue()-t4));
		for(int i1 = 0; i1 < cmp1nsub_; ++i1){
			Vecf n1 = rotation_matrix_degrees(cmp1axs_,a1 + 360.0/(Real)cmp1nsub_*(Real)i1) * n1_0;
			Vecf c1 = rotation_matrix_degrees(cmp1axs_,a1 + 360.0/(Real)cmp1nsub_*(Real)i1) * c1_0;
			for(int i2 = 0; i2 < cmp2nsub_; ++i2){
				Vecf n2 = rotation_matrix_degrees(cmp2axs_,a2 + 360.0/(Real)cmp2nsub_*(Real)i2) * n2_0;
				Vecf c2 = rotation_matrix_degrees(cmp2axs_,a2 + 360.0/(Real)cmp2nsub_*(Real)i2) * c2_0;
				if(usen1) td = min( projperp(sicaxis,n1-c2).length(), td );
				if(usec1) td = min( projperp(sicaxis,c1-n2).length(), td );
			}
		}
		//std::cerr << td << std::endl;
		return td;
	}
	Real termini_score(Real d1, Real a1, Real d2, Real a2){
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		Real td = min_termini_dis(d1,a1,d2,a2);
		return option[tcdock::termini_weight]() * max(0.0,option[tcdock::termini_cutoff]()-max(option[tcdock::termini_cutoff_short](),td));
	}
	void justone(int icmp1, int icmp2, int iori) {
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		int ilm=0;
		LMAX const h(0,0,icmp1,icmp2,iori);
		Real d,dcmp2,dcmp1,icbc,cmp1cbc,cmp2cbc;
		int N = option[tcdock::peak_grid_size]();
		ObjexxFCL::FArray3D<Real> grid(2*N+1,2*N+1,2*N+1,0.0);
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for(int di = -N; di <= N; ++di) {
			for(int dj = -N; dj <= N; ++dj) {
				for(int dk = -N; dk <= N; ++dk) {
					if( Vecf(di,dj,dk).length() > (Real)N+0.5 ) {
						grid(di+N+1,dj+N+1,dk+N+1) = -9e9;
					} else {
						int i = (h.icmp2+di+gscore.size1())%gscore.size1();
						int j = (h.icmp1+dj+gscore.size2())%gscore.size2();
						int k = (h.iori +dk+gscore.size3())%gscore.size3();
						dock_no_score(i,j,k);
						grid(di+N+1,dj+N+1,dk+N+1) = gradii(i+1,j+1,k+1);
					}
				}
			}
		}
		vector1<Real> ffhist(25,0);
		// #ifdef USE_OPENMP
		// #pragma omp parallel for schedule(dynamic,1)
		// #endif
		for(int ifh = 0; ifh < (int)ffhist.size(); ++ifh) {
			protocols::sic_dock::flood_fill3D(N+1,N+1,N+1,grid, grid(N+1,N+1,N+1)-0.000001 - 0.2);
			Real count = 0;
			int Nedge = option[tcdock::peak_grid_smooth]();
			ObjexxFCL::FArray3D<Real> grid2(grid);
			for(int i = 1+Nedge; i <= (int)grid.size1()-Nedge; ++i) {
				for(int j = 1+Nedge; j <= (int)grid.size2()-Nedge; ++j) {
					for(int k = 1+Nedge; k <= (int)grid.size3()-Nedge; ++k) {
						if( grid(i,j,k)!=grid(N+1,N+1,N+1) ) continue;
						int ninside = 0;
						for(int di = -Nedge; di <= Nedge; ++di) {
							for(int dj = -Nedge; dj <= Nedge; ++dj) {
								for(int dk = -Nedge; dk <= Nedge; ++dk) {
									ninside += grid(i+di,j+dj,k+dk)==grid(N+1,N+1,N+1);
								}
							}
						}
						Real w = max(0.0, 1.0 - Vecf(N+1-i,N+1-j,N+1-k).length() / (Real)N );
						count += w*((Real)ninside/(Real)((2*Nedge+1)*(2*Nedge+1)*(2*Nedge+1)));
						// cerr << F(7,3,((Real)ninside/(Real)((2*Nedge+1)*(2*Nedge+1)*(2*Nedge+1)))) << " " << F(5,3,w) << " " << I(2,ninside) << endl;
						// grid2(i,j,k) += allgood;
					}
				}
			}
			ffhist[ifh+1] = pow(count,1.0/3.0);
			vector1<int> dumpg = option[tcdock::dump_peak_grids]();
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			if( std::find(dumpg.begin(),dumpg.end(),ilm)!=dumpg.end() ) {
				utility::io::ozstream o(("out/grid_"+ObjexxFCL::string_of(ilm)+"_"+ObjexxFCL::string_of(ifh)+".dat.gz"));
				for(int i = 1; i <= (int)grid.size1(); ++i) {
					for(int j = 1; j <= (int)grid.size2(); ++j) {
						for(int k = 1; k <= (int)grid.size3(); ++k) {
							o << grid2(i,j,k) << endl;
						}
					}
				}
				o.close();
				dump_pdb(h.icmp2,h.icmp1,h.iori,"test_"+symtype_+"_"+ObjexxFCL::string_of(ilm)+".pdb",0,true);
			}
		}
		Real score = dock_get_geom(h.icmp2,h.icmp1,h.iori,d,dcmp2,dcmp1,icbc,cmp2cbc,cmp1cbc);
		Real da1=0,da2=0;
		get_best_sub1_contact_delta_rotations(h.icmp1,dcmp1,h.icmp2,dcmp2,da1,da2);

		string fn = cmp1name_+"_"+cmp2name_+"_"+symtype_+(option[tcdock::reverse]()?"R":"F")+"_"+ObjexxFCL::string_of(ilm);
		cout << "| " << fn << ((ilm<10)?"  ":" ")
             << F(8,3,score) << " "
		     << F(6,2,2*max(fabs(dcmp1)+(dcmp1>0?cmp1diapos_:cmp1dianeg_),fabs(dcmp2)+(dcmp2>0?cmp2diapos_:cmp2dianeg_))) << " "
		     << F(6,2, min(999.99,min_termini_dis(dcmp1,h.icmp1,dcmp2,h.icmp2) )) << " "
             << F(7,3,icbc) << " "
		     << F(6,2,cmp1cbc) << " "
		     << F(6,2,cmp2cbc) << " "
		     << I(4,cmp1in_.n_residue()) << " "
		     << I(3,h.icmp1+(int)da1) << " "
		     << F(8,3,dcmp1) << " "
		     << I(4,cmp2in_.n_residue()) << " "
		     << I(3,h.icmp2+(int)da2) << " "
		     << F(8,3,dcmp2) << " "
			 << I(3,h.iori);
		for(int ifh = 1; ifh <= (int)ffhist.size(); ++ifh) {
			cout << " " << F(5,2,ffhist[ifh]);
		}
		cout << endl;
		if(option[tcdock::dump_pdb]() ) dump_pdb(h.icmp2+int(da2),h.icmp1+(int)da1,h.iori,fn+".pdb"+(option[tcdock::dump_gz]()?".gz":""),0,true);

	}
	void run() {
		if(abort_) {
			cout << "abort run" << endl;
			return;
		}
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using namespace core::id;
		using numeric::conversions::radians;
		// Real const CONTACT_D  = basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]();
		// Real const CLASH_D    = basic::options::option[basic::options::OptionKeys::sicdock::	clash_dis]();
		// Real const CONTACT_D2 = sqr(CONTACT_D);
		// Real const CLASH_D2   = sqr(CLASH_D);
		Pose const cmp1init(cmp1in_);
		Pose const cmp2init(cmp2in_);
		{                                                         // 3deg loop
			// Real max1=0;
			// dump_onecomp();
			precompute_intra();

			// start_time_ = time_highres();
			// conf_count_ = 0;

			cout << "main loop 1 over icmp2, icmp1, iori every 3 degrees" << endl;
			Real max_score = 0;
			int max_i1=0,max_i2=0,max_io=0;
			if(option[tcdock::fast_stage_one]()){
				for(int icmp1 = 0; icmp1 < cmp1nangle_; icmp1+=3) {
					if(icmp1%15==0 && icmp1!=0){
						// Real rate = Real(conf_count_) / Real(time_highres()-start_time_);
						cout<<" lowres dock "
						    <<cmp1name_<<" "
						    <<cmp2name_<<" "
						    << (option[tcdock::reverse]()?"R ":"F ")
						    <<I(2,100*icmp1/cmp1nangle_)
						    <<"% done, max_score: "
						    <<F(10,6,max_score)
						    // <<" rate: "<<rate<<" confs/sec"
						    // <<" "<<rate/Real(num_threads())<<" confs/sec/thread"
						    <<endl;
					}
					#ifdef USE_OPENMP
					#pragma omp parallel for schedule(dynamic,1)
					#endif
					for(int icmp2 = 0; icmp2 < cmp2nangle_; icmp2+=3) {
						vector1<Vecf> pb,cbb; get_cmp2(icmp2,pb,cbb);
						int iori = -1, stg = 1;	bool newstage = true;
						while(stg < 5) {
							if(newstage) {
								iori = (int)(((stg>2)?270.0:90.0)+angle_degrees(cmp1axs_,Vecf(0,0,0),cmp2axs_)/2.0);
								iori = (iori / 3) * 3; // round to closest multiple of angle incr
								if(stg==2||stg==4) iori -= 3;
								newstage = false;
							} else {
								iori += (stg%2==0) ? -3 : 3;
							}
							Real const score = dock_score(icmp2,icmp1,iori);
							if( -12345==score ) { stg++; newstage=true; continue; }
							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							{
								if(score > max_score) max_score = score;
							}
						}
					}
				}
			} else {
				vector1<numeric::xyzVector<int> > tasks;
				for(int icmp1 = 0; icmp1 < cmp1nangle_; icmp1+=3) {
					for(int icmp2 = 0; icmp2 < cmp2nangle_; icmp2+=3) {
						for(int iori = 0; iori < 360; iori+=3) {
							tasks.push_back( numeric::xyzVector<int>(icmp1,icmp2,iori) );
						}
					}
				}
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int i = 1; i <= (int)tasks.size(); ++i){
					int icmp1 = tasks[i].x();
					int icmp2 = tasks[i].y();
					int iori  = tasks[i].z();
					if(i%28800==0){
						// Real rate = Real(conf_count_) / Real(time_highres()-start_time_);
						cout<<" lowres dock "
						    <<cmp1name_<<" "
						    <<I(2,100*icmp1/cmp1nangle_)
						    <<"% done, max_score: "
						    <<F(10,6,max_score)
						    // <<" rate: "<<rate<<" confs/sec"
						    // <<" "<<rate/num_threads()<<" confs/sec/thread"<<" " <<num_threads()
						    <<endl;
					}

					Real const score = dock_score(icmp2,icmp1,iori);
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					{
						++conf_count_;
						if(score > max_score){
							max_score = score;
							// max_i1 = icmp1;
							// max_i2 = icmp2;
							// max_io = iori;
						}
					}
				}
			}
			if(max_score<0.00001) utility_exit_with_message("cmp1 or cmp2 too large, no contacts!");
			cout << "MAX3 " << max_score << " " << max_i1 << " " << max_i2 << " " << max_io << endl;
		}
		utility::vector1<vector1<int> > cmp2lmx,cmp1lmx,orilmx; { // set up work for main loop 2
			Real topX_3 = 0;
			vector1<Real> cbtmp;
			for(Size i = 0; i < gscore.size(); ++i) if(gscore[i] > 0) cbtmp.push_back(gscore[i]);
			std::sort(cbtmp.begin(),cbtmp.end());
			topX_3 = cbtmp[max(1,(int)cbtmp.size()-option[tcdock::nsamp1]()+1)];
			cout << "scanning top "<<min(option[tcdock::nsamp1](),(int)cbtmp.size())<<" with score3 >= " << topX_3 << endl;
			for(int icmp2 = 0; icmp2 < cmp2nangle_; icmp2+=3) {
				for(int icmp1 = 0; icmp1 < cmp1nangle_; icmp1+=3) {
					for(int iori = 0; iori < 360; iori+=3) {
						if( gscore(icmp2+1,icmp1+1,iori+1) >= topX_3) {
							vector1<int> cmp2,cmp1,ori;
							for(int i = -1; i <= 1; ++i) cmp2.push_back( (icmp2+i+ cmp2nangle_)% cmp2nangle_ );
							for(int j = -1; j <= 1; ++j) cmp1.push_back( (icmp1+j+cmp1nangle_)%cmp1nangle_ );
							for(int k = -1; k <= 1; ++k) ori.push_back( (iori+k+360)%360 );
							cmp2lmx.push_back(cmp2);
							cmp1lmx.push_back(cmp1);
							orilmx.push_back(ori);
						}
					}
				}
			}
		}
		{                                                         //main loop 2
			#ifdef USE_OPENMP
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for(int ilmx = 1; ilmx <= (int)cmp2lmx.size(); ++ilmx)  {       //  MAIN LOOP 2
				if( (ilmx-1)%(option[tcdock::nsamp1]()/10)==0 && ilmx!=1){
					cout<<" highres dock "<<cmp1name_<<" "<<(Real(ilmx-1)/(option[tcdock::nsamp1]()/100))<<"% done"<<endl;
				}
				for(vector1<int>::const_iterator picmp2 = cmp2lmx[ilmx].begin(); picmp2 != cmp2lmx[ilmx].end(); ++picmp2) {
					int icmp2 = *picmp2;
					for(vector1<int>::const_iterator picmp1 = cmp1lmx[ilmx].begin(); picmp1 != cmp1lmx[ilmx].end(); ++picmp1) {
						int icmp1 = *picmp1;
						for(vector1<int>::const_iterator piori = orilmx[ilmx].begin(); piori != orilmx[ilmx].end(); ++piori) {
							int iori = *piori;
							dock_score(icmp2,icmp1,iori);
						}
					}
				}
			}
		}
		cout << "gather local minima" << endl;
		int npad = (cmp1name_+"_"+cmp2name_).size()+6;
		vector1<LMAX> local_maxima;	{                             // get local radial disp maxima (minima, really)
			Real highscore = -9e9;
			for(int i = 1; i <=(int) gradii.size1(); ++i){
				for(int j = 1; j <= (int)gradii.size2(); ++j){
					for(int k = 1; k <= (int)gradii.size3(); ++k){
						Real const val = gradii(i,j,k);
						if( val < -9e8 ) continue;

						// Real nbmax = -9e9, scoremax = -9e9;
						// if( option[tcdock::geometric_minima_only]() ){
						// 	// int nedge = 0;
						// 	for(int di = -1; di <= 1; ++di){
						// 		for(int dj = -1; dj <= 1; ++dj){
						// 			for(int dk = -1; dk <= 1; ++dk){
						// 				if(di==0 && dj==0 && dk==0) continue;
						// 				int i2 = (i+di+gradii.size1()-1)%gradii.size1()+1;
						// 				int j2 = (j+dj+gradii.size2()-1)%gradii.size2()+1;
						// 				int k2 = (k+dk+gradii.size3()-1)%gradii.size3()+1;
						// 				Real const nbval = gradii(i2,j2,k2);
						// 				nbmax = max(nbmax,nbval);
						// 			}
						// 		}
						// 	}
						// }
						// if( option[tcdock::score_minima_only]() ){
						// 	// int nedge = 0;
						// 	for(int di = -1; di <= 1; ++di){
						// 		for(int dj = -1; dj <= 1; ++dj){
						// 			for(int dk = -1; dk <= 1; ++dk){
						// 				if(di==0 && dj==0 && dk==0) continue;
						// 				int i2 = (i+di+gradii.size1()-1)%gradii.size1()+1;
						// 				int j2 = (j+dj+gradii.size2()-1)%gradii.size2()+1;
						// 				int k2 = (k+dk+gradii.size3()-1)%gradii.size3()+1;
						// 				Real const nbval = gradii(i2,j2,k2);
						// 				nbmax = max(nbmax,nbval);
						// 			}
						// 		}
						// 	}
						// }

						// if( /*nbmax != -9e9 &&*/ val >= nbmax ) {
							Real score = gscore(i,j,k); //dock_score(i-1,j-1,k-1);
							if(score <= 0) continue;
							local_maxima.push_back( LMAX(score,gradii(i,j,k),i-1,j-1,k-1) );
							highscore = max(score,highscore);
						// }
					}
				}
			}
			std::sort(local_maxima.begin(),local_maxima.end(),compareLMAX);
			cout << "N maxima: " << local_maxima.size() << ", best score: " << highscore << endl;
			string nc1=cmp1type_.substr(1,1), nc2=cmp2type_.substr(1,1);
			string pad = ""; for(int i = 1; i <= npad; ++i) pad += " ";
			cout << "  tag" << pad << "    score   diam   tdis   inter    ";
			cout << "sc"+nc1+"    sc"+nc2+"  nr"+nc1+"  a"+nc1+"       r"+nc1+"  nr"+nc2+"  a"+nc2+"       r"+nc2+" ori";
			cout << "  v0.2  v0.4  v0.6  v0.8  v1.0  v1.2  v1.4  v1.6  v1.8  v2.0";
			cout << "  v2.2  v2.4  v2.6  v2.8  v3.0  v3.2  v3.4  v3.6  v3.8  v4.0  v4.2  v4.4  v4.6  v4.8  v5.0" << endl;
		}
		vector1<Vec> dumpedit;
		int ilm = 0, nout = 0;
		Real redundant_thresh2 = option[tcdock::redundant_angle]() * option[tcdock::redundant_angle]();
		while( ++ilm <= (int)local_maxima.size() && nout < option[tcdock::topx]() ){
		// for(Size ilm = 1; ilm <= min(local_maxima.size(),(Size)option[tcdock::topx]()); ++ilm) { // dump top hit info
			LMAX const & h(local_maxima[ilm]);

			// redundency check
			bool toosimilar = false;
			for(vector1<Vec>::const_iterator i = dumpedit.begin(); i != dumpedit.end(); ++i){
				Real x = cmp1scale_*(Real)h.icmp1;
				Real y = cmp2scale_*(Real)h.icmp2;
				Real z = (cmp1scale_+cmp2scale_)/2.0*h.iori;
				if( i->distance_squared( Vec(x,y,z) ) < redundant_thresh2 ){
					toosimilar = true;
				}
			}
			if(toosimilar) continue;

			Real d,dcmp2,dcmp1,icbc,cmp1cbc,cmp2cbc;
			int N = option[tcdock::peak_grid_size]();
			ObjexxFCL::FArray3D<Real> grid(2*N+1,2*N+1,2*N+1,0.0);
			#ifdef USE_OPENMP
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for(int di = -N; di <= N; ++di) {
				for(int dj = -N; dj <= N; ++dj) {
					for(int dk = -N; dk <= N; ++dk) {
						if( Vecf(di,dj,dk).length() > (Real)N+0.5 ) {
							grid(di+N+1,dj+N+1,dk+N+1) = -9e9;
						} else {
							int i = (h.icmp2+di+gscore.size1())%gscore.size1();
							int j = (h.icmp1+dj+gscore.size2())%gscore.size2();
							int k = (h.iori +dk+gscore.size3())%gscore.size3();
							dock_no_score(i,j,k);
							grid(di+N+1,dj+N+1,dk+N+1) = gradii(i+1,j+1,k+1);
						}
					}
				}
			}
			vector1<Real> ffhist(25,0);
			int const Nedge = option[tcdock::peak_grid_smooth]();
			// #ifdef USE_OPENMP
			// #pragma omp parallel for schedule(dynamic,1)
			// #endif
			for(int ifh = 0; ifh < (int)ffhist.size(); ++ifh) {
				protocols::sic_dock::flood_fill3D(N+1,N+1,N+1, grid, grid(N+1,N+1,N+1)-0.000001 - 0.2);
				Real count = 0;
				// ObjexxFCL::FArray3D<Real> grid2(grid);
				for(int i = 1+Nedge; i <= (int)grid.size1()-Nedge; ++i) {
					for(int j = 1+Nedge; j <= (int)grid.size2()-Nedge; ++j) {
						for(int k = 1+Nedge; k <= (int)grid.size3()-Nedge; ++k) {
							if( grid(i,j,k)!=grid(N+1,N+1,N+1) ) continue;
							int ninside = 0;
							for(int di = -Nedge; di <= Nedge; ++di) {
								for(int dj = -Nedge; dj <= Nedge; ++dj) {
									for(int dk = -Nedge; dk <= Nedge; ++dk) {
										ninside += grid(i+di,j+dj,k+dk)==grid(N+1,N+1,N+1);
									}
								}
							}
							Real w = max(0.0, 1.0 - Vecf(N+1-i,N+1-j,N+1-k).length() / (Real)N );
							count += w*((Real)ninside/(Real)((2*Nedge+1)*(2*Nedge+1)*(2*Nedge+1)));
							// cerr << F(7,3,((Real)ninside/(Real)((2*Nedge+1)*(2*Nedge+1)*(2*Nedge+1)))) << " " << F(5,3,w) << " " << I(2,ninside) << endl;
							// grid2(i,j,k) += allgood;
						}
					}
				}
				ffhist[ifh+1] = pow(count,1.0/3.0);
				// vector1<int> dumpg = option[tcdock::dump_peak_grids]();
				// #ifdef USE_OPENMP
				// #pragma omp critical
				// #endif
				// if( std::find(dumpg.begin(),dumpg.end(),ilm)!=dumpg.end() ) {
				// 	utility::io::ozstream o(("out/grid_"+ObjexxFCL::string_of(ilm)+"_"+ObjexxFCL::string_of(ifh)+".dat.gz"));
				// 	for(int i = 1; i <= (int)grid.size1(); ++i) {
				// 		for(int j = 1; j <= (int)grid.size2(); ++j) {
				// 			for(int k = 1; k <= (int)grid.size3(); ++k) {
				// 				o << grid2(i,j,k) << endl;
				// 			}
				// 		}
				// 	}
				// 	o.close();
				// }
			}
			Real score = dock_get_geom(h.icmp2,h.icmp1,h.iori,d,dcmp2,dcmp1,icbc,cmp2cbc,cmp1cbc);
			Real da1=0,da2=0;
			get_best_sub1_contact_delta_rotations(h.icmp1,dcmp1,h.icmp2,dcmp2,da1,da2);

			string fn = cmp1name_+"_"+cmp2name_+"_"+symtype_+(option[tcdock::reverse]()?"R":"F")+"_"+ObjexxFCL::string_of(ilm);
			cout << "| " << ObjexxFCL::format::LJ(npad+3,fn) << " "
                 << F(8,3,score) << " "
			     << F(6,2,2*max(fabs(dcmp1)+(dcmp1>0?cmp1diapos_:cmp1dianeg_),fabs(dcmp2)+(dcmp2>0?cmp2diapos_:cmp2dianeg_))) << " "
			     << F(6,2, min_termini_dis(dcmp1,h.icmp1,dcmp2,h.icmp2) ) << " "
                 << F(7,3,icbc) << " "
			     << F(6,2,cmp1cbc) << " "
			     << F(6,2,cmp2cbc) << " "
			     << I(4,cmp1in_.n_residue()) << " "
			     << I(3,h.icmp1+(int)da1) << " "
			     << F(8,3,dcmp1) << " "
			     << I(4,cmp2in_.n_residue()) << " "
			     << I(3,h.icmp2+(int)da2) << " "
			     << F(8,3,dcmp2) << " "
				 << I(3,h.iori);
			if(lnscore_){
				core::kinematics::Stub x1(rotation_matrix_degrees(cmp1axs_,(Real)h.icmp1),dcmp1*cmp1axs_);
				core::kinematics::Stub x2(rotation_matrix_degrees(cmp2axs_,(Real)h.icmp2),dcmp2*cmp2axs_);
				cout << "   " << F(6,2,lnscore_->score(toXform(x1),toXform(x2))) << "  ";
			}
			for(int ifh = 1; ifh <= (int)ffhist.size(); ++ifh) {
				cout << " " << F(5,2,ffhist[ifh]);
			}
			cout << endl;
			vector1<int> dumpg = option[tcdock::dump_pdb_grid]();
			if(option[tcdock::dump_pdb_primary_sub]() ) {
				dump_pdb(h.icmp2+(int)da2,h.icmp1+(int)da1,h.iori,fn+"_sub1.pdb"+(option[tcdock::dump_gz]()?".gz":""),ilm, false, true );
			}
			if(option[tcdock::dump_pdb]() || std::find(dumpg.begin(),dumpg.end(),ilm)!=dumpg.end() ) {
				dump_pdb(h.icmp2+(int)da2,h.icmp1+(int)da1,h.iori,fn+"_full.pdb"+(option[tcdock::dump_gz]()?".gz":""),ilm, true, option[tcdock::separate_components]() );
			}
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			{
				nout++;
				Real x = cmp1scale_*(Real)h.icmp1;
				Real y = cmp2scale_*(Real)h.icmp2;
				Real z = (cmp1scale_+cmp2scale_)/2.0*h.iori;
				// this is kinda dumb... to handle periodicity
				dumpedit.push_back(Vec(x            ,y            ,z));
				dumpedit.push_back(Vec(x+cmp2nangle_,y            ,z));
				dumpedit.push_back(Vec(x            ,y+cmp2nangle_,z));
				dumpedit.push_back(Vec(x+cmp2nangle_,y+cmp2nangle_,z));
				dumpedit.push_back(Vec(x-cmp2nangle_,y            ,z));
				dumpedit.push_back(Vec(x            ,y-cmp2nangle_,z));
				dumpedit.push_back(Vec(x-cmp2nangle_,y-cmp2nangle_,z));
				dumpedit.push_back(Vec(x+cmp2nangle_,y-cmp2nangle_,z));
				dumpedit.push_back(Vec(x-cmp2nangle_,y+cmp2nangle_,z));
			}

		}
		cout << "DONE symdock_enum_3_1 " << cmp1name_+" "+cmp2name_+" "+symtype_+(option[tcdock::reverse]()?" R":" F") << endl;
	}

};
int main (int argc, char *argv[]) {

	try {

	register_options();
	devel::init(argc,argv);
	using basic::options::option;
	using namespace basic::options::OptionKeys;

// 	for(Size i = 1; i <= option[tcdock::I5]().size(); ++i) {
// 		Pose pose1;
// 		std::cout << option[tcdock::I5]()[i];
// 		core::import_pose::pose_from_file(pose1,option[tcdock::I5]()[i], core::import_pose::PDB_file);
// 		std::cout << " " << pose1.n_residue() << " DONE" << std::endl;
// 		continue;
// 	}
// 	return 0;

	vector1<string> compkind;
	vector1<vector1<string> > compfiles;
	vector1<Size> compnfold;
	if( option[tcdock::I5].user() ) { compnfold.push_back(5); compkind.push_back("I5"); compfiles.push_back(option[tcdock::I5]()); }
	if( option[tcdock::I3].user() ) { compnfold.push_back(3); compkind.push_back("I3"); compfiles.push_back(option[tcdock::I3]()); }
	if( option[tcdock::I2].user() ) { compnfold.push_back(2); compkind.push_back("I2"); compfiles.push_back(option[tcdock::I2]()); }
	if( option[tcdock::O4].user() ) { compnfold.push_back(4); compkind.push_back("O4"); compfiles.push_back(option[tcdock::O4]()); }
	if( option[tcdock::O3].user() ) { compnfold.push_back(3); compkind.push_back("O3"); compfiles.push_back(option[tcdock::O3]()); }
	if( option[tcdock::O2].user() ) { compnfold.push_back(2); compkind.push_back("O2"); compfiles.push_back(option[tcdock::O2]()); }
	if( option[tcdock::T3].user() ) { compnfold.push_back(3); compkind.push_back("T3"); compfiles.push_back(option[tcdock::T3]()); }
	if( option[tcdock::T2].user() ) { compnfold.push_back(2); compkind.push_back("T2"); compfiles.push_back(option[tcdock::T2]()); }
	if(compkind.size()!=2) utility_exit_with_message("must specify two components!");
	if(compkind[1][0]!=compkind[2][0]) utility_exit_with_message("components must be of same sym icos, tetra or octa");
	Real ttrim_cut = option[tcdock::trim_floppy_termini]();
	for(Size i = 1; i <= compfiles[1].size(); ++i) {
		Pose pose1;
		core::import_pose::pose_from_file(pose1,compfiles[1][i], core::import_pose::PDB_file);
		if(ttrim_cut > 0.001){
			TR << "triming tails on " << compfiles[1][i] << std::endl;
			protocols::sic_dock::auto_trim_floppy_termini(pose1,ttrim_cut,compnfold[1]);
		}
		if( pose1.n_residue() > (Size)option[tcdock::max_res]() ){
			TR << "skip " << compfiles[1][i] << " " << pose1.n_residue() << std::endl;
			continue;
		}
		for(Size j = 1; j <= compfiles[2].size(); ++j) {
			Pose pose2;
			core::import_pose::pose_from_file(pose2,compfiles[2][j], core::import_pose::PDB_file);
			if(ttrim_cut > 0.001){
				TR << "triming tails on " << compfiles[2][j] << std::endl;
				protocols::sic_dock::auto_trim_floppy_termini(pose2,ttrim_cut,compnfold[2]);
			}
			if( pose2.n_residue() > (Size)option[tcdock::max_res]() ){
				TR << "skip " << compfiles[2][j] << " " << pose2.n_residue() << std::endl;
				continue;
			}
			TCDock tcd(pose1,pose2,compfiles[1][i],compfiles[2][j],compkind[1],compkind[2]);
			if(option[tcdock::justone].user()) {
				int icm1 = option[tcdock::justone]()[1];
				int icm2 = option[tcdock::justone]()[2];
				int iori = option[tcdock::justone]()[3];
				tcd.justone(icm2,icm1,iori);
			} else {
				try{
					cout << "RUN " << i << " " << j << endl;
					tcd.run();
				} catch(int e) {
					continue;
				}
			}
		}
	}
	cout << "DONE symdock_enum_3_1" << endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


