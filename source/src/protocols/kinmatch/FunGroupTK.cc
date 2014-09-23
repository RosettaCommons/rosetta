// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <protocols/kinmatch/FunGroupTK.hh>

#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/AtomTree.hh>



static thread_local basic::Tracer TR( "protocols.kinmatchFunGroupTK" );

namespace protocols {
namespace kinmatch {

/// @details Auto-generated virtual destructor
FunGroupTK::~FunGroupTK() {}

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using core::Real;
using core::Size;
using core::pose::Pose;
using core::pose::PoseOP;
using core::pose::PoseCOP;
using core::id::AtomID;
using core::kinematics::Stub;
using core::conformation::Residue;
using core::conformation::ResidueOP;
using protocols::scoring::ImplicitFastClashCheck;
using std::string;
using utility::vector1;
using numeric::min;
using numeric::min;
using numeric::max;
using numeric::constants::d::pi;
using numeric::conversions::degrees;
using utility::io::ozstream;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;

typedef utility::vector1<core::Real> Reals;
typedef utility::vector1<core::Size> Sizes;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef utility::vector1<Vec> Vecs;

template< typename T >
inline
std::string
str( T const & t )
{
	return ObjexxFCL::string_of(t);
}

template< typename T >
inline
std::string
lzs(
	T const & t,
	int const w // Minimum width
)
{
	return ObjexxFCL::lead_zero_string_of(t,w);
}




		inline Real sqr(Real const r) { return r*r; }

		inline Vec xyz(Pose const & p, Size const & ia, Size const & ir) {
    return p.xyz(AtomID(ia,ir));
		}

		PoseOP alapose(Pose const & pose_in) {
				PoseOP rpose( new Pose(pose_in) );
    Pose & pose(*rpose);
    for(Size i=1; i<=pose.n_residue(); ++i) {
						core::pose::replace_pose_residue_copying_existing_coordinates(pose,i,pose.residue(i).residue_type_set().name_map("ALA"));
    }
    return rpose;
		}

		vector1<Size> allifnone(vector1<Size> v, Size n) {
    if( v.size()==0 ) {
						v.resize(n);
						for(Size i = 1; i <=n; ++i) v[i] = i;
    }
    return(v);
		}


		void dumpcgo(Vec v, string l) {
    std::cerr << "cmd.load_cgo(Vec( " << v.x() << "," << v.y() << "," << v.z() << ").cgo(), '"+l+"')"  << std::endl;
		}

		void xform_rsd_gl2(Stub const & s, Residue & rsd) {
    for(Size i = 1; i <= rsd.natoms(); ++i) rsd.set_xyz(i, s.global2local(rsd.xyz(i)) );
		}


KRSQuery::KRSQuery(KRSQueryType typ                                                   ) : type(typ),cen(0.0,0.0,0.0),axs(1,0,0),ori(0.0,1.0,0.0),disth(0.5),angth(0.0175),clash( 2.8) {}
KRSQuery::KRSQuery(KRSQueryType typ, Vec c, Vec a, Vec o, Real dt, Real at, Real clsh ) : type(typ),cen(  c  ),axs(  a  ),ori(  o  ),disth( dt),angth(  at  ),clash(clsh) {}
KRSQuery::KRSQuery(KRSQueryType typ, Vec c, Vec a,        Real dt, Real at, Real clsh ) : type(typ),cen(  c  ),axs(  a  ),ori(0.0,0.0,0.0),disth( dt),angth(  at  ),clash(clsh) {}


FunGroupTK::FunGroupTK(
				Pose & p_in,
				vector1<Size> & pos
				)
				: pose_(alapose(p_in)), pos_(allifnone(pos,pose_->n_residue()))
{
		ifc_ = protocols::scoring::ImplicitFastClashCheckCOP( new ImplicitFastClashCheck(*pose_,2.2) );
		frs_ = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueTypeSetCOP frs( frs_ );
		vector1<string> res_types;
		res_types.push_back("ASP");
		res_types.push_back("CYS");
		res_types.push_back("HIS");
		for(vector1<string>::const_iterator it = res_types.begin(); it != res_types.end(); ++it) {
				rsd_[*it].resize(pose_->n_residue());
				for(vector1<Size>::const_iterator i=pos_.begin(); i!=pos_.end(); ++i) {
						ResidueOP rsd = core::conformation::ResidueFactory::create_residue(frs->name_map(*it),pose_->residue(*i),pose_->conformation());
						stb_[*it].push_back(Stub(rsd->xyz(5),rsd->xyz(2),rsd->xyz(1)));
						rsd_[*it][*i] = rsd;
				}
		}
}

ImplicitFastClashCheck const &
FunGroupTK::ifc() const {
    return *ifc_;
}

BruteFunGroupTK::BruteFunGroupTK( Pose & p_in, vector1<Size> pos ) : FunGroupTK(p_in,pos) {}

void
BruteFunGroupTK::place_c(
    KRSQuery const & q,
    Residue const & qrsd,
    vector1<core::conformation::ResidueOP> & hits
    ) const
{
		vector1<Size>::const_iterator ip = pos_.begin();
		vector1<Stub> const & stb(stb_.find("CYS")->second);
		for(vector1<Stub>::const_iterator is = stb.begin(); is!=stb.end(); ++is,++ip) {
				Real cbd2 = q.cen.distance_squared(is->v);
				if(*ip==qrsd.seqpos()) continue;
				if( cbd2 > 14.0 ) continue;
				ResidueOP bestrsd;
				Real bestsc = 9e9;
				ResidueOP rsd = rsd_.find("CYS")->second[*ip];
				rsd->set_xyz(11, rsd->xyz("SG") + 2.1*(rsd->xyz("HG")-rsd->xyz("SG")).normalized() );
				for(Real ch1 = 0; ch1 <= 360; ch1+=3.0) {
						rsd->set_chi(1,ch1);
						Vec sg = rsd->xyz("SG");
						if( sg.distance_squared(q.cen) > 8.0 ) continue;
						for(Real ch2 = 0; ch2 <= 360; ch2+=3.0) {
								rsd->set_chi(2,ch2);
								Vec hg = rsd->xyz("HG");
								Vec cen = sg+2.1*((hg-sg).normalized());
								Real const dsq = cen.distance_squared(q.cen);
								if( dsq > q.disth*q.disth) continue;
								Vec axs = ((sg-q.cen).normalized());
								Real const da = numeric::angle_radians( axs, Vec(0,0,0), q.axs);

								// Vec qc = is->global2local(q.cen);
								// qc = qc * 2.911171 / qc.length();
								// ozstream out("test_"+str(qc.x())+"_"+lzs(qrsd.seqpos(),3)+"_"+lzs(*ip,3)+"_"+lzs((Size)ch1,3)+"_"+lzs((Size)ch2,3)+".pdb");
								// Size ano=0;
								// core::io::pdb::dump_pdb_residue(*rsd,ano,out);
								// out.close();

								if( fabs(da-q.ori.x()) > q.angth ) continue;
								bool clash = false;
								for(Size ia = 6; ia <= rsd->nheavyatoms(); ++ia) {
										if(!ifc().clash_check(rsd->xyz(ia),*ip) ) { clash = true; break; }
										for(Size ja = 1; ja <= qrsd.nheavyatoms(); ++ja) if( qrsd.xyz(ja).distance_squared(rsd->xyz(ia)) < q.clash ) { clash = true; break; }
										if(clash) break;
								}
								if(clash) continue;
								if( da < bestsc ) {
										//TR << da << " " << q.ori.x() << " " << q.angth << std::endl;
										bestrsd = rsd->clone();
										bestsc = da;
								}
						}
				}
				if(bestsc < 9e8) {
						hits.push_back(bestrsd);
				}
		}
}

void
BruteFunGroupTK::place_h(
    KRSQuery const & q,
    Residue const & qrsd,
    vector1<core::conformation::ResidueOP> & hits
    ) const
{
		vector1<Size>::const_iterator ip = pos_.begin();
		vector1<Stub> const & stb(stb_.find("HIS")->second);
		for(vector1<Stub>::const_iterator is = stb.begin(); is!=stb.end(); ++is,++ip) {
				Real cbd2 = q.cen.distance_squared(is->v);
				if(*ip==qrsd.seqpos()) continue;
				if( cbd2 > 14.0 ) continue;
				ResidueOP bestrsd;
				Real bestsc = 9e9;
				ResidueOP rsd = rsd_.find("HIS")->second[*ip];
				rsd->set_xyz(11, rsd->xyz("SG") + 2.1*(rsd->xyz("HG")-rsd->xyz("SG")).normalized() );
				for(Real ch1 = 0; ch1 <= 360; ch1+=3.0) {
						rsd->set_chi(1,ch1);
						Vec sg = rsd->xyz("SG");
						if( sg.distance_squared(q.cen) > 8.0 ) continue;
						for(Real ch2 = 0; ch2 <= 360; ch2+=3.0) {
								rsd->set_chi(2,ch2);
								Vec hg = rsd->xyz("HG");
								Vec cen = sg+2.1*((hg-sg).normalized());
								Real const dsq = cen.distance_squared(q.cen);
								if( dsq > q.disth*q.disth) continue;
								Vec axs = ((sg-q.cen).normalized());
								Real const da = numeric::angle_radians( axs, Vec(0,0,0), q.axs);

								// Vec qc = is->global2local(q.cen);
								// qc = qc * 2.911171 / qc.length();
								// ozstream out("test_"+str(qc.x())+"_"+lzs(qrsd.seqpos(),3)+"_"+lzs(*ip,3)+"_"+lzs((Size)ch1,3)+"_"+lzs((Size)ch2,3)+".pdb");
								// Size ano=0;
								// core::io::pdb::dump_pdb_residue(*rsd,ano,out);
								// out.close();

								if( fabs(da-q.ori.x()) > q.angth ) continue;
								bool clash = false;
								for(Size ia = 6; ia <= rsd->nheavyatoms(); ++ia) {
										if(!ifc().clash_check(rsd->xyz(ia),*ip) ) { clash = true; break; }
										for(Size ja = 1; ja <= qrsd.nheavyatoms(); ++ja) if( qrsd.xyz(ja).distance_squared(rsd->xyz(ia)) < q.clash ) { clash = true; break; }
										if(clash) break;
								}
								if(clash) continue;
								if( da < bestsc ) {
										//TR << da << " " << q.ori.x() << " " << q.angth << std::endl;
										bestrsd = rsd->clone();
										bestsc = da;
								}
						}
				}
				if(bestsc < 9e8) {
						hits.push_back(bestrsd);
				}
		}
}

void
BruteFunGroupTK::place_d(
    KRSQuery const & q,
    Residue const & qrsd,
    vector1<core::conformation::ResidueOP> & hits
    ) const
{
	vector1<Size>::const_iterator ip = pos_.begin();
	vector1<Stub> const & stb(stb_.find("ASP")->second);
	for(vector1<Stub>::const_iterator is = stb.begin(); is!=stb.end(); ++is,++ip) {
		Real cbd2 = q.cen.distance_squared(is->v);
		if(*ip==qrsd.seqpos()) continue;
		if( cbd2 > 36.0 ) continue;
		ResidueOP bestrsd;
		Real bestsc = 9e9;
		ResidueOP rsd = rsd_.find("ASP")->second[*ip];
		Vec cb = rsd->xyz("CB");
		for(Real ch1 = 0; ch1 <= 360; ch1+=3.0) {
			rsd->set_chi(1,ch1);
			Vec cg = rsd->xyz("CG");
			if( cg.distance_squared(q.cen) > 25.0 ) continue;
			for(Real ch2 = 0; ch2 <= 360; ch2+=3.0) {
				rsd->set_chi(2,ch2);
				Vec od1 = rsd->xyz("OD1");
				Vec od2 = rsd->xyz("OD2");
				for(int od12 = 0; od12 <= 1; ++od12) {
				//int od12 = 0;
					for(int prpp = 0; prpp <= 1; ++prpp) {
						Vec  const oda(od12?od1:od2);
						Vec  const odb(od12?od2:od1);
						Vec  const axs( ((prpp?(cg-cb):(oda-odb)).normalized()) );
						Vec  const cen = oda + 1.8 * axs;
						Real const dsq = cen.distance_squared(q.cen);
						if( dsq > q.disth*q.disth) continue;
						Real const dt( axs.dot(q.axs) );
						if( dt < q.angth ) continue; // dot, not ang
						bool clash = false;
						for(Size ia = 6; ia <= rsd->nheavyatoms()-2; ++ia) {
							if(!ifc().clash_check(rsd->xyz(ia),*ip) ) { clash = true; break; }
							for(Size ja = 1; ja <= qrsd.nheavyatoms(); ++ja) if( qrsd.xyz(ja).distance_squared(rsd->xyz(ia)) < q.clash ) { clash = true; break; }
							if(clash) break;
						}
						if(clash) continue;
						Real const da = dsq - 5*dt;
						if( da < bestsc ) {
							//TR << da << " " << q.ori.x() << " " << q.angth << std::endl;
							bestrsd = rsd->clone();
							bestsc = da;
						}
					}
				}
			}
		}
		if(bestsc < 9e8) {
			hits.push_back(bestrsd);
		}
	}
}


KinFunGroupTK::KinFunGroupTK(
    Pose & p_in,
    vector1<Size> pos )
    : FunGroupTK(p_in,pos) {}

void
KinFunGroupTK::place_c(
    KRSQuery const & q,
    Residue const & qrsd,
    vector1<core::conformation::ResidueOP> & hits_out
    ) const
{
#define CYS_CB_HG_DIS      2.911171 // 3.388379 // 2.771698
#define CYS_SG_CB_H        0.7494478 // height of SG "above" CB
#define CYS_HG_SG_PRJLEN   2.088496  //1.818653 // sin(CB-SG-HG)*len(HG-SG)
#define CYS_1oSIN_CB_SG_HG 1.092026   //
#define CYS_HG_SG_X_DROP   0.0998    // h drop due measured, 2.1*tan((90-84.011803)/180*pi)*tan((90-65.522644)/180*pi)
#define CYS_CB_SG_PERP     1.646237
		if( fabs(q.axs.length()-1.0) > 0.0000001 ) utility_exit_with_message("place_c query axs not kormalized()!");
		Real const r3o2 = sqrt(3.0)/2.0;
		Real const dis2ub = sqr(        CYS_CB_HG_DIS+q.disth );
		Real const dis2lb = sqr(max(0.0,CYS_CB_HG_DIS-q.disth));
		vector1<Size>::const_iterator ip = pos_.begin();
		vector1<Stub> const & stb(stb_.find("CYS")->second);
		vector1<ResidueOP> const & rsdlst(rsd_.find("CYS")->second);
		for(vector1<Stub>::const_iterator is = stb.begin(); is!=stb.end(); ++is,++ip) {
				if(*ip==qrsd.seqpos()) continue;
				core::conformation::ResidueOP rtmp = rsdlst[*ip];;
				Real const cbd2 = q.cen.distance_squared(is->v);
				if( dis2ub < cbd2 || cbd2 < dis2lb ) continue;
				Vec  const qcen0(is->global2local(q.cen));
				Real const pdis = sqrt(sqr(q.disth)-sqr(sqrt(cbd2)-CYS_CB_HG_DIS)) * qcen0.length() / CYS_CB_HG_DIS;
				Size const NPOS = (pdis > 0.2) ? 7 : 1;
				Vec  const Y = Vec(0,1,0).cross(qcen0).normalized();
				Vec  const Z =          Y.cross(qcen0).normalized();
				vector1<core::conformation::ResidueOP> hits;
				for(int flp = -1; flp <= 1; flp+=2) {
						Real best_angerr = 9e9, best_ch1=0.0, best_ch2=0.0, mn_ang=9e9, mx_ang=-9e9;
						for(Size ipos = 0; ipos < NPOS; ++ipos) {
								Real y=0,z=0;
								if     (ipos==1) { y =  1.0; z =   0.0; }
								else if(ipos==2) { y =  0.5; z =  r3o2; }
								else if(ipos==3) { y = -0.5; z =  r3o2; }
								else if(ipos==4) { y = -1.0; z =   0.0; }
								else if(ipos==5) { y = -0.5; z = -r3o2; }
								else if(ipos==6) { y =  0.5; z = -r3o2; }
								Vec  const qcen = qcen0 + y*Y*pdis + z*Z*pdis;
								// calc chi2
								Real const lqcen = qcen.length();
								Real const qcenxsphere = qcen.x() * CYS_CB_HG_DIS / lqcen;
								Real const h = ( qcenxsphere - CYS_SG_CB_H) * CYS_1oSIN_CB_SG_HG - CYS_HG_SG_X_DROP;
								if( CYS_HG_SG_PRJLEN < fabs(h) ) continue;
								Real const ch2_0 = acos( h/CYS_HG_SG_PRJLEN ); // chi2 = pi +- ch2_0
								//calc chi1
								Real const lqcenpx = sqrt(qcen.y()*qcen.y()+qcen.z()*qcen.z());
								Real const ch1_0 = (qcen.y()>0)? asin(qcen.z()/lqcenpx) : pi-asin(qcen.z()/lqcenpx);
								Real const hgz = sin(pi-ch2_0)*CYS_HG_SG_PRJLEN;
								Real const hgdzy = lqcenpx * CYS_CB_HG_DIS / lqcen;
								Real const ch1_ofst = asin( hgz/hgdzy );
								Real const ch1 = ch1_0 + flp*ch1_ofst;
								Real const ch2 = pi    + flp*ch2_0   ;
								// check angle
								Vec  const sg  = is->local2global(Vec(CYS_SG_CB_H, cos(ch1)*CYS_CB_SG_PERP, sin(ch1)*CYS_CB_SG_PERP));
								if(!ifc().clash_check(sg,*ip) ) continue;

								Real const bang = acos( (sg-q.cen).normalized().dot(q.axs) );
								// // make rsd
								// core::conformation::ResidueOP rtmp = core::conformation::ResidueFactory::create_residue(frs_->name_map("CYS"),pose_.residue(*ip),pose_.conformation());
								// rtmp->set_xyz(11, rtmp->xyz(6) + 2.1*(rtmp->xyz(11)-rtmp->xyz(6)).normalized() );
								// rtmp->set_chi(1,degrees(ch1));
								// rtmp->set_chi(2,degrees(ch2));
								// //        ResidueOP hrsd = qrsd.clone();
								// //        xform_rsd_gl2(*is,*rtmp);
								// //        xform_rsd_gl2(*is,*hrsd);
								// //hrsd->set_xyz("HE2", hrsd->xyz("HE2") * CYS_CB_HG_DIS / hrsd->xyz("HE2").length() );
								// Size tmp = 1;
								// ozstream out1("cys_"+lzs(*ip,3)+"_"+str(ipos)+"_"+str(flp+1)+".pdb");
								// core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
								// //core::io::pdb::dump_pdb_residue(*hrsd,tmp,out1);
								// out1.close();
								// //        utility_exit_with_message(str(degrees(bang)));
								mn_ang = min(bang,mn_ang);
								mx_ang = max(bang,mx_ang);
								Real const angerr = fabs( bang - q.ori.x() );
								if( angerr < best_angerr ) {
										best_angerr = angerr;
										best_ch1 = ch1;
										best_ch2 = ch2;
								}
						}
						if( q.ori.x() < mn_ang-q.angth || mx_ang+q.angth < q.ori.x() ) continue;

						// position rsd
						rtmp->set_chi(1,degrees(best_ch1));
						rtmp->set_chi(2,degrees(best_ch2));

						// if( hits.size() == 1 ) {
						//   Real tmp = fabs( acos((hits[1]->xyz("SG")-q.cen).normalized().dot(q.axs)) - q.ori.x() );
						//   if( tmp < best_angerr ) hits[1] = rtmp;
						// } else {
						hits.push_back(rtmp->clone());
						//        }
						//if(hits.size()) TR << *ip << " " << hits.size() << std::endl;

						// ResidueOP hrsd = qrsd.clone();
						// xform_rsd_gl2(*is,*rtmp);
						// xform_rsd_gl2(*is,*hrsd);
						// hrsd->set_xyz("HE2", hrsd->xyz("HE2") * CYS_CB_HG_DIS / hrsd->xyz("HE2").length() );
						// Size tmp = 1;
						// ozstream out1("cys"+str(ipos)+"_"+str(flp+1)+".pdb");
						// core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
						// core::io::pdb::dump_pdb_residue(*hrsd,tmp,out1);
						// out1.close();

				}
				hits_out.insert(hits_out.end(),hits.begin(),hits.end());
		}
}

void
KinFunGroupTK::place_h(
    KRSQuery const & q,
    Residue const & qrsd,
    vector1<core::conformation::ResidueOP> & hits_out
    ) const
{
#define HIS_CB_HG_DIS       5.6570235
#define HIS_CB_HG_PAR       3.946279
#define HIS_CG_CB_H         0.6019999
#define HIS_HE2_CG_PRJLEN   1.5392682574
#define HIS_1oSIN_CA_CB_CG  1.092026358  // 1/sina(ca-cb-cg)
#define HIS_HE2_CG_X_DROP   (2.4821480185 - 0.087) // 5.6570235*math.tan((90-66.309441 )*math.pi/180)
#define HIS_CB_CG_PERP      1.371000

		if( fabs(q.axs.length()-1.0) > 0.0000001 ) utility_exit_with_message("place_h query axs not kormalized()!");
		Real const r3o2 = sqrt(3.0)/2.0;
		Real const dis2ub = sqr(        HIS_CB_HG_DIS+q.disth );
		Real const dis2lb = sqr(max(0.0,HIS_CB_HG_DIS-q.disth));
		vector1<Size>::const_iterator ip = pos_.begin();
		vector1<Stub> const & stb(stb_.find("HIS")->second);
		vector1<ResidueOP> const & rsdlst(rsd_.find("HIS")->second);
		for(vector1<Stub>::const_iterator is = stb.begin(); is!=stb.end(); ++is,++ip) {
				if(*ip==qrsd.seqpos()) continue;
				core::conformation::ResidueOP rtmp = rsdlst[*ip];

				// xform_rsd_gl2(*is,*rtmp);
				// rtmp->set_xyz("HE2", rtmp->xyz("NE2") + 2.1*(rtmp->xyz("HE2")-rtmp->xyz("NE2")).normalized() );
				// Size tmp = 1;
				// ozstream out1("his_"+ObjexxFCL::string_of(*ip,3)+".pdb");
				// core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
				// out1.close();
				// utility_exit_with_message("his test");

				Real const cbd2 = q.cen.distance_squared(is->v);
				if( dis2ub < cbd2 || cbd2 < dis2lb ) continue;
				Vec  const qcen0(is->global2local(q.cen));
				Real const pdis = sqrt(sqr(q.disth)-sqr(sqrt(cbd2)-HIS_CB_HG_DIS)) * qcen0.length() / HIS_CB_HG_DIS;
				Size const NPOS = (pdis > 0.2) ? 7 : 1;
				Vec  const Y = Vec(0,1,0).cross(qcen0).normalized();
				Vec  const Z =          Y.cross(qcen0).normalized();
				vector1<core::conformation::ResidueOP> hits;
				for(int flp = -1; flp <= 1; flp+=2) {
						Real best_angerr = 9e9, best_ch1=0.0, best_ch2=0.0, mn_ang=9e9, mx_ang=-9e9;
						for(Size ipos = 0; ipos < NPOS; ++ipos) {
								Real y=0,z=0;
								if     (ipos==1) { y =  1.0; z =   0.0; }
								else if(ipos==2) { y =  0.5; z =  r3o2; }
								else if(ipos==3) { y = -0.5; z =  r3o2; }
								else if(ipos==4) { y = -1.0; z =   0.0; }
								else if(ipos==5) { y = -0.5; z = -r3o2; }
								else if(ipos==6) { y =  0.5; z = -r3o2; }
								Vec  const qcen = qcen0 + y*Y*pdis + z*Z*pdis;
								// calc chi2
								Real const lqcen = qcen.length();
								Real const qcenxsphere = qcen.x() * HIS_CB_HG_DIS / lqcen;
								Real const h = ( qcenxsphere ) * HIS_1oSIN_CA_CB_CG - HIS_HE2_CG_X_DROP;
								if( h < -HIS_HE2_CG_PRJLEN || HIS_HE2_CG_PRJLEN < h ) continue;
								Real const ch2_0 = (h>0) ? acos( -h/HIS_HE2_CG_PRJLEN ) : 2*pi-acos( -h/HIS_HE2_CG_PRJLEN ); // chi2 = pi +- ch2_0
								//calc chi1
								Real const lqcenpx = sqrt(qcen.y()*qcen.y()+qcen.z()*qcen.z());
								Real const ch1_0 = (qcen.y()>0)? asin(qcen.z()/lqcenpx) : pi-asin(qcen.z()/lqcenpx);
								Real const hgz = sin(pi-ch2_0)*HIS_HE2_CG_PRJLEN;
								Real const hgdzy = lqcenpx * HIS_CB_HG_DIS / lqcen;
								Real const ch1_ofst = asin( hgz/hgdzy );
								Real const ch1 = ch1_0 - flp*ch1_ofst;
								Real const ch2 = pi    + flp*ch2_0   ;
								// check angle
								Vec const ne2 = is->local2global(rotation_matrix(Vec(1,0,0),ch1)*rotation_matrix(Vec(0.602000, 1.371000, 0.000000),ch2) * Vec(1.884000, 3.152000, -0.001000));
								//Vec  const ne2  = is->local2global(Vec(HIS_CG_CB_H, cos(ch1)*HIS_CB_CG_PERP, sin(ch1)*HIS_CB_CG_PERP));
								if(!ifc().clash_check(ne2,*ip) ) continue;

								Real const bang = acos( (ne2-q.cen).normalized().dot(q.axs) );

								// // make rsd
								//  core::conformation::ResidueOP rtmp = core::conformation::ResidueFactory::create_residue(frs_->name_map("HIS"),pose_->residue(*ip),pose_->conformation());
								//  rtmp->set_xyz("HE2", rtmp->xyz("NE2") + 2.1*(rtmp->xyz("HE2")-rtmp->xyz("NE2")).normalized() );
								//  rtmp->set_chi(1,degrees(ch1));
								//  rtmp->set_chi(2,degrees(ch2));

								// // dumpcgo(ne2,"rast");
								//  TR << (rtmp->xyz("NE2") - ne2 ).length() << " " << degrees(bang) << " " << angle_degrees(rtmp->xyz("NE2"),q.cen,q.cen+q.axs) << std::endl;

								// ResidueOP hrsd = qrsd.clone();
								// xform_rsd_gl2(*is,*rtmp);
								// xform_rsd_gl2(*is,*hrsd);
								// hrsd->set_xyz("HE2", hrsd->xyz("HE2") * HIS_CB_HG_DIS / hrsd->xyz("HE2").length() );
								// Size tmp = 1;
								// ozstream out1("his_"+lzs(qrsd.seqpos(),3)+"_"+lzs(q.ori.z(),3)+"_"+lzs(*ip,3)+"_"+str(ipos)+"_"+str(flp+1)+".pdb");
								// core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
								// core::io::pdb::dump_pdb_residue(*hrsd,tmp,out1);
								// out1.close();
								// TR << h << std::endl;
								// TR << lqcenpx << std::endl;
								// TR << qcen.z() << std::endl;
								// TR << degrees(ch1_ofst) << std::endl;
								// TR << degrees(ch1_0) << std::endl;
								// TR << degrees(ch2_0) << std::endl;
								// TR << degrees(ch2_0) << std::endl;
								// utility_exit_with_message("aoirstne");
								mn_ang = min(bang,mn_ang);
								mx_ang = max(bang,mx_ang);
								Real const angerr = fabs( bang - q.ori.x() );
								if( angerr < best_angerr ) {
										best_angerr = angerr;
										best_ch1 = ch1;
										best_ch2 = ch2;
								}
						}
						if( q.ori.x() < mn_ang-q.angth || mx_ang+q.angth < q.ori.x() ) continue;

						// position rsd
						rtmp->set_chi(1,degrees(best_ch1));
						rtmp->set_chi(2,degrees(best_ch2));

						// if( hits.size() == 1 ) {
						//   Real tmp = fabs( acos((hits[1]->xyz("NE2")-q.cen).normalized().dot(q.axs)) - q.ori.x() );
						//   if( tmp < best_angerr ) hits[1] = rtmp;
						// } else {
						hits.push_back(rtmp->clone());
						//        }
						//if(hits.size()) TR << *ip << " " << hits.size() << std::endl;

						// ResidueOP hrsd = qrsd.clone();
						// xform_rsd_gl2(*is,*rtmp);
						// xform_rsd_gl2(*is,*hrsd);
						// hrsd->set_xyz("HE2", hrsd->xyz("HE2") * HIS_CB_HG_DIS / hrsd->xyz("HE2").length() );
						// Size tmp = 1;
						// ozstream out1("his"+str(ipos)+"_"+str(flp+1)+".pdb");
						// core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
						// core::io::pdb::dump_pdb_residue(*hrsd,tmp,out1);
						// out1.close();

				}
				hits_out.insert(hits_out.end(),hits.begin(),hits.end());
				//				if(hits.size()) utility_exit_with_message("aoirstne");
		}
}


void
KinFunGroupTK::place_d(
    KRSQuery const & q,
    Residue const & qrsd,
    vector1<core::conformation::ResidueOP> & hits_out
    ) const
{
#define HIS_CB_HG_DIS       5.6570235
#define HIS_CB_HG_PAR       3.946279
#define HIS_CG_CB_H         0.6019999
#define HIS_HE2_CG_PRJLEN   1.5392682574
#define HIS_1oSIN_CA_CB_CG  1.092026358  // 1/sina(ca-cb-cg)
#define HIS_HE2_CG_X_DROP   (2.4821480185 - 0.087) // 5.6570235*math.tan((90-66.309441 )*math.pi/180)
#define HIS_CB_CG_PERP      1.371000

		if( fabs(q.axs.length()-1.0) > 0.0000001 ) utility_exit_with_message("place_h query axs not kormalized()!");
		Real const r3o2 = sqrt(3.0)/2.0;
		Real const dis2ub = sqr(        HIS_CB_HG_DIS+q.disth );
		Real const dis2lb = sqr(max(0.0,HIS_CB_HG_DIS-q.disth));
		vector1<Size>::const_iterator ip = pos_.begin();
		vector1<Stub> const & stb(stb_.find("HIS")->second);
		vector1<ResidueOP> const & rsdlst(rsd_.find("HIS")->second);
		for(vector1<Stub>::const_iterator is = stb.begin(); is!=stb.end(); ++is,++ip) {
				if(*ip==qrsd.seqpos()) continue;
				core::conformation::ResidueOP rtmp = rsdlst[*ip];

				// xform_rsd_gl2(*is,*rtmp);
				// rtmp->set_xyz("HE2", rtmp->xyz("NE2") + 2.1*(rtmp->xyz("HE2")-rtmp->xyz("NE2")).normalized() );
				// Size tmp = 1;
				// ozstream out1("his_"+ObjexxFCL::string_of(*ip,3)+".pdb");
				// core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
				// out1.close();
				// utility_exit_with_message("his test");

				Real const cbd2 = q.cen.distance_squared(is->v);
				if( dis2ub < cbd2 || cbd2 < dis2lb ) continue;
				Vec  const qcen0(is->global2local(q.cen));
				Real const pdis = sqrt(sqr(q.disth)-sqr(sqrt(cbd2)-HIS_CB_HG_DIS)) * qcen0.length() / HIS_CB_HG_DIS;
				Size const NPOS = (pdis > 0.2) ? 7 : 1;
				Vec  const Y = Vec(0,1,0).cross(qcen0).normalized();
				Vec  const Z =          Y.cross(qcen0).normalized();
				vector1<core::conformation::ResidueOP> hits;
				for(int flp = -1; flp <= 1; flp+=2) {
						Real best_angerr = 9e9, best_ch1=0.0, best_ch2=0.0, mn_ang=9e9, mx_ang=-9e9;
						for(Size ipos = 0; ipos < NPOS; ++ipos) {
								Real y=0,z=0;
								if     (ipos==1) { y =  1.0; z =   0.0; }
								else if(ipos==2) { y =  0.5; z =  r3o2; }
								else if(ipos==3) { y = -0.5; z =  r3o2; }
								else if(ipos==4) { y = -1.0; z =   0.0; }
								else if(ipos==5) { y = -0.5; z = -r3o2; }
								else if(ipos==6) { y =  0.5; z = -r3o2; }
								Vec  const qcen = qcen0 + y*Y*pdis + z*Z*pdis;
								// calc chi2
								Real const lqcen = qcen.length();
								Real const qcenxsphere = qcen.x() * HIS_CB_HG_DIS / lqcen;
								Real const h = ( qcenxsphere ) * HIS_1oSIN_CA_CB_CG - HIS_HE2_CG_X_DROP;
								if( h < -HIS_HE2_CG_PRJLEN || HIS_HE2_CG_PRJLEN < h ) continue;
								Real const ch2_0 = (h>0) ? acos( -h/HIS_HE2_CG_PRJLEN ) : 2*pi-acos( -h/HIS_HE2_CG_PRJLEN ); // chi2 = pi +- ch2_0
								//calc chi1
								Real const lqcenpx = sqrt(qcen.y()*qcen.y()+qcen.z()*qcen.z());
								Real const ch1_0 = (qcen.y()>0)? asin(qcen.z()/lqcenpx) : pi-asin(qcen.z()/lqcenpx);
								Real const hgz = sin(pi-ch2_0)*HIS_HE2_CG_PRJLEN;
								Real const hgdzy = lqcenpx * HIS_CB_HG_DIS / lqcen;
								Real const ch1_ofst = asin( hgz/hgdzy );
								Real const ch1 = ch1_0 - flp*ch1_ofst;
								Real const ch2 = pi    + flp*ch2_0   ;
								// check angle
								Vec const ne2 = is->local2global(rotation_matrix(Vec(1,0,0),ch1)*rotation_matrix(Vec(0.602000, 1.371000, 0.000000),ch2) * Vec(1.884000, 3.152000, -0.001000));
								//Vec  const ne2  = is->local2global(Vec(HIS_CG_CB_H, cos(ch1)*HIS_CB_CG_PERP, sin(ch1)*HIS_CB_CG_PERP));
								if(!ifc().clash_check(ne2,*ip) ) continue;

								Real const bang = acos( (ne2-q.cen).normalized().dot(q.axs) );

								// // make rsd
								//  core::conformation::ResidueOP rtmp = core::conformation::ResidueFactory::create_residue(frs_->name_map("HIS"),pose_->residue(*ip),pose_->conformation());
								//  rtmp->set_xyz("HE2", rtmp->xyz("NE2") + 2.1*(rtmp->xyz("HE2")-rtmp->xyz("NE2")).normalized() );
								//  rtmp->set_chi(1,degrees(ch1));
								//  rtmp->set_chi(2,degrees(ch2));

								// // dumpcgo(ne2,"rast");
								//  TR << (rtmp->xyz("NE2") - ne2 ).length() << " " << degrees(bang) << " " << angle_degrees(rtmp->xyz("NE2"),q.cen,q.cen+q.axs) << std::endl;

								// ResidueOP hrsd = qrsd.clone();
								// xform_rsd_gl2(*is,*rtmp);
								// xform_rsd_gl2(*is,*hrsd);
								// hrsd->set_xyz("HE2", hrsd->xyz("HE2") * HIS_CB_HG_DIS / hrsd->xyz("HE2").length() );
								// Size tmp = 1;
								// ozstream out1("his_"+lzs(qrsd.seqpos(),3)+"_"+lzs(q.ori.z(),3)+"_"+lzs(*ip,3)+"_"+str(ipos)+"_"+str(flp+1)+".pdb");
								// core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
								// core::io::pdb::dump_pdb_residue(*hrsd,tmp,out1);
								// out1.close();
								// TR << h << std::endl;
								// TR << lqcenpx << std::endl;
								// TR << qcen.z() << std::endl;
								// TR << degrees(ch1_ofst) << std::endl;
								// TR << degrees(ch1_0) << std::endl;
								// TR << degrees(ch2_0) << std::endl;
								// TR << degrees(ch2_0) << std::endl;
								// utility_exit_with_message("aoirstne");
								mn_ang = min(bang,mn_ang);
								mx_ang = max(bang,mx_ang);
								Real const angerr = fabs( bang - q.ori.x() );
								if( angerr < best_angerr ) {
										best_angerr = angerr;
										best_ch1 = ch1;
										best_ch2 = ch2;
								}
						}
						if( q.ori.x() < mn_ang-q.angth || mx_ang+q.angth < q.ori.x() ) continue;

						// position rsd
						rtmp->set_chi(1,degrees(best_ch1));
						rtmp->set_chi(2,degrees(best_ch2));

						// if( hits.size() == 1 ) {
						//   Real tmp = fabs( acos((hits[1]->xyz("NE2")-q.cen).normalized().dot(q.axs)) - q.ori.x() );
						//   if( tmp < best_angerr ) hits[1] = rtmp;
						// } else {
						hits.push_back(rtmp->clone());
						//        }
						//if(hits.size()) TR << *ip << " " << hits.size() << std::endl;

						// ResidueOP hrsd = qrsd.clone();
						// xform_rsd_gl2(*is,*rtmp);
						// xform_rsd_gl2(*is,*hrsd);
						// hrsd->set_xyz("HE2", hrsd->xyz("HE2") * HIS_CB_HG_DIS / hrsd->xyz("HE2").length() );
						// Size tmp = 1;
						// ozstream out1("his"+str(ipos)+"_"+str(flp+1)+".pdb");
						// core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
						// core::io::pdb::dump_pdb_residue(*hrsd,tmp,out1);
						// out1.close();

				}
				hits_out.insert(hits_out.end(),hits.begin(),hits.end());
				//				if(hits.size()) utility_exit_with_message("aoirstne");
		}
}


} // namespace kinmatch
} // namespace protocols


