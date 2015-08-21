// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/sic_dock/designability_score.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

namespace protocols {
namespace sic_dock {

using platform::Real;
using numeric::min;
using numeric::max;
typedef numeric::xyzVector<platform::Real> Vec;
typedef numeric::xyzMatrix<platform::Real> Mat;

void
get_xform_stats(
	Xform const & sir,
	Xform const & sjr,
	Real& dx, Real& dy, Real& dz,
	Real& ex, Real& ey, Real& ez
){
	Vec d = ~sir * sjr.t;
	dx = d.x();
	dy = d.y();
	dz = d.z();
	Mat R = sir.R.transposed()*sjr.R;
	// Real ang;
	// Vec axis = rotation_axis(R,ang);
	// ang = numeric::conversions::degrees(ang);
	// std::cout << axis << " " << ang << std::endl;

	Real phi,psi,theta;
	if ( fabs(R.zx())-1.0 < 0.00001 ) {
		theta = -asin(R.zx());
		Real const ctheta = cos(theta);
		psi = atan2(R.zy()/ctheta,R.zz()/ctheta);
		phi = atan2(R.yx()/ctheta,R.xx()/ctheta);
	} else {
		if ( R.zx() < 0.0 ) {
			theta = numeric::constants::d::pi / 2.0;
			psi = atan2(R.xy(),R.xz());
		} else {
			theta = -numeric::constants::d::pi / 2.0;
			psi = atan2(-R.xy(),-R.xz());
		}
		phi = 0;
	}
	ex = psi;
	ey = theta;
	ez = phi;
}


XfoxmScore::XfoxmScore(
	std::string datadir
){
	hh = new char[16*16*16*24*12*24];
	he = new char[16*16*16*24*12*24];
	hl = new char[16*16*16*24*12*24];
	ee = new char[16*16*16*24*12*24];
	el = new char[16*16*16*24*12*24];
	ll = new char[16*16*16*24*12*24];
	fillarray(hh,datadir+"/hhpb.dat.gz.hist6.dat.gz.bin");
	fillarray(he,datadir+"/hepb.dat.gz.hist6.dat.gz.bin");
	fillarray(hl,datadir+"/hlpb.dat.gz.hist6.dat.gz.bin");
	fillarray(ee,datadir+"/eepb.dat.gz.hist6.dat.gz.bin");
	fillarray(el,datadir+"/elpb.dat.gz.hist6.dat.gz.bin");
	fillarray(ll,datadir+"/llpb.dat.gz.hist6.dat.gz.bin");
	// makebinary(hh,datadir+"/hhpb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
	// makebinary(he,datadir+"/hepb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
	// makebinary(hl,datadir+"/hlpb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
	// makebinary(ee,datadir+"/eepb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
	// makebinary(el,datadir+"/elpb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
	// makebinary(ll,datadir+"/llpb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
	// utility_exit_with_message("made binary file");
}
void
XfoxmScore::fillarray(
	char *a, std::string fname
){
	std::cout << "reading " << fname << std::endl;
	utility::io::izstream in(fname,std::ios::binary);
	if ( !in.good() ) utility_exit_with_message("bad file");
	in.read(a,16*16*16*24*12*24);
	in.close();
}
void
XfoxmScore::makebinary(
	char *a, std::string fname
){
	std::cout << "reading " << fname << std::endl;
	utility::io::izstream in(fname);
	if ( !in.good() ) utility_exit_with_message("bad file");
	for ( int i = 0; i < 16*16*16*24*12*24; ++i ) {
		float f;
		in >> f;
		int tmp = (int)((log(f)+2.5)*12.0);
		a[i] =  (f==0.0f) ? ((char)-127) : ((char)max(-126,min(128,tmp)));
		// std::cout << f << " " << tmp << " " << (int)a[i] << std::endl;
		// if( i > 100) utility_exit_with_message("FOO");
	}
	in.close();
	std::cout << "writing: " << fname+".bin.gz" << std::endl;
	utility::io::ozstream out(fname+".bin",std::ios::out | std::ios::binary);
	out.write(a,16*16*16*24*12*24);
	out.close();
	// std::cout << "writing: " << fname+".txt.gz for testing" << std::endl;
	// utility::io::ozstream out2(fname+".txt.gz",std::ios::out);
	// for(int i = 0; i < 16*16*16*24*12*24; ++i) out2 << (int)a[i] << std::endl;
	// out2.close();
}

float
XfoxmScore::score(
	Xform const & s1,
	Xform const & s2,
	char ss1, char ss2
) const {
	using numeric::constants::d::pi_2;
	if ( ss1=='L' || ss2=='L' ) return 0.0;
	if ( (~s1*(s2.t)).length_squared() > 64.0 ) return 0.0;
	char *a = hh;
	if ( ss1=='E' && ss2=='E' ) a = ee;
	else if ( ss2=='E' || ss1=='E' ) a = he;
	// if( ss1=='E' && ss2=='E' ) a = ee;
	// // if( ss1=='L' && ss2=='L' ) a = ll;
	// else if( ss1=='H' && ss2=='E' || ss1=='E' && ss2=='H' ) a = he;
	// // if( ss1=='H' && ss2=='L' || ss1=='L' && ss2=='H' ) a = hl;
	// // if( ss1=='E' && ss2=='L' || ss1=='L' && ss2=='E' ) a = el;
	Real dx,dy,dz,ex,ey,ez;
	if ( (ss1=='E' && ss2=='H') || (ss1=='L' && ss2=='H') || (ss1=='L' && ss2=='E') ) {
		get_xform_stats(s2,s1,dx,dy,dz,ex,ey,ez); // reverse
	} else get_xform_stats(s1,s2,dx,dy,dz,ex,ey,ez);
	int idx = static_cast<int>(dx+8.0); // 0.0-0.999 -> 10
	int idy = static_cast<int>(dy+8.0);
	int idz = static_cast<int>(dz+8.0);
	int iex = static_cast<int>(ex/pi_2*24.0 + 12.0);
	int iey = static_cast<int>(ey/pi_2*24.0 +  6.0);
	int iez = static_cast<int>(ez/pi_2*24.0 + 12.0);
	int index = idx + 16*idy + 16*16*idz + 16*16*16*iex + 16*16*16*24*iey + 16*16*16*24*12*iez;
	if ( 0 > index || index >= 16*16*16*24*12*24 ) utility_exit_with_message("FOO");
	// expensive memory lookup
	char val = a[index];
	//
	// std::cout << (int)val << "'" << val << "'"<< std::endl;
	return val == -127 ? 0.0 : min(10.0,exp(((float)val)/12.0-2.5));
	// return val == -127 ? -15.0 : ((float)val)/12.0;
}

float
XfoxmScore::score(
	core::pose::Pose const & pose,
	platform::Size rsd1,
	platform::Size rsd2
) const {
	if ( !pose.residue(rsd1).is_protein() ) return 0.0;
	if ( !pose.residue(rsd2).is_protein() ) return 0.0;
	if ( !pose.residue(rsd1).has("CB") ) return 0.0;
	if ( !pose.residue(rsd2).has("CB") ) return 0.0;
	Vec CBi = pose.residue(rsd1).xyz("CB");
	Vec CAi = pose.residue(rsd1).xyz("CA");
	Vec  Ni = pose.residue(rsd1).xyz( "N");
	Xform sir(CBi,CAi,Ni);
	Vec CBj = pose.residue(rsd2).xyz("CB");
	Vec CAj = pose.residue(rsd2).xyz("CA");
	Vec  Nj = pose.residue(rsd2).xyz( "N");
	if ( CBi.distance_squared(CBj) > 64.0 ) return 0.0;
	Xform sjr(CBj,CAj,Nj);
	return score(sir,sjr,pose.secstruct(rsd1),pose.secstruct(rsd2));
}

float
XfoxmScore::score(
	core::pose::Pose & pose,
	bool compute_ss
) const {
	if ( compute_ss ) {
		core::scoring::dssp::Dssp dssp(pose);
		dssp.insert_ss_into_pose(pose);
	}
	float tot_score = 0.0;
	for ( platform::Size ir = 1; ir <= pose.n_residue(); ++ir ) {
		for ( platform::Size jr = ir+1; jr <= pose.n_residue(); ++jr ) {
			float s1 = score(pose,ir,jr);
			float s2 = score(pose,jr,ir);
			// if( s1 < 0.0f ) std::cout << s1 << std::endl;
			// if( s2 < 0.0f ) std::cout << s2 << std::endl;
			tot_score += s1;
			tot_score += s2;
		}
	}
	return tot_score;
}
float
XfoxmScore::score(
	core::pose::Pose const & pose
) const {
	float tot_score = 0.0;
	for ( platform::Size ir = 1; ir <= pose.n_residue(); ++ir ) {
		for ( platform::Size jr = ir+1; jr <= pose.n_residue(); ++jr ) {
			float s1 = score(pose,ir,jr);
			float s2 = score(pose,jr,ir);
			// if( s1 < 0.0f ) std::cout << s1 << std::endl;
			// if( s2 < 0.0f ) std::cout << s2 << std::endl;
			tot_score += s1;
			tot_score += s2;
		}
	}
	return tot_score;
}


} // end namespace sic_dock
} // end namespace protocols
