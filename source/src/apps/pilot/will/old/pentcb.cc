// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include <apps/pilot/will/will_util.ihh>
#include <apps/pilot/will/mynamespaces.ihh>

#include <sys/stat.h>


#define CONTACT_D2 20.25
#define CONTACT_TH 5
#define NSS 8192

typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

inline Real const sqr(Real const r) { return r*r; }
inline Real sigmoidish_neighbor( Real const & sqdist ) {
	Real ub = CONTACT_D2+5.0; Real lb = CONTACT_D2-5.0;
	if     ( sqdist > ub ) return 0.0;
	else if ( sqdist < lb ) return 1.0;
	else return sqr(1.0  - sqr( (sqrt(sqdist) - lb) / (ub - lb) ) );
}

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using core::pose::Pose;
using core::conformation::ResidueOP;

static basic::Tracer TR( "pentcb" );
static core::io::silent::SilentFileData sfd;

double
sicfast(
	vector1<Vec> & pa , vector1<Vec> & pb ,
	vector1<Vec> & cba, vector1<Vec> & cbb,
	Real & cbcount,      double BIN = 1.8
){
	// get bounds for plane hashes
	double xmx1=-9e9,xmn1=9e9,ymx1=-9e9,ymn1=9e9,xmx=-9e9,xmn=9e9,ymx=-9e9,ymn=9e9;
	for ( vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia ) {
		xmx1 = max(xmx1,ia->x()); xmn1 = min(xmn1,ia->x());
		ymx1 = max(ymx1,ia->y()); ymn1 = min(ymn1,ia->y());
	}
	for ( vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib ) {
		xmx = max(xmx,ib->x()); xmn = min(xmn,ib->x());
		ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
	}
	xmx = min(xmx,xmx1); xmn = max(xmn,xmn1);
	ymx = min(ymx,ymx1); ymn = max(ymn,ymn1);
	int xlb = (int)floor(xmn/BIN)-2; int xub = (int)ceil(xmx/BIN)+2; // one extra on each side for correctness,
	int ylb = (int)floor(ymn/BIN)-2; int yub = (int)ceil(ymx/BIN)+2; // and one extra for outside atoms
	// insert points into hashes
	int const xsize = xub-xlb+1;
	int const ysize = yub-ylb+1;
	ObjexxFCL::FArray2D<Vec> ha(xsize,ysize,Vec(0,0,-9e9)),hb(xsize,ysize,Vec(0,0,9e9));
	for ( vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia ) {
		int const ix = (int)ceil(ia->x()/BIN)-xlb;
		int const iy = (int)ceil(ia->y()/BIN)-ylb;
		if ( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
		if ( ha(ix,iy).z() < ia->z() ) ha(ix,iy) = *ia;
	}
	for ( vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib ) {
		int const ix = (int)ceil(ib->x()/BIN)-xlb;
		int const iy = (int)ceil(ib->y()/BIN)-ylb;
		if ( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
		if ( hb(ix,iy).z() > ib->z() ) hb(ix,iy) = *ib;
	}
	// check hashes for min dis
	int imna=0,jmna=0,imnb=0,jmnb=0;
	double mindis = 9e9;
	for ( int i = 1; i <= xsize; ++i ) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
		for ( int j = 1; j <= ysize; ++j ) {
			for ( int k = -2; k <= 2; ++k ) {
				if ( i+k < 1 || i+k > xsize ) continue;
				for ( int l = -2; l <= 2; ++l ) {
					if ( j+l < 1 || j+l > ysize ) continue;
					double const xa = ha(i  ,j  ).x();
					double const ya = ha(i  ,j  ).y();
					double const xb = hb(i+k,j+l).x();
					double const yb = hb(i+k,j+l).y();
					double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);
					if ( d2 < BIN*BIN*4.0 ) {
						double dz = hb(i+k,j+l).z() - ha(i,j).z() - sqrt(BIN*BIN*4.0-d2);
						if ( dz < mindis ) {
							mindis = dz;
						}
					}
				}
			}
		}
	}
	cbcount = 0.0;
	for ( vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia ) {
		for ( vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib ) {
			cbcount += sigmoidish_neighbor( ib->distance_squared( (*ia) + (mindis*Vec(0,0,1)) ) );
		}
	}
	return mindis;
}

struct Hit {
	int iss,irt,sym;
	Real cbc;
	core::kinematics::Stub s1,s2;
	Hit(int is, int ir, Real cb, int sm) : iss(is),irt(ir),cbc(cb),sym(sm) {}
};
bool cmp(Hit i,Hit j) { return i.cbc > j.cbc; }


void TEST(Pose const & init, Hit const & h, Real cbcount) {
	static int ntest = 0;
	++ntest;
	TR << "passed: " << ntest << endl;
	Pose p(init),q(init);
	core::kinematics::Stub s(init.xyz(AtomID(1,1)),init.xyz(AtomID(2,1)),init.xyz(AtomID(3,1)));
	xform_pose_rev(p,s); xform_pose(p,h.s1);
	xform_pose_rev(q,s); xform_pose(q,h.s2);
	Real mycbc = 0.0;
	Real mxclash = 9e9;
	for ( Size ir = 1; ir <= p.size(); ++ir ) {
		for ( Size ia = 1; ia <= ((p.residue(ir).has("CB"))?5:4); ++ia ) {
			Vec ip = p.xyz(AtomID(ia,ir));
			for ( Size jr = 1; jr <= q.size(); ++jr ) {
				for ( Size ja = 1; ja <= ((q.residue(jr).has("CB"))?5:4); ++ja ) {
					Vec jp = q.xyz(AtomID(ja,jr));
					if ( ia==5 && ja==5 ) mycbc += sigmoidish_neighbor(ip.distance_squared(jp));
					mxclash = min(mxclash,ip.distance_squared(jp));
				}
			}
		}
	}
	if ( mxclash < 3.2*3.2 ) {
		p.dump_pdb("fail1.pdb");
		q.dump_pdb("fail2.pdb");
		utility_exit_with_message("fail: clash! "+string_of(sqrt(mxclash)));
	}
	if ( fabs(cbcount-mycbc) > 0.1 ) {
		utility_exit_with_message("cbcount mismatch "+str(cbcount)+" "+str(mycbc));
	}
}


void dock(Pose const init, std::string const & fn, vector1<xyzVector<double> > const & ssamp) {
	vector1<double> asamp; for ( Real i = 0; i < 180; ++i ) asamp.push_back(i);

	vector1<Vec> bb0,cb0;
	for ( int ir = 1; ir <= init.size(); ++ir ) {
		if ( !init.residue(ir).is_protein() ) continue;
		for ( int ia = 1; ia <= ((init.residue(ir).has("CB"))?5:4); ++ia ) bb0.push_back(init.xyz(AtomID(ia,ir)));
		if ( init.secstruct(ir)=='H' && init.residue(ir).has("CB") ) cb0.push_back(init.xyz(AtomID(5,ir)));
	}

	vector1<Mat> Rsym(6);
	Rsym[2] = rotation_matrix_degrees(Vec(1,0,0),180.0);
	Rsym[3] = rotation_matrix_degrees(Vec(1,0,0),120.0);
	Rsym[4] = rotation_matrix_degrees(Vec(1,0,0), 90.0);
	Rsym[5] = rotation_matrix_degrees(Vec(1,0,0), 72.0);
	Rsym[6] = rotation_matrix_degrees(Vec(1,0,0), 60.0);

	vector1<vector1<Hit> > hits(6);

	for ( int iss = 1; iss <= ssamp.size(); ++iss ) {
		if ( iss%1000==0 ) TR << iss << " of " << NSS << std::endl;
		Vec axs = ssamp[iss];
		for ( int irt = 1; irt <= asamp.size(); ++irt ) {
			//if(axs.x() < 0.0) continue;
			//Real const u = uniform();
			Real const rot = asamp[irt];
			Mat const R = rotation_matrix_degrees(axs,rot);

			vector1<Vec> bb1 = bb0;
			vector1<Vec> cb1 = cb0;
			for ( vector1<Vec>::iterator i = bb1.begin(); i != bb1.end(); ++i ) *i = R*(*i);
			for ( vector1<Vec>::iterator i = cb1.begin(); i != cb1.end(); ++i ) *i = R*(*i);

			for ( int ic = 2; ic < 7; ic++ ) {
				vector1<Vec> bb2 = bb1;
				vector1<Vec> cb2 = cb1;
				for ( vector1<Vec>::iterator i = bb2.begin(); i != bb2.end(); ++i ) *i = Rsym[ic]*(*i);
				for ( vector1<Vec>::iterator i = cb2.begin(); i != cb2.end(); ++i ) *i = Rsym[ic]*(*i);
				Real cbc;
				Real t = sicfast(bb1,bb2,cb1,cb2,cbc);
				if ( cbc >= CONTACT_TH ) {
					Hit h(iss,irt,cbc,ic);
					h.s1.from_four_points(bb1[1],bb1[1],bb1[2],bb1[3]);
					h.s2.from_four_points(bb2[1],bb2[1],bb2[2],bb2[3]);
					h.s1.v += t*Vec(0,0,1);
					hits[ic].push_back(h);

					TEST(init,h,h.cbc);

				}
			}
		}
	}

	for ( int ic = 2; ic < 7; ic++ ) {
		std::sort(hits[ic].begin(),hits[ic].end(),cmp);
		for ( int i = 1; i <= min((Size)10,hits[ic].size()); ++i ) {
			Hit & h(hits[ic][i]);
			cout << "RESULT " << fn << " " << h.sym << " " << h.iss << " " << NSS  << " " << h.irt << " " << h.cbc << endl;
		}
	}

}


int main(int argc, char *argv[]) {

	try {

		devel::init(argc,argv);
		using namespace basic::options::OptionKeys;

		vector1<xyzVector<double> > ssamp(NSS);
		{
			izstream is;
			basic::database::open(is,"geometry/sphere_"+str(NSS)+".dat");
			for ( int i = 1; i <= NSS; ++i ) {
				double x,y,z;
				is >> x >> y >> z;
				ssamp[i] = xyzVector<double>(x,y,z);
			}
			is.close();
		}
		for ( Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn ) {
			string fn = option[in::file::s]()[ifn];
			Pose pnat;
			TR << "searching " << fn << std::endl;
			core::import_pose::pose_from_file(pnat,fn, core::import_pose::PDB_file);
			trans_pose(pnat,-center_of_geom(pnat,1,pnat.size()));
			core::scoring::dssp::Dssp dssp(pnat);
			dssp.insert_ss_into_pose(pnat);
			//if( pnat.size() > 150 ) continue;
			Size cyscnt=0, nhelix=0;
			for ( Size ir = 2; ir <= pnat.size()-1; ++ir ) {
				if ( pnat.secstruct(ir) == 'H' ) nhelix++;
				//if(!pnat.residue(ir).is_protein()) goto cont1;
				if ( pnat.residue(ir).is_lower_terminus() ) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
				if ( pnat.residue(ir).is_upper_terminus() ) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
				if ( pnat.residue(ir).name3()=="CYS" ) { if ( ++cyscnt > 3 ) goto cont1; }
			} goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
			if ( nhelix < 20 ) continue;
			Pose pala(pnat);
			dock(pala,fn,ssamp);
		}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

