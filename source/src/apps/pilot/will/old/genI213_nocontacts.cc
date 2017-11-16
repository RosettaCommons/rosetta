// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
//#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/id/SequenceMapping.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_constants.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <sys/stat.h>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>


#define CONTACT_D2 36.0
#define CONTACT_TH1 12
#define CONTACT_TH2 15

#define ATET 54.735610317245360079 // asin(sr2/sr3)
#define AOCT 35.264389682754668343 // asin(sr1/sr3)

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using core::pose::Pose;
using core::conformation::ResidueOP;

static basic::Tracer TR( "genI213" );
static core::io::silent::SilentFileData sfd;


//////////////////////////////////////////////////////////
class SameResidueTypeConstraint : public core::scoring::constraints::Constraint
{
public:

	SameResidueTypeConstraint(
		Size seqpos1,
		Size seqpos2,
		Real bonus
	) : core::scoring::constraints::Constraint(core::scoring::res_type_constraint),seqpos1_(seqpos1),seqpos2_(seqpos2),bonus_(bonus) {
		atom_ids_.push_back(AtomID(1,seqpos1_));
		atom_ids_.push_back(AtomID(1,seqpos2_));
	}

	virtual ~SameResidueTypeConstraint() {}

	virtual
	Size
	natoms() const { return 2; }

	void
	show( std::ostream & out ) const { out << "SameResidueTypeConstraint " << seqpos1_ << " " << seqpos2_ << " " << bonus_; }

	virtual
	core::scoring::constraints::ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const {
		core::Size newseqpos1 = seqmap[ seqpos1_ ];
		core::Size newseqpos2 = seqmap[ seqpos2_ ];
		if ( newseqpos1 != 0 && newseqpos2 != 0 ) {
			return core::scoring::constraints::ConstraintOP( new SameResidueTypeConstraint( newseqpos1, newseqpos2, bonus_ ) );
		} else {
			return NULL;
		}
	}

	virtual
	AtomID const &
	atom( Size const index ) const {
		if ( index==1 ) return atom_ids_[1];
		if ( index==2 ) return atom_ids_[2];
		return core::id::GLOBAL_BOGUS_ATOM_ID;
	}

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	bool operator == ( core::scoring::constraints::Constraint const & otherin ) const {
		if ( !dynamic_cast< SameResidueTypeConstraint const * > ( &otherin ) ) return false;
		SameResidueTypeConstraint const & other( static_cast< SameResidueTypeConstraint const & > (otherin) );
		return seqpos1_==other.seqpos1_ && seqpos2_==other.seqpos2_ && seqpos1_==other.bonus_;
	}

	// virtual ConstraintOP remapped_clone(
	//  pose::Pose const & src,
	//  pose::Pose const & dest,
	//  id::SequenceMappingCOP map = NULL
	// ) const;

	virtual
	void
	score( core::scoring::constraints::XYZ_Func const & xyz_func, core::scoring::EnergyMap const & , core::scoring::EnergyMap & emap ) const {
		// Real const weight(weights[ this->score_type() ] );
		if ( xyz_func.residue(seqpos1_).name3() == xyz_func.residue(seqpos2_).name3() ) emap[ this->score_type() ] -= bonus_;
	}

	virtual
	void
	fill_f1_f2(
		AtomID const & ,
		core::scoring::constraints::XYZ_Func const & ,
		Vec & ,
		Vec & ,
		core::scoring::EnergyMap const &
	) const {}

	virtual core::scoring::constraints::ConstraintOP
	clone() const {
		return new SameResidueTypeConstraint( *this );
	}

private:
	Size seqpos1_,seqpos2_;
	core::Real bonus_;
	utility::vector1< AtomID > atom_ids_;
};


//////////////////////////////////////////////////////////


inline Real sqr(Real const r) { return r*r; }
inline Real sigmoidish_neighbor( Real const & sqdist ) {
	if ( sqdist > 9.*9. ) {
		return 0.0;
	} else if ( sqdist < 6.*6. ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0  - sqr( (dist - 6.) / (9. - 6.) ) );
	}
}

vector1<Size> get_scanres(Pose const & pose) {
	vector1<Size> scanres;
	//if(basic::options::option[basic::options::OptionKeys::willmatch::residues].user()) {
	//TR << "input scanres!!!!!!" << std::endl;
	//scanres = basic::options::option[basic::options::OptionKeys::willmatch::residues]();
	//} else {
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( !pose.residue(i).has("N" ) ) { continue; }
		if ( !pose.residue(i).has("CA") ) { continue; }
		if ( !pose.residue(i).has("C" ) ) { continue; }
		if ( !pose.residue(i).has("O" ) ) { continue; }
		if ( !pose.residue(i).has("CB") ) { continue; }
		if ( pose.residue(i).name3()=="PRO" ) { continue; }
		scanres.push_back(i);
	}
	//}
	return scanres;
}

bool file_exists(string strFilename) {
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;
	// Attempt to get the file attributes
	intStat = stat(strFilename.c_str(),&stFileInfo);
	if ( intStat == 0 ) {
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	} else {
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}
	return(blnReturn);
}


void dumpsym(Pose const & pose, Mat R2, Mat R3a, Mat R3b, Vec cen2, string fname) {
	vector1<Vec> seenit;
	Mat R3[3];
	R3[0] = Mat::identity();
	R3[1] = R3a;
	R3[2] = R3b;
	TR << "output" << std::endl;
	vector1<string> ANAME(3);
	ANAME[1] = " N  ";
	ANAME[2] = " CA ";
	ANAME[3] = " C  ";
	string const & CHAIN( utility::UPPERCASE_LETTERS );
	ozstream out( fname );
	Size acount=0,rcount=0,ccount=0;
	for ( Size i3a = 0; i3a < 3; i3a++ ) {
		for ( Size i2a = 0; i2a < 2; i2a++ ) {
			for ( Size j3a = 0; j3a < 3; j3a++ ) {
				for ( Size j2a = 0; j2a < 2; j2a++ ) {
					for ( Size k3a = 0; k3a < 3; k3a++ ) {
						for ( Size k2a = 0; k2a < 2; k2a++ ) {
							for ( Size l3a = 0; l3a < 3; l3a++ ) {
								for ( Size l2a = 0; l2a < 2; l2a++ ) {
									for ( Size m3a = 0; m3a < 2; m3a++ ) {
										for ( Size m2a = 0; m2a < 2; m2a++ ) {
											for ( Size n3a = 0; n3a < 2; n3a++ ) {
												for ( Size n2a = 0; n2a < 2; n2a++ ) {

													Vec chk( pose.xyz(AtomID(1,1)) );
													chk = R3[i3a]*chk; if ( i2a ) chk = R2*(chk-cen2)+cen2;
													chk = R3[j3a]*chk; if ( j2a ) chk = R2*(chk-cen2)+cen2;
													chk = R3[k3a]*chk; if ( k2a ) chk = R2*(chk-cen2)+cen2;
													chk = R3[l3a]*chk; if ( l2a ) chk = R2*(chk-cen2)+cen2;
													chk = R3[m3a]*chk; if ( m2a ) chk = R2*(chk-cen2)+cen2;
													chk = R3[n3a]*chk; if ( n2a ) chk = R2*(chk-cen2)+cen2;
													for ( vector1<Vec>::const_iterator i = seenit.begin(); i != seenit.end(); ++i ) {
														if ( i->distance_squared(chk) < 1.0 ) goto cont2;
													}
													goto done2; cont2: continue; done2:
													seenit.push_back(chk);

													char chain = CHAIN[ccount];
													if ( (i2a+j2a+k2a+l2a+m2a+n2a) % 2 == 1 ) chain = 'B';
													for ( Size ir = 1; ir <= pose.size(); ++ir ) {
														if ( rcount >= 9999 ) {
															rcount = 0;
															ccount++;
														}
														Size rn = ++rcount;
														for ( Size ia = 1; ia <= 3; ia++ ) {
															Vec tmp(pose.residue(ir).xyz(ia));
															tmp = R3[i3a]*tmp; if ( i2a ) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[j3a]*tmp; if ( j2a ) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[k3a]*tmp; if ( k2a ) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[l3a]*tmp; if ( l2a ) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[m3a]*tmp; if ( m2a ) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[n3a]*tmp; if ( n2a ) tmp = R2*(tmp-cen2)+cen2;
															string X = F(8,3,tmp.x());
															string Y = F(8,3,tmp.y());
															string Z = F(8,3,tmp.z());
															out<<"ATOM  "<<I(5,++acount)<<' '<<ANAME[ia]<<' '<<"ALA"<<' '<<chain<<I(4,rn)<<"    "<<X<<Y<<Z<<F(6,2,1.0)<<F(6,2,0.0)<<'\n';
														}
													}
													out << "TER" << std::endl;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	out.close();
}


int pose_cbcount(Pose const & a, Pose const & b) {
	int count = 0;
	for ( Size i = 1; i <= a.size(); ++i ) {
		for ( Size j = 1; j <= b.size(); ++j ) {
			if ( a.residue(i).xyz(2).distance_squared(b.residue(j).xyz(2)) < 100.0 ) {
				count++;
			}
		}
	}
	return count;
}

double
sicfast(
	vector1<Vec>  pa,
	vector1<Vec>  pb,
	vector1<Vec> & cba,
	vector1<Vec> & cbb,
	Vec ori,
	int & cbcount,
	bool debug = false
) {
	double BIN = 2.0;

	// get points, rotated ro ori is 0,0,1, might already be done
	Mat rot = Mat::identity();
	if     ( ori.dot(Vec(0,0,1)) < -0.99999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	else if ( ori.dot(Vec(0,0,1)) <  0.99999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	if ( rot != Mat::identity() ) {
		for ( vector1<Vec>::iterator ia = pa.begin(); ia != pa.end(); ++ia ) *ia = rot*(*ia);
		for ( vector1<Vec>::iterator ib = pb.begin(); ib != pb.end(); ++ib ) *ib = rot*(*ib);
	}

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

	// TR << "BOUNDS " << xmn << " " << xmx << " " << ymn << " " << ymx << std::endl;
	// TR << "BOUNDS " << xlb << " " << xub << " " << ylb << " " << yub << std::endl;

	// insert points into hashes
	int const xsize = xub-xlb+1;
	int const ysize = yub-ylb+1;
	ObjexxFCL::FArray2D<Vec> ha(xsize,ysize,Vec(0,0,-9e9)),hb(xsize,ysize,Vec(0,0,9e9));
	for ( vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia ) {
		// int const ix = min(xsize,max(1,(int)ceil(ia->x()/BIN)-xlb));
		// int const iy = min(ysize,max(1,(int)ceil(ia->y()/BIN)-ylb));
		int const ix = (int)ceil(ia->x()/BIN)-xlb;
		int const iy = (int)ceil(ia->y()/BIN)-ylb;
		if ( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
		if ( ha(ix,iy).z() < ia->z() ) ha(ix,iy) = *ia;
	}
	for ( vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib ) {
		// int const ix = min(xsize,max(1,(int)ceil(ib->x()/BIN)-xlb));
		// int const iy = min(ysize,max(1,(int)ceil(ib->y()/BIN)-ylb));
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

					if ( d2 < 16.0 ) {
						double dz = hb(i+k,j+l).z() - ha(i,j).z() - sqrt(16.0-d2);
						if ( dz < mindis ) {
							mindis = dz;
							imna = i;
							jmna = j;
							imnb = i+k;
							jmnb = j+l;
						}
					}
				}
			}
		}
	}

	// {
	//  utility::io::ozstream out("cba.pdb");
	//  for(vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
	//    Vec viz = (*ia) + (mindis*ori);
	//    out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"A"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	//  }
	//  out.close();
	// }
	// {
	//  utility::io::ozstream out("cbb.pdb");
	//  for(vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
	//    Vec viz = (*ib);
	//    out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	//  }
	//  out.close();
	// }

	cbcount = 0;
	// utility::io::ozstream out("cb8.pdb");
	// TR << "CB0 " << cbcount << std::endl;
	for ( vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia ) {
		for ( vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib ) {
			if ( ib->distance_squared( (*ia) + (mindis*ori) ) < CONTACT_D2 ) {
				cbcount++;
				// Vec viz = (*ia) + (mindis*ori);
				// out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				// viz = *ib;
				// out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"B"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			}
		}
	}
	// out.close();
	// TR << "CB1 " << cbcount << std::endl;

	// // rotate points back -- needed iff pa/pb come by reference
	rot = rot.transposed();
	// if( rot != Mat::identity() ) {
	//  for(vector1<Vec>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
	//  for(vector1<Vec>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
	// }

	// uncomment this to get hashes in local space
	// rot = Mat::identity();
	// ori = Vec(0,0,1);

	if ( debug ) {
		{
			utility::io::ozstream out("hasha.pdb");
			for ( int i = 2; i <= xsize-1; ++i ) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
				for ( int j = 2; j <= ysize-1; ++j ) {
					Vec viz = rot*ha(i,j) + mindis*ori;
					if ( viz.z() < -9e8 || 9e8 < viz.z() ) continue;
					out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100+j)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				}
			}
			Vec viz = rot*ha(imna,jmna) + mindis*ori;
			out<<"HETATM"<<I(5,1000+imna)<<' '<<"MIN "<<' ' <<  "MIN"<<' '<<"B"<<I(4,100+jmna)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			out.close();
		}
		{
			utility::io::ozstream out("hashb.pdb");
			for ( int i = 2; i <= xsize-1; ++i ) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
				for ( int j = 2; j <= ysize-1; ++j ) {
					Vec viz = rot*hb(i,j);
					if ( viz.z() < -9e8 || 9e8 < viz.z() ) continue;
					out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"C"<<I(4,100+j)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				}
			}
			Vec viz = rot*hb(imnb,jmnb);
			out<<"HETATM"<<I(5,1000+imnb)<<' '<<"MIN "<<' ' <<  "MIN"<<' '<<"C"<<I(4,100+jmnb)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			out.close();
		}
	}

	return mindis;
}

double sicfast(
	Pose const & a,
	Pose const & b,
	Vec ori_in,
	int & cbcount
) {
	// get points, rotated ro ori is 0,0,1
	vector1<Vec> pa,pb;
	vector1<Vec> cba,cbb;
	Vec ori = ori_in.normalized();
	Mat rot = Mat::identity();
	if     ( ori.dot(Vec(0,0,1)) < -0.999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	else if ( ori.dot(Vec(0,0,1)) <  0.999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	for ( int i = 1; i <= (int)a.size(); ++i ) {
		if ( a.residue(i).name3()=="GLY" ) {
			cba.push_back(rot*Vec(a.residue(i).xyz(2)));
			for ( int j = 1; j <= 4; ++j ) pa.push_back(rot*Vec(a.residue(i).xyz(j)));
		} else {
			cba.push_back(rot*Vec(a.residue(i).xyz(5)));
			for ( int j = 1; j <= 5; ++j ) pa.push_back(rot*Vec(a.residue(i).xyz(j)));
		}
	}
	for ( int i = 1; i <= (int)b.size(); ++i ) {
		if ( b.residue(i).name3()=="GLY" ) {
			cbb.push_back(rot*Vec(b.residue(i).xyz(2)));
			for ( int j = 1; j <= 4; ++j ) pb.push_back(rot*Vec(b.residue(i).xyz(j)));
		} else {
			cbb.push_back(rot*Vec(b.residue(i).xyz(5)));
			for ( int j = 1; j <= 5; ++j ) pb.push_back(rot*Vec(b.residue(i).xyz(j)));
		}
	}
	return sicfast( pa, pb, cba, cbb, Vec(0,0,1), cbcount );
}


struct DsfHit {
	Real chi11,chi12,chi21,chi22;
	Size rsd1,rsd2;
};

vector1<DsfHit> find_dsf( Pose & p3a, Pose & p3b, ImplicitFastClashCheck const & ifc, vector1<vector1<Real> > cyschi1,
	Mat const & R2, Vec const & C2, Real , Real  ) {
	vector1<DsfHit> hits;
	core::chemical::ResidueTypeSetCAP  rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	Pose cys;
	make_pose_from_sequence(cys,"C",core::chemical::FA_STANDARD,false);
	remove_lower_terminus_type_from_pose_residue(cys,1);
	remove_upper_terminus_type_from_pose_residue(cys,1);
	Size Nres = p3a.size()/3;
	Mat R3a = rotation_matrix_degrees(Vec(0,0,1),120.0);
	Mat R3b = rotation_matrix_degrees(Vec(0,0,1),240.0);
	for ( Size ic = 2; ic < p3a.size(); ++ic ) {
		Vec cb1 = p3a.xyz(AtomID(5,ic));
		if ( cb1.distance_squared( R2*(cb1-C2)+C2 ) < 225.0 ) continue;
		Size rsd1 = (ic-1)%Nres+1;
		for ( Size jc = ic+1; jc < p3b.size(); ++jc ) {
			Vec cb2 = p3b.xyz(AtomID(5,jc));
			if ( cb1.distance_squared(cb2) > 5*5 ) continue;
			Size rsd2 = (jc-1)%Nres+1;
			core::conformation::Residue const rtmp1 = p3a.residue(ic);
			core::conformation::Residue const rtmp2 = p3b.residue(jc);
			p3a.replace_residue(ic,cys.residue(1),true);
			p3b.replace_residue(jc,cys.residue(1),true);
			for ( vector1<Real>::const_iterator i = cyschi1[rsd1].begin(); i != cyschi1[rsd1].end(); ++i ) {
				Real chi11 = *i;
				p3a.set_chi(1,ic,chi11);
				Vec sg1 = p3a.xyz(AtomID(6,ic));
				if ( !ifc.clash_check(    (R2*(sg1-C2)+C2),rsd1) ) continue;
				if ( !ifc.clash_check(R3a*(R2*(sg1-C2)+C2),rsd1) ) continue;
				if ( !ifc.clash_check(R3b*(R2*(sg1-C2)+C2),rsd1) ) continue;
				for ( vector1<Real>::const_iterator j = cyschi1[rsd2].begin(); j != cyschi1[rsd2].end(); ++j ) {
					Real chi12 = *j;
					p3b.set_chi(1,jc,chi12);
					Vec sg2 = p3b.xyz(AtomID(6,jc));
					if ( fabs(sg1.distance(sg2)-2.2) > 0.4 ) continue;
					if ( fabs(        angle_degrees(cb1,sg1,sg2    )-102.0) > 15.0 ) continue;
					if ( fabs(        angle_degrees(    sg1,sg2,cb2)-102.0) > 15.0 ) continue;
					if ( fabs(fabs(dihedral_degrees(cb1,sg1,sg2,cb2))-90.0) > 15.0 ) continue;
					if ( !ifc.clash_check(    sg2,rsd2) ) continue;
					if ( !ifc.clash_check(R3a*sg2,rsd2) ) continue;
					if ( !ifc.clash_check(R3b*sg2,rsd2) ) continue;
					DsfHit h;
					h.chi11 = chi11;
					h.chi12 = chi12;
					h.chi21 = dihedral_degrees( p3a.xyz(AtomID(2,ic)),cb1,sg1,sg2 );
					h.chi22 = dihedral_degrees( p3b.xyz(AtomID(2,jc)),cb2,sg2,sg1 );
					h.rsd1 = ic;
					h.rsd2 = jc;
					hits.push_back(h);
					goto done1;
					// p3a.set_chi(2,ic,h.chi21);
					// p3b.set_chi(2,jc,h.chi22);
					// p3a.dump_pdb("test1.pdb");
					// p3b.dump_pdb("test2.pdb");
					// utility_exit_with_message("aroistn");
				}
			}
			done1:
			p3a.replace_residue(ic,rtmp1,false);
			p3b.replace_residue(jc,rtmp2,false);
		}
	}

	// #define DSF_CB_SG2 3.1
	//   Real const r = 1.6454076;
	//   for( Size i = 1; i <= p3a.size(); ++i) {
	//     // SG1 pos where chi1 SG1 circle iscts chi1(2)/chi2(2) SG1 sphere
	//     Vec const bn1(p3a.xyz(AtomID(1,i)));
	//     Vec const ca1(p3a.xyz(AtomID(2,i)));
	//     Vec const cb1(p3a.xyz(AtomID(5,i)));
	//     if( cb1.distance_squared( p3b.xyz(AtomID(5,i)) ) < 100.0 ) continue;
	//     Vec const a1 = (cb1-ca1).normalized();
	//     Vec const c1 = ca1+a1*2.27887; // cen of SG circle 1
	//     for( Size j = i+1; j <= p3b.size(); ++j) {
	//       Vec const cb2(p3b.xyz(AtomID(5,j)));
	//       if( cb1.distance_squared(cb2) > 4.8*4.8 ) continue;
	//       Vec const bn2(p3b.xyz(AtomID(1,j)));
	//       Vec const ca2(p3b.xyz(AtomID(2,j)));
	//       Real const d1 = a1.dot(c1-cb2); // dist to cen of plane/sphere isct
	//       if(d1>=DSF_CB_SG2) continue;
	//       Real const r1 = sqrt(DSF_CB_SG2*DSF_CB_SG2-d1*d1); //radius of isct circle
	//       Vec  const e1 = cb2 + d1*a1; // cen of plane/sphere isct
	//       Real const f1 = c1.distance(e1);
	//       Real const g1 = (r*r+f1*f1-r1*r1)/2.0/f1; // cosA*r = dst to isct along c1 to e1
	//       Vec  const m1 = (e1-c1)/f1; // c1 to e1 unit vec
	//       Vec  const b1 = g1*m1 + c1; // midpt to two isct pts
	//       if(g1>=r) continue;
	//       Real const h1  = sqrt(r*r-g1*g1); // half dist between isct pts
	//       Vec  const n1 = a1.cross(m1); // unit vec pointing from midpt to isct pts
	//       Vec  const a2 = (cb2-ca2).normalized();
	//       Vec  const c2 = ca2+a2*2.27887;
	//       Real const d2 = a2.dot(c2-cb1);
	//       if(d2>=DSF_CB_SG2) continue;
	//       Real const r2 = sqrt(DSF_CB_SG2*DSF_CB_SG2-d2*d2);
	//       Vec  const e2 = cb1 + d2*a2;
	//       Real const f2 = c2.distance(e2);
	//       Real const g2 = (r*r+f2*f2-r1*r1)/2.0/f2;
	//       Vec  const m2 = (e2-c2)/f2;
	//       Vec  const b2 = g2*m2 + c2;
	//       if(g2>=r) continue;
	//       Real const h2  = sqrt(r*r-g2*g2);
	//       Vec  const n2 = a2.cross(m2);
	//       for( Size idr1 = 0; idr1 <= 1; ++idr1 ) {
	//         Vec const sg1 = b1+(idr1?h1:-h1)*n1;
	//         if( !ifc.clash_check(    (    sg1       ),(i-1)%Nres+1) ) continue;
	//         if( !ifc.clash_check(    (R2*(sg1-C2)+C2),(i-1)%Nres+1) ) continue;
	//         if( !ifc.clash_check(R3a*(    sg1       ),(i-1)%Nres+1) ) continue;
	//         if( !ifc.clash_check(R3a*(R2*(sg1-C2)+C2),(i-1)%Nres+1) ) continue;
	//         if( !ifc.clash_check(R3b*(    sg1       ),(i-1)%Nres+1) ) continue;
	//         if( !ifc.clash_check(R3b*(R2*(sg1-C2)+C2),(i-1)%Nres+1) ) continue;
	//         for( Size idr2 = 0; idr2 <= 1; ++idr2 ) {
	//           Vec const sg2 = b2+(idr2?h2:-h2)*n2;
	//           if( fabs(sg1.distance(sg2)-2.05) > 0.5 ) continue;
	//           if( fabs(fabs(dihedral_degrees(cb1,sg1,sg2,cb2))-90.0) > 15.0 ) continue; // dihedral ~ 90
	//           if( fabs(angle_degrees(sg1,cb1,ca1)-105.0) > 15.0 ) continue;
	//           if( fabs(angle_degrees(sg2,cb2,ca2)-105.0) > 15.0 ) continue;
	//           if( !ifc.clash_check(    (    sg2       ),(j-1)%Nres+1) ) continue;
	//           if( !ifc.clash_check(    (R2*(sg2-C2)+C2),(j-1)%Nres+1) ) continue;
	//           if( !ifc.clash_check(R3a*(    sg2       ),(j-1)%Nres+1) ) continue;
	//           if( !ifc.clash_check(R3a*(R2*(sg2-C2)+C2),(j-1)%Nres+1) ) continue;
	//           if( !ifc.clash_check(R3b*(    sg2       ),(j-1)%Nres+1) ) continue;
	//           if( !ifc.clash_check(R3b*(R2*(sg2-C2)+C2),(j-1)%Nres+1) ) continue;
	//           DsfHit h;
	//           h.rsd1 = (i-1)%Nres+1;
	//           h.rsd2 = (j-1)%Nres+1;
	//           h.chi11 = dihedral_degrees(bn1,ca1,cb1,sg1);
	//           h.chi12 = dihedral_degrees(bn2,ca2,cb2,sg2);
	//           // p3a.replace_residue(i,cys.residue(1),true);
	//           // p3b.replace_residue(i,cys.residue(1),true);
	//           // p3a.replace_residue(j,cys.residue(1),true);
	//           // p3b.replace_residue(j,cys.residue(1),true);
	//           // p3a.set_chi(1,i,h.chi11);
	//           // p3a.set_chi(1,j,h.chi12);
	//           // p3b.set_chi(1,i,h.chi11);
	//           // p3b.set_chi(1,j,h.chi12);
	//           // p3a.dump_pdb(lzs(i,3)+" "+lzs(j,3)+"_dsf1.pdb");
	//           // p3b.dump_pdb(lzs(i,3)+" "+lzs(j,3)+"_dsf2.pdb");
	//           //          utility_exit_with_message("arst");
	//           hits.push_back(h);
	//         }
	//       }
	//     }
	//   }
	return(hits);

	// def dsf(CA1,CB1,CA2,CB2):
	//    c = CA1+(CB1-CA1).normalized()*2.27887
	//    a = (CB1-CA1).normalized()
	//    r = 1.6454076
	//    drawringcar(c,a,r,[1,1,0],'cr')
	//    d = a.dot(c-CB2)
	//    r2 = sqrt( 3.4*3.4 - d*d )
	//    c2 = CB2 + d*a
	//    a2 = a
	//    drawringcar(c2,a2,r2,[1,1,1],'cd')
	//    d = (c-c2).length()
	//    d2 = (r*r+d*d-r2*r2)/2/d
	//    x = d2 * (c2-c).normalized() + c
	//    h = sqrt(r*r - d2*d2)
	//    a3 = a.cross(c2-c).normalized()
	//    (x+h*a3).show('p1')
	//    (x-h*a3).show('p2')

}

struct Hit {
	Size rsd,cbc;
	Real chi1,chi2;
	Vec axs,cen;
	bool sym;
};

void design_1comp(Pose & pose, ScoreFunctionOP sf, Size Ntri ){
	using namespace core::pack::task;
	//TR << "design" << std::endl;
	// sasa
	Pose psasa(pose);
	for ( Size ir = 1; ir <= psasa.size(); ++ir ) {
		if ( psasa.residue(ir).name3()!="GLY" ) {
			if ( !psasa.residue(ir).is_protein() ) continue;
			Vec CA = psasa.residue(ir).xyz("CA");
			Vec CB = psasa.residue(ir).xyz("CB");
			psasa.set_xyz( AtomID(5,ir), CA + 2.3*(CB-CA).normalized() );
		}
	}
	core::id::AtomID_Map< bool > atom_map;
	core::pose::initialize_atomid_map( atom_map, psasa, false );
	for ( Size ir = 1; ir <= psasa.size(); ++ir ) {
		if ( !psasa.residue(ir).is_protein() ) continue;
		atom_map.set(AtomID(1,ir),true);
		atom_map.set(AtomID(2,ir),true);
		atom_map.set(AtomID(3,ir),true);
		atom_map.set(AtomID(4,ir),true);
		if ( psasa.residue(ir).name3()!="GLY" )  atom_map.set(AtomID(5,ir), true );
		//else                                  atom_map.set(AtomID(2,ir), true );
		//for( Size ia = 1; ia <= min(Size(5),psasa.residue(ir).nheavyatoms()); ++ia ) {
		//atom_map.set(AtomID(ia,ir) , true );
		//}
	}
	core::id::AtomID_Map<Real> atom_sasa; utility::vector1<Real> sasa;
	core::scoring::calc_per_atom_sasa( psasa, atom_sasa, sasa, 3.0, false, atom_map );
	for ( Size ir = 1; ir <= Ntri; ++ir ) {
		if ( psasa.residue(ir).name3()=="GLY" ) sasa[ir] = 0;
		else                                 sasa[ir] = atom_sasa[AtomID(5,ir)];
	}

	// std::cout << "select bur=resi ";
	// for(Size i = 1; i <= Ntri; ++i) if(sasa[i] < 10.0) std::cout << i << "+";
	// std::cout << std::endl;
	// std::cout << "color white, chain A; color grey, chain B; color red, bur" << std::endl;

	// pose.dump_pdb("test.pdb");
	// utility_exit_with_message("dbg sasa");

	// Size nres = pose.size();
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	// task->initialize_extra_rotamer_flags_from_command_line();
	vector1< bool > aac(20,false);
	aac[core::chemical::aa_ala] = true;
	//aac[core::chemical::aa_cys] = true;
	aac[core::chemical::aa_asp] = true;
	//aac[core::chemical::aa_glu] = true;
	aac[core::chemical::aa_phe] = true;
	//aac[core::chemical::aa_gly] = true;
	//aac[core::chemical::aa_his] = true;
	aac[core::chemical::aa_ile] = true;
	//aac[core::chemical::aa_lys] = true;
	aac[core::chemical::aa_leu] = true;
	aac[core::chemical::aa_met] = true;
	aac[core::chemical::aa_asn] = true;
	//aac[core::chemical::aa_pro] = true;
	//aac[core::chemical::aa_gln] = true;
	//aac[core::chemical::aa_arg] = true;
	aac[core::chemical::aa_ser] = true;
	aac[core::chemical::aa_thr] = true;
	aac[core::chemical::aa_val] = true;
	//aac[core::chemical::aa_trp] = true;
	//aac[core::chemical::aa_tyr] = true;

	vector1< bool > aap(20,false);
	//aap[core::chemical::aa_ala] = true;
	//aap[core::chemical::aa_cys] = true;
	aap[core::chemical::aa_asp] = true;
	aap[core::chemical::aa_glu] = true;
	//aap[core::chemical::aa_phe] = true;
	//aap[core::chemical::aa_gly] = true;
	//aap[core::chemical::aa_his] = true;
	//aap[core::chemical::aa_ile] = true;
	aap[core::chemical::aa_lys] = true;
	//aap[core::chemical::aa_leu] = true;
	//aap[core::chemical::aa_met] = true;
	aap[core::chemical::aa_asn] = true;
	//aap[core::chemical::aa_pro] = true;
	aap[core::chemical::aa_gln] = true;
	aap[core::chemical::aa_arg] = true;
	aap[core::chemical::aa_ser] = true;
	aap[core::chemical::aa_thr] = true;
	//aap[core::chemical::aa_val] = true;
	//aap[core::chemical::aa_trp] = true;
	//aap[core::chemical::aa_tyr] = true;

	vector1<Size> iface(Ntri,0);
	for ( Size ir = 1; ir <= Ntri; ++ir ) {
		if ( pose.residue(ir).name3()=="CYS" || pose.residue(ir).name3()=="GLY" || pose.residue(ir).name3()=="PRO" ) {
			iface[ir] = 3;
			continue;
		}
		Real closestcb = 9e9;
		for ( Size jr = Ntri+1; jr <= 2*Ntri; ++jr ) {
			if ( pose.residue(jr).name3()=="GLY" ) continue;
			Real d = pose.xyz(AtomID(5,ir)).distance( pose.xyz(AtomID(5,jr)) );
			closestcb = min(closestcb,d);
		}
		if ( closestcb < 9.0 ) {
			//if     ( sasa[ir] > 10.0 ) iface[ir] = 1; // if exposed, is priphery
			//else if(closestcb <  7.0 ) iface[ir] = 2;
			//else                       iface[ir] = 0;
			iface[ir] = 2;
		}
	}
	for ( Size i=1; i<=Ntri; i++ ) iface.push_back(iface[i]);

	// std::cout << "select priph=resi ";
	// for(Size i = 1; i <= Ntri; ++i) if(iface[i]==1) std::cout << i << "+";
	// std::cout << std::endl;
	// std::cout << "select iface=resi ";
	// for(Size i = 1; i <= Ntri; ++i) if(iface[i]==2) std::cout << i << "+";
	// std::cout << std::endl;
	// std::cout << "color white, chain A; color grey, chain B; color blue, priph; color red, iface" << std::endl;

	// pose.dump_pdb("test.pdb");
	// utility_exit_with_message("dbg iface sel");

	for ( Size i = 1; i <= 2*Ntri; ++i ) {
		if       ( iface[i] == 3 ) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if ( iface[i] == 2 ) {
			bool tmp = aac[pose.residue(i).aa()];
			aac[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aac);
			aac[pose.residue(i).aa()] = tmp;
			//task->nonconst_residue_task(i).or_ex1_sample_level(core::pack::task::EX_TWO_HALF_STEP_STDDEVS);
			//task->nonconst_residue_task(i).or_ex2_sample_level(core::pack::task::EX_TWO_HALF_STEP_STDDEVS);
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else if ( iface[i] == 1 ) {
			bool tmp = aap[pose.residue(i).aa()];
			aap[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aap);
			aap[pose.residue(i).aa()] = tmp;
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else if ( iface[i] == 0 ) {
			task->nonconst_residue_task(i).prevent_repacking();
			//task->nonconst_residue_task(i).restrict_to_repacking();
			//task->nonconst_residue_task(i).or_include_current(true);
			//task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		}

	}

	//Real rorig = sf->get_weight(core::scoring::fa_rep);
	//sf->set_weight(core::scoring::fa_rep,rorig/4.0);
	Real worig = sf->get_weight(core::scoring::res_type_constraint);
	if ( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);
	utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);

	utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst2;
	if ( !core::pose::symmetry::is_symmetric(pose) ) {
		for ( Size i = 1; i <= Ntri; ++i ) {
			if ( iface[i]==2||iface[i]==1 ) {
				res_cst2.push_back(new SameResidueTypeConstraint(i,Ntri+i,-3.0));
			}
		}
	}
	pose.add_constraints( res_cst );
	pose.add_constraints( res_cst2 );

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::pack::make_symmetric_PackerTask_by_truncation(pose,task);
		protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
		repack.apply(pose);
	} else {
		protocols::simple_moves::PackRotamersMover repack( sf, task );
		repack.apply(pose);
	}

	// REMOVE SER NOT HBONDED ACROSS IFACE

	// cleanup 2
	pose.remove_constraints( res_cst );
	pose.remove_constraints( res_cst2 );
	sf->set_weight(core::scoring::res_type_constraint,worig);

	//sf->set_weight(core::scoring::fa_rep,rorig);

	// //TR << "sc min" << std::endl;
	// core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	// movemap->set_jump(false);
	// movemap->set_bb(false);
	// movemap->set_chi(true);

	// if( core::pose::symmetry::is_symmetric(pose) ) {
	//  core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	//  protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );
	//  m.apply(pose);
	// } else {
	//  protocols::simple_moves::MinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );
	//  m.apply(pose);
	// }
	// //TR << "done" << std::endl;


}

void repack_iface(Pose & p, ScoreFunctionOP sf, Size Ntri, vector1<bool> & iface_io) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(p);
	bool useiface_io = iface_io.size() == p.size();
	for ( Size i = 1; i <= p.size(); ++i ) {
		bool iface = false;
		if ( i <= Ntri ) {
			for ( Size j = Ntri+1; j <= 2*Ntri; ++j ) if ( p.xyz(AtomID(5,i)).distance_squared(p.xyz(AtomID(5,j))) < 81.0 ) iface = true;
		} else if ( i <= 2*Ntri ) {
			for ( Size j =      1; j <=   Ntri; ++j ) if ( p.xyz(AtomID(5,i)).distance_squared(p.xyz(AtomID(5,j))) < 81.0 ) iface = true;
		}
		if ( useiface_io ) iface = iface_io[i];
		else            iface_io.push_back( iface );
		if ( !iface ) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_include_current(true);
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();

			task->nonconst_residue_task(i).and_extrachi_cutoff(1);
		}
	}
	//Real rorig = sf->get_weight(core::scoring::fa_rep);
	//sf->set_weight(core::scoring::fa_rep,rorig/4.0);

	if ( core::pose::symmetry::is_symmetric(p) ) {
		core::pack::make_symmetric_PackerTask_by_truncation(p,task);
		protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
		repack.apply(p);
	} else {
		protocols::simple_moves::PackRotamersMover repack( sf, task );
		repack.apply(p);
	}
	//sf->set_weight(core::scoring::fa_rep,rorig);
}


Real ddg(Pose const & p_in, ScoreFunctionOP sf, Size Ntri, Real & rholes, Real & , Real & , Real & ) {
	TR << "ddg" << std::endl;
	ScoreFunctionOP sfhsp,sfssp,sfsht;
	if ( core::pose::symmetry::is_symmetric(p_in) ) {
		sfhsp = new core::scoring::symmetry::SymmetricScoreFunction;
		sfssp = new core::scoring::symmetry::SymmetricScoreFunction;
		sfsht = new core::scoring::symmetry::SymmetricScoreFunction;
	} else {
		sfhsp = new core::scoring::ScoreFunction;
		sfssp = new core::scoring::ScoreFunction;
		sfsht = new core::scoring::ScoreFunction;
	}
	sfhsp->set_weight(core::scoring::hs_pair,1.0);
	sfssp->set_weight(core::scoring::ss_pair,1.0);
	sfsht->set_weight(core::scoring::sheet  ,1.0);

	Pose p(p_in);

	vector1<bool> iface;
	repack_iface(p,sf,Ntri,iface);

	Real s0 = sf->score(p);
	Real p0 = core::scoring::packing::compute_dec15_score(p);
	// Real hsp0 = sfhsp->score(p);
	// Real ssp0 = sfssp->score(p);
	// Real sht0 = sfsht->score(p);
	//p.dump_pdb("test1.pdb");
	for ( Size ir = 1; ir <= Ntri; ++ir ) for ( Size ia = 1; ia <= p.residue_type(ir).natoms(); ++ia ) p.set_xyz(AtomID(ia,ir),p.xyz(AtomID(ia,ir))+Vec(9999,9999,9999));

	repack_iface(p,sf,Ntri,iface);

	Real s1 = sf->score(p);
	// Real hsp1 = sfhsp->score(p);
	// Real ssp1 = sfssp->score(p);
	// Real sht1 = sfsht->score(p);

	for ( Size ir = 1; ir <= Ntri; ++ir ) for ( Size ia = 1; ia <= p.residue_type(ir).natoms(); ++ia ) p.set_xyz(AtomID(ia,ir),p.xyz(AtomID(ia,ir))-Vec(9999,9999,9999));
	//  p.dump_pdb("test2.pdb");
	//  utility_exit_with_message("ddg "+str(s0-s1));

	Pose tmp(p);
	for ( Size i = 1; i <= Ntri; ++i ) {
		if ( p.residue(i).is_lower_terminus() ) tmp.append_residue_by_jump(p.residue(i),1);
		else tmp.append_residue_by_bond(p.residue(i));
	}
	Real p1 = core::scoring::packing::compute_dec15_score(tmp);

	rholes = p0-p1;
	return s0 - s1;
}

void replace_nat_seq(Pose & p, Pose const & nat) {
	for ( Size i = 1; i <= p.size(); ++i ) {
		if ( p.residue(i).name3() == "ALA" ) {
			p.replace_residue(i, nat.residue((i-1)%nat.size()+1), true );
		}
	}
}

vector1<Hit> dock(Pose & init, string fname) {

	using namespace basic::options::OptionKeys;

	core::chemical::ResidueTypeSetCAP  rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	Pose cys,ala;
	make_pose_from_sequence(cys,"C",core::chemical::FA_STANDARD,false);
	make_pose_from_sequence(ala,"A",core::chemical::FA_STANDARD,false);
	remove_lower_terminus_type_from_pose_residue(cys,1);
	remove_upper_terminus_type_from_pose_residue(cys,1);
	remove_lower_terminus_type_from_pose_residue(ala,1);
	remove_upper_terminus_type_from_pose_residue(ala,1);
	//  add_variant_type_to_pose_residue(cys,"DISULF_PARTNER",1);

	core::scoring::dssp::Dssp dssp(init);
	dssp.insert_ss_into_pose(init);

	Pose pnat(init);

	//Real rholes0 = core::scoring::packing::compute_dec15_score(pnat);
	//utility_exit_with_message("aoriestn");


	Real ofst = option[matdes::design::grid_size_radius]();
	string D = option[out::file::o]();

	Vec com(0,0,0);
	for ( Size ir = 1; ir <= init.size(); ++ir ) {
		init.replace_residue(ir,ala.residue(1),true);
		com += init.xyz(AtomID(2,ir));
	}
	com /= init.size();
	Real mxd = 0;
	for ( Size ir = 1; ir <= init.size(); ++ir ) {
		if ( init.xyz(AtomID(5,ir)).distance(com) > mxd ) mxd = init.xyz(AtomID(5,ir)).distance(com);
	}
	mxd = (2*mxd+4.0)*(2*mxd+4.0);
	//TR << mxd << std::endl;

	protocols::scoring::ImplicitFastClashCheck ifc3(init,3.0);
	protocols::scoring::ImplicitFastClashCheck ifc2(init,2.5);

	vector1<vector1<Real> > cyschi1(init.size());
	for ( Size ir = 2; ir < init.size(); ++ir ) {
		init.replace_residue(ir,cys.residue(1),true);
		for ( Real chi = 0; chi < 360.0; chi += 5.0 ) {
			init.set_chi(1,ir,chi);
			if ( ifc3.clash_check(init.residue(ir).xyz("SG"),ir) ) {
				cyschi1[ir].push_back( chi );
			}
		}
		init.replace_residue(ir,ala.residue(1),true);
	}

	// Size nres = init.size();
	// ScoreFunctionOP sf = core::scoring::get_score_function();
	ScoreFunctionOP sf = new core::scoring::symmetry::SymmetricScoreFunction(core::scoring::get_score_function());

	Mat R3a = rotation_matrix_degrees(Vec(0,0,1), 120.0);
	Mat R3b = rotation_matrix_degrees(Vec(0,0,1),-120.0);
	Mat R3[3];
	R3[0] = Mat::identity();
	R3[1] = R3a;
	R3[2] = R3b;

	Pose p3a;
	for ( Size i3a = 0; i3a < 3; i3a++ ) {
		Pose tmp(init);
		rot_pose(tmp,R3[i3a]);
		p3a.append_residue_by_jump(tmp.residue(1),1);
		for ( Size ir = 2; ir <= tmp.size(); ++ir ) p3a.append_residue_by_bond(tmp.residue(ir));
	}

	vector1<Vec> pa;
	vector1<Vec> cba;
	for ( int i = 1; i <= (int)p3a.size(); ++i ) {
		if ( p3a.residue(i).name3()=="GLY" ) {
			cba.push_back(Vec(p3a.residue(i).xyz(2)));
			for ( int j = 1; j <= 4; ++j ) pa.push_back(Vec(p3a.residue(i).xyz(j)));
		} else {
			cba.push_back(Vec(p3a.residue(i).xyz(5)));
			for ( int j = 1; j <= 5; ++j ) pa.push_back(Vec(p3a.residue(i).xyz(j)));
		}
	}
	vector1<Hit> hits;
	for ( Size iaxs = 0; iaxs < 120; ++iaxs ) {
		if ( iaxs%12==0 ) TR << fname << " " << iaxs << std::endl;
		for ( Size isg1 = 0; isg1 <= 0; ++isg1 ) {
			for ( Size isg2 = 0; isg2 <= 1; ++isg2 ) {
				Vec axs(0,0,1);
				Real ang = isg1 ? ATET : AOCT;
				if ( isg2 ) ang = 180.0-ang;
				axs = rotation_matrix_degrees(Vec(1,0,0),      ang     ) * axs;
				axs = rotation_matrix_degrees(Vec(0,0,1),Real(iaxs)/2.0) * axs;
				Mat R2 = rotation_matrix_degrees(axs,180.0);
				Pose p3b(p3a);
				rot_pose(p3b,R2);
				vector1<Vec> pb;
				vector1<Vec> cbb;
				for ( int i = 1; i <= (int)p3b.size(); ++i ) {
					if ( p3b.residue(i).name3()=="GLY" ) {
						cbb.push_back(Vec(p3b.residue(i).xyz(2)));
						for ( int j = 1; j <= 4; ++j ) pb.push_back(Vec(p3b.residue(i).xyz(j)));
					} else {
						cbb.push_back(Vec(p3b.residue(i).xyz(5)));
						for ( int j = 1; j <= 5; ++j ) pb.push_back(Vec(p3b.residue(i).xyz(j)));
					}
				}
				for ( Size iori = 0; iori <= 360; ++iori ) {
					Vec ori = Vec(0,0,1).cross(axs).normalized();
					ori = rotation_matrix_degrees(axs,(Real)iori)*ori;
					int cbc=0;
					//Real t = sicfast(p3a,p3b,ori,cbc);
					Real t = sicfast( pa, pb, cba, cbb, ori, cbc ) + ofst; // !!!!!!!!
					Vec C2 = -t*ori/2.0;

					//correct cbc (ofst)
					cbc = 0;
					for ( Size ir = 1; ir <= p3a.size(); ++ir ) {
						Vec x = R2*(p3a.xyz(AtomID(5,ir))-C2)+C2;
						for ( Size jr = 1; jr <= p3a.size(); ++jr ) {
							Vec y = p3a.xyz(AtomID(5,jr));
							if ( x.distance_squared(y) < 36.0 ) cbc++;
						}
					}

					if ( cbc < CONTACT_TH1 ) continue;

					//dumpsym(init,R2,R3a,R3b,C2,"xtal.pdb");

					//ozstream xout("lattice.pdb");
					//TR << "BEGIN XTAL CHECK " << fname << " " << iaxs << " " << isg1 << " " << isg2 << " " << iori << " " << cbc << std::endl;
					vector1<Pose> toadd;
					vector1<Vec> olap,ori1,ori2,toaddcom;
					{ // now check xtal symm
						{
							Vec chk( com );
							olap.push_back(        chk);
							olap.push_back(    R3a*chk);
							olap.push_back(    R3b*chk);
							olap.push_back(R2*(    chk-C2)+C2);
							olap.push_back(R2*(R3a*chk-C2)+C2);
							olap.push_back(R2*(R3b*chk-C2)+C2);
							chk = Vec(0,0,1);
							ori1.push_back(       chk);
							ori1.push_back(   R3a*chk);
							ori1.push_back(   R3b*chk);
							ori1.push_back(R2*    chk);
							ori1.push_back(R2*R3a*chk);
							ori1.push_back(R2*R3b*chk);
							chk = Vec(0,1,0);
							ori2.push_back(       chk);
							ori2.push_back(   R3a*chk);
							ori2.push_back(   R3b*chk);
							ori2.push_back(R2*    chk);
							ori2.push_back(R2*R3a*chk);
							ori2.push_back(R2*R3b*chk);
						}
						for ( Size i3a = 0; i3a < 3; i3a++ ) {
							for ( Size i2a = 0; i2a < 2; i2a++ ) {
								for ( Size j3a = 0; j3a < 3; j3a++ ) {
									for ( Size j2a = 0; j2a < 2; j2a++ ) {
										for ( Size k3a = 0; k3a < 3; k3a++ ) {
											for ( Size k2a = 0; k2a < 2; k2a++ ) {
												for ( Size l3a = 0; l3a < 3; l3a++ ) {
													for ( Size l2a = 0; l2a < 2; l2a++ ) {
														for ( Size m3a = 0; m3a < 3; m3a++ ) {
															for ( Size m2a = 0; m2a < 2; m2a++ ) {
																for ( Size n3a = 0; n3a < 3; n3a++ ) {
																	for ( Size n2a = 0; n2a < 2; n2a++ ) {
																		Vec chk(com),or1(0,0,1),or2(0,1,0);
																		chk = R3[i3a] * chk; if ( i2a ) chk = R2*(chk-C2)+C2;
																		chk = R3[j3a] * chk; if ( j2a ) chk = R2*(chk-C2)+C2;
																		chk = R3[k3a] * chk; if ( k2a ) chk = R2*(chk-C2)+C2;
																		chk = R3[l3a] * chk; if ( l2a ) chk = R2*(chk-C2)+C2;
																		chk = R3[m3a] * chk; if ( m2a ) chk = R2*(chk-C2)+C2;
																		chk = R3[n3a] * chk; if ( n2a ) chk = R2*(chk-C2)+C2;
																		or1 = R3[i3a] * or1; if ( i2a ) or1 = R2*or1;
																		or1 = R3[j3a] * or1; if ( j2a ) or1 = R2*or1;
																		or1 = R3[k3a] * or1; if ( k2a ) or1 = R2*or1;
																		or1 = R3[l3a] * or1; if ( l2a ) or1 = R2*or1;
																		or1 = R3[m3a] * or1; if ( m2a ) or1 = R2*or1;
																		or1 = R3[n3a] * or1; if ( n2a ) or1 = R2*or1;
																		or2 = R3[i3a] * or2; if ( i2a ) or2 = R2*or2;
																		or2 = R3[j3a] * or2; if ( j2a ) or2 = R2*or2;
																		or2 = R3[k3a] * or2; if ( k2a ) or2 = R2*or2;
																		or2 = R3[l3a] * or2; if ( l2a ) or2 = R2*or2;
																		or2 = R3[m3a] * or2; if ( m2a ) or2 = R2*or2;
																		or2 = R3[n3a] * or2; if ( n2a ) or2 = R2*or2;
																		for ( vector1<Vec>::const_iterator i0=olap.begin(),i1=ori1.begin(),i2=ori2.begin(); i0 != olap.end(); ++i0,++i1,++i2 ) {
																			if ( i0->distance_squared(chk) < 1.0 && or1.dot( *i1 ) > 0.95 && or2.dot( *i2 ) > 0.95 ) goto cont6;
																		} goto done6; cont6: continue; done6:
																		olap.push_back(chk);
																		ori1.push_back(or1);
																		ori2.push_back(or2);
																		string X = F(8,3,chk.x());
																		string Y = F(8,3,chk.y());
																		string Z = F(8,3,chk.z());
																		//xout<<"ATOM  "<<I(5,1)<<' '<<" CA "<<' '<<"ALA"<<' '<<"Z"<<I(4,1)<<"    "<<X<<Y<<Z<<F(6,2,1.0)<<F(6,2,0.0)<<'\n';
																		if ( chk.distance_squared(com) > mxd ) continue;
																		//TR << chk << std::endl;
																		bool contact = false;
																		for ( Size ir = 1; ir <= init.size(); ++ir ) {
																			for ( Size ia = 1; ia <= 5; ++ia ) {
																				Vec X( init.xyz(AtomID(ia,ir)) );
																				X = R3[i3a] * X; if ( i2a ) X = R2*(X-C2)+C2;
																				X = R3[j3a] * X; if ( j2a ) X = R2*(X-C2)+C2;
																				X = R3[k3a] * X; if ( k2a ) X = R2*(X-C2)+C2;
																				X = R3[l3a] * X; if ( l2a ) X = R2*(X-C2)+C2;
																				X = R3[m3a] * X; if ( m2a ) X = R2*(X-C2)+C2;
																				X = R3[n3a] * X; if ( n2a ) X = R2*(X-C2)+C2;
																				if ( ! ifc3.clash_check(X) ) goto cont7;
																				if ( ia == 5 ) {
																					for ( vector1<Vec>::const_iterator jcb = cba.begin(); jcb != cba.end(); ++jcb ) {
																						if ( jcb->distance_squared(X) < CONTACT_D2 ) {
																							contact = true;
																						}
																					}
																				}
																			}
																		}
																		if ( contact ) {
																			Pose tmp(init);
																			//core::util::switch_to_residue_type_set(tmp,"centroid");
																			rot_pose(tmp,R3[i3a]); if ( i2a ) rot_pose(tmp,R2,C2);
																			rot_pose(tmp,R3[j3a]); if ( j2a ) rot_pose(tmp,R2,C2);
																			rot_pose(tmp,R3[k3a]); if ( k2a ) rot_pose(tmp,R2,C2);
																			rot_pose(tmp,R3[l3a]); if ( l2a ) rot_pose(tmp,R2,C2);
																			rot_pose(tmp,R3[m3a]); if ( m2a ) rot_pose(tmp,R2,C2);
																			rot_pose(tmp,R3[n3a]); if ( n2a ) rot_pose(tmp,R2,C2);
																			toadd.push_back(tmp);
																			toaddcom.push_back(chk);
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
						goto done7; cont7: continue; done7:
						if ( olap.size() < 25 ) continue;
					}

					// //         xout.close();
					// TR << "DONE XTAL SCAN " << olap.size() << std::endl;
					// add symm units which are relevant to iface
					// core::util::switch_to_residue_type_set(p3a,"centroid");
					// p3a.dump_pdb("testA.pdb");
					// trans_pose(p3b,-t*ori);
					// core::util::switch_to_residue_type_set(p3b,"centroid");
					// p3b.dump_pdb("testB.pdb");

					vector1<bool> keepme(toadd.size(),true);
					for ( Size iadd = 1; iadd <= toadd.size(); ++iadd ) {
						Vec const X( R2*(toaddcom[iadd]-C2)+C2 );
						for ( Size jadd = iadd+1; jadd <= toadd.size(); ++jadd ) {
							if ( X.distance_squared(toaddcom[jadd]) < 0.01 ) {
								if ( toaddcom[iadd].length() < toaddcom[jadd].length() ) {
									keepme[jadd] = false;
								} else {
									keepme[iadd] = false;
								}
							}
						}
					}
					Pose todes(p3a);
					int xcbc = 0;
					for ( Size iadd = 1; iadd <= toadd.size(); ++iadd ) {
						if ( !keepme[iadd] ) continue;
						Pose & pi(toadd[iadd]);
						//pi.dump_pdb("test"+str(iadd)+".pdb");
						// bool contact = false;
						for ( Size ir = 1; ir <= pi.size(); ++ir ) {
							Vec const X(  R2*(pi.xyz(AtomID(5,ir))-C2)+C2  );
							for ( vector1<Vec>::const_iterator jcb = cba.begin(); jcb != cba.end(); ++jcb ) {
								if ( X.distance_squared( *jcb ) < CONTACT_D2 ) {
									xcbc++;
									// contact = true;
								}
							}
							for ( Size jadd = 1; jadd <= toadd.size(); ++jadd ) {
								if ( !keepme[jadd] ) continue;
								Pose & pj(toadd[jadd]);
								for ( Size jr = 1; jr <= pj.size(); ++jr ) {
									if ( X.distance_squared( pj.xyz(AtomID(5,jr)) ) < CONTACT_D2 ) {
										xcbc++;
										// contact = true;
									}
								}
							}
						}
						// if(contact) {
						//   todes.append_residue_by_jump(pi.residue(1),1,"","",true);
						//   for(Size ir = 2; ir <= pi.size(); ++ir) {
						//     todes.append_residue_by_bond(pi.residue(ir));
						//   }
						// }
					}
					//if( xcbc > 0 ) continue;
					//xcbc += cbc;


					if ( cbc+xcbc < CONTACT_TH1 ) continue;

					trans_pose(p3b,-t*ori);
					vector1<DsfHit> dhits = find_dsf(p3a,p3b,ifc3,cyschi1,R2,C2,0,0);
					trans_pose(p3b,t*ori);

					if ( cbc+xcbc < CONTACT_TH2 && dhits.size() == 0 ) continue;

					//TR << "CBC " << iaxs << " " << isg1 << " " << isg2 << " " << iori << " " << cbc << std::endl;

					string tag = utility::file_basename(fname)+"_"+lzs(iaxs,3)+"_"+lzs(isg1,1)+"_"+lzs(isg2,1)+"_"+lzs(iori,3)+"_"+F(4,2,ofst)+"_"+lzs(cbc,3)+"_"+lzs(xcbc,3);

					trans_pose(todes,-C2);
					alignaxis(todes,Vec(0,0,1),axs,Vec(0,0,0));

					// Size realxcbc = 0;
					// for(Size ir = 1; ir <= todes.size(); ++ir) {
					//   Vec x = rotation_matrix_degrees(Vec(0,0,1),180.0)*todes.xyz(AtomID(5,ir));
					//   for(Size jr = 1; jr <= todes.size(); ++jr) {
					//     Vec y = todes.xyz(AtomID(5,jr));
					//     if( x.distance_squared(y) < 36.0 ) realxcbc++;
					//   }
					// }
					// Size realcbc = 0;
					// for(Size ir = 1; ir <= p3a.size(); ++ir) {
					//   Vec x = rotation_matrix_degrees(Vec(0,0,1),180.0)*todes.xyz(AtomID(5,ir));
					//   for(Size jr = 1; jr <= p3a.size(); ++jr) {
					//     Vec y = todes.xyz(AtomID(5,jr));
					//     if( x.distance_squared(y) < 36.0 ) realcbc++;
					//   }
					// }
					// TR << cbc << " " << xcbc << " " << realcbc << " " << realxcbc << std::endl;
					// todes.dump_pdb("test0.pdb");
					// rot_pose(todes,Vec(0,0,1),180.0);
					// todes.dump_pdb("test1.pdb");
					// utility_exit_with_message("raost");

					replace_nat_seq(todes,pnat);
					Real rholes0 = core::scoring::packing::compute_dec15_score(todes);

					for ( Size id = 1; id <= dhits.size(); id++ ) {
						Pose sym(todes);
						core::pose::symmetry::make_symmetric_pose(sym);
						//TR << dhits[id].rsd1 << " " << dhits[id].rsd2 << std::endl;
						core::conformation::Residue rtmp1 = sym.residue(dhits[id].rsd1);
						core::conformation::Residue rtmp2 = sym.residue(dhits[id].rsd2);
						for ( int i = 0; i < (int)todes.size()/(int)init.size(); ++i ) {
							sym.replace_residue(dhits[id].rsd1+init.size()*i,cys.residue(1),true);
							sym.replace_residue(dhits[id].rsd2+init.size()*i,cys.residue(1),true);
							sym.set_chi(1,dhits[id].rsd1+init.size()*i,dhits[id].chi11);
							sym.set_chi(1,dhits[id].rsd2+init.size()*i,dhits[id].chi12);
							sym.set_chi(2,dhits[id].rsd1+init.size()*i,dhits[id].chi21);
							sym.set_chi(2,dhits[id].rsd2+init.size()*i,dhits[id].chi22);
						}

						ScoreFunctionOP sf = core::scoring::get_score_function();
						TR << "design w/ dsf" << std::endl;
						design_1comp(sym,sf,todes.size());
						Pose ssym(sym);
						for ( int i = 0; i < (int)todes.size()/(int)init.size(); ++i ) {
							ssym.replace_residue(dhits[id].rsd1+init.size()*i,ala.residue(1),true);
							ssym.replace_residue(dhits[id].rsd2+init.size()*i,ala.residue(1),true);
						}

						Real rholes=0,hsp=0,ssp=0,sht=0;
						core::io::silent::SilentStructOP ss( new core::io::silent::ScoreFileSilentStruct );
						sf->score(ssym);
						ss->fill_struct(ssym,option[out::file::o]+"/"+tag+"_dsf"+lzs(id,3)+".pdb.gz");
						ss->add_energy( "xcbc", xcbc );
						ss->add_energy( "cbc", cbc );
						ss->add_energy( "nres", init.size() );
						ss->add_energy( "rholes", rholes0 - core::scoring::packing::compute_dec15_score(ssym) );
						ss->add_energy( "ddg", ddg(ssym,sf,todes.size(),rholes,hsp,ssp,sht) );
						ss->add_energy( "ddgrh", rholes/todes.size() );
						ss->add_energy( "ddghsp", hsp );
						ss->add_energy( "ddgssp", ssp );
						ss->add_energy( "ddgsht", sht );
						sfd.write_silent_struct( *ss, D+"/disulf.sc" );

						core::pose::symmetry::make_asymmetric_pose(sym);
						alignaxis(sym,axs,Vec(0,0,1),Vec(0,0,0));
						trans_pose(sym,C2);
						sym.dump_pdb(option[out::file::o]+"/"+tag+"_dsf"+lzs(id,3)+".pdb.gz");

					}
					if ( cbc+xcbc >= CONTACT_TH2 ) {
						Pose sym(todes);
						replace_nat_seq(sym,pnat);
						core::pose::symmetry::make_symmetric_pose(sym);
						ScoreFunctionOP sf = core::scoring::get_score_function();
						TR << "design 1 component" << std::endl;
						design_1comp(sym,sf,todes.size());
						{
							Real rholes=0,hsp=0,ssp=0,sht=0;
							core::io::silent::SilentStructOP ss( new core::io::silent::ScoreFileSilentStruct );
							sf->score(sym);
							ss->fill_struct(sym,option[out::file::o]+"/"+tag+"_1comp.pdb.gz");
							ss->add_energy( "xcbc", xcbc );
							ss->add_energy( "cbc", cbc );
							ss->add_energy( "nres", init.size() );
							ss->add_energy( "rholes", rholes0 - core::scoring::packing::compute_dec15_score(sym) );
							ss->add_energy( "ddg", ddg(sym,sf,todes.size(),rholes, hsp, ssp, sht) );
							ss->add_energy( "ddgrh", rholes/todes.size() );
							ss->add_energy( "ddghsp", hsp );
							ss->add_energy( "ddgssp", ssp );
							ss->add_energy( "ddgsht", sht );
							sfd.write_silent_struct( *ss, D+"/1comp.sc" );
						}
						core::pose::symmetry::make_asymmetric_pose(sym);
						alignaxis(sym,axs,Vec(0,0,1),Vec(0,0,0));
						trans_pose(sym,C2);
						sym.dump_pdb(option[out::file::o]+"/"+tag+"_1comp.pdb.gz");

						string tmp = option[symmetry::symmetry_definition]();
						option[symmetry::symmetry_definition].clear();

						ScoreFunctionOP sf2 = core::scoring::get_score_function();
						replace_nat_seq(sym,pnat);
						TR << "design 2 component" << std::endl;
						design_1comp(sym,sf2,todes.size());

						Real ddg2comp = 0;
						{
							trans_pose(sym,-C2);
							alignaxis(sym,Vec(0,0,1),axs,Vec(0,0,0));
							option[symmetry::symmetry_definition](tmp);
							Pose comp1,comp2;
							for ( Size ir = 1; ir <= todes.size(); ++ir ) {
								if ( sym.residue(ir                  ).is_lower_terminus() ) comp1.append_residue_by_jump(sym.residue(ir                  ),1);
								else                                                        comp1.append_residue_by_bond(sym.residue(ir                  )  );
								if ( sym.residue(ir+todes.size()).is_lower_terminus() ) comp2.append_residue_by_jump(sym.residue(ir+todes.size()),1);
								else                                                        comp2.append_residue_by_bond(sym.residue(ir+todes.size())  );
							}
							core::pose::symmetry::make_symmetric_pose(comp1);
							core::pose::symmetry::make_symmetric_pose(comp2);
							core::pose::symmetry::make_asymmetric_pose(comp1);
							core::pose::symmetry::make_asymmetric_pose(comp2);

							// sym.dump_pdb("sym.pdb");
							// comp1.dump_pdb("comp1.pdb");
							// comp2.dump_pdb("comp2.pdb");
							vector1<bool> iface;
							repack_iface(sym  ,sf2,todes.size(),iface); iface.clear();
							repack_iface(comp1,sf2,todes.size(),iface); iface.clear();
							repack_iface(comp2,sf2,todes.size(),iface); iface.clear();
							ddg2comp = max(sf2->score(sym)-sf2->score(comp1), sf2->score(sym)-sf2->score(comp2));
							// TR << sf2->score(sym) << " " << sf2->score(comp1) << " " << sf2->score(comp2) << " " << ddg2comp << std::endl;
							// utility_exit_with_message("roistdn");
							option[symmetry::symmetry_definition].clear();
							alignaxis(sym,axs,Vec(0,0,1),Vec(0,0,0));
							trans_pose(sym,C2);
						}

						{
							Real rholes=0,hsp=0,ssp=0,sht=0;
							core::io::silent::SilentStructOP ss( new core::io::silent::ScoreFileSilentStruct );
							sf2 ->score(sym);
							ss->fill_struct(sym,option[out::file::o]+"/"+tag+"_2comp.pdb.gz");
							ss->add_energy( "xcbc", xcbc );
							ss->add_energy( "cbc", cbc );
							ss->add_energy( "nres", init.size() );
							ss->add_energy( "rholes", rholes0 - core::scoring::packing::compute_dec15_score(sym) );
							ss->add_energy( "ddg", ddg(sym,sf2,todes.size(),rholes, hsp, ssp, sht) );
							ss->add_energy( "ddgrh", rholes/todes.size() );
							ss->add_energy( "ddghsp", hsp );
							ss->add_energy( "ddgssp", ssp );
							ss->add_energy( "ddgsht", sht );
							ss->add_energy( "ddg2comp", ddg2comp );
							sfd.write_silent_struct( *ss, D+"/2comp.sc" );
						}
						sym.dump_pdb(option[out::file::o]+"/"+tag+"_2comp.pdb.gz");
						option[symmetry::symmetry_definition](tmp);

					}

				}
			}
		}
	}

	return hits;
}


int main (int argc, char *argv[]) {

	try {

		devel::init(argc,argv);
		using namespace basic::options::OptionKeys;

		string D = option[out::file::o]();
		//if( file_exists("/scratch/USERS/sheffler/xtal/"+D) ) {
		//D = "/scratch/USERS/sheffler/xtal/"+D;
		//}
		option[out::file::o](D);
		std::cout << "output to " << option[out::file::o]() << std::endl;


		for ( Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn ) {
			string fn = option[in::file::s]()[ifn];
			Pose pnat;
			TR << "searching " << fn << std::endl;
			core::import_pose::pose_from_file(pnat,fn, core::import_pose::PDB_file);
			if ( pnat.size() > 150 ) continue;
			Size cyscnt = 0;
			for ( Size ir = 2; ir <= pnat.size()-1; ++ir ) {
				//if(!pnat.residue(ir).is_protein()) goto cont1;
				if ( pnat.residue(ir).is_lower_terminus() ) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
				if ( pnat.residue(ir).is_upper_terminus() ) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
				if ( pnat.residue(ir).name3()=="CYS" ) { if ( ++cyscnt > 3 ) goto cont1; }
			} goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
			//continue;
			Pose pala(pnat);
			vector1<Hit> hts = dock(pala,fn);
			//design_hits(pnat,fn,hts.first,hts.second);
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

