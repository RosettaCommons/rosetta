// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/atomic_depth/AtomicDepth.cc
/// @brief  Calculates depth of atoms from the edge of the Sasa surface
/// @author Dong Xu and Yang Zhang
/// @author Brian Coventry (port to Rosetta)

// This class is a port of "EDTSurf: Quick and accurate construction of macromolecular surfaces" by Dong Xu and Yang Zhang
//   https://zhanglab.ccmb.med.umich.edu/EDTSurf/
// References:
//   D. Xu, Y. Zhang (2009) Generating Triangulated Macromolecular Surfaces by Euclidean Distance Transform. PLoS ONE 4(12): e8140.
//   D. Xu, H. Li, Y. Zhang (2013) Protein Depth Calculation and the Use for Improving Accuracy of Protein Fold Recognition.
//         Journal of Computational Biology 20(10):805-816.

// The following is the license originally packaged with EDTSurf
//
/*////////////////////////////////////////////////////////////////
Permission to use, copy, modify, and distribute this program for
any purpose, with or without fee, is hereby granted, provided that
the notices on the head, the reference information, and this
copyright notice appear in all copies or substantial portions of
the Software. It is provided "as is" without express or implied
warranty.
*////////////////////////////////////////////////////////////////

// EDTSurf has been modified here to remove everything besides what is needed for the atomic depth calculation


#include <core/scoring/atomic_depth/AtomicDepth.hh>

#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <basic/Tracer.hh>

#include <utility/io/ozstream.hh>

#include <boost/format.hpp>


static basic::Tracer TR("core.scoring.atomic_depth.AtomicDepth");

namespace core {
namespace scoring {
namespace atomic_depth {

/// @brief Constructor for AtomicDepth object
/// @detail If poly_leu_depth is true, mutate pose to poly_leu before depth calculations.
/// @detail   The intent is to make this aa independent to avoid surface lysines from adding noise.
AtomicDepth::AtomicDepth( pose::Pose const & pose, Real probe_radius /*=1.4*/, bool poly_leu_depth /*=false*/, Real resolution_scale/*=4.00*/ )
{
	if ( pose.size() == 0 ) {
		utility_exit_with_message("core.scoring.atomic_depth.AtomicDepth passed empty pose!");
	}
	boxlength_=128;
	flagradius_=false;
	scalefactor_=1;
	proberadius_=probe_radius;
	widxz_.resize(13);
	depty_.resize(13);
	// for(i=0;i<13;i++)
	//  depty_[i]=NULL;
	vp_.clear();
	pheight_=0;
	pwidth_=0;
	plength_=0;

	fixsf_=resolution_scale;

	pose::Pose poly_leu_pose;
	if ( poly_leu_depth ) {
		poly_leu_pose = pose;
		conformation::Residue leu = *conformation::get_residue_from_name1( 'L' );
		leu.set_chi(1, 240); // This leu rotamer fills the space where phe likes to sit
		leu.set_chi(2, 120);
		for ( Size i = 1; i <= poly_leu_pose.size(); i++ ) poly_leu_pose.replace_residue( i, leu, true );
	}
	pose::Pose const & use_pose = poly_leu_depth ? poly_leu_pose : pose;

	initpara( use_pose );
	TR << boost::str(boost::format("actual boxlength %3d, box[%3d*%3d*%3d], scale factor %6.3f\n")%
		boxlength_%plength_%pwidth_%pheight_%scalefactor_) << std::endl;
	fillvoxels( use_pose );
	buildboundary();
	fastdistancemap(1);
}

utility::vector1<Real> AtomicDepth::calcdepth( utility::vector1<conformation::Atom> const & atoms, chemical::AtomTypeSet const & type_set ) const
{
	int ox,oy,oz;
	point3d cp;
	Real radius;
	utility::vector1<Real> depval( atoms.size() );

	for ( Size i = 1; i <= atoms.size(); i++ ) {
		conformation::Atom const & atom = atoms[i];
		cp.x=atom.xyz().x()+ptran_.x;
		cp.y=atom.xyz().y()+ptran_.y;
		cp.z=atom.xyz().z()+ptran_.z;
		cp.x*=scalefactor_;
		cp.y*=scalefactor_;
		cp.z*=scalefactor_;
		ox=int(cp.x+0.5);
		oy=int(cp.y+0.5);
		oz=int(cp.z+0.5);
		if ( ox >= 0 && oy >= 0 && oz >= 0 && ox < plength_ && oy < pwidth_ && oz < pheight_ ) {
			depval[i]=vp_[ox][oy][oz].distance/scalefactor_-proberadius_;
			radius = type_set[atom.type()].lj_radius();
			if ( depval[i] < radius ) depval[i]=radius;
			// depval[i]+=proberadius_;                 // criticial change for Rosetta
			//                                          //  removing final proberadius addition so that atomic depth does not include
			//                                          //   radius of sphere
		} else {
			depval[i] = type_set[atom.type()].lj_radius();
		}
	}

	return depval;
}

void AtomicDepth::visualize_at_depth( Real depth, std::string const & fname, Real fraction/*=1.0*/ ) const {
	int anum=1, rnum=1, i, j, k;

	utility::io::ozstream out( fname );

	Size dump_every = Size( 1 / fraction );
	Size count = 0;

	for ( i=0; i<plength_; i++ ) {
		for ( j=0; j<pwidth_; j++ ) {
			for ( k=0; k<pheight_; k++ ) {
				Real this_depth = vp_[i][j][k].distance/scalefactor_-proberadius_;
				if ( this_depth < depth ) continue;
				count += 1;
				if ( count % dump_every != 0 ) continue;

				Real x = i / scalefactor_ - ptran_.x;
				Real y = j / scalefactor_ - ptran_.y;
				Real z = k / scalefactor_ - ptran_.z;

				char buf[128];
				snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
					"HETATM",
					anum++,
					"BURR",
					"BUR",
					'B',
					rnum++,
					x,y,z,
					1.0,
					1.0,
					"B"
				);
				out << buf;

				rnum %= 10000;
				anum %= 100000;
			}
		}
	}

	out.close();
}

void AtomicDepth::boundbox( pose::Pose const & pose,
	point3d & minp,point3d & maxp)
{
	minp.x=100000;minp.y=100000;minp.z=100000;
	maxp.x=-100000;maxp.y=-100000;maxp.z=-100000;

	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		conformation::Residue const & res = pose.residue(seqpos);
		for ( Size atno = 1; atno <= res.nheavyatoms(); atno++ ) {
			if ( res.atom_type(atno).is_virtual() ) continue;
			numeric::xyzVector<Real> xyz = res.xyz(atno);
			if ( xyz.x()<minp.x ) {
				minp.x=xyz.x();
			}
			if ( xyz.y()<minp.y ) {
				minp.y=xyz.y();
			}
			if ( xyz.z()<minp.z ) {
				minp.z=xyz.z();
			}
			if ( xyz.x()>maxp.x ) {
				maxp.x=xyz.x();
			}
			if ( xyz.y()>maxp.y ) {
				maxp.y=xyz.y();
			}
			if ( xyz.z()>maxp.z ) {
				maxp.z=xyz.z();
			}
		}
	}

}
void AtomicDepth::buildboundary()
{
	int i,j,k;
	int ii;
	bool flagbound;
	for ( i=0; i<plength_; i++ ) {
		for ( j=0; j<pheight_; j++ ) {
			for ( k=0; k<pwidth_; k++ ) {
				if ( vp_[i][k][j].inout ) {
					//6 neighbors
					//                  if(( k-1>-1 && !vp_[i][k-1][j].inout) || ( k+1<pwidth_ &&!vp_[i][k+1][j].inout)
					//                  || ( j-1>-1 && !vp_[i][k][j-1].inout) || ( j+1<pheight_ &&!vp_[i][k][j+1].inout)
					//                  || ( i-1>-1 && !vp_[i-1][k][j].inout) || ( i+1<plength_ &&!vp_[i+1][k][j].inout))
					//                      vp_[i][k][j].isbound=true;
					//  /*
					//26 neighbors
					flagbound=false;
					ii=0;
					while ( !flagbound && ii<26 )
							{
						if ( i+nb[ii][0]>-1 && i+nb[ii][0]<plength_
								&& k+nb[ii][1]>-1 && k+nb[ii][1]<pwidth_
								&& j+nb[ii][2]>-1 && j+nb[ii][2]<pheight_
								&& !vp_[i+nb[ii][0]][k+nb[ii][1]][j+nb[ii][2]].inout ) {
							vp_[i][k][j].isbound=true;
							flagbound=true;
						} else ii++;
					}
					//      */
				}
			}

		}
	}
}

void AtomicDepth::boundingatom( pose::Pose const & pose )
{

	chemical::AtomTypeSet const & type_set = pose.residue(1).type().atom_type_set();

	// first figure out which atoms are being used
	std::vector<bool> atom_list( type_set.n_atomtypes()+1, false );
	// std::vector<int> radii( type_set.n_atomtypes()+1, 0 );

	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		conformation::Residue const & res = pose.residue(seqpos);

		for ( Size atno = 1; atno <= res.nheavyatoms(); atno++ ) {
			if ( res.atom_type(atno).is_virtual() ) continue;
			atom_list.at(res.atom_type_index(atno)) = true;

		}
	}



	widxz_.clear();
	depty_.clear();

	widxz_.resize( atom_list.size() );
	depty_.resize( atom_list.size() );


	for ( Size atype = 1; atype < atom_list.size(); atype++ ) {
		depty_[atype].clear();
		if ( ! atom_list[atype] ) continue;

		Real tradius = (type_set[ atype ].lj_radius() + proberadius_)*scalefactor_+0.5;
		// Real tradius = (radii[ atype ] + proberadius_)*scalefactor_+0.5;
		Real sradius = tradius*tradius;
		widxz_[atype] = int(tradius)+1;

		depty_[atype].resize(widxz_[atype]*widxz_[atype]);

		Size indx = 0;
		for ( int j = 0; j < widxz_[atype]; j++ ) {
			for ( int k = 0; k < widxz_[atype]; k++ ) {
				Real txz = j*j+k*k;
				if ( txz > sradius ) {
					depty_[atype][indx]=-1;
				} else {
					Real tdept = sqrt(sradius-txz);
					depty_[atype][indx] = int(tdept);
				}
				indx++;
			}
		}
	}

}


void AtomicDepth::initpara(pose::Pose const & pose)
{
	double fmargin=2.5;
	vp_.clear();
	boundbox(pose,pmin_,pmax_);


	pmin_.x-=proberadius_+fmargin;
	pmin_.y-=proberadius_+fmargin;
	pmin_.z-=proberadius_+fmargin;
	pmax_.x+=proberadius_+fmargin;
	pmax_.y+=proberadius_+fmargin;
	pmax_.z+=proberadius_+fmargin;


	ptran_.x=-pmin_.x;
	ptran_.y=-pmin_.y;
	ptran_.z=-pmin_.z;
	scalefactor_=pmax_.x-pmin_.x;
	if ( (pmax_.y-pmin_.y)>scalefactor_ ) {
		scalefactor_=pmax_.y-pmin_.y;
	}
	if ( (pmax_.z-pmin_.z)>scalefactor_ ) {
		scalefactor_=pmax_.z-pmin_.z;
	}
	scalefactor_=(boxlength_-1.0)/double(scalefactor_);
	///////////////////////////add this automatically first fix sf then fix boxlength_
	//  /*
	boxlength_=int(boxlength_*fixsf_/scalefactor_);
	scalefactor_=fixsf_;
	double threshbox=300;
	if ( boxlength_>threshbox ) {
		double sfthresh=threshbox/double(boxlength_);
		boxlength_=int(threshbox);
		scalefactor_=scalefactor_*sfthresh;
	}
	//  */

	plength_=int(ceil(scalefactor_*(pmax_.x-pmin_.x))+1);
	pwidth_=int(ceil(scalefactor_*(pmax_.y-pmin_.y))+1);
	pheight_=int(ceil(scalefactor_*(pmax_.z-pmin_.z))+1);
	if ( plength_>boxlength_ ) {
		plength_=boxlength_;
	}
	if ( pwidth_>boxlength_ ) {
		pwidth_=boxlength_;
	}
	if ( pheight_>boxlength_ ) {
		pheight_=boxlength_;
	}

	boundingatom( pose );
	cutradis_=proberadius_*scalefactor_;
}
void AtomicDepth::fillatom( conformation::Atom const & atom )
{
	int cx,cy,cz;
	point3d cp;
	cp.x=atom.xyz().x()+ptran_.x;
	cp.y=atom.xyz().y()+ptran_.y;
	cp.z=atom.xyz().z()+ptran_.z;
	cp.x*=scalefactor_;
	cp.y*=scalefactor_;
	cp.z*=scalefactor_;
	cx=int(cp.x+0.5);
	cy=int(cp.y+0.5);
	cz=int(cp.z+0.5);
	int at=atom.type();
	int i,j,k;
	int ii,jj,kk;
	int mi,mj,mk;
	int si,sj,sk;
	int nind=0;
	for ( i=0; i<widxz_[at]; i++ ) {
		for ( j=0; j<widxz_[at]; j++ ) {
			if ( depty_[at][nind]!=-1 ) {

				for ( ii=-1; ii<2; ii++ ) {
					for ( jj=-1; jj<2; jj++ ) {
						for ( kk=-1; kk<2; kk++ ) {
							if ( ii!=0 && jj!=0 && kk!=0 ) {
								mi=ii*i;
								mk=kk*j;
								for ( k=0; k<=depty_[at][nind]; k++ ) {
									mj=k*jj;
									si=cx+mi;
									sj=cy+mj;
									sk=cz+mk;
									if ( si<0 || sj<0 || sk<0 || si>=plength_ || sj>=pwidth_ || sk>=pheight_ ) {
										continue;
									}

									if ( vp_[si][sj][sk].inout==false ) {
										vp_[si][sj][sk].inout=true;
									}
								}//k

							}//if
						}//kk
					}//jj
				}//ii


			}//if
			nind++;
		}//j
	}//i
}
void AtomicDepth::fill_vp() {
	int i,j;
	vp_.resize( plength_ );
	for ( i=0; i<plength_; i++ ) {
		vp_[i].resize(pwidth_);
	}
	for ( i=0; i<plength_; i++ ) {
		for ( j=0; j<pwidth_; j++ ) {
			vp_[i][j].resize(pheight_);
		}
	}
}

//sas use inout
void AtomicDepth::fillvoxels( pose::Pose const & pose )
{

	int i,j,k;
	if ( vp_.size() == 0 ) {
		fill_vp();
	}

	for ( i=0; i<plength_; i++ ) {
		for ( j=0; j<pwidth_; j++ ) {
			for ( k=0; k<pheight_; k++ ) {
				vp_[i][j][k].inout=false;
				vp_[i][j][k].isdone=false;
				vp_[i][j][k].isbound=false;
				vp_[i][j][k].distance=-1;
			}
		}
	}

	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		conformation::Residue const & res = pose.residue(seqpos);

		for ( Size atno = 1; atno <= res.nheavyatoms(); atno++ ) {
			if ( res.atom_type(atno).is_virtual() ) continue;

			fillatom( res.atom(atno) );

		}
	}

	for ( i=0; i<plength_; i++ ) {
		for ( j=0; j<pwidth_; j++ ) {
			for ( k=0; k<pheight_; k++ ) {
				if ( vp_[i][j][k].inout ) {
					vp_[i][j][k].isdone=true;
				}
			}
		}
	}
}

void AtomicDepth::fastdistancemap(int type)
{
	int i,j,k;
	int positin,positout,eliminate;
	totalsurfacevox_=0;
	totalinnervox_=0;

	std::vector<std::vector<std::vector<voxel2> > > boundpoint(plength_);
	for ( i=0; i<plength_; i++ ) {
		boundpoint[i].resize(pwidth_);
	}
	for ( i=0; i<plength_; i++ ) {
		for ( j=0; j<pwidth_; j++ ) {
			boundpoint[i][j].resize(pheight_);
		}
	}
	for ( i=0; i<plength_; i++ ) {
		for ( j=0; j<pwidth_; j++ ) {
			for ( k=0; k<pheight_; k++ ) {
				vp_[i][j][k].isdone=false;
				if ( vp_[i][j][k].inout ) {
					if ( vp_[i][j][k].isbound ) {
						totalsurfacevox_++;
						boundpoint[i][j][k].ix=i;
						boundpoint[i][j][k].iy=j;
						boundpoint[i][j][k].iz=k;
						vp_[i][j][k].distance=0;
						vp_[i][j][k].isdone=true;
					} else {
						totalinnervox_++;
					}
				}
			}
		}
	}
	int allocin=int(1.2*totalsurfacevox_);
	int allocout=int(1.2*totalsurfacevox_);
	if ( allocin>totalinnervox_ ) {
		allocin=totalinnervox_;
	}
	if ( allocin<totalsurfacevox_ ) {
		allocin=totalsurfacevox_;
	}
	if ( allocout>totalinnervox_ ) {
		allocout=totalinnervox_;
	}

	inarray_=std::make_shared<std::vector<voxel2> >( allocin );
	outarray_=std::make_shared<std::vector<voxel2> >( allocout );
	positin=0;positout=0;

	for ( i=0; i<plength_; i++ ) {
		for ( j=0; j<pwidth_; j++ ) {
			for ( k=0; k<pheight_; k++ ) {
				if ( vp_[i][j][k].isbound ) {
					(*inarray_)[positin].ix=i;
					(*inarray_)[positin].iy=j;
					(*inarray_)[positin].iz=k;
					positin++;
					vp_[i][j][k].isbound=false;//as flag of outarray_
				}
			}
		}
	}
	///////////////////////////////////////////////////
	if ( type==0 ) { //do part
		do {
			fastoneshell(positin, allocout, boundpoint, positout,eliminate);

			positin=0;
			for ( i=0; i<positout; i++ ) {
				vp_[(*outarray_)[i].ix][(*outarray_)[i].iy][(*outarray_)[i].iz].isbound=false;
				if ( vp_[(*outarray_)[i].ix][(*outarray_)[i].iy][(*outarray_)[i].iz].distance<=1.02*cutradis_ ) {
					(*inarray_)[positin].ix=(*outarray_)[i].ix;
					(*inarray_)[positin].iy=(*outarray_)[i].iy;
					(*inarray_)[positin].iz=(*outarray_)[i].iz;
					positin++;
				}
				if ( positin>=allocin ) {
					allocin*=2;
					if ( allocin>totalinnervox_ ) allocin=totalinnervox_;

					inarray_->resize(allocin);
				}
			}
		}
		while(positin!=0);
	} else if ( type==1 ) { //do all
		std::shared_ptr<std::vector<voxel2> > tpoint;
		do {

			fastoneshell( positin, allocout, boundpoint, positout,eliminate);//inarray_, outarray_,

			tpoint=inarray_;
			inarray_=outarray_;
			outarray_=tpoint;
			positin=positout;
			int alloctmp;
			alloctmp=allocin;
			allocin=allocout;
			allocout=alloctmp;
			for ( i=0; i<positin; i++ ) {
				vp_[(*inarray_)[i].ix][(*inarray_)[i].iy][(*inarray_)[i].iz].isbound=false;
			}


		}
		while(positout!=0);
	}

	inarray_ = nullptr;
	outarray_ = nullptr;

	double cutsf=scalefactor_-0.5;
	if ( cutsf<0 ) cutsf=0;
	//   cutsf=100000000;
	for ( i=0; i<plength_; i++ ) {
		for ( j=0; j<pwidth_; j++ ) {
			for ( k=0; k<pheight_; k++ ) {
				vp_[i][j][k].isbound=false;
				//ses solid
				if ( vp_[i][j][k].inout ) {
					if ( !vp_[i][j][k].isdone
							|| (vp_[i][j][k].isdone && vp_[i][j][k].distance>=cutradis_-0.50/(0.1+cutsf))//0.33  0.75/scalefactor_
							) {
						vp_[i][j][k].isbound=true;
					}
				}
			}
		}
	}
}

void AtomicDepth::fastoneshell(int innum,int & allocout,std::vector<std::vector<std::vector<voxel2> > > & boundpoint, int & outnum, int & elimi)
{
	int i, number,positout;
	int tx,ty,tz;
	int dx,dy,dz;
	int eliminate=0;
	float squre;
	positout=0;
	number=innum;
	if ( number==0 ) return;
	//new code
	int j,q;
	voxel tnv;

	for ( q=0; q<3; q++ ) {

		for ( i=0; i<number; i++ ) {
			if ( positout>=allocout-fast_one_shell_cut[q] ) {
				allocout=int(1.2*allocout);
				if ( allocout>totalinnervox_ ) allocout=totalinnervox_;
				// outarray_=(voxel2 *)realloc(outarray_,(*allocout)*sizeof(voxel2));
				outarray_->resize(allocout);
			}
			tx=(*inarray_)[i].ix;
			ty=(*inarray_)[i].iy;
			tz=(*inarray_)[i].iz;
			for ( j=fast_one_shell_lowb[q]; j<fast_one_shell_highb[q]; j++ ) {
				tnv.ix=tx+nb[j][0];
				tnv.iy=ty+nb[j][1];
				tnv.iz=tz+nb[j][2];
				if ( tnv.ix<plength_ && tnv.ix>-1 &&
						tnv.iy<pwidth_ && tnv.iy>-1 &&
						tnv.iz<pheight_ && tnv.iz>-1 &&
						vp_[tnv.ix][tnv.iy][tnv.iz].inout &&
						!vp_[tnv.ix][tnv.iy][tnv.iz].isdone ) {
					boundpoint[tnv.ix][tnv.iy][tz+nb[j][2]].ix=boundpoint[tx][ty][tz].ix;
					boundpoint[tnv.ix][tnv.iy][tz+nb[j][2]].iy=boundpoint[tx][ty][tz].iy;
					boundpoint[tnv.ix][tnv.iy][tz+nb[j][2]].iz=boundpoint[tx][ty][tz].iz;
					dx=tnv.ix-boundpoint[tx][ty][tz].ix;
					dy=tnv.iy-boundpoint[tx][ty][tz].iy;
					dz=tnv.iz-boundpoint[tx][ty][tz].iz;
					squre=float(dx*dx+dy*dy+dz*dz);
					vp_[tnv.ix][tnv.iy][tnv.iz].distance=float(sqrt(squre));
					vp_[tnv.ix][tnv.iy][tnv.iz].isdone=true;
					vp_[tnv.ix][tnv.iy][tnv.iz].isbound=true;
					(*outarray_)[positout].ix=tnv.ix;
					(*outarray_)[positout].iy=tnv.iy;
					(*outarray_)[positout].iz=tnv.iz;
					positout++;eliminate++;
				} else if ( tnv.ix<plength_ && tnv.ix>-1 &&
						tnv.iy<pwidth_ && tnv.iy>-1 &&
						tnv.iz<pheight_ && tnv.iz>-1 &&
						vp_[tnv.ix][tnv.iy][tnv.iz].inout &&
						vp_[tnv.ix][tnv.iy][tnv.iz].isdone ) {

					dx=tnv.ix-boundpoint[tx][ty][tz].ix;
					dy=tnv.iy-boundpoint[tx][ty][tz].iy;
					dz=tnv.iz-boundpoint[tx][ty][tz].iz;
					squre=float(dx*dx+dy*dy+dz*dz);
					squre=float(sqrt(squre));
					if ( squre<vp_[tnv.ix][tnv.iy][tnv.iz].distance ) {
						boundpoint[tnv.ix][tnv.iy][tnv.iz].ix=boundpoint[tx][ty][tz].ix;
						boundpoint[tnv.ix][tnv.iy][tnv.iz].iy=boundpoint[tx][ty][tz].iy;
						boundpoint[tnv.ix][tnv.iy][tnv.iz].iz=boundpoint[tx][ty][tz].iz;
						vp_[tnv.ix][tnv.iy][tnv.iz].distance=float(squre);
						if ( !vp_[tnv.ix][tnv.iy][tnv.iz].isbound ) {
							vp_[tnv.ix][tnv.iy][tnv.iz].isbound=true;
							(*outarray_)[positout].ix=tnv.ix;
							(*outarray_)[positout].iy=tnv.iy;
							(*outarray_)[positout].iz=tnv.iz;
							positout++;
						}
					}

				}
			}
		}
	}

	outnum=positout;
	elimi=eliminate;

}



} // namespace atomic_depth
} // namespace scoring
} // namespace core

