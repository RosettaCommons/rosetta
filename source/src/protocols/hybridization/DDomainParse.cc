// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief A reimplementation of DDomain domain parsing algorithm, Protein Sci.16,947--955 (2007)
/// @author Frank DiMaio

#include <protocols/hybridization/DDomainParse.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>
#include <core/types.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/pose/PDBInfo.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/task/PackerTask.hh>

#include <list>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

// parameters
//const core::Real DDomainParse::pcut_ = 0.81;
//const core::Real DDomainParse::hcut_ = 0.18;
//const core::Real DDomainParse::length_ = 38;

void
DDomainParse::pulldomain( int isize0, int ist, int ilast, utility::vector1< core::Real > &resect, int &idip) {
	idip = 1;

	utility::vector1< int> maxe1(mseq_), maxe2(mseq_);
	core::Real /*tmp0 = -100.0,*/ tmp10 = 100.0, tmp20 = 100.0;
	//int ikp = 1;
	//int i0 = 0, idip0 = 1, kk1 = 0, kk2 = 0;

	for ( int i=isize0; i<=ilast-isize0; ++i ) {
		core::Real tmp = resect[i];

		maxe1[i] = i;
		core::Real tmp1 = tmp;
		for ( int j=ist; j<=i; ++j ) {
			core::Real tmp2 = resect[j];
			if ( tmp2 > tmp1 ) {
				tmp1 = tmp2;
				maxe1[i] = j;
			}
		}
		maxe2[i] = i;
		tmp1 = tmp;
		for ( int j=i; j<=ilast-1; ++j ) {
			core::Real tmp2 = resect[j];
			if ( tmp2 > tmp1 ) {
				tmp1 = tmp2;
				maxe2[i] = j;
			}
		}

		core::Size k1 = 0;
		for ( int k=maxe1[i]; k<=i; ++k ) {
			core::Real de1 = resect[k] - tmp;
			if ( de1 > hcut_ ) {
				k1++;
			} else {
				break;
			}
		}
		for ( int k=maxe1[i]; k>=ist; --k ) {
			core::Real de1 = resect[k] - tmp;
			if ( de1 > hcut_ ) {
				k1++;
			} else {
				break;
			}
		}

		core::Size k2 = 0;
		for ( int k=maxe2[i]; k<=ilast; ++k ) {
			core::Real de1 = resect[k] - tmp;
			if ( de1 > hcut_ ) {
				k2++;
			} else {
				break;
			}

		}
		for ( int k=maxe2[i]; k>=i; --k ) {
			core::Real de1 = resect[k] - tmp;
			if ( de1 > hcut_ ) {
				k2++;
			} else {
				break;
			}
		}


		core::Real de1 = resect[maxe1[i]] - tmp;
		core::Real de2 = resect[maxe2[i]] - tmp;

		if ( tmp < tmp20 && de1 > hcut_ && de2 > hcut_ ) {
			//ikp = i;  // set but never used ~Labonte
			tmp20 = tmp;
		}
		if ( (k1 > length_ && k2 > length_) ) {
			if ( tmp < tmp10 ) {
				tmp10 = tmp;
				idip = i;
				//kk1 = k1;  // set but never used ~Labonte
				//kk2 = k2;  // set but never used ~Labonte
			}
		}
	}


	if ( resect[idip] > pcut_ ) {
		idip = ist;
	}

	int i1 = idip;
	int i2 = ilast - idip;
	if ( i1 > isize0 && i2 > isize0 ) {
		// do nothing
	} else {
		idip = ist;
	}
}


void
DDomainParse::findpos(
	int, int ndom, int id,
	utility::vector1<int> &ipdom, int &ist, int &ilast)
{
	if ( ndom > 1 ) {
		if ( id == 1 ) {
			ist = 1;
			ilast = ipdom[id];
		} else if ( id == ndom ) {
			ist = ipdom[id - 1] + 1;
			ilast = nseq_;
		} else {
			ist = ipdom[id - 1] + 1;
			ilast = ipdom[id];
		}
	}
}

void
DDomainParse::findpos2(
	int ist, int ilast, int imid, int i,
	int &ista, int &ilasta)
{
	if ( i == 1 ) {
		ista = ist;
		ilasta = imid;
	} else {
		ista = imid + 1;
		ilasta = ilast;
	}
}

void
DDomainParse::small_big(
	int const& mdom,
	int const& ipdd,
	utility::vector1<int> &ipdom)
{
	int /*imin = 0, */ip = 0;
	utility::vector1<int> ifg_min(mdom, 0);
	utility::vector1<int> itmp(mdom, 0);

	for ( int id=1; id<=ipdd; ++id ) {
		int imin = 100000;
		for ( int i=1; i<=ipdd; ++i ) {
			if ( ifg_min[i] == 0 ) {
				if ( ipdom[i] < imin ) {
					imin = ipdom[i];
					ip = i;
				}
			}
		}
		ifg_min[ip] = 1;
		itmp[id] = imin;
	}
	for ( int i=1; i<=ipdd; ++i ) {
		ipdom[i] = itmp[i];
	}
}


void
DDomainParse::distance(
	int ist, int ilast,
	utility::vector1< numeric::xyzVector<core::Real> > const &xs,
	utility::vector1< utility::vector1< core::Real > > &dij)
{
	dij.resize(mseq_, utility::vector1< core::Real >(mseq_,0));

	for ( int i=ist; i<=ilast-1; ++i ) {
		for ( int j=i+1; j<=ilast; ++j ) {
			dij[i][j] = dij[j][i] = xs[i].distance( xs[j] );
		}
	}
}

void
DDomainParse::ddomain_pot(
	core::Real pw, int ist, int ilast,
	utility::vector1< utility::vector1< core::Real > > const &dij,
	utility::vector1< core::Real > &resect)
{
	resect.resize(mseq_);

	core::Real ave = 0.0, rcut = 6.5, /*ect = 0,*/ dr = 0;
	//int ipt1 = 0, ipt2 = 0;
	for ( int k=ist; k<=ilast - 1; ++k ) {
		core::Real ect = 0.0;
		for ( int i=ist; i<=k; ++i ) {
			for ( int j=k+1; j<=ilast; ++j ) {
				dr = dij[i][j];
				if ( dr < rcut ) {
					ect -= 1.0;
				}
			}
		}
		int ipt1 = k - ist + 1;
		int ipt2 = (ilast - 1) - k + 1;

		ect /= std::pow((ipt1 * ipt2), pw);
		resect[k] = ect;
		ave += ect;
	}
	core::Real tmp = ave / ((ilast - ist + 1) - 1);
	for ( int k=ist; k<=ilast - 1; ++k ) {
		resect[k] /= tmp;
	}
}

utility::vector1< loops::Loops >
DDomainParse::split( core::pose::Pose const &templ ) {
	//runtime_assert( nres >= templ.size() );

	utility::vector1< loops::Loops > retval;

	//std::string fname;
	mseq_ = templ.size();
	while ( mseq_>0 && !templ.residue(mseq_).is_protein() ) mseq_--;
	nseq_ = mseq_;

	if ( mseq_ == 0 ) { return retval; }

	utility::vector1<int> ipdom(mseq_,0);
	utility::vector1< numeric::xyzVector<core::Real> > xs(mseq_, numeric::xyzVector<core::Real>(0,0,0));
	utility::vector1< utility::vector1< core::Real > > dij;
	utility::vector1<core::Real> resect;

	// grab "centroid" from pose
	// use same definition as DDomain
	for ( core::Size i=1; i<=mseq_; ++i ) {
		core::conformation::Residue const &rsd_i = templ.residue(i);
		core::Size nsc=0;
		if ( rsd_i.aa() == core::chemical::aa_gly ) {
			xs[i] = rsd_i.atom(2).xyz();
			nsc++;
		}
		for ( int j=rsd_i.first_sidechain_atom(); j<=(int)rsd_i.nheavyatoms(); ++j ) {
			if ( rsd_i.aa() == core::chemical::aa_pro && j==(int)rsd_i.nheavyatoms() ) continue;
			xs[i] += rsd_i.atom(j).xyz();
			nsc++;
		}
		xs[i] /= nsc;
	}

	const int mdom = 20;  // maximum domains

	core::Real pw = 0.43f;
	int isize0 = 40;
	int ist=1, ilast=nseq_, ndom=1, idom=1, ipdd=0, ifg_dom=0;
	int imid=0, istx=0, ilastx=0, imidx=0;

	bool done=false;

	while ( !done ) {
		done = true;
		for ( int id=1; id <= ndom; ++id ) {
			findpos(mdom, ndom, id, ipdom, ist, ilast);
			distance(ist, ilast, xs, dij);
			ddomain_pot(pw, ist, ilast, dij, resect);
			pulldomain(isize0, ist, ilast, resect, imid);

			if ( imid > ist && imid < ilast ) {
				idom++;
				if ( idom > mdom ) {
					std::cerr << "DDomainParse::split(): idom > mdom [1]";
					ifg_dom = 0;
					break;
				}
				ipdd++;
				ipdom[ipdd] = imid;
				ifg_dom = 1;

				for ( int i=1; i<=2; ++i ) {
					findpos2(ist, ilast, imid, i, istx, ilastx);
					if ( ilastx - istx >= isize0 ) {
						distance(istx, ilastx, xs, dij);
						ddomain_pot(pw, istx, ilastx, dij, resect);
						pulldomain(isize0, istx, ilastx, resect, imidx);
						if ( imidx > istx && imidx < ilastx ) {
							idom++;
							if ( idom > mdom ) {
								std::cerr << "DDomainParse::split(): idom > mdom [2]";
								ifg_dom = 0;
								break;
							}
							ipdd++;
							ipdom[ipdd] = imidx;
							ifg_dom = 1;
						}
					}
				}
			}
		}

		if ( ifg_dom == 1 ) {
			ndom = idom;
			ifg_dom = 0;
			small_big(mdom, ipdd, ipdom);
			done=false;
		}
	}

	for ( int i=1; i<=ndom; ++i ) {
		if ( ndom == 1 ) {
			loops::Loops newloops;
			newloops.add_loop( 1, nseq_ );
			retval.push_back( newloops );
		} else {
			if ( i == 1 ) {
				loops::Loops newloops;
				newloops.add_loop( 1, ipdom[i] );
				retval.push_back( newloops );
			} else if ( i == ndom ) {
				loops::Loops newloops;
				newloops.add_loop( ipdom[i-1]+1, nseq_ );
				retval.push_back( newloops );
			} else {
				loops::Loops newloops;
				newloops.add_loop( ipdom[i-1]+1, ipdom[i] );
				retval.push_back( newloops );
			}
		}
	}

	return retval;
}

} // hybridize
//} // comparative_modeling
} // protocols
