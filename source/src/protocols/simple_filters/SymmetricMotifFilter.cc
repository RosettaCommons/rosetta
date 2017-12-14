// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/SymmetricMotifFilter.cc
/// @brief  position-independent RMS filter evaluating how close a set of interfaces is to symmetric
/// @author Frank DiMaio

#include <protocols/simple_filters/SymmetricMotifFilter.hh>
#include <protocols/simple_filters/SymmetricMotifFilterCreator.hh>
#include <protocols/filters/Filter.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/conformation/Conformation.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/PDBInfo.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>

#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

// tracer
static basic::Tracer TR( "protocols.simple_filters.SymmetricMotifFilter" );

// creator
// XRW TEMP protocols::filters::FilterOP
// XRW TEMP SymmetricMotifFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SymmetricMotifFilter ); }

// XRW TEMP std::string
// XRW TEMP SymmetricMotifFilterCreator::keyname() const { return "SymmetricMotif"; }

// helper functions
void R2quat( numeric::xyzMatrix< core::Real > R, Quat &Q ) {
	core::Real S;
	if ( R.xx() > R.yy() && R.xx() > R.zz() )  {
		S  = sqrt( 1.0 + R.xx() - R.yy() - R.zz() ) * 2;
		Q.x = 0.25 * S;
		Q.y = (R.xy() + R.yx() ) / S;
		Q.z = (R.zx() + R.xz() ) / S;
		Q.w = (R.zy() - R.yz() ) / S;
	} else if ( R.yy() > R.zz() ) {
		S  = sqrt( 1.0 + R.yy() - R.xx() - R.zz() ) * 2;
		Q.x = (R.yx() + R.xy() ) / S;
		Q.y = 0.25 * S;
		Q.z = (R.zy() + R.yz() ) / S;
		Q.w = (R.xz() - R.zx() ) / S;
	} else {
		S  = sqrt( 1.0 + R.zz() - R.xx() - R.yy() ) * 2;
		Q.x = (R.xz() + R.zx() ) / S;
		Q.y = (R.zy() + R.yz() ) / S;
		Q.z = 0.25 * S;
		Q.w = (R.yx() - R.xy()) / S;
	}
}

void quat2R( Quat &Q ,numeric::xyzMatrix< core::Real > R ) {
	core::Real xx = Q.x*Q.x; core::Real xy = Q.x*Q.y; core::Real xz = Q.x*Q.z;
	core::Real xw = Q.x*Q.w; core::Real yy = Q.y*Q.y; core::Real yz = Q.y*Q.z;
	core::Real yw = Q.y*Q.w; core::Real zz = Q.z*Q.z; core::Real zw = Q.z*Q.w;

	R.xx(1-2*(yy+zz)); R.xy(  2*(xy-zw)); R.xz(  2*(xz+yw));
	R.yx(  2*(xy+zw)); R.yy(1-2*(xx+zz)); R.yz(  2*(yz-xw));
	R.zx(  2*(xz-yw)); R.zy(  2*(yz+xw)); R.zz(1-2*(xx+yy));
}

//rms wrapper
core::Real RMSwrapper( utility::vector1 <numeric::xyzVector< core::Real > > chainA,
	utility::vector1 <numeric::xyzVector< core::Real > > chainB,
	numeric::xyzMatrix< core::Real > &R,
	numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT) {
	core::Size motiflen = std::min( chainA.size(), chainB.size() );
	ObjexxFCL::FArray2D< core::Real > final_coords( 3, motiflen );
	ObjexxFCL::FArray2D< core::Real > init_coords( 3, motiflen );
	ObjexxFCL::FArray1D< numeric::Real > ww( motiflen, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );

	numeric::Real ctx;
	preT = postT = numeric::xyzVector< core::Real >(0,0,0);
	for ( int j=1; j<=(int)motiflen; ++j ) {
		for ( int k=0; k<3; ++k ) {
			init_coords(k+1,j)  = chainA[j][k];
			final_coords(k+1,j) = chainB[j][k];
		}
		preT  += chainA[j];
		postT += chainB[j];
	}
	preT /= motiflen;
	postT /= motiflen;
	for ( int j=1; j<=(int)motiflen; ++j ) {
		for ( int k=0; k<3; ++k ) {
			init_coords(k+1,j)  -= preT[k];
			final_coords(k+1,j) -= postT[k];
		}
	}

	float rms;
	numeric::model_quality::findUU( final_coords, init_coords, ww, motiflen, uu, ctx );
	numeric::model_quality::calc_rms_fast( rms, final_coords, init_coords, ww, motiflen, ctx );

	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

	return ((core::Real)rms);
}


// filter definition
SymmetricMotifFilter::SymmetricMotifFilter() :
	protocols::filters::Filter( "SymmetricMotif" )
{
	set_defaults();
}

SymmetricMotifFilter::SymmetricMotifFilter( utility::vector1<core::pose::PoseOP> const & reference_motifs, std::string const & symm_type_in)
: protocols::filters::Filter( "SymmetricMotif" ), symm_type_(symm_type_in)
{
	ref_motifs_ = reference_motifs;
	process_motifs();
	set_defaults();
}

SymmetricMotifFilter::~SymmetricMotifFilter() = default;

void
SymmetricMotifFilter::set_defaults() {
	angle_thresh_ = 5.0;
	trans_thresh_ = 4.0;
	rmsd_thresh_ = 2.0;
	clash_thresh_ = 0;

	angle_wt_ = 1.0;
	trans_wt_ = 5.0;
	rmsd_wt_ = 10.0;
	clash_wt_ = 100.0;
}


protocols::filters::FilterOP
SymmetricMotifFilter::clone() const {
	return protocols::filters::FilterOP( new SymmetricMotifFilter( *this ) );
}

bool
SymmetricMotifFilter::apply( core::pose::Pose const & pose ) const {
	core::Real bestscore;
	std::string hitLocation;
	bool found_motif = compute( pose, bestscore, hitLocation );
	if ( found_motif ) {
		TR.Debug << "Motif hit at " << hitLocation << std::endl;
	}
	return( found_motif );
}

void SymmetricMotifFilter::add_motif( core::pose::PoseOP motif )
{
	ref_motifs_.push_back( motif );
}

void
SymmetricMotifFilter::report( std::ostream & /*out*/, core::pose::Pose const & pose ) const {
	core::Real bestscore;
	std::string hitLocation;
	bool found_motif = compute( pose, bestscore, hitLocation );
	if ( found_motif ) {
		TR << "Motif hit at " << hitLocation << std::endl;
	} else {
		TR << "No motif hits found." << std::endl;
	}
	return;
}

core::Real
SymmetricMotifFilter::report_sm( core::pose::Pose const & /*pose*/ ) const {
	core::Real bestscore = 0.0;  // arbitrarily set to zero to silence compiler warning ~Labonte
	//std::string hitLocation;
	//bool found_motif = compute( pose, bestscore, hitLocation );
	return bestscore;
}


// compute
bool
SymmetricMotifFilter::compute( core::pose::Pose const & pose, core::Real &best_score, std::string &motifhit ) const {
	if ( symm_type_ == "D2" ) {
		return compute_d2( pose, best_score, motifhit );
	}
	return false;  // added to remove compiler warning; I assume false is warranted here? ~ Labonte
}


// compute
bool
SymmetricMotifFilter::compute_d2( core::pose::Pose const & pose, core::Real &best_score, std::string &motifhit ) const {
	using numeric::constants::d::pi;

	// three-stage approach
	//  1) find backbone segments satisfying each motif (fast)
	//  2) [symm dependent!] for each backbone hit pair, measure angle error (fast)
	//  3) for each hit passing these filters, do CA clash check (slow)
	motifhit = "";
	bool foundOne = false;
	best_score = 999;
	bool noforce = (forced_pos_.size() == 0);

	// ca trace of pose
	utility::vector1< numeric::xyzVector< core::Real > > cas_pose;
	core::Size nres = pose.size();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres = symm_info->num_independent_residues();
	}
	for ( Size i=1; i<=nres; ++i ) {
		if ( !pose.residue_type(i).is_protein() ) continue;
		cas_pose.push_back( pose.residue(i).atom(" CA ").xyz() );
		cas_pose.push_back( pose.residue(i).atom(" C  ").xyz() );
		cas_pose.push_back( pose.residue(i).atom(" N  ").xyz() );
		cas_pose.push_back( pose.residue(i).atom(" O  ").xyz() );
	}
	core::Size nres_prot = cas_pose.size() / 4;

	//  1) find backbone segments satisfying each motif
	utility::vector1< utility::vector1< int > > motif_hits1(1), motif_hits2(1);
	utility::vector1< core::Real > motif_rms1, motif_rms2;
	utility::vector1< numeric::xyzMatrix< core::Real > > Rs1, Rs2;
	utility::vector1< numeric::xyzVector< core::Real > > preTs1, preTs2, postTs1, postTs2;

	core::Size segmentCounter = 0;

	for ( int i=1; i<=2; ++i ) {
		utility::vector1< utility::vector1< int > > &motif_hits_i = (i==1 ? motif_hits1 : motif_hits2);
		utility::vector1< core::Real > &motif_rms_i = (i==1 ? motif_rms1 : motif_rms2);
		utility::vector1< numeric::xyzMatrix< core::Real > > &Rs_i = (i==1 ? Rs1 : Rs2);
		utility::vector1< numeric::xyzVector< core::Real > > &preTs_i = (i==1 ? preTs1 : preTs2);
		utility::vector1< numeric::xyzVector< core::Real > > &postTs_i = (i==1 ? postTs1 : postTs2);
		utility::vector1< numeric::xyzVector< core::Real > > const &cas_tgt_i = cas_chainA[i];
		utility::vector1< core::Size > const &motif_cuts_i = motif_cuts[i];

		utility::vector1< utility::vector1< int > > prev_Is; // intermediate hits
		motif_hits_i[1].push_back( -1 );


		// build up multi-segment motifs one segment at a time
		for ( Size j=2; j<=motif_cuts_i.size(); ++j ) {
			prev_Is = motif_hits_i;
			motif_hits_i.clear();
			segmentCounter++;

			for ( Size k=1; k<=prev_Is.size(); ++k ) {
				utility::vector1< numeric::xyzVector< core::Real > > prevCAs;
				utility::vector1< bool > elim(nres_prot, false);

				// load prev hits' CA coords
				for ( Size j_prev = 2; j_prev<j; ++j_prev ) {
					core::Size segstart_x = prev_Is[k][j_prev-1];
					core::Size seglen_x   = motif_cuts_i[j_prev] - motif_cuts_i[j_prev-1] - 1;
					for ( Size l=segstart_x; l<=segstart_x+seglen_x; ++l ) {
						prevCAs.push_back( cas_pose[4*l-3] );
						prevCAs.push_back( cas_pose[4*l-2] );
						prevCAs.push_back( cas_pose[4*l-1] );
						prevCAs.push_back( cas_pose[4*l] );
						elim[l] = true;
					}
				}

				// scan through all possible placements of this motif
				int lstart=1, lstop=nres_prot;
				if ( !noforce && forced_pos_[segmentCounter]!=-1 ) {
					lstart=lstop=forced_pos_[segmentCounter];
				}

				for ( int l=lstart; l<=lstop; ++l ) {
					core::Size segstart_l = l;
					core::Size seglen_l   = motif_cuts_i[j] - motif_cuts_i[j-1] - 1;
					bool overlap=false;
					if ( segstart_l+seglen_l > nres_prot ) continue;
					for ( Size m=segstart_l; m<=segstart_l+seglen_l && !overlap; ++m ) {
						overlap |= elim[m];
					}
					if ( overlap ) continue;

					utility::vector1< numeric::xyzVector< core::Real > > ca_chunk_pose = prevCAs;
					for ( Size m=segstart_l; m<=segstart_l+seglen_l; ++m ) {
						ca_chunk_pose.push_back( cas_pose[4*m-3] );
						ca_chunk_pose.push_back( cas_pose[4*m-2] );
						ca_chunk_pose.push_back( cas_pose[4*m-1] );
						ca_chunk_pose.push_back( cas_pose[4*m  ] );
					}

					// align pose CAs to tgt
					numeric::xyzVector< core::Real > preT, postT;
					numeric::xyzMatrix< core::Real > R;

					// this will only align the subset of CAs present in 'ca_chunk_pose'
					core::Real rms = RMSwrapper( cas_tgt_i, ca_chunk_pose, R, preT, postT);

					if ( rms < rmsd_thresh_ ) {
						utility::vector1< int > newI;
						if ( j != 2 ) newI = prev_Is[k];
						newI.push_back(l);
						motif_hits_i.push_back( newI );

						TR.Debug << "Motif " << i << " hit at ";
						for ( Size z=1; z<=newI.size(); ++z ) TR.Debug << newI[z] << " ";
						TR.Debug << " rms = " << rms << std::endl;

						if ( j == motif_cuts_i.size() ) {   // last segment has been placed
							Rs_i.push_back( R );
							motif_rms_i.push_back( rms );
							preTs_i.push_back( preT );
							postTs_i.push_back( postT );
						}
					}
				}
			}
		}
	}

	//  2) for each backbone hit pair combination, measure angle error
	for ( Size i=1; i<=motif_hits1.size(); ++i ) {
		for ( Size j=1; j<=motif_hits2.size(); ++j ) {
			// a) check for overlap
			utility::vector1< bool > elim(nres_prot, false);
			for ( Size x=1; x<=motif_hits1[i].size(); ++x ) {
				core::Size segstart_x = motif_hits1[i][x];
				core::Size seglen_x = motif_cuts[1][x+1] - motif_cuts[1][x] - 1;
				for ( Size k=segstart_x; k<=segstart_x+seglen_x; ++k ) {
					elim[k] = true;
				}
			}
			bool overlap = false;
			for ( Size x=1; x<=motif_hits2[j].size() && !overlap; ++x ) {
				core::Size segstart_x = motif_hits2[j][x];
				core::Size seglen_x = motif_cuts[2][x+1] - motif_cuts[2][x] - 1;
				for ( Size k=segstart_x; k<=segstart_x+seglen_x && !overlap; ++k ) {
					overlap |= elim[k];
				}
			}
			if ( overlap ) continue;

			// symm axes in global frame
			numeric::xyzVector< core::Real > x1 = numeric::inverse(Rs1[i])*symm_axes[1];
			numeric::xyzVector< core::Real > x2 = numeric::inverse(Rs2[j])*symm_axes[2];

			core::Real angle12 = angle_of( x1, x2 ) * 180/pi;
			core::Real angle_i = std::fabs(90-angle12);
			TR.Debug << "Hit " << i << "," << j << ": angle = " << angle_i << std::endl;
			if ( angle_i > angle_thresh_ ) continue;

			// centers of mass of motifs in all subunits ... somewhat tricky
			numeric::xyzMatrix< core::Real > R2 = (numeric::inverse(Rs1[i])*Rdimers[1])*Rs1[i];
			numeric::xyzMatrix< core::Real > R3 = (numeric::inverse(Rs2[j])*Rdimers[2])*Rs2[j];
			numeric::xyzMatrix< core::Real > R4a = R2*R3;
			numeric::xyzMatrix< core::Real > R4b = R3*R2;

			// motif centers in chain A
			numeric::xyzVector< core::Real > m1a = postTs1[i];
			numeric::xyzVector< core::Real > m2a = postTs2[j];
			numeric::xyzVector< core::Real > m1to2 = m2a-m1a;

			// ... in chain B
			numeric::xyzVector< core::Real > m1b = m1a + (numeric::inverse(Rs1[i])*Rdimers[1])*delta_coms[1];
			numeric::xyzVector< core::Real > m2b = m1b + (R2)*m1to2;

			// ... in chain C
			numeric::xyzVector< core::Real > m2c = m2a + (numeric::inverse(Rs2[j])*Rdimers[2])*delta_coms[2];
			numeric::xyzVector< core::Real > m1c = m2c - (R3)*m1to2;

			// ... in chain D
			numeric::xyzVector< core::Real > m1d = m1c + R4a*numeric::inverse(Rs1[i])*delta_coms[1];
			numeric::xyzVector< core::Real > m2d = m2b + R4b*numeric::inverse(Rs2[j])*delta_coms[2];
			numeric::xyzVector< core::Real > m2d_alt = m1d + R4a*m1to2;

			// translation error is the distance between m2_d and m2_d_alt
			core::Real trans_i = (m2d).distance(m2d_alt);

			TR.Debug << "Hit " << i << "," << j << ": trans = " << trans_i << std::endl;
			if ( trans_i > trans_thresh_ ) continue;

			core::Real rms_i = std::max( motif_rms1[i], motif_rms2[j] );

			// i,j has passed RMS and angle filters

			//  3) for each hit passing these filters, do CA clash check (to do: CB??)
			numeric::xyzMatrix <core::Real> &R1A = Rs1[i];
			numeric::xyzMatrix <core::Real> R1B = (numeric::inverse(Rdimers[1])*Rs1[i]);
			numeric::xyzMatrix <core::Real> &R2A = Rs2[j];
			numeric::xyzMatrix <core::Real> R2B = (numeric::inverse(Rdimers[2])*Rs2[j]);

			numeric::xyzVector< core::Real > com1A = preTs1[i] - postTs1[i];
			numeric::xyzVector< core::Real > com1B = preTs1[i] - postTs1[i] + Rdimers[1]*delta_coms[1];
			numeric::xyzVector< core::Real > com2A = preTs2[j] - postTs2[j];
			numeric::xyzVector< core::Real > com2B = preTs2[j] - postTs2[j] + Rdimers[2]*delta_coms[2];

			bool clashcheck = false;
			Size nclashes = 0;
			int CUTOFF2 = 3*3;
			for ( Size x=1; x<=nres_prot && !clashcheck; ++x ) {
				for ( Size y=1; y<=nres_prot && !clashcheck; ++y ) {
					numeric::xyzVector< core::Real > x_x = R1A*(cas_pose[4*x-3]-postTs1[i]) + com1A;  // just check CA
					numeric::xyzVector< core::Real > x_y = R1B*(cas_pose[4*y-3]-postTs1[i]) + com1B;  // just check CA
					if ( x_x.distance_squared(x_y) < CUTOFF2 ) {
						nclashes++;
						clashcheck = (nclashes>clash_thresh_);
					}
				}
			}
			for ( Size x=1; x<=cas_pose.size() && !clashcheck; ++x ) {
				for ( Size y=1; y<=cas_pose.size() && !clashcheck; ++y ) {
					numeric::xyzVector< core::Real > x_x = R2A*(cas_pose[4*x]-postTs2[j]) + com2A;
					numeric::xyzVector< core::Real > x_y = R2B*(cas_pose[4*y]-postTs2[j]) + com2B;
					if ( x_x.distance_squared(x_y) < CUTOFF2 ) {
						nclashes++;
						clashcheck = (nclashes>clash_thresh_);
					}
				}
			}

			if ( clashcheck ) {
				TR.Debug << "Hit " << i << "," << j << ": clash > " << clash_thresh_ << std::endl;
				continue;
			}

			// passed! report rms, angle and clash
			foundOne = true;
			core::Real score_i = score_d2(rms_i,angle_i,trans_i,nclashes);
			if ( score_i < best_score ) {
				// dump ID to a string for reporting purposes
				std::ostringstream oss;
				oss << "( ";
				for ( Size k=1; k<=motif_hits1[i].size(); ++k ) oss << motif_hits1[i][k] << " ";
				oss << ") ( ";
				for ( Size k=1; k<=motif_hits2[j].size(); ++k ) oss << motif_hits2[j][k] << " ";
				oss << ")";
				motifhit = oss.str();
				best_score = score_i;
			}
		}
	}
	return foundOne;
}


// process_motifs
//fpd currently only D2 is supported!
//fpd at some point more spacegroups will be supported
void
SymmetricMotifFilter::process_motifs() {
	using numeric::constants::d::pi;

	// process motifs
	core::Size nmotifs = ref_motifs_.size();
	cas_chainA.resize(nmotifs);
	cas_chainB.resize(nmotifs);
	motif_cuts.resize(nmotifs);

	// parse motifs
	nsegs_ = 0;
	for ( Size i=1; i<=nmotifs; ++i ) {
		motif_cuts[i].push_back( 0 );
		core::pose::PoseOP motif_i = ref_motifs_[i];
		for ( Size j=1; j<=motif_i->size(); ++j ) {
			if ( motif_i->pdb_info()->chain(j) == 'A' ) {
				cas_chainA[i].push_back( motif_i->residue(j).atom(" CA ").xyz() );
				cas_chainA[i].push_back( motif_i->residue(j).atom(" C  ").xyz() );
				cas_chainA[i].push_back( motif_i->residue(j).atom(" N  ").xyz() );
				cas_chainA[i].push_back( motif_i->residue(j).atom(" O  ").xyz() );

				if ( j>1 && (motif_i->pdb_info()->number(j) != motif_i->pdb_info()->number(j-1)+1) ) {
					motif_cuts[i].push_back( j-1 );
					TR.Debug << i << ": add cut " << j-1 << std::endl;
				}
			} else {
				cas_chainB[i].push_back( motif_i->residue(j).atom(" CA ").xyz() );
				cas_chainB[i].push_back( motif_i->residue(j).atom(" C  ").xyz() );
				cas_chainB[i].push_back( motif_i->residue(j).atom(" N  ").xyz() );
				cas_chainB[i].push_back( motif_i->residue(j).atom(" O  ").xyz() );
			}
		}
		runtime_assert( cas_chainA[i].size() == cas_chainB[i].size() ); // both chains should be same length
		motif_cuts[i].push_back( cas_chainA[i].size()/4 );
		nsegs_ += (motif_cuts[i].size() - 1);
	}

	// get superposition
	Qs.resize(nmotifs);
	Rdimers.resize(nmotifs);
	delta_coms.resize(nmotifs);
	symm_orders.resize(nmotifs);
	symm_axes.resize(nmotifs);
	for ( Size i=1; i<=nmotifs; ++i ) {
		numeric::xyzVector< core::Real > preT(0,0,0), postT(0,0,0);
		core::Real rms = RMSwrapper( cas_chainB[i], cas_chainA[i], Rdimers[i], preT, postT);

		// check if near identity
		core::Real residual =
			std::fabs(Rdimers[i].xx()-1) + std::fabs(Rdimers[i].yy()-1) + std::fabs(Rdimers[i].zz()-1) +
			std::fabs(Rdimers[i].xy()) + std::fabs(Rdimers[i].yx()) + std::fabs(Rdimers[i].zy()) +
			std::fabs(Rdimers[i].xz()) + std::fabs(Rdimers[i].yz()) + std::fabs(Rdimers[i].zx());

		if ( residual < 1e-6 ) {
			utility_exit_with_message( "Chains related by transformation only!" );
		}

		if ( rms > 1.0 ) {
			TR << "RMS = " << rms << std::endl;
			utility_exit_with_message( "Interchain RMS > 1 ... aborting" );
		}

		// make transformation perfectly symmetrical
		R2quat( Rdimers[i], Qs[i]);
		core::Real Wmult = 1;
		if ( Qs[i].w < 0 ) { Qs[i].w = -Qs[i].w; Wmult = -1; }
		core::Real omega = acos( Qs[i].w );
		symm_orders[i] = (core::Size)floor(pi/omega + 0.5);

		core::Real newW = -Wmult * cos( pi/symm_orders[i] );
		core::Real newS = sqrt ( (1-newW*newW)/(Qs[i].x*Qs[i].x+Qs[i].y*Qs[i].y+Qs[i].z*Qs[i].z) );
		Qs[i].x *= newS; Qs[i].y *= newS; Qs[i].z *= newS;
		Qs[i].w = newW;

		symm_axes[i] = numeric::xyzVector< core::Real >( Qs[i].x, Qs[i].y, Qs[i].z );
		symm_axes[i].normalize();
		quat2R( Qs[i], Rdimers[i] );
		delta_coms[i] = postT - preT;
		delta_coms[i].project_normal( symm_axes[i] );
	}

	// symmetry / motif count agreement
	if ( symm_type_ == "D2" ) {
		if ( nmotifs != 2 ) {
			TR << "nmotifs = " << nmotifs << std::endl;
			utility_exit_with_message( "Symmetry group D2: 2 motifs expected!" );
		}
		if ( symm_orders[1] != 2 || symm_orders[2] != 2 ) {
			TR << "Symm_orders = " << symm_orders[1] << "," << symm_orders[2] << std::endl;
			utility_exit_with_message( "Symmetry group D2: 2 homodimeric motifs expected!" );
		}
	} else {
		utility_exit_with_message( "Symmetry type unknown!" );
	}
}


// parse_my_tag
void
SymmetricMotifFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*reference_pose*/ ) {
	symm_type_ = tag->getOption<std::string>( "symm_type", "D2" );

	utility::vector1<std::string> motif_files( utility::string_split( tag->getOption< std::string >("motifs"), ',') );
	for ( Size i=1; i<=motif_files.size(); ++i ) {
		core::pose::PoseOP motif( new core::pose::Pose() );
		core::import_pose::pose_from_file( *motif, motif_files[i] , core::import_pose::PDB_file);
		ref_motifs_.push_back( motif );
	}
	core::Size nmotifs = ref_motifs_.size();
	TR << "Read " << nmotifs << " motifs" << std::endl;
	//to do: save on datamap?

	// override default options
	if ( tag->hasOption("angle_thresh") ) {
		angle_thresh_ = tag->getOption<core::Real>( "angle_thresh" );
	}
	if ( tag->hasOption("trans_thresh") ) {
		trans_thresh_ = tag->getOption<core::Real>( "trans_thresh" );
	}
	if ( tag->hasOption("rmsd_thresh") ) {
		rmsd_thresh_ = tag->getOption<core::Real>( "rmsd_thresh" );
	}
	if ( tag->hasOption("clash_thresh") ) {
		clash_thresh_ = tag->getOption<core::Size>( "clash_thresh" );
	}
	if ( tag->hasOption("angle_wt") ) {
		angle_wt_ = tag->getOption<core::Real>( "angle_wt" );
	}
	if ( tag->hasOption("trans_wt") ) {
		trans_wt_ = tag->getOption<core::Real>( "trans_wt" );
	}
	if ( tag->hasOption("rmsd_wt") ) {
		rmsd_wt_ = tag->getOption<core::Real>( "rmsd_wt" );
	}
	if ( tag->hasOption("clash_wt") ) {
		clash_wt_ = tag->getOption<core::Real>( "clash_wt" );
	}

	if ( tag->hasOption("force_pos") ) {
		utility::vector1<std::string> forced( utility::string_split( tag->getOption< std::string >("force_pos"), ',') );
		for ( Size i=1; i<=forced.size(); ++i ) {
			forced_pos_.push_back( std::atoi( forced[i].c_str() ) );
		}
	}

	process_motifs();

	// make sure that if we force motifs we force all of them
	runtime_assert( forced_pos_.size() == 0 || forced_pos_.size() == nsegs_ );
}

std::string SymmetricMotifFilter::name() const {
	return class_name();
}

std::string SymmetricMotifFilter::class_name() {
	return "SymmetricMotif";
}

void SymmetricMotifFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "symm_type" , xs_string , "Symmetric arrangement of motifs, i.e. Symmetry group D2: 2 motifs expected! Currently only D2 is supported." , "D2" )
		+ XMLSchemaAttribute( "motifs" , xs_string , "PDB file containing the symmetrized motif." )
		+ XMLSchemaAttribute( "angle_thresh" , xsct_real , "Angle threshold (degrees) for accepting a motif hit." )
		+ XMLSchemaAttribute( "trans_thresh" , xsct_real , "Translation threshold (Angstrom) for accepting a motif hit." )
		+ XMLSchemaAttribute( "rmsd_thresh" , xsct_real , "RMSD threshold for accepting a motif hit." )
		+ XMLSchemaAttribute( "clash_thresh" , xsct_non_negative_integer , "Most clashes allowed for accepting a motif hit." )
		+ XMLSchemaAttribute( "angle_wt" , xsct_real , "The angle weight for ranking all the passing hits." )
		+ XMLSchemaAttribute( "trans_wt" , xsct_real , "The translation weight for ranking all the passing hits." )
		+ XMLSchemaAttribute( "rmsd_wt" , xsct_real , "The RMSD for ranking all the passing hits." )
		+ XMLSchemaAttribute( "clash_wt" , xsct_real , "The clash weight for ranking all the passing hits." )
		+ XMLSchemaAttribute( "force_pos" , xs_string , "Force the motif to match at a particular position." ) ;

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Position-independent RMS filter evaluating how close a set of interfaces is to symmetric.", attlist );
}

std::string SymmetricMotifFilterCreator::keyname() const {
	return SymmetricMotifFilter::class_name();
}

protocols::filters::FilterOP
SymmetricMotifFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SymmetricMotifFilter );
}

void SymmetricMotifFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymmetricMotifFilter::provide_xml_schema( xsd );
}



} // filters
}
