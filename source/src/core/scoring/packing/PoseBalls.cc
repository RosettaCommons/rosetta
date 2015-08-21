// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packing/PoseBalls.cc
/// @brief
/// @author

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/packing/PoseBalls.hh>
#include <core/scoring/packstat/AtomRadiusMap.hh>

#include <basic/Tracer.hh>

#include <map>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

namespace core {
namespace scoring {
namespace packing {

/// @details Auto-generated virtual destructor
PoseBalls::~PoseBalls() {}

using namespace ObjexxFCL;

static thread_local basic::Tracer TR( "core.scoring.packing.PoseBalls" );

inline core::Real sqr ( core::Real x ) {
	return x*x;
}

// smoothed neighbor is between 9 and 11 A
inline core::Real sigmoidish_neighbor( core::Real sqdist ) {
	if ( sqdist > 121.0 ) {
		return 0.0;
	} else if ( sqdist < 81.0 ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0 - sqr( (dist - 9.0) / (11.0 - 9.0) ) );
	}
}

/// @details Remove spaces from given string.
inline std::string strip_whitespace( std::string const & name )
{
	std::string trimmed_name( name );
	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to dothis?
	return trimmed_name;
}


PoseBalls::PoseBalls(
	core::pose::Pose const & pose,
	core::Size Hmode,
	bool ignore_water
)
{
	using namespace numeric;

	// std::cerr << "PoseBalls.cc:67 (" << pose.total_residue() << ")" << std::endl;
	ObjexxFCL::FArray1D_char dssp_reduced_secstruct(pose.n_residue());
	core::scoring::dssp::Dssp(pose).dssp_reduced(dssp_reduced_secstruct);

	// initialize index and vars
	core::Size index = 0;
	Size num_unrec = 0;
	if ( pose.pdb_info() ) num_unrec = pose.pdb_info()->get_num_unrecognized_atoms();
	balls_      .reserve( pose.total_residue()*5 + num_unrec );
	index_to_id_.reserve( pose.total_residue()*5 + num_unrec );
	atom_name_  .reserve( pose.total_residue()*5 + num_unrec );
	atom_type_  .reserve( pose.total_residue()*5 + num_unrec );
	atom_parent_.reserve( pose.total_residue()*5 + num_unrec );
	secstruct_  .reserve( pose.total_residue()*5 + num_unrec );
	res_name_   .reserve( pose.total_residue()*5 + num_unrec );
	is_heavy_   .reserve( pose.total_residue()*5 + num_unrec );
	bfac_       .reserve( pose.total_residue()*5 + num_unrec );
	core::pose::initialize_atomid_map( id_to_index_, pose );
	id_to_index_.resize( pose.total_residue() + num_unrec );

	// std::cerr << "PoseBalls.cc:85 (" << ")" << std::endl;

	// add atoms in pose
	core::Size skippedH = 0;
	for ( core::Size ir = 1; ir <= pose.total_residue(); ++ir ) {
		for ( core::Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia ) {
			// std::cerr << "PoseBalls.cc:91 (" << ")" << std::endl;
			if ( pose.residue(ir).is_virtual(ia) ) continue;
			if ( Hmode == 0 ) { // no H's
				if ( pose.residue(ir).atom_type(ia).is_hydrogen() ) {
					skippedH++;
					continue;
				}
			} else if ( Hmode == 1 ) { // polar H's only
				if ( pose.residue(ir).atom_type(ia).is_hydrogen() &&
						!pose.residue(ir).atom_type(ia).is_polar_hydrogen() ) {
					skippedH++;
					continue;
				}
			} else if ( Hmode == 2 ) {
				; // all Hs
			}
			// std::cerr << "PoseBalls.cc:107 (" << ")" << std::endl;
			core::id::AtomID aid(ia,ir);
			id_to_index_[ aid ] = ++index;
			index_to_id_.push_back( aid );
			balls_.push_back( Ball( pose.xyz(aid), pose.residue(ir).atom_type(ia).lj_radius() ) );
			res_name_.push_back(pose.residue(ir).name3());
			atom_num_.push_back(ia);
			// if( pose.residue(ir).atom_type(ia).is_polar_hydrogen() ) {
			//  atom_type_.push_back(100+pose.residue(ir).atom_type_index(pose.residue(ir).atom_base(ia)));
			// } else {
			atom_type_.push_back(pose.residue(ir).atom_type_index(ia));
			// }
			secstruct_.push_back( dssp_reduced_secstruct(ir) );
			atom_name_.push_back(pose.residue(ir).atom_type(ia).name());
			if ( pose.pdb_info() ) bfac_.push_back( pose.pdb_info()->temperature(ir,ia) );
			else                  bfac_.push_back( 0.0 );
			is_heavy_.push_back( pose.residue(ir).atom_type(ia).is_heavyatom() );
			if ( pose.residue(ir).atom_type(ia).is_hydrogen() /*&&
					!pose.residue(ir).atom_type(ia).is_polar_hydrogen()*/ ) {
				// if is apolar H, set add surf to base heavy atom rather than H itself
				core::Size iabase = pose.residue(ir).atom_base(ia);
				atom_parent_.push_back( id_to_index(core::id::AtomID(iabase,ir)) );
			} else {
				atom_parent_.push_back(index);
			}
			res_num_.push_back(ir);
		}
	}
	// std::cerr << "PoseBalls.cc:134 (" << ")" << std::endl;

	// add hetero atoms not in pose
	// res num is *NOT* necessarially the same as PDBInfo res num!!!!!
	core::scoring::packstat::AtomRadiusMap arm;
	std::map<core::Size,core::Size> atomcount;
	std::map<core::Size,core::Size> resnum;
	core::Size rescount = 0, skippedlig = 0, skippedligH = 0;
	for ( core::Size i = 1; i <= num_unrec; ++i ) {
		core::pose::UnrecognizedAtomRecord const & a( pose.pdb_info()->get_unrecognized_atoms()[i] );
		core::Real radius = arm.get_radius(a.atom_name(),a.res_name());
		if ( Hmode < 2 && radius <= 1.2 ) {
			skippedligH++;
			continue;
		}
		if ( radius <= 0.1 ) {
			skippedlig++;
			continue;
		}
		if ( ignore_water ) { // TODO other names for water?
			if ( "HOH" == (a.res_name()) || "DOD" == (a.res_name()) ) {
				skippedlig++;
				continue;
			}
		}
		if ( resnum.find(a.res_num()) == resnum.end() ) {
			resnum[a.res_num()] = ++rescount + pose.total_residue();
			atomcount[a.res_num()] = 0;
		}
		id_to_index_.resize( resnum[a.res_num()], pose.pdb_info()->get_unrecognized_res_size(a.res_num()) );
		core::id::AtomID aid( ++atomcount[a.res_num()], resnum[a.res_num()] );
		id_to_index_[ aid ] = ++index;
		index_to_id_.push_back( aid );
		atom_name_.push_back(a.atom_name());
		atom_type_.push_back(0);
		secstruct_.push_back('U');
		atom_num_.push_back(atomcount[a.res_num()]);
		atom_parent_.push_back(index);
		is_heavy_.push_back( radius > 1.3 );
		res_name_.push_back(a.res_name());
		res_num_.push_back(resnum[a.res_num()]);
		bfac_.push_back( a.temp() );
		balls_.push_back( Ball( a.coords(), radius ) );

	}
	nballs_ = index;
	reset_surf();
	compute_smooth_nb();
}

PoseBalls::PoseBalls(
	core::pose::Pose const & pose,
	core::id::AtomID_Mask const & whichatoms
){
	using namespace numeric;

	// std::cerr << "PoseBalls.cc:67 (" << pose.total_residue() << ")" << std::endl;
	ObjexxFCL::FArray1D_char dssp_reduced_secstruct(pose.n_residue());
	core::scoring::dssp::Dssp(pose).dssp_reduced(dssp_reduced_secstruct);

	// initialize index and vars
	core::Size index = 0;
	Size num_unrec = 0;
	if ( pose.pdb_info() ) num_unrec = pose.pdb_info()->get_num_unrecognized_atoms();
	balls_      .reserve( pose.total_residue()*5 + num_unrec );
	index_to_id_.reserve( pose.total_residue()*5 + num_unrec );
	atom_name_  .reserve( pose.total_residue()*5 + num_unrec );
	atom_type_  .reserve( pose.total_residue()*5 + num_unrec );
	atom_parent_.reserve( pose.total_residue()*5 + num_unrec );
	secstruct_  .reserve( pose.total_residue()*5 + num_unrec );
	res_name_   .reserve( pose.total_residue()*5 + num_unrec );
	is_heavy_   .reserve( pose.total_residue()*5 + num_unrec );
	bfac_       .reserve( pose.total_residue()*5 + num_unrec );
	core::pose::initialize_atomid_map( id_to_index_, pose );
	id_to_index_.resize( pose.total_residue() + num_unrec );

	// std::cerr << "PoseBalls.cc:85 (" << ")" << std::endl;

	// add atoms in pose
	for ( core::Size ir = 1; ir <= pose.total_residue(); ++ir ) {
		for ( core::Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia ) {
			core::id::AtomID aid(ia,ir);
			if ( ! whichatoms[aid] ) continue;
			id_to_index_[ aid ] = ++index;
			index_to_id_.push_back( aid );
			balls_.push_back( Ball( pose.xyz(aid), pose.residue(ir).atom_type(ia).lj_radius() ) );
			res_name_.push_back(pose.residue(ir).name3());
			atom_num_.push_back(ia);
			// if( pose.residue(ir).atom_type(ia).is_polar_hydrogen() ) {
			//  atom_type_.push_back(100+pose.residue(ir).atom_type_index(pose.residue(ir).atom_base(ia)));
			// } else {
			atom_type_.push_back(pose.residue(ir).atom_type_index(ia));
			// }
			secstruct_.push_back( dssp_reduced_secstruct(ir) );
			atom_name_.push_back(pose.residue(ir).atom_type(ia).name());
			if ( pose.pdb_info() ) bfac_.push_back( pose.pdb_info()->temperature(ir,ia) );
			else                  bfac_.push_back( 0.0 );
			is_heavy_.push_back( pose.residue(ir).atom_type(ia).is_heavyatom() );
			if ( pose.residue(ir).atom_type(ia).is_hydrogen() /*&&
					!pose.residue(ir).atom_type(ia).is_polar_hydrogen()*/ ) {
				// if is apolar H, set add surf to base heavy atom rather than H itself
				core::Size iabase = pose.residue(ir).atom_base(ia);
				atom_parent_.push_back( id_to_index(core::id::AtomID(iabase,ir)) );
			} else {
				atom_parent_.push_back(index);
			}
			res_num_.push_back(ir);
		}
	}
	// std::cerr << "PoseBalls.cc:134 (" << ")" << std::endl;

	// add hetero atoms not in pose
	// res num is *NOT* necessarially the same as PDBInfo res num!!!!!
	core::scoring::packstat::AtomRadiusMap arm;
	std::map<core::Size,core::Size> atomcount;
	std::map<core::Size,core::Size> resnum;
	core::Size rescount = 0, skippedlig = 0;
	for ( core::Size i = 1; i <= num_unrec; ++i ) {
		core::pose::UnrecognizedAtomRecord const & a( pose.pdb_info()->get_unrecognized_atoms()[i] );
		core::Real radius = arm.get_radius(a.atom_name(),a.res_name());
		if ( radius <= 0.1 ) {
			skippedlig++;
			continue;
		}
		if ( resnum.find(a.res_num()) == resnum.end() ) {
			resnum[a.res_num()] = ++rescount + pose.total_residue();
			atomcount[a.res_num()] = 0;
		}
		id_to_index_.resize( resnum[a.res_num()], pose.pdb_info()->get_unrecognized_res_size(a.res_num()) );
		core::id::AtomID aid( ++atomcount[a.res_num()], resnum[a.res_num()] );
		id_to_index_[ aid ] = ++index;
		index_to_id_.push_back( aid );
		atom_name_.push_back(a.atom_name());
		atom_type_.push_back(0);
		secstruct_.push_back('U');
		atom_num_.push_back(atomcount[a.res_num()]);
		atom_parent_.push_back(index);
		is_heavy_.push_back( radius > 1.3 );
		res_name_.push_back(a.res_name());
		res_num_.push_back(resnum[a.res_num()]);
		bfac_.push_back( a.temp() );
		balls_.push_back( Ball( a.coords(), radius ) );

	}
	nballs_ = index;
	reset_surf();
	compute_smooth_nb();
}


void PoseBalls::compute_smooth_nb() {

	using namespace numeric;

	// compute NB counts -- no Hs for neighbors
	smooth_nb_.resize(nballs_);
	// nbhist_.resize(29);
	// for( core::Size i = 1; i <= 29; ++i ) nbhist_[i].resize(100,0);
	for ( core::Size i = 1; i <= nballs_; i++ ) smooth_nb_[i] = 0.0;
	for ( core::Size i = 1; i <= nballs_; i++ ) {
		// if( is_heavy_[i] ) smooth_nb_[i] += 1.0;
		// else continue;
		xyzVector<Real> const & ixyz( balls_[i].xyz() );
		core::Size atype = atom_type_[i];
		if ( atype >= 18 ) {
			if ( 'E'==secstruct_[i] ) atype = 18+3*(atype-18)+0;
			else if ( 'H'==secstruct_[i] ) atype = 18+3*(atype-18)+1;
			else if ( 'L'==secstruct_[i] ) atype = 18+3*(atype-18)+2;
			else atype = 999;//std::cerr << "secstruct other than H E or L!!!! " << secstruct_[i] << std::endl;
		}
		for ( core::Size j = 1; j <= nballs_; j++ ) {
			if ( !is_heavy_[j] ) continue;
			xyzVector<Real> const & jxyz( balls_[j].xyz() );
			Real const d2( ixyz.distance_squared(jxyz) );
			if ( d2 < 121.0 ) {
				Real const sn( sigmoidish_neighbor(d2) );
				smooth_nb_[i] += sn;
				// smooth_nb_[j] += sn;
				//     if( atype > 29 || atype < 1 ) continue;
				//     if( i==j ) continue;
				//     Real d = sqrt(d2);
				//     for( core::Size idth = 1; idth <= 100; idth++ ) {
				//      Real dth = (Real)(idth) / 10.0;
				//      // std::cerr << "DIST " << d << " " << dth << " " << idth << std::endl;
				//      if( d < dth ) {
				//       nbhist_[atype][idth]++;
				//       break;
				//      }
				//     }
			}
		}
	}
}


void
PoseBalls::output_pdb( std::ostream & out ) const {

	using namespace ObjexxFCL::format;

	for ( Size i = 1; i <= nballs(); i++ ) {

		out << std::string("ATOM  " + I( 5, i ) + " "+LJ(4,atom_name(i))+" "+LJ(3,res_name(i))+" " + " "
			+ I( 4, res_num(i) ) + "    "
			+ F( 8, 3, ball(i).x() ) + F( 8, 3, ball(i).y() ) + F( 8, 3, ball(i).z() )
			+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, ball(i).r() ) ) + "\n";
	}

}


} // namespace packing
} // namespace scoring
} // namespace core
