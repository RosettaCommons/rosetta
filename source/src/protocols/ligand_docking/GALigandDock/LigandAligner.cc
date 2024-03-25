// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/LigandAligner.cc
///
/// @brief  Align a ligand to a reference ligand or an inferred reference ligand
/// @author Hahnbeom Park and Frank DiMaio


#include <protocols/ligand_docking/GALigandDock/LigandAligner.hh>
#include <protocols/ligand_docking/GALigandDock/GriddedAtomTreeMultifunc.hh>
#include <protocols/ligand_docking/GALigandDock/util.hh>

#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/func/SumFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>
#include <core/scoring/func/ConstantFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/LinearPenaltyFunction.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/residue_datacache.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/Bond.fwd.hh>
#include <core/id/AtomID.hh>


#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random_xyz.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// density aligner headers
#include <ObjexxFCL/FArray3D.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <numeric/fourier/SHT.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>


// C++ headers
#include <fstream>
#include <string>
//boost header
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <core/optimization/MinimizerMap.hh> // AUTO IWYU For MinimizerMap
#include <core/scoring/hbonds/hbonds.hh> // AUTO IWYU For fill_hbond_set
#include <core/scoring/lkball/LK_BallEnergy.hh> // AUTO IWYU For LK_BallEnergy

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

static basic::Tracer TR( "ligand_align" );

// property resolution: does a pair of atoms match?  What is the scale?
// halogen_specific treats halogens separately, only matching to other halogens
core::Real
AtomProperties::match( AtomProperties const &other, core::Real polar_scale, bool halogen_specific ) const {
	// nonpolar heavy <-> nonpolar heavy
	bool nonpolarHeavy = !is_donor && !is_acceptor && !is_H && !is_halogen;
	bool other_nonpolarHeavy = !other.is_donor && !other.is_acceptor && !other.is_H && !other.is_halogen;

	if ( other.is_generic ) return 1.0;

	if ( nonpolarHeavy != other_nonpolarHeavy ) return 0.0;

	// polar heavy <-> polar heavy (w/ donor/acceptor agreement)
	if ( is_donor != other.is_donor && is_acceptor != other.is_acceptor ) return 0.0;

	// halogen
	if ( halogen_specific && (is_halogen != other.is_halogen) ) return 0.0;

	// H <-> H //???
	if ( is_H != other.is_H || is_polarH != other.is_polarH ) return 0.0;

	// we have a match!  are we polar?  if so, upweight
	if ( is_donor || is_acceptor || is_polarH ) return polar_scale;

	return 1.0;
}

std::string
AtomProperties::show() const
{ return std::to_string(is_donor)+
	std::to_string(is_acceptor)+
	std::to_string(is_H)+
	std::to_string(is_polarH)+
	std::to_string(is_halogen);
}

/////////////////////////////////////////////
Pharmacophore::Pharmacophore()
{
	atms_.resize( 0 );
	props_.resize( 0 );
	dist_.resize( 0 );
	com_ = numeric::xyzVector< core::Real >( 0.0 );
	index_ = "";
}

// general constructor
Pharmacophore::Pharmacophore( utility::vector1< core::Size > atmno,
	utility::vector1< AtomProperties > const & p,
	utility::vector1< utility::vector1< core::Real > > const & Dmtrx,
	utility::vector1< numeric::xyzVector< core::Real > > const & coords
)
{
	core::Size const n( atmno.size() );
	dist_ = utility::vector1< utility::vector1< core::Real > >
		( n, utility::vector1< core::Real >( n, 0.0 ) );

	com_ = numeric::xyzVector< core::Real >( 0.0 );
	for ( core::Size i = 1; i <= n; ++i ) {
		atms_.push_back( atmno[i] );
		props_.push_back( p[atmno[i]] );
		com_ += coords[atmno[i]];
	}
	com_ /= n;

	for ( core::Size i = 1; i < n; ++i ) {
		for ( core::Size j = i+1; j < n; ++j ) {
			dist_[i][j] = dist_[j][i] = Dmtrx[atmno[i]][atmno[j]];
		}
	}

	index_ = "";
	for ( core::Size i = 1; i < n; ++i ) index_ += std::to_string(atmno[i])+"/";
	index_ += std::to_string(atmno[n]);
}

void
Pharmacophore::update_com( utility::vector1< numeric::xyzVector< core::Real > > const & coords )
{
	com_ = numeric::xyzVector< core::Real >( 0.0 );
	for ( core::Size i = 1; i <= atms_.size(); ++i ) com_ += coords[atms_[i]];
	com_ /= atms_.size();
}

std::string
Pharmacophore::show() const{
	std::string tag("");
	for ( core::Size iatm = 1; iatm < atms_.size(); ++iatm ) {
		tag += props_[iatm].tag()+"/";
	}
	tag += props_[atms_.size()].tag();
	return tag;
}

// close to all member
bool
Pharmacophore::is_close( utility::vector1< core::Real > const &Dcol,
	core::Real dcut ) const
{
	bool close( true );
	for ( core::Size i = 1; i <= atms_.size(); ++i ) {
		if ( Dcol[atms_[i]] > dcut ) {
			close = false;
			break;
		}
	}
	return close;
}

// works within Ligand using receptor pharmacophore
core::Real
Pharmacophore::find_type_matches( Pharmacophore const &other,
	bool const just_count ) const
{
	core::Size nH_lig( 0 ), nA_lig( 0 ), nH_rec( 0 ), nA_rec( 0 );
	core::Real ligH_scores( 0.0 ), ligA_scores( 0.0 ), recH_scores( 0.0 ), recA_scores( 0.0 );

	for ( core::Size iatm = 1; iatm <= size(); ++iatm ) {
		if ( is_polarH( iatm ) ) {
			//ligH_scores.push_back( properties_[iatm].score() );
			ligH_scores += props_[iatm].score();
			nH_lig++;
		}
		if ( is_acceptor( iatm ) ) {
			//ligA_scores.push_back( properties_[iatm].score() );
			ligA_scores += props_[iatm].score();
			nA_lig++;
		}
	}

	// make sure not add Zscore > 0.0....
	utility::vector1< AtomProperties > other_prop = other.props();
	for ( core::Size jatm = 1; jatm <= other.size(); ++jatm ) {
		if ( other.is_polarH( jatm ) ) {
			//recH_scores.push_back( other_prop[jatm].score() );
			recH_scores += other_prop[jatm].score();
			nH_rec++;
		}
		if ( other.is_acceptor( jatm ) ) {
			//recH_scores.push_back( other_prop[jatm].score() );
			recA_scores += other_prop[jatm].score();
			nA_rec++;
		}
	}

	//core::Real denom = (core::Real)(std::max(nH_rec+nA_rec,nH_lig+nA_lig)); // maximum possible
	core::Size nH_match( std::min(nH_rec,nH_lig) );
	core::Size nA_match( std::min(nA_rec,nA_lig) );

	//core::Real norm = (core::Real)(nH_match+nA_match); // overlap

	if ( nH_rec > 0 ) recH_scores /= nH_rec;
	if ( nA_rec > 0 ) recA_scores /= nA_rec;
	if ( nH_lig > 0 ) ligH_scores /= nH_lig;
	if ( nA_lig > 0 ) ligA_scores /= nA_lig;

	//core::Real score = norm; // more the better
	// take average score weighted by num matches
	core::Real score( 0.0 );
	if ( just_count ) {
		score = nH_match+nA_match;
	} else {
		if ( nH_match > 0 ) score += nH_match*(recH_scores+ligH_scores);
		if ( nA_match > 0 ) score += nA_match*(recA_scores+ligA_scores);
	}

	return score;
}

/*
// works within Ligand using receptor pharmacophore
bool
Pharmacophore::find_matching_permutation( Pharmacophore const &other,
utility::vector1< core::Size > &map_index,
core::Real &best_score ) const
{
// make sure combinatorinal order is calculated
bool mapped( false );

// make combinations; note that ligand contains up to 3-body phores only
utility::vector1< utility::vector1< core::Size > > combs;
if ( size() == 3 ) {
combs.resize( 6 );
combs[1].push_back( atms_[1] ); combs[1].push_back( atms_[2] ); combs[1].push_back( atms_[3] );
combs[2].push_back( atms_[1] ); combs[2].push_back( atms_[3] ); combs[2].push_back( atms_[2] );
combs[3].push_back( atms_[2] ); combs[3].push_back( atms_[1] ); combs[3].push_back( atms_[3] );
combs[4].push_back( atms_[2] ); combs[4].push_back( atms_[3] ); combs[4].push_back( atms_[1] );
combs[5].push_back( atms_[3] ); combs[5].push_back( atms_[1] ); combs[5].push_back( atms_[2] );
combs[6].push_back( atms_[3] ); combs[6].push_back( atms_[2] ); combs[6].push_back( atms_[1] );
} else if ( size() == 2 ) {
combs.resize( 2 );
combs[1].push_back( atms_[1] ); combs[1].push_back( atms_[2] );
combs[2].push_back( atms_[2] ); combs[1].push_back( atms_[1] );
} else {
combs.resize( 1 );
combs[1].push_back( atms_[1] );
}

// iter through all possible combination
for ( core::Size i = 1; i <= combs.size(); ++i ) {
utility::vector1< core::Size > const &comb_i = combs[i];
core::Real score = match( other, comb_i, best_score );
if ( score < best_score ) {
map_index = comb_i;
mapped = true;
best_score = score;
}
}
return mapped;
}
*/

core::Real
Pharmacophore::match( Pharmacophore const &other,
	utility::vector1< core::Size > const &map_index,
	core::Real const unmatch_score ) const
{
	int n( other.size() ), m( size() );

	// first check chemical property
	// this still works even n == 1 or m == 1
	for ( int i = 1; i <= (int)map_index.size(); ++i ) {
		if ( i > n ) continue;
		core::Size my_i = map_index[i];
		bool is_complementary =
			(is_polarH(my_i) && other.is_acceptor(i)) || (is_acceptor(my_i) && other.is_polarH(i));
		if ( !is_complementary ) return unmatch_score+1.0;
	}

	// second, check distance relation
	core::Real const score_missing( 3.0 ); // treat as 3 Ang apart for missing
	core::Real score( std::abs(n-m)*score_missing );

	for ( int i = 1; i <= (int)map_index.size(); ++i ) {
		core::Size my_i = map_index[i];
		for ( int j = i+1; j <= (int)map_index.size(); ++j ) {
			core::Size my_j = map_index[j];
			if ( i > n || j > n ) continue;
			score += std::abs( dist(my_i,my_j) - other.dist(i,j) );
		}
	}
	score /= core::Real(std::min(n,m));
	return score;
}

void
ConstraintInfo::set_default()
{
	limit_ = 4.0;
	weight_ = 3.0;
	polar_scale_ = 5.0;
	phore_dcut_ = 3.0; // use tighter for ligand
	n_match_per_ligand_motif_ = 15;
	debug_ = TR.Debug.visible();

	// hard cutoff is sometimes quite risky... (was 11/6 for 8Ang cut)
	NNEIGH_CUT_FOR_EXPOSED_SC = 18;
	NNEIGH_CUT_FOR_EXPOSED_BB = 15;
	RECEPTOR_PHORE_MIN = 5;
	RECEPTOR_PHORE_MAX = 15;

	//RECEPTOR_VSITE_MIN = 25;
	//RECEPTOR_VSITE_MAX = 40;
}

/////////////////////////////////////////////
//  initialize from a (possibly different) ligand conformation
void
ConstraintInfo::init_from_ligand(
	core::pose::Pose const &pose,
	utility::vector1< core::Size > const &ligids
) {
	set_default();
	is_ligand_ = true;

	// this handles multi-res reference poses, as long as reference pose has FEWER residues than target
	// TO DO: handle the opposite case
	for ( auto ligid : ligids ) {
		if ( ligid > pose.total_residue() ) {
			TR << "Warn: ligid " << ligid << " not found in reference pose.  Continuing!" << std::endl;
			continue;
		}
		core::conformation::Residue const &mylig = pose.residue(ligid);
		for ( core::Size iatm_src=1; iatm_src<=mylig.natoms(); ++iatm_src ) {
			// not including aliphatic hydrogen... too messy
			if ( mylig.atom_type(iatm_src).is_hydrogen() &&
					!mylig.atom_type(iatm_src).is_polar_hydrogen() ) continue;

			coords_.push_back( mylig.xyz(iatm_src) );
			atmids_.push_back( core::id::AtomID( iatm_src, ligid ) ); // resno 0 is ligand

			bool is_halogen = (
				mylig.atom_type(iatm_src).element() == "F" ||
				mylig.atom_type(iatm_src).element() == "Cl" ||
				mylig.atom_type(iatm_src).element() == "Br" ||
				mylig.atom_type(iatm_src).element() == "I" );

			properties_.push_back( AtomProperties(
				mylig.atom_type(iatm_src).is_donor(),
				mylig.atom_type(iatm_src).is_acceptor(),
				iatm_src > mylig.nheavyatoms(),    // isH
				mylig.atom_type(iatm_src).is_polar_hydrogen(),
				is_halogen,
				false, //is_generic
				0.0, // non-ambiguous
				1.0, // score weight; default 1
				+"L."+utility::strip(mylig.atom_name(iatm_src)) // tag
				));
		}
	}
}

void
ConstraintInfo::init_from_map( utility::vector1< numeric::xyzVector< core::Real > > density_points ) {
	//limit_ = 4.0;
	//weight_ = 3.0;
	set_default();

	for ( core::Size ipoint = 1; ipoint <= density_points.size(); ++ipoint ) {
		coords_.push_back( density_points[ipoint] );

		properties_.push_back( AtomProperties(
			false,
			false,
			false,
			false,
			false,
			true,
			0.0,
			1.0,
			+"L."+ipoint
			));
	}
}

//  initialize from a (possibly different) ligand conformation
void
ConstraintInfo::init_from_ligand_pharmacophore(
	core::pose::Pose const &pose,
	utility::vector1< core::Size > const &ligids,
	bool report_phore_info
) {
	set_default();
	is_ligand_ = true;

	for ( auto resid : ligids ) {
		core::conformation::Residue const &mylig = pose.residue(resid);
		for ( core::Size iatm_src=1; iatm_src<=mylig.natoms(); ++iatm_src ) {
			bool is_acceptor( mylig.atom_type(iatm_src).is_acceptor() );
			// hard-coded generalization logic to append missing "pseudo-acceptors"...
			// N or O; No attached polarH; large negative charge
			if ( iatm_src <= mylig.nheavyatoms() ) {
				if ( !mylig.heavyatom_has_polar_hydrogens(iatm_src) &&
						mylig.atomic_charge(iatm_src) < -0.5 &&
						( mylig.atom_type(iatm_src).element() == "N" || mylig.atom_type(iatm_src).element() == "O" )
						) {
					is_acceptor = true;
				}
			}

			// just to get estimation of numbers... should be changed if more rules being added
			bool is_using( mylig.atom_type(iatm_src).is_donor() ||
				is_acceptor ||
				mylig.atom_type(iatm_src).is_polar_hydrogen() );
			if ( !is_using ) continue;

			bool is_halogen = (
				mylig.atom_type(iatm_src).element() == "F" ||
				mylig.atom_type(iatm_src).element() == "Cl" ||
				mylig.atom_type(iatm_src).element() == "Br" ||
				mylig.atom_type(iatm_src).element() == "I" );

			atmids_.push_back( core::id::AtomID( iatm_src, resid ) ); // resno 0 is ligand
			coords_.push_back( mylig.xyz(iatm_src) );
			properties_.push_back( AtomProperties(
				mylig.atom_type(iatm_src).is_donor(),
				is_acceptor,
				iatm_src > mylig.nheavyatoms(),    // isH
				mylig.atom_type(iatm_src).is_polar_hydrogen(),
				is_halogen,
				false, //is_generic (used for density)
				0.0, // everything non-ambiguous... maybe assign differently for hxl?
				1.0,
				+"L."+utility::strip(mylig.atom_name(iatm_src)) // tag
				));
		}

	}

	define_all_ligand_phores( 1, report_phore_info );
}

void
ConstraintInfo::init_from_receptor(
	core::pose::Pose const &pose,
	GridScorerOP gridscore,
	bool const simple
)
{
	set_default();
	phore_dcut_ = 4.0; // default 3.0; use looser for receptor
	is_ligand_ = false;

	// used for seed selection
	utility::vector1< std::pair< core::Real, core::Size > > Vdonor_sort, Vacceptor_sort;

	// first try with default setup of exposure params
	define_active_virtual_sites( pose, gridscore, Vdonor_sort, Vacceptor_sort,
		NNEIGH_CUT_FOR_EXPOSED_SC, NNEIGH_CUT_FOR_EXPOSED_BB );

	// PHdock: adjust solvent exposure params so that #receptor phores be 10~15
	// do not add 1-site phore yet
	define_receptor_phores( Vdonor_sort, Vacceptor_sort,  0, false );

	if ( simple ) return;

	// get best optimal numbers
	int shift( 0 );
	if ( phores_.size() < RECEPTOR_PHORE_MIN ) {
		for ( shift = -1; shift >= -5; --shift ) {
			TR << "Too few sites, redefining with receptor sites with NNEIGH_CUT_EXPOSED_BB/SC "

				<< NNEIGH_CUT_FOR_EXPOSED_SC << " -> " << NNEIGH_CUT_FOR_EXPOSED_SC+shift << std::endl;
			define_active_virtual_sites( pose, gridscore, Vdonor_sort, Vacceptor_sort,
				NNEIGH_CUT_FOR_EXPOSED_SC+shift, NNEIGH_CUT_FOR_EXPOSED_BB );
			define_receptor_phores( Vdonor_sort, Vacceptor_sort, 0, false );
			if ( phores_.size() >= RECEPTOR_PHORE_MIN ) break;
		}

	} else if ( phores_.size() > RECEPTOR_PHORE_MAX ) {
		// this is dangerous for exposed bb; only sc
		for ( shift = 1; shift <= 5; ++shift ) {
			TR << "Too many sites, redefining with receptor sites with NNEIGH_CUT_EXPOSED_BB/SC "
				//<< NNEIGH_CUT_FOR_EXPOSED_BB << " -> " << NNEIGH_CUT_FOR_EXPOSED_BB+shift << " "
				<< NNEIGH_CUT_FOR_EXPOSED_SC << " -> " << NNEIGH_CUT_FOR_EXPOSED_SC+shift << std::endl;
			define_active_virtual_sites( pose, gridscore, Vdonor_sort, Vacceptor_sort,
				NNEIGH_CUT_FOR_EXPOSED_SC+shift, NNEIGH_CUT_FOR_EXPOSED_BB );

			define_receptor_phores( Vdonor_sort, Vacceptor_sort, 0, false );
			if ( phores_.size() <= RECEPTOR_PHORE_MAX ) break;
		}
	}

	// report final status
	if ( debug_ ) {
		TR << "Get receptor phores from final guess of NNEIGH_CUT_EXPOSED_BB/SC "
			<< NNEIGH_CUT_FOR_EXPOSED_BB << " " << NNEIGH_CUT_FOR_EXPOSED_SC+shift << std::endl;
		define_active_virtual_sites( pose, gridscore, Vdonor_sort, Vacceptor_sort,
			NNEIGH_CUT_FOR_EXPOSED_SC+shift, NNEIGH_CUT_FOR_EXPOSED_BB, true );
		// add more 1-site phores if < 10 defined
		define_receptor_phores( Vdonor_sort, Vacceptor_sort, 10, true );
	}

}

void
ConstraintInfo::define_active_virtual_sites(
	core::pose::Pose const &pose,
	GridScorerOP gridscore,
	utility::vector1< std::pair< core::Real, core::Size > > &Vdonor_sort,
	utility::vector1< std::pair< core::Real, core::Size > > &Vacceptor_sort,
	core::Size const nneigh_cut_for_exposed_sc,
	core::Size const nneigh_cut_for_exposed_bb,
	bool const report
)
{
	atmids_.resize( 0 );
	coords_.resize( 0 );
	properties_.resize( 0 );
	Vdonor_sort.resize( 0 );
	Vacceptor_sort.resize( 0 );

	core::scoring::methods::EnergyMethodOptions e_opts;
	core::scoring::lkball::LK_BallEnergy lkb(e_opts);

	core::scoring::hbonds::HBondSet hbset;
	core::scoring::hbonds::fill_hbond_set( pose, false, hbset );

	// Exposure check: 10 Angstrom!!
	utility::vector1< core::Size > nneigh_res_bb_naive = count_neighbors( pose, "CA", 10.0 );
	utility::vector1< core::Size > nneigh_res_sc_naive = count_neighbors( pose, "nbr", 10.0 );

	// smoothening
	utility::vector1< core::Real > nneigh_res_bb( pose.size(), 0 ), nneigh_res_sc( pose.size(), 0 );
	int const nres( pose.size() );

	for ( int i=1; i<=nres; ++i ) {
		// non-smoothing version
		nneigh_res_bb[i] = (core::Real)(nneigh_res_bb_naive[i]);
		nneigh_res_sc[i] = (core::Real)(nneigh_res_sc_naive[i]);
	}

	std::string rsd_exposed("");
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( nneigh_res_sc[i] <= nneigh_cut_for_exposed_sc ) rsd_exposed += std::to_string(i)+"+";
	}

	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		core::conformation::Residue const &rsd = pose.residue(i);

		// logic to skip already h-bonding backbones
		utility::vector1< core::scoring::hbonds::HBondCOP > const & resHBs = hbset.residue_hbonds(i);
		bool H_has_hbond( false ), O_has_hbond( false );
		for ( core::Size ihb = 1; ihb <= resHBs.size(); ++ihb ) {
			core::scoring::hbonds::HBondCOP hb = resHBs[ihb];
			if ( hb->don_res() == i && hb->don_hatm_is_backbone() && (hb->energy() < -0.5) ) H_has_hbond = true;
			if ( hb->acc_res() == i && hb->acc_atm_is_backbone() && (hb->energy() < -0.5) ) O_has_hbond = true;
		}

		utility::vector1< core::Vector > rsd_waters;
		utility::vector1< core::Size > rsd_n_attached_waters;

		core::scoring::lkball::LKB_ResidueInfo const &lkbinfo = static_cast< core::scoring::lkball::LKB_ResidueInfo const & > (
			rsd.data_ptr()->get( core::conformation::residue_datacache::LK_BALL_INFO ));
		rsd_n_attached_waters = lkbinfo.n_attached_waters();
		rsd_waters = ( lkbinfo.waters() );
		utility::vector1< core::Size > const &offset = lkbinfo.water_offset_for_atom();

		for ( core::Size j=1; j<=rsd.natoms(); ++j ) {
			if ( rsd_n_attached_waters[j]==0 ) continue;
			bool is_bb( j <= pose.residue(i).last_backbone_atom() );

			//bool is_exposed_sc( !is_bb && (nneigh_res_sc[i] <= nneigh_cut_for_exposed_sc) && pose.residue(i).is_protein() );
			//bool is_exposed_bb( is_bb  && (nneigh_res_bb[i] <= nneigh_cut_for_exposed_bb) && pose.residue(i).is_protein() );
			bool is_acc = rsd.atom_type(j).is_acceptor();
			bool is_don = rsd.atom_type(j).is_donor();
			bool is_OCbb = (rsd.atom_type(j).atom_type_name() == "OCbb");
			bool is_ONH2 = (rsd.atom_type(j).atom_type_name() == "ONH2");
			bool is_Hbb = (rsd.atom_type(j).atom_type_name() == "HNbb");
			bool is_hxl = (rsd.atom_type(j).atom_type_name() == "OH");

			bool is_hbonding( (is_OCbb && O_has_hbond) || (is_Hbb && H_has_hbond) );
			//bool accepted_in_stage1( !is_exposed_sc && !is_exposed_bb && !is_hbonding );
			bool accepted_in_stage1( !is_hbonding );

			// logic for Vsite definitions:
			//   if we are a donor, put acceptor atoms at water positions
			//   if we are an acceptor, put donor atoms at water positions, and H pointing at atom
			//   if both, ignore current H position, and treat as both at all water positions

			utility::vector1< numeric::xyzVector< core::Real > > virtual_sites;

			std::string atype = is_acc? "Hpol":"OH"; //reference grids in GridScorer
			std::string aname = is_acc?(is_bb?"H  ":"HA "):(is_bb?"O  ":"OXT");

			// take average position for OCbb instead of taking both...

			if ( (is_OCbb || is_hxl || is_ONH2 ) && rsd_n_attached_waters[j] > 1 ) {
				numeric::xyzVector< core::Real > avrg_water_xyz( 0.0 );

				for ( core::Size iwat = 1; iwat <= rsd_n_attached_waters[j]; ++iwat ) {
					avrg_water_xyz += rsd_waters[offset[j]+iwat];
				}
				avrg_water_xyz /= rsd_n_attached_waters[j];

				numeric::xyzVector< core::Real > delD = avrg_water_xyz - rsd.xyz(j);
				avrg_water_xyz = rsd.xyz(j) + 1.7 * delD / delD.length(); // hydrogen

				virtual_sites.push_back( avrg_water_xyz );

			} else {
				// type differs?
				if ( is_acc ) { // virtual is donor
					for ( core::Size iwat = 1; iwat <= rsd_n_attached_waters[j]; ++iwat ) {
						numeric::xyzVector< core::Real > delH = rsd_waters[offset[j]+iwat] - rsd.xyz(j);
						numeric::xyzVector< core::Real > H = rsd.xyz(j) + 1.7*delH / delH.length();
						virtual_sites.push_back( H );
					}
				} else {
					for ( core::Size iwat = 1; iwat <= rsd_n_attached_waters[j]; ++iwat ) virtual_sites.push_back( rsd_waters[offset[j]+iwat] );
				}
			}

			if ( is_don || is_acc ) {
				for ( core::Size k=1; k<=virtual_sites.size(); ++k ) {

					//core::Size nneigh_res = count_neighbors_on_coord( pose, virtual_sites[k], "nbr", 10.0 );
					core::Size nneigh_res = is_bb?nneigh_res_bb_naive[i] : nneigh_res_sc_naive[i];

					bool is_exposed_sc( !is_bb && (nneigh_res <= nneigh_cut_for_exposed_sc) && pose.residue(i).is_protein() );
					bool is_exposed_bb( is_bb  && (nneigh_res <= nneigh_cut_for_exposed_bb) && pose.residue(i).is_protein() );

					bool accepted_in_stage2( !is_exposed_sc && !is_exposed_bb && accepted_in_stage1 );

					// 1) ensure water binding site is w/i grid
					// Position check2. Vsite it self
					if ( !gridscore->is_point_in_grid(virtual_sites[k]) ) continue;

					// Position check3. Also make a projection of base atom position and see if still possible
					numeric::xyzVector< core::Real > delV = virtual_sites[k] - rsd.xyz(j);
					numeric::xyzVector< core::Real > Vbase = rsd.xyz(j) + 6.0*delV / delV.length(); // is 6 Angstrom projection safe?
					if ( !gridscore->is_point_in_grid(Vbase) ) continue;

					// 2) ensure water binding site does not clash with background
					//core::Real solvscore = gridscore->point_solvation_energy( virtual_sites[k] );
					core::Real score_soft = gridscore->point_energy( virtual_sites[k], atype, 3.0 ); //upweight elec
					core::Real clashscore = gridscore->point_clash_energy( virtual_sites[k], atype, false );

					if ( clashscore > 5.0 ) {
						accepted_in_stage2 = false;
						clashscore = 9.9;
					}

					// just for visualization
					if ( !accepted_in_stage2 ) {
						if ( report ) {
							//core::Real nneigh_res = is_bb?nneigh_res_bb[i]:nneigh_res_sc[i];
							printf ("HETATM %4d %3s  OFF  %4d     %7.3f %7.3f %7.3f %4.1f %6.2f %3d\n",
								int(properties_.size()), aname.c_str(), int(i),
								virtual_sites[k][0], virtual_sites[k][1], virtual_sites[k][2], clashscore, score_soft, int(nneigh_res) );
						}
						continue;
					}

					core::Real ambiguity = (is_OCbb||is_hxl) ? 1.0 : 0.0; // ambiguity of virtual site location
					if ( is_acc ) { // virtual is hydrogen
						atmids_.push_back( core::id::AtomID( j, i ) );
						coords_.push_back( virtual_sites[k] );
						properties_.push_back( AtomProperties(
							false,   // is_don (referring to the virt position)
							false,   // is_acc (referring to the virt position)
							true,    // isH
							true,    // isPolarH
							false,
							false,
							ambiguity,
							1.0,  // score: place holder for now...
							std::to_string(i)+"."+utility::strip(pose.residue(i).atom_name(j)) // take original root atomname
							));

						if ( report ) {
							//core::Real nneigh_res = is_bb?nneigh_res_bb[i]:nneigh_res_sc[i];
							printf ("HETATM %4d %3s  UNK  %4d     %7.3f %7.3f %7.3f %4.1f %6.2f %3d %3.1f\n",
								int(properties_.size()), aname.c_str(), int(i),
								virtual_sites[k][0], virtual_sites[k][1], virtual_sites[k][2], clashscore, score_soft, int(nneigh_res),
								ambiguity);
						}

						Vdonor_sort.push_back( std::make_pair( score_soft, coords_.size() ) );

					} else {  // is_don; skip defining twice with both-donor-and-acceptor such as hxls; those go to is_acc above only
						atmids_.push_back( core::id::AtomID( j, i ) );
						coords_.push_back( virtual_sites[k] );

						properties_.push_back( AtomProperties(
							is_acc,   // is_don (referring to the virt position)
							is_don,   // is_acc (referring to the virt position)
							false,    // isH
							false,    // isPolarH
							false,
							false,
							ambiguity,
							1.0,  // score: place holder for now...
							std::to_string(i)+"."+utility::strip(pose.residue(i).atom_name(j)) // take original root atomname
							));
						if ( report ) {
							core::Real nneigh_res = is_bb?nneigh_res_bb[i]:nneigh_res_sc[i];
							printf ("HETATM %4d %3s  UNK  %4d     %7.3f %7.3f %7.3f %4.1f %6.2f %3d %3.1f\n",
								int(properties_.size()), aname.c_str(), int(i),
								virtual_sites[k][0], virtual_sites[k][1], virtual_sites[k][2], clashscore, score_soft, int(nneigh_res),
								ambiguity);
						}

						Vacceptor_sort.push_back( std::make_pair( score_soft, coords_.size() ) );
					}
				} // virtual sites
			} //
		} // nwaters

		// VSITES from Metal
		if ( pose.residue(i).is_metal() ) append_metal_vsites( pose, i, Vacceptor_sort, report );
	} //total residue
	TR << "Defined " << properties_.size() << " receptor Vsites." << std::endl;
	TR.Debug << "Excluded sidechain positions for the residues exposed: " << rsd_exposed << std::endl;
}

void
ConstraintInfo::append_metal_vsites( core::pose::Pose const &pose,
	core::Size const &ires,
	utility::vector1< std::pair< core::Real, core::Size > > &Vacceptor_sort,
	bool const report
)
{
	// This is a coarse implementation of getting vsites for the uncoordinated site of metal
	// given this low-resolution problem this should be enough
	TR.Debug << "Adding metal vsites from: " << pose.residue(ires).name() << " " << ires << " " << std::endl;

	// first get list of metal-binding residues in pose
	core::conformation::Residue const &metal = pose.residue(ires);
	core::Vector const &xyz_i = metal.xyz(1);

	core::Real const D2CUT( 3.0*3.0 );

	utility::vector1< core::Vector > coordinating_xyzs;
	for ( core::Size jres = 1; jres <= pose.size(); ++jres ) {
		if ( !pose.residue(jres).is_protein() ) continue;

		core::Vector const &xyz_jnbr = pose.residue(jres).nbr_atom_xyz();
		if ( xyz_i.distance_squared( xyz_jnbr ) > 100.0 ) continue;

		// do not use residue.is_metalbinding; it's all aa
		core::chemical::AA const &aa_j = pose.residue(jres).aa();
		if ( aa_j == core::chemical::aa_his || aa_j == core::chemical::aa_cys ||
				aa_j == core::chemical::aa_asp || aa_j == core::chemical::aa_glu ) {

			utility::vector1< core::Size > metal_binding_atoms;
			pose.residue(jres).get_metal_binding_atoms( metal_binding_atoms );
			for ( core::Size j = 1; j <= metal_binding_atoms.size(); ++j ) {
				core::Size const &jatm = metal_binding_atoms[j];
				if ( jatm < pose.residue(jres).first_sidechain_atom() ) continue; //skip backbones

				core::Vector const &xyz_j = pose.residue(jres).xyz(jatm);
				core::Vector dv = xyz_i - xyz_j; //sidechain-to-metal
				if ( dv.length_squared() < D2CUT ) {
					TR.Debug << "Coordination from " << jres << " " << pose.residue(jres).atom_name(jatm) << " " << dv.length_squared() << std::endl;
					coordinating_xyzs.push_back( dv.normalize() ); // push normalized ones
				}
			}
		}
	}

	// Get vector sum
	core::Vector dvsum( 0.0 );
	for ( core::Size icrd = 1; icrd <= coordinating_xyzs.size(); ++icrd ) {
		//TR.Debug << " " << dvsum[0] << " " << dvsum[1] << " " << dvsum[2] << std::endl;
		dvsum += coordinating_xyzs[icrd];
	}

	TR.Debug << "Calculated Net dv: " << coordinating_xyzs.size() << " " << dvsum[0] << " " << dvsum[1] << " " << dvsum[2] << std::endl;

	// works only if summed vector has certain magnitude
	if ( dvsum.length() < 0.5 ) return;

	core::Real ambiguity( 0.0 );
	if ( coordinating_xyzs.size() < 3 ) ambiguity = 2.0; // put ambiguity if not fully coordinated elsewhere

	TR.Debug << "Adding vsite metal with index " << properties_.size()+1 << std::endl;

	core::Vector vcoord = xyz_i + dvsum.normalize()*2.0; // 2.0 Ang apart from metal
	atmids_.push_back( core::id::AtomID( 1, ires ) );
	coords_.push_back( vcoord );

	// metal-coordination site is acceptor!!
	properties_.push_back( AtomProperties( false, //donor
		true, //acceptor
		false, //isH
		false, //polarH
		false, // halogen
		false, //generic
		ambiguity, // ambiguity
		1.0, // score weight; default 1
		+"METAL"+std::to_string(ires) ));

	core::Real const score_soft( -99.0 ); // metals having highest-priority
	Vacceptor_sort.push_back( std::make_pair( score_soft, coords_.size() ) );
	if ( report ) {
		printf ("HETATM %4d %3s  UNK  %4d     %7.3f %7.3f %7.3f %4.1f %6.2f %3d %3.1f\n",
			int(properties_.size()), ("V"+metal.name()).c_str(), int(ires),
			vcoord[0], vcoord[1], vcoord[2], 0.0, score_soft, int(coordinating_xyzs.size()), ambiguity);
	}
}

void
ConstraintInfo::build_Dmtrx()
{
	core::Size const ncrds( coords_.size() );
	Dmtrx_ = utility::vector1< utility::vector1< core::Real > >( ncrds, utility::vector1< core::Real >( ncrds, 0.0 ) );

	for ( core::Size icrd = 1; icrd < ncrds; ++icrd ) {
		for ( core::Size jcrd = icrd+1; jcrd <= ncrds; ++jcrd ) {
			core::Real d = coords_[icrd].distance( coords_[jcrd] );
			d -= std::max(properties_[icrd].ambiguity(), properties_[jcrd].ambiguity() );
			Dmtrx_[icrd][jcrd] = d; Dmtrx_[jcrd][icrd] = d;
		}
	}
}

void
ConstraintInfo::define_ligand_multibody_phores(bool /*report_phore_info*/){
	utility::vector1< Pharmacophore > multibodies;
	utility::vector1< bool > atom_in_multibody( coords_.size(), false );
	// find all cores
	std::map< std::pair<core::Size, core::Size>, utility::vector1< core::Size > > core_maps;
	typedef
		boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS> Graph;
	Graph g(properties_.size());
	for ( core::Size i=1; i <= properties_.size(); ++i ) {
		if ( !properties_[i].used_for_phore() ) continue;
		for ( core::Size j=i+1; j <= properties_.size(); ++j ) {
			if ( !properties_[j].used_for_phore() ) continue;
			if ( Dmtrx_[i][j] < phore_dcut_ ) boost::add_edge(i-1,j-1,g);
		}
	}
	std::vector<int> component(boost::num_vertices(g));
	size_t num = boost::connected_components(g, &component[0]);
	if ( TR.Debug.visible() ) TR.Debug << "Number of core components: " << num << std::endl;
	for ( size_t icom=0; icom<num; icom++ ) {
		utility::vector1<core::Size> atomids;
		for ( size_t ii=0; ii<component.size(); ii++ ) {
			if ( component[ii] != (int)icom ) continue;
			runtime_assert( !atom_in_multibody[ii+1] ); // assert that the atom is not in use
			atomids.push_back(ii+1);
			atom_in_multibody[ii+1]=true;
		}
		if ( atomids.size() >= 3 ) multibodies.push_back( Pharmacophore( atomids, properties_, Dmtrx_, coords_ ) );
	}

	for ( core::Size i=1; i<=multibodies.size(); ++i ) {
		multibody_phores_.push_back( multibodies[i] );
	}
}

void
ConstraintInfo::define_all_ligand_phores( core::Size const minsize,
	bool report_phore_info
) {
	debug_assert( is_ligand_ );
	build_Dmtrx();

	phore_dcut_ = 4.0; // Default is 3.0 A for ligands, 4.0 for receptor; Using 4.0 for both now.
	define_ligand_multibody_phores(report_phore_info);

	utility::vector1< Pharmacophore > twobodies;

	// First make a list of two-body phores
	for ( core::Size i = 1; i < coords_.size(); ++i ) {
		if ( !properties_[i].used_for_phore() ) continue;
		for ( core::Size j = i+1; j <= coords_.size(); ++j ) {
			if ( !properties_[j].used_for_phore() ) continue;
			if ( Dmtrx_[i][j] < phore_dcut_ ) {
				utility::vector1< core::Size > ij; ij.push_back( i ); ij.push_back( j );
				twobodies.push_back( Pharmacophore( ij, properties_, Dmtrx_, coords_ ) );
			}
		}
	}

	// three-body; add if any two-body is expandable, that is, allow multiple 3-bodies sharing same two-body...
	// is this okay though?
	utility::vector1< bool > atm_defined( coords_.size(), false );
	utility::vector1< utility::vector1< bool > > pair_defined( coords_.size(), atm_defined );
	utility::vector1< utility::vector1< utility::vector1< bool > > > triple_defined( coords_.size(), pair_defined );

	// Should we think of better logic for 3body? (clustering based on closest members...)
	for ( core::Size itwo = 1; itwo <= twobodies.size(); ++itwo ) {
		Pharmacophore const &phore_2body = twobodies[itwo];

		for ( core::Size j = 1; j <= coords_.size(); ++j ) {
			if ( phore_2body.has(j) ) continue;
			if ( !properties_[j].used_for_phore() ) continue;

			utility::vector1< core::Size > atms = phore_2body.atms();
			core::Size k( atms[1] ), l( atms[2] );

			// skip if any two are already part of other phore
			if ( triple_defined[k][l][j] || pair_defined[k][l] || pair_defined[k][j] || pair_defined[l][j] ) continue;

			atms.push_back( j );
			bool expandable = phore_2body.is_close( Dmtrx_[j], phore_dcut_ );
			if ( expandable ) {
				atm_defined[k] = atm_defined[l] = atm_defined[j] = true;
				pair_defined[k][l] = pair_defined[l][k] = true;
				pair_defined[j][k] = pair_defined[k][j] = true;
				pair_defined[j][l] = pair_defined[l][j] = true;

				triple_defined[k][l][j] = triple_defined[k][j][l] = true;
				triple_defined[l][k][j] = triple_defined[l][j][k] = true;
				triple_defined[j][k][l] = triple_defined[j][l][k] = true;

				Pharmacophore phore_3body( atms, properties_, Dmtrx_, coords_ );
				/*
				// Finally skip if already close to appended ones
				core::Size jphore( 0 );
				if( phore_overlaps_with_existing( phore_3body, jphore ) ){
				TR.Debug << "EXCLUDE: Ligand phore " << phore_2body.show() << ": close to " << jphore << std::endl;
				}
				*/
				phores_.push_back( phore_3body );
			}
		}
	}

	// add orignal 2body if no 3body was derived from it
	if ( minsize <= 2 ) {
		for ( core::Size itwo = 1; itwo <= twobodies.size(); ++itwo ) {
			Pharmacophore const &phore_2body = twobodies[itwo];
			utility::vector1< core::Size > atms = phore_2body.atms();

			// distance check aginst existing ones before adding
			core::Size jphore( 0 );
			if ( phore_overlaps_with_existing( phore_2body, jphore ) ) continue;

			if ( /*!(atm_defined[atms[1]] || atm_defined[atms[2]] ) &&*/ // check COM instead
					!pair_defined[atms[1]][atms[2]] ) {
				atm_defined[atms[1]] = true; atm_defined[atms[2]] = true;
				phores_.push_back( phore_2body );
			}
		}
	}

	// onebody
	if ( minsize <= 1 ) {
		for ( core::Size i = 1; i <= coords_.size(); ++i ) {
			if ( !properties_[i].used_for_phore() ) continue;
			if ( !atm_defined[i] ) phores_.push_back( Pharmacophore( utility::vector1< core::Size >( 1, i ), properties_, Dmtrx_, coords_ ) );
		}
	}

	// report
	if ( report_phore_info ) {
		TR << "Total " << phores_.size() << " ligand phores defined:" << std::endl;
		for ( core::Size i = 1; i <= phores_.size(); ++i ) {
			TR << "Ligand phore " << i << ": " << phores_[i].show() << std::endl;
		}
	}

	// debuging
	if ( debug_ ) {
		std::string resname = "LIG";

		for ( core::Size iphore = 1; iphore <= phores_.size(); ++iphore ) {
			for ( core::Size imem = 1; imem<= phores_[iphore].atms().size(); ++imem ) {
				core::Size idx( phores_[iphore].atms()[imem] );
				numeric::xyzVector< core::Real > const &xyz = coords_[idx];
				core::Real score = 0.0; // score of seed
				printf ("HETATM %4d VRT  %3s  %4d     %7.3f %7.3f %7.3f      %6.2f\n",
					int(imem), resname.c_str(), int(iphore), xyz[0], xyz[1], xyz[2], score );
			}
		}
	}
}

bool
ConstraintInfo::phore_overlaps_with_existing( Pharmacophore const &phore_i,
	core::Size &jphore ) const
{
	for ( jphore = 1; jphore <= phores_.size(); ++jphore ) {
		core::Real dist = phores_[jphore].com().distance( phore_i.com() );
		if ( dist < 3.0 ) {
			core::Real matchscore = phore_i.find_type_matches( phores_[jphore], true ); // just count
			if ( matchscore == (core::Real)(std::min(phore_i.size(),phores_[jphore].size())) ) {
				return true;
			}
		}
	}
	return false;
}

void
ConstraintInfo::define_receptor_phores( utility::vector1< std::pair< core::Real, core::Size > > &Vdonor_sort,
	utility::vector1< std::pair< core::Real, core::Size > > &Vacceptor_sort,
	core::Size const nmin,
	bool report
)
{
	debug_assert( !is_ligand_ );
	phores_.resize(0);

	// make distance matrix b/w all (optional: drop if apart to anything)
	build_Dmtrx();

	// sort based on energy; treat separately because score scales are different
	std::sort( Vdonor_sort.begin(), Vdonor_sort.end(),
		[&](std::pair<core::Real,core::Size> const &a,std::pair<core::Real,core::Size> const &b)
		{ return a.first < b.first; } );

	std::sort( Vacceptor_sort.begin(), Vacceptor_sort.end(),
		[&](std::pair<core::Real,core::Size> const &a,std::pair<core::Real,core::Size> const &b)
		{ return a.first < b.first; } );

	// None of these values are being used
	//// calculate Zscore
	//core::Real meanD( 0.0 ), sdevD( 1.0e-6 ), meanA( 0.0 ), sdevA( 1.0e-6 );
	//if ( Vdonor_sort.size() > 1 ) {
	// for ( auto a = Vdonor_sort.begin(); a != Vdonor_sort.end(); ++a ) meanD += a->first;
	// meanD /= Vdonor_sort.size();
	// for ( auto a = Vdonor_sort.begin(); a != Vdonor_sort.end(); ++a ) sdevD += (a->first-meanD)*(a->first-meanD);
	// sdevD /= Vdonor_sort.size(); sdevD = std::sqrt(sdevD);
	//}

	//if ( Vacceptor_sort.size() > 1 ) {
	// for ( auto a = Vacceptor_sort.begin(); a != Vacceptor_sort.end(); ++a ) meanA += a->first;
	// meanA /= Vacceptor_sort.size();
	// for ( auto a = Vacceptor_sort.begin(); a != Vacceptor_sort.end(); ++a ) sdevA += (a->first-meanA)*(a->first-meanA);
	// sdevA /= Vacceptor_sort.size(); sdevA = std::sqrt(sdevA);
	//}

	// put into a shared store; set inverse Zscore so that it ranges +0.5~+1.0 & more positive the better
	utility::vector1< std::pair< core::Real, core::Size > > Vtotal_sort;
	core::Size i( 0 );
	for ( auto a = Vdonor_sort.begin(); a != Vdonor_sort.end(); ++a, ++i ) {
		core::Real score = 1.0 - core::Real(i)/Vdonor_sort.size();
		Vtotal_sort.push_back( std::make_pair( score, a->second ) );
	}
	i = 0;
	for ( auto a = Vacceptor_sort.begin(); a != Vacceptor_sort.end(); ++a, ++i ) {
		core::Real score = 1.0 - core::Real(i)/Vacceptor_sort.size();
		Vtotal_sort.push_back( std::make_pair( score, a->second ) );
	}

	// Descending order as scores are positive
	std::sort( Vtotal_sort.begin(), Vtotal_sort.end(),
		[&](std::pair<core::Real,core::Size> const &a,std::pair<core::Real,core::Size> const &b)
		{ return a.first > b.first; } );

	core::Size const ncrds( coords_.size() );

	// update AtomProperty scores by Zscore
	// property score is used whichever method below
	for ( core::Size i = 1; i <= ncrds; ++i ) { // sorted
		core::Real score = std::max(0.5,Vtotal_sort[i].first);
		core::Size icrd = Vtotal_sort[i].second;
		properties_[icrd].score( score );
	}

	//bool by_cluster( true );
	//if ( !by_cluster ) {
	//define_all_phores( 2, report ); // CAREFUL! This function has been changed to only run in ligand context

	//} else {
	utility::vector1< utility::vector1< core::Size > > clusters;
	utility::vector1< core::Size > nused( ncrds, 0 );

	for ( core::Size i = 1; i <= ncrds; ++i ) { // priority
		core::Size icrd = Vtotal_sort[i].second;
		if ( nused[icrd] > 3 ) continue; // do not use as a seed if taken multiple times by other higher priority clusters

		// First, sort others by distance
		utility::vector1< core::Size > cluster_i; cluster_i.push_back( icrd ); //start with seed
		utility::vector1< std::pair< core::Real, core::Size > > distance_sort;
		for ( core::Size jcrd = 1; jcrd <= ncrds; ++jcrd ) {
			if ( icrd == jcrd ) continue;
			//if( Dmtrx_[icrd][jcrd] < phore_dcut_ ) cluster_i.push_back( jcrd );
			distance_sort.push_back( std::make_pair( Dmtrx_[icrd][jcrd], jcrd ) );
		}

		// Start from closest, add until self-consistent
		std::sort( distance_sort.begin(), distance_sort.end(),
			[&](std::pair<core::Real,core::Size> const &a,std::pair<core::Real,core::Size> const &b)
			{ return a.first < b.first; } );

		for ( auto it = distance_sort.begin(); it != distance_sort.end(); ++it ) {
			if ( it->first > phore_dcut_ ) break; // distance to seed

			// distance to all other members
			bool close_to_all( true );
			for ( core::Size imem = 1; imem <= cluster_i.size(); ++imem ) {
				if ( Dmtrx_[cluster_i[imem]][it->second] > phore_dcut_ ) {
					close_to_all = false;
					break;
				}
			}
			if ( close_to_all ) cluster_i.push_back(it->second);
		}

		// minimum size as 2
		if ( cluster_i.size() < 2 ) continue;
		//if( cluster_overlaps_any( clusters, cluster_i, 3, 1.0 ) ) continue; // more than three or half the same from either side

		clusters.push_back( cluster_i );
		for ( core::Size imem = 1; imem <= cluster_i.size(); ++imem ) {
			nused[ cluster_i[imem] ] ++;
		}
	}

	// add more clusters with singletons if number not enough
	if ( clusters.size() < nmin ) {
		for ( core::Size i = 1; i <= ncrds; ++i ) { // priority
			core::Size icrd = Vtotal_sort[i].second;
			//core::Real score = Vtotal_sort[i].first; // score of seed

			if ( nused[icrd] > 0 ) continue;
			utility::vector1< core::Size > cluster_i; cluster_i.push_back( icrd ); //start with seed
			clusters.push_back( cluster_i );
			if ( clusters.size() >= nmin ) break;
		}
	}

	// remove if too close to other higher priority clusters
	for ( core::Size icl = 1; icl <= clusters.size(); ++icl ) {
		utility::vector1< core::Size > const &cluster_i = clusters[icl];
		Pharmacophore phore_i( cluster_i, properties_, Dmtrx_, coords_ );

		core::Size jcl( 0 );
		bool overlaps = phore_overlaps_with_existing( phore_i, jcl );

		std::string resname( "CL " );
		if ( overlaps ) {
			resname = "OVE";
		} else {
			phores_.push_back( phore_i );
		}

		// debuging
		if ( debug_ && report ) {
			if ( overlaps ) { TR.Debug << "Excluded " << icl << ", close to " << jcl << ": "; }
			else { TR.Debug << "Cluster  " << icl << ":"; }

			for ( core::Size imem = 1; imem <= cluster_i.size(); ++imem ) TR.Debug << " " << cluster_i[imem];
			TR.Debug << std::endl;

			for ( core::Size imem = 1; imem <= cluster_i.size(); ++imem ) {
				core::Size idx( cluster_i[imem] );
				numeric::xyzVector< core::Real > const &xyz = coords_[idx];
				core::Real score = Vtotal_sort[icl].first; // score of seed
				printf ("HETATM %4d VRT  %3s  %4d     %7.3f %7.3f %7.3f      %6.2f\n",
					int(imem), resname.c_str(), int(icl), xyz[0], xyz[1], xyz[2], score );
			}
		}
	}
	TR << "Made total " << phores_.size() << " phores from " << clusters.size() << " clusters." << std::endl;
}

// SHOULD be called from ligand side
void
ConstraintInfo::map_phores( utility::vector1< Pharmacophore > const &receptor_phores,
	bool const update,
	core::Size const nmax
)
{
	if ( phore_match_.size() > 0 && !update ) return;

	debug_assert( is_ligand_ );
	phore_match_.resize( 0 );

	utility::vector1< std::pair< core::Real, std::pair< core::Size, Pharmacophore > > > phores_sort_rest; // total
	utility::vector1< bool > added_from_ligand_phore( phores_.size(), false );

	// iter through motifs in ligand
	for ( core::Size i = 1; i <= phores_.size(); ++i ) {
		utility::vector1< std::pair< core::Real, Pharmacophore > > phores_sort; // sort by each ligand phore
		Pharmacophore const &phore_lig = phores_[i];

		// iter over receptor phores
		for ( core::Size j = 1; j <= receptor_phores.size(); ++j ) {
			Pharmacophore phore_rec( receptor_phores[j] ); // make as a copy instead
			utility::vector1< core::Size > map_index;
			core::Real score( 1.0e6 );

			// just by type counts
			score = phore_lig.find_type_matches( phore_rec );
			if ( score == 0.0 ) continue; // no type match
			phores_sort.push_back( std::make_pair( -score, phore_rec ) ); // lower the better

			// more complicated; use exact atom mapping
			/*
			if( phore_lig.find_matching_permutation( phore_rec, map_index, score ) ){
			phore_rec.set_map_index( map_index ); // store info
			phores_sort.push_back( std::make_pair( score, phore_rec ) );
			}
			*/
		}

		if ( phores_sort.size() == 0 ) continue;

		// take top N from each ligand site; lower the better
		std::sort( phores_sort.begin(), phores_sort.end(),
			[&](std::pair<core::Real,Pharmacophore> const &a,std::pair<core::Real,Pharmacophore> const &b)
			{ return a.first < b.first; } );
		if ( phores_sort.size() > n_match_per_ligand_motif_ ) phores_sort.resize( n_match_per_ligand_motif_ );

		for ( auto it = phores_sort.begin(); it != phores_sort.end(); ++it ) {
			if ( !added_from_ligand_phore[i] ) { // add at least on phore match per ligand-phore
				added_from_ligand_phore[i] = true;
				phore_match_.push_back( std::make_pair( i, it->second ) );
				if ( update && debug_ ) {
					TR << "Match on ligand phore " << phores_[i].show() << " (" << i << "/ " << phores_.size() << "): "
						<<  it->second.show() << ", score " << it->first
						<< ", total " << phore_match_.size() << std::endl;
				}
			} else {
				phores_sort_rest.push_back( std::make_pair( it->first, std::make_pair( i, it->second ) ) );
			}
		}
	}

	// Then add rest by ligand score
	std::sort( phores_sort_rest.begin(), phores_sort_rest.end(),
		[&](std::pair< core::Real, std::pair<core::Size,Pharmacophore> > const &a,
		std::pair< core::Real, std::pair<core::Size,Pharmacophore> > const &b)
		{ return a.first < b.first; } );

	// add phores based on score
	for ( auto it = phores_sort_rest.begin(); it != phores_sort_rest.end(); ++it ) {
		phore_match_.push_back( it->second );

		if ( update && debug_ ) {
			std::pair<core::Size,Pharmacophore> const &match = it->second;

			core::Size ilig = match.first;
			Pharmacophore const &phore_rec = match.second;
			TR << "Match on ligand phore " << phores_[ilig].show() << " (" << ilig << "/ " << phores_.size() << "): "
				<<  phore_rec.show() << ", score " << it->first
				<< ", total " << phore_match_.size() << std::endl;
		}

		if ( phore_match_.size() > nmax ) break;
	}


	if ( update ) {
		TR << "Found " << phore_match_.size() << " ligand-receptor phore matches." << std::endl;
	}
}

void
ConstraintInfo::select_phore_match( core::Size const run_index )
{
	debug_assert( run_index > 0 );
	core::Size imatch = (run_index-1)%phore_match_.size()+1; // enumerate through matches
	current_phore_match_ = phore_match_[imatch];

	Pharmacophore lig_phore = phores_[ current_phore_match_.first ]; // saved in index
	Pharmacophore const &rec_phore = current_phore_match_.second;

	// reporting purpose
	TR << "Istruct " << run_index << ", apply pharmacophore match "
		<< imatch << " (" << lig_phore.index() << "<->" << rec_phore.index() << "): "
		<< lig_phore.show() << " <-> " << rec_phore.show() << std::endl;
}

// SHOULD be called from ligand side
void
ConstraintInfo::align_to_current_phore_match( core::pose::Pose &pose,
	ConstraintInfo const &, //recinfo,
	utility::vector1< core::Size > const &ligids,
	utility::vector1< std::pair< core::Size, core::Size > > &, //marked_pairs,
	utility::vector1< core::Size > &SrcPriorIDs,
	utility::vector1< core::Size > &TgtPriorIDs
)
{
	debug_assert( is_ligand_ );

	// update coord in case not synced with pose (like randomize)
	update_ligand_coord( pose );

	Pharmacophore lig_phore = phores_[ current_phore_match_.first ]; // saved in index
	Pharmacophore const &rec_phore = current_phore_match_.second;

	lig_phore.update_com( coords_ );

	// Put bonus b/w these atoms at stage1 matching
	SrcPriorIDs = lig_phore.atms();
	TgtPriorIDs = rec_phore.atms();

	numeric::xyzVector< core::Real > dcom = rec_phore.com() - lig_phore.com();

	// add slight noise
	dcom += numeric::random::random_translation<core::Real>( 0.5, numeric::random::rg() );

	// just translate to match phores keeping rotation
	for ( auto ligid : ligids ) {
		for ( core::Size iatm=1; iatm<=pose.residue(ligid).natoms(); ++iatm ) {
			pose.set_xyz( core::id::AtomID(iatm,ligid),
				pose.xyz(core::id::AtomID(iatm,ligid)) + dcom );
		}
	}
}

void
ConstraintInfo::update_ligand_coord( core::pose::Pose const & pose )
{
	coords_.resize(0);
	for ( core::Size i = 1; i <= atmids_.size(); ++i ) {
		coords_.push_back( pose.xyz( atmids_[i] ) );
	}
}

//////////////////////////////////////////////////////////////////
LigandAligner::LigandAligner() = default;

LigandAligner::LigandAligner( bool use_pharmacophore,
	utility::vector1< core::Size > const &movable_scs,
	bool faster ) :
	movable_scs_( movable_scs ),
	use_pharmacophore_( use_pharmacophore ),
	faster_( faster )
{}

void
LigandAligner::set_pharmacophore_reference( core::pose::Pose const &pose )
{
	// this part is sometimes confusing because of multiple way of constructions
	// but let's keep it for now...
	target_ = ConstraintInfo( pose, sf_, faster_ );
}

void
LigandAligner::set_constraints(
	core::pose::Pose & pose,
	utility::vector1<core::Size> ligids,
	utility::vector1< std::pair< core::Size, core::Size > > &marked_pairs,
	core::Real const w_prior,
	utility::vector1< core::Size > const &SrcPriorIDs,
	utility::vector1< core::Size > const &TgtPriorIDs
)
{
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	core::Real MAXSCORE=1e6;
	ConstraintSetOP cst_set(new ConstraintSet);
	core::Size ftroot = pose.fold_tree().root();

	ConstraintInfo cst_source( pose, ligids, use_pharmacophore() );

	core::Size nSrcAtms = cst_source.natoms();//=pose.residue(ligid).natoms();
	core::Size nTgtAtms = target_.natoms();

	if ( nSrcAtms == 0 || nTgtAtms == 0 ) return;

	// populate a full M x N matrix of tgt - src distances
	utility::vector1< utility::vector1< core::Real > > score_ij(
		nSrcAtms, utility::vector1< core::Real >(nTgtAtms, MAXSCORE));

	for ( core::Size i=1; i<=nSrcAtms; ++i ) {
		for ( core::Size j=1; j<=nTgtAtms; ++j ) {
			core::Real weight_ij = cst_source.properties(i).match(
				target_.properties(j), target_.polar_scale(), true );
			if ( weight_ij == 0.0 ) continue;

			// sqrt won't be the time determining calculation
			core::Real dist = (cst_source.coord(i) - target_.coord(j) ).length();
			core::Real ambiguity = std::max(cst_source.properties(i).ambiguity(),target_.properties(j).ambiguity());
			dist = (dist > ambiguity)? dist-ambiguity : 0.0;

			score_ij[i][j] = dist;
		}
	}

	// now greedily find the best unassigned target atom / unassigned source atom
	utility::vector1< bool > src_marked( nSrcAtms, false );
	utility::vector1< bool > tgt_marked( nTgtAtms, false );
	marked_pairs.resize( 0 );

	bool done = false;
	while ( !done ) {
		done = true;

		core::Real mindist=MAXSCORE;
		numeric::xyzVector<core::Real> tgtAtm;
		core::Size min_idx_i=0, min_idx_j=0;

		for ( core::Size i=1; i<=nSrcAtms; ++i ) {
			if ( src_marked[i] ) continue;
			for ( core::Size j=1; j<=nTgtAtms; ++j ) {
				if ( tgt_marked[j] ) continue;
				if ( score_ij[i][j] < mindist ) {
					mindist = score_ij[i][j];
					min_idx_i = i;
					min_idx_j = j;
				}
			}
		}
		bool is_prior_pair = (SrcPriorIDs.contains(min_idx_i) && TgtPriorIDs.contains(min_idx_j));

		// we have our best correspondence
		if ( mindist < MAXSCORE ) {
			core::Real weight_i = target_.weight();
			core::Real limit_i = target_.limit();
			core::Real atomtype_scale = cst_source.properties(min_idx_i).match(
				target_.properties(min_idx_j), target_.polar_scale(), true );
			weight_i *= atomtype_scale;
			if ( is_prior_pair ) weight_i *= w_prior;

			SumFuncOP cstfunc(new SumFunc);
			cstfunc->add_func( FuncOP( new ConstantFunc(-weight_i*limit_i*limit_i) ) );
			cstfunc->add_func( FuncOP( new TopOutFunc(weight_i,0.0,limit_i) ) );

			ConstraintCOP newCST(new CoordinateConstraint(
				cst_source.atomid(min_idx_i),   //core::id::AtomID(min_idx_i, ligid),
				core::id::AtomID(1, ftroot), // safe! we're guaranteed to have jump upstream of ligand
				target_.coord(min_idx_j),
				cstfunc ) );
			//newCST->show_def( TR, pose );
			cst_set->add_constraint( newCST );

			src_marked[min_idx_i] = true;
			tgt_marked[min_idx_j] = true;
			// order of src/tgt
			marked_pairs.push_back( std::make_pair( min_idx_i, min_idx_j ) );
			done = false;
		}
		//if( marked_pairs.size() >= nmax_match ) break;
	}

	pose.constraint_set( cst_set );

}

void
LigandAligner::set_hard_constraint_on_marked(
	core::pose::Pose & pose,
	utility::vector1< core::Size > ligids,
	utility::vector1< std::pair< core::Size, core::Size > > const &marked_pairs
) const
{
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	core::Size ftroot = pose.fold_tree().root();
	ConstraintInfo cst_source( pose, ligids, use_pharmacophore() );

	ConstraintSetOP cst_set(new ConstraintSet);
	for ( core::Size ipair = 1; ipair <= marked_pairs.size(); ++ipair ) {
		core::Size i_src = marked_pairs[ipair].first;
		core::Size i_tgt = marked_pairs[ipair].second;

		FuncOP cstfunc( new HarmonicFunc( 0.0, 1.0 ) );
		ConstraintCOP newCST(new CoordinateConstraint(
			cst_source.atomid(i_src),
			core::id::AtomID(1, ftroot), // safe! we're guaranteed to have jump upstream of ligand
			target_.coord(i_tgt),
			cstfunc ) );
		//newCST->show_def( TR, pose );
		cst_set->add_constraint( newCST );
	}

	pose.constraint_set( cst_set );
}

// randomize
void
LigandAligner::randomize_lig(
	core::pose::Pose & pose,
	utility::vector1< core::Size > ligids,
	numeric::xyzVector<core::Real> const &T
) {

	for ( auto ligid : ligids ) {
		core::Size nchi = pose.residue_type(ligid).nchi();
		for ( core::Size ichi_src=1; ichi_src<=nchi; ++ichi_src ) {
			core::Real angle_i = 360.0 * numeric::random::rg().uniform();

			core::chemical::AtomIndices const & chiatoms = pose.residue_type(ligid).chi_atoms( ichi_src );
			core::chemical::BondName bondtype( pose.residue_type(ligid).bond_type( chiatoms[2], chiatoms[3] ) );

			if ( bondtype == core::chemical::DoubleBond ) {
				if ( angle_i < 90.0 || angle_i > 270.0 ) {
					angle_i = 0.0;
				} else {
					angle_i = 180.0;
				}
			}
			pose.set_chi( ichi_src, ligid, angle_i );
		}
		if ( pose.residue_type(ligid).is_protein() ) {
			core::Real phi_i = 360.0 * numeric::random::rg().uniform();
			core::Real psi_i = 360.0 * numeric::random::rg().uniform();
			core::Real omega_i = numeric::random::rg().uniform() > 0.5 ? 0.0 : 180.0;
			pose.set_phi( ligid, phi_i );
			pose.set_psi( ligid, psi_i );
			pose.set_omega( ligid, omega_i );
		}
	}

	core::Vector Tcom( 0.0,0.0,0.0 );
	core::Size natm = 0;
	for ( auto ligid : ligids ) {
		for ( core::Size iatm_src=1; iatm_src<=pose.residue_type(ligid).natoms(); ++iatm_src ) {
			Tcom += pose.xyz(core::id::AtomID(iatm_src,ligid));
		}
		natm += pose.residue_type(ligid).natoms();
	}
	Tcom /= natm;

	core::Real len = trans_step_ * numeric::random::rg().uniform();
	core::Vector Toff( len*numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );

	numeric::xyzMatrix< core::Real > R=numeric::random::random_rotation();
	for ( auto ligid : ligids ) {
		for ( core::Size iatm_src=1; iatm_src<=pose.residue_type(ligid).natoms(); ++iatm_src ) {
			pose.set_xyz(
				core::id::AtomID(iatm_src,ligid),
				R*( pose.xyz(core::id::AtomID(iatm_src,ligid)) - Tcom) + T + Toff
			);
		}
	}
}


// perturb
void
LigandAligner::perturb_lig(
	core::pose::Pose & pose,
	utility::vector1< core::Size > ligids
) {
	core::Vector Raxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
	core::Real angle = rot_step_ * numeric::NumericTraits<core::Real>::deg2rad() * numeric::random::rg().gaussian();
	numeric::xyzMatrix< core::Real > R = numeric::rotation_matrix( Raxis, angle );

	core::Vector Tcom( 0.0,0.0,0.0 );
	core::Size natm = 0;
	for ( auto ligid : ligids ) {
		for ( core::Size iatm_src=1; iatm_src<=pose.residue_type(ligid).natoms(); ++iatm_src ) {
			Tcom += pose.xyz(core::id::AtomID(iatm_src,ligid));
		}
		natm += pose.residue_type(ligid).natoms();
	}
	Tcom /= natm;

	core::Real len = trans_step_ * numeric::random::rg().uniform();
	core::Vector T( len*numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
	for ( auto ligid : ligids ) {
		for ( core::Size iatm_src=1; iatm_src<=pose.residue_type(ligid).natoms(); ++iatm_src ) {
			pose.set_xyz(
				core::id::AtomID(iatm_src,ligid),
				R*( pose.xyz(core::id::AtomID(iatm_src,ligid)) - Tcom) + T + Tcom
			);
		}
	}

	// internal torsions
	for ( auto ligid : ligids ) {
		core::Size nchi = pose.residue_type(ligid).nchi();
		for ( core::Size ichi_src=1; ichi_src<=nchi; ++ichi_src ) {
			core::chemical::AtomIndices const & chiatoms = pose.residue_type(ligid).chi_atoms( ichi_src );
			core::chemical::BondName bondtype( pose.residue_type(ligid).bond_type( chiatoms[2], chiatoms[3] ) );

			core::Real angle_i = pose.chi( ichi_src, ligid );
			if ( bondtype == core::chemical::DoubleBond ) {
				// keep amide bonds at 0 or 180
				if ( numeric::random::rg().uniform() < (chi_step_ / 180.0) ) {
					angle_i = std::fmod( angle_i + 180.0, 360.0);
				}
			} else {
				angle_i = std::fmod( angle_i + chi_step_ * numeric::random::rg().gaussian(), 360.0);
			}
			pose.set_chi( ichi_src, ligid, angle_i );
		}
		if ( pose.residue_type(ligid).is_protein() ) {
			core::Real phi_i = std::fmod( pose.phi(ligid) + chi_step_ * numeric::random::rg().gaussian(), 360.0);
			core::Real psi_i = std::fmod( pose.psi(ligid) + chi_step_ * numeric::random::rg().gaussian(), 360.0);
			core::Real omega_i = pose.omega(ligid);
			if ( numeric::random::rg().uniform() < (chi_step_ / 180.0) ) {
				omega_i = std::fmod( omega_i + 0.1 * chi_step_ * numeric::random::rg().gaussian(), 360.0);
			}
			pose.set_phi( ligid, phi_i );
			pose.set_psi( ligid, psi_i );
			pose.set_omega( ligid, omega_i );
		}
	}
}


core::Size
LigandAligner::estimate_nstruct_sample(
	core::pose::Pose const &pose,
	utility::vector1< core::Size > const &ligids,
	core::Size const ntotal
) {
	debug_assert( use_pharmacophore_ ); // make sure this called only with PHdock

	// Initialize and report ligand phores here
	// needs to call this function everytime ligand type gets updated

	// TODO
	// Is this part time consuming? cst_source info is temporary here and needs to be recalculated in apply function
	// if then, move this to apply function...
	ConstraintInfo cst_source( pose, ligids, true, true );

	// initialize phore matching
	cst_source.map_phores( target_.phores(), true ); //force update

	core::Size n_phore_matches = cst_source.n_phore_match();
	core::Size nsample( 0 );

	// Rule:
	// Nsites_receptor: 25~40
	// Nsites_lig: 2~10
	// Nmatches ranges within 30~100; nsample should reasonably cover this number to enumerate most of the matches
	// Can think of reducing Nsites in rec or in lig, or speeding up per-match speed

	// Do not exceed 70% of rest of the pool or n=70... too slow with 1sec-per-struct
	if ( n_phore_matches == 0 ) {
		return 0;
	} else if ( n_phore_matches <= 50 ) {
		nsample = 50; // 50 struct takes ~1 minute
		TR << "Small num. of n_phore_matches (<50); Nautogen set by ";
	} else {
		nsample = std::min( int(n_phore_matches), 70 );
		TR << "Large num. of n_phore_matches (" << n_phore_matches << "); Nautogen set by ";
	}
	if ( nsample > (core::Size)(0.7*ntotal) ) nsample = (core::Size)(0.7*ntotal);
	TR << nsample << std::endl;
	return nsample;
}

//skeletonizes the density and stores the skeleton in points_to_search_
void
LigandAligner::select_points( core::pose::Pose const & pose, core::Size const ligid, core::Real skeleton_threshold_const, core::Size neighborhood_size ) {

	ObjexxFCL::FArray3D< float > const & densdata = core::scoring::electron_density::getDensityMap().get_data();

	//points_to_search_.clear();

	//calculate center of mass for ligand; probably a Rosetta function that does this
	core::conformation::Residue const &lig( pose.residue( ligid ) );
	core::Vector lig_com( 0.0 );
	for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
		lig_com += lig.xyz( iatm );
	}
	lig_com /= lig.nheavyatoms();
	TR.Debug << "Ligand center of mass: " << lig_com[0] << " " << lig_com[1] << " " << lig_com[2] << std::endl;

	// sort points and store if within 10 units of ligand center of mass
	core::Real distance_to_com = 0;
	utility::vector1< std::pair< numeric::xyzVector< core::Real >, core::Real > > point_score_pairs; //voxel location and density value
	for ( int z=1; z<=(int)densdata.u3(); z+=gridStep_ ) {
		for ( int y=1; y<=(int)densdata.u2(); y+=gridStep_ ) {
			for ( int x=1; x<=(int)densdata.u1(); x+=gridStep_ ) {
				numeric::xyzVector< core::Real > x_idx(x,y,z);
				numeric::xyzVector< core::Real > x_cart( x_idx[0], x_idx[1], x_idx[2] );
				core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
				core::Real dens_value = -densdata(x,y,z);
				std::pair< numeric::xyzVector< core::Real >, core::Real > point_score_pair = std::make_pair(x_idx, dens_value);
				distance_to_com = sqrt( pow( x_cart[0] - lig_com[0], 2 ) + pow( x_cart[1] - lig_com[1], 2 ) + pow( x_cart[2] - lig_com[2], 2 ) );

				if ( distance_to_com <= 10.0 ) { //default 10.0
					point_score_pairs.push_back(point_score_pair);
				}
			}
		}
	}

	std::sort(point_score_pairs.begin(), point_score_pairs.end(), PointScoreComparator()); //sorts by density value
	utility::vector1< std::pair< numeric::xyzVector< core::Real >, core::Real > > points_to_erode;
	double voxel_volume = core::scoring::electron_density::getDensityMap().get_voxel_volume();

	TR.Debug << "Number of heavy atoms: " << lig.nheavyatoms() << std::endl;
	TR.Debug << "Voxel volume: " << voxel_volume << std::endl;

	for ( core::Size i=1; i<=point_score_pairs.size(); i++ ) {
		//bool hasneighbor = false;
		numeric::xyzVector< core::Real > x_idx = point_score_pairs[i].first;
		numeric::xyzVector< core::Real > x_cart(x_idx[0],x_idx[1],x_idx[2]);
		core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );

		//visualize erosion process; outputs all points in PDB format
		TR.Debug << "HETATM" << std::setw(5) << utility::to_string(i) << " " <<  std::setw(4) << std::left << "O1" << " " <<  std::setw(3) << std::right << "ALA" << " A" << std::setw(4) << utility::to_string(i) << "    " << std::setw(8) << x_cart[0] << std::setw(8) << x_cart[1] << std::setw(8) << x_cart[2] << std::setw(6) << "1.00" << std::setw(6) << std::fixed << std::setprecision(2) << point_score_pairs[i].second << std::endl;

		//finds points that are not occupied and saves them to search
		if ( !check_voxel_distance_to_receptor( x_cart, pose, ligid ) && points_to_erode.size() < skeleton_threshold_const * lig.nheavyatoms() / voxel_volume ) {
			points_to_erode.push_back( point_score_pairs[i] );
		}

		//if ( points_to_search_.size() >= topNtrans_ ) break; //I have no idea what this does. A remnant from the code I ripped off
	}

	//prints points before erosion could go inside debug condition
	for ( core::Size i = 1; i <= points_to_erode.size(); i++ ) {
		numeric::xyzVector< core::Real > x_idx = points_to_erode[i].first;
		numeric::xyzVector< core::Real > x_cart(x_idx[0],x_idx[1],x_idx[2]);
		core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );

		TR.Debug << "HETATM" << std::setw(5) << utility::to_string(i) << " " <<  std::setw(4) << std::left << "O1" << " " <<  std::setw(3) << std::right << "ALA" << " A" << std::setw(4) << utility::to_string(i) << "    " << std::setw(8) << x_cart[0] << std::setw(8) << x_cart[1] << std::setw(8) << x_cart[2] << std::setw(6) << "1.00" << std::setw(6) << std::fixed << std::setprecision(2) << "20.00" << std::endl;

	}

	//first pass erosion with small neighborhood
	utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > eroded_points = erode_points( points_to_erode, neighborhood_size );
	//gets largest connection of points
	eroded_points = find_biggest_skeleton( eroded_points );
	//erodes again with larger neighborhood
	eroded_points = erode_points( eroded_points, 19 );

	TR.Debug << "Done eroding" << std::endl;
	for ( core::Size i = 1; i <= eroded_points.size(); i++ ) {
		numeric::xyzVector< core::Real > x_idx = eroded_points[i].first;
		numeric::xyzVector< core::Real > x_cart(x_idx[0],x_idx[1],x_idx[2]);
		core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );

		points_to_search_.push_back( x_cart );

	}

}

//erosion algorithm; if a voxel contains any non-empty neighbors, it is removed;
//neighbors are added back if after all points are removed, they have an empty neighborhood
utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > >
LigandAligner::erode_points( utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > selected_points, core::Size neighborhood_size ) {
	utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > removed_points;
	core::Size begin_size = selected_points.size();

	for ( core::Size i = 1; i <= selected_points.size(); i++ ) {
		numeric::xyzVector< core::Real > voxel = selected_points[i].first;
		std::map< core::Size, numeric::xyzVector< core::Real > > neighbor_points = assign_neighbors( voxel, neighborhood_size );
		for ( core::Size j = 1; j <= neighborhood_size; j++ ) {

			//remove point from skeleton if one of its neighbors is empty
			if ( !is_point_in_search_group( selected_points, neighbor_points[j] ) ) {
				removed_points.push_back( selected_points[i] );
				selected_points.erase( selected_points.begin() + (i-1) );
				j = 28;
			}

		}

	}

	for ( core::Size ipoint = 1; ipoint <= removed_points.size(); ipoint++ ) {
		bool is_in_selected_points = false;
		numeric::xyzVector< core::Real > voxel = removed_points[ipoint].first;
		std::map< core::Size, numeric::xyzVector< core::Real > > neighbor_points = assign_neighbors( voxel, neighborhood_size );

		for ( core::Size j = 1; j <= neighborhood_size; j++ ) {
			if ( is_point_in_search_group ( selected_points, neighbor_points[j] ) ) { //if neighboring point is in the skeleton
				is_in_selected_points = true;
			}
		}

		if ( !is_in_selected_points ) { //if a voxels neighbors are empty, add it back
			selected_points.push_back( removed_points[ipoint] );
		}
	}

	//prints skeleton after each round of erosion
	TR.Debug << "Printing select_points" << std::endl;
	for ( core::Size i = 1; i <= selected_points.size(); ++i ) {
		numeric::xyzVector< core::Real > voxel = selected_points[i].first;
		numeric::xyzVector< core::Real > x_cart(voxel[0],voxel[1],voxel[2]);
		core::scoring::electron_density::getDensityMap().idx2cart( voxel, x_cart );

		TR.Debug << "HETATM" << std::setw(5) << utility::to_string(i) << " " <<  std::setw(4) << std::left << "O1" << " " <<  std::setw(3) << std::right << "ALA" << " A" << std::setw(4) << utility::to_string(i) << "    " << std::setw(8) << x_cart[0] << std::setw(8) << x_cart[1] << std::setw(8) << x_cart[2] << std::setw(6) << "1.00" << std::setw(6) << std::fixed << std::setprecision(2) << "20.00" << std::endl;

	}

	//if no erosion occurs, return the fully eroded skeleton
	if ( begin_size == selected_points.size() ) {
		return selected_points;
	} else { //else erode again
		return erode_points( selected_points, neighborhood_size );
	}
}

//creates networks of points and selects the largest one
utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > >
LigandAligner::find_biggest_skeleton( utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > eroded_points ) {

	std::sort( eroded_points.begin(), eroded_points.end(), PointScoreComparator());
	utility::vector1 < utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > > networks;

	for ( core::Size ipoint = 1; ipoint <= eroded_points.size(); ++ipoint ) {
		//if a point is not in a network, make a new network and find all the points in the network
		if ( !(is_point_in_network( eroded_points[ipoint].first, networks ) ) ) {
			utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > network;
			network.push_back( eroded_points[ipoint] );

			network = find_network( eroded_points[ipoint].first, network, eroded_points );
			networks.push_back( network );
		}
	}

	utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > largest_network;
	core::Size largest_size = 0;

	for ( core::Size inetwork = 1; inetwork <= networks.size(); ++inetwork ) {
		TR.Debug << "Network " << inetwork << ":" << std::endl;
		for ( core ::Size ipoint = 1; ipoint <= networks[inetwork].size(); ++ipoint ) {
			numeric::xyzVector < core::Real > point_to_print = networks[inetwork][ipoint].first;
			numeric::xyzVector< core::Real > x_cart(point_to_print[0],point_to_print[1],point_to_print[2]);
			core::scoring::electron_density::getDensityMap().idx2cart( point_to_print, x_cart );

			//prints out networks
			TR.Debug << "HETATM" << std::setw(5) << utility::to_string(ipoint) << " " <<  std::setw(4) << std::left << "O1" << " " <<  std::setw(3) << std::right << "ALA" << " A" << std::setw(4) << utility::to_string(ipoint) << "    " << std::setw(8) << x_cart[0] << std::setw(8) << x_cart[1] << std::setw(8) << x_cart[2] << std::setw(6) << "1.00" << std::setw(6) << std::fixed << std::setprecision(2) << "20.00" << std::endl;
		}

		if ( networks[inetwork].size() > largest_size ) {
			largest_network = networks[inetwork];
			largest_size = largest_network.size();
		}
	}

	return largest_network;

}

//goes through all points in all networks and sees if a point belongs to a network
bool LigandAligner::is_point_in_network ( numeric::xyzVector< core::Real > point, utility::vector1 < utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > > networks ) {
	for ( core::Size inetwork = 1; inetwork <= networks.size(); ++inetwork ) {
		for ( core::Size ipoint = 1; ipoint <= networks[inetwork].size(); ++ipoint ) {
			numeric::xyzVector < core::Real > point_from_network = networks[inetwork][ipoint].first;
			if ( point_from_network[0] - point[0] == 0 && point_from_network[1] - point[1] == 0 && point_from_network[2] - point[2] == 0 ) {
				return true;
			}
		}
	}

	return false;
}

//creates the network of points from a single head point
utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > >
LigandAligner::find_network ( numeric::xyzVector< core::Real > start_point, utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > network, utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > eroded_points ) {

	for ( core::Size ipoint = 1; ipoint <= eroded_points.size(); ++ipoint ) {
		numeric::xyzVector < core::Real > point_to_compare = eroded_points[ipoint].first;
		bool is_in_network = false;

		for ( core::Size inetwork_point = 1; inetwork_point <= network.size(); ++inetwork_point ) {
			numeric::xyzVector < core::Real > point_from_network = network[inetwork_point].first;
			if ( point_from_network[0] - point_to_compare[0] == 0 && point_from_network[1] - point_to_compare[1] == 0 && point_from_network[2] - point_to_compare[2] == 0 ) {
				is_in_network = true;
			}

		}

		//if a point is not in a network check its distance
		//if the point is less than two single C bonds away, include it in the network
		//expand the network with that point
		if ( !is_in_network ) {
			core::Real distance = sqrt( pow( start_point[0] - point_to_compare[0], 2 ) + pow( start_point[1] - point_to_compare[1], 2 ) + pow( start_point[2] - point_to_compare[2], 2 ) );

			if ( distance < 3.1 ) {
				network.push_back(eroded_points[ipoint]);
				network = find_network ( point_to_compare, network, eroded_points );
			}
		}
	}

	return network;
}

//gets the neighbors for a point
std::map< core::Size, numeric::xyzVector< core::Real > >
LigandAligner::assign_neighbors(numeric::xyzVector< core::Real > point, core::Size neighborhood_size){
	std::map< core::Size, numeric::xyzVector< core::Real > > neighbor_map;
	if ( neighborhood_size != 27 && neighborhood_size != 7 && neighborhood_size != 19 ) {
		TR << "Not a valid neighborhood option. Using default of 27." << std::endl;
		neighborhood_size = 27;
	}

	if ( neighborhood_size == 27 ) {
		//z = -1
		neighbor_map[1] = numeric::xyzVector< core::Real >(point[0]-1,point[1]+1,point[2]-1);
		neighbor_map[2] = numeric::xyzVector< core::Real >(point[0],point[1]+1,point[2]-1);
		neighbor_map[3] = numeric::xyzVector< core::Real >(point[0]+1,point[1]+1,point[2]-1);

		neighbor_map[4] = numeric::xyzVector< core::Real >(point[0]-1,point[1],point[2]-1);
		neighbor_map[5] = numeric::xyzVector< core::Real >(point[0],point[1],point[2]-1);
		neighbor_map[6] = numeric::xyzVector< core::Real >(point[0]+1,point[1],point[2]-1);

		neighbor_map[7] = numeric::xyzVector< core::Real >(point[0]-1,point[1]-1,point[2]-1);
		neighbor_map[8] = numeric::xyzVector< core::Real >(point[0],point[1]-1,point[2]-1);
		neighbor_map[9] = numeric::xyzVector< core::Real >(point[0]+1,point[1]-1,point[2]-1);

		//z = 0
		neighbor_map[10] = numeric::xyzVector< core::Real >(point[0]-1,point[1]+1,point[2]);
		neighbor_map[11] = numeric::xyzVector< core::Real >(point[0],point[1]+1,point[2]);
		neighbor_map[12] = numeric::xyzVector< core::Real >(point[0]+1,point[1]+1,point[2]);

		neighbor_map[13] = numeric::xyzVector< core::Real >(point[0]+1,point[1],point[2]);
		neighbor_map[14] = numeric::xyzVector< core::Real >(point[0],point[1],point[2]);
		neighbor_map[15] = numeric::xyzVector< core::Real >(point[0]-1,point[1],point[2]);

		neighbor_map[16] = numeric::xyzVector< core::Real >(point[0]+1,point[1]-1,point[2]);
		neighbor_map[17] = numeric::xyzVector< core::Real >(point[0],point[1]-1,point[2]);
		neighbor_map[18] = numeric::xyzVector< core::Real >(point[0]-1,point[1]-1,point[2]);

		//z = +1
		neighbor_map[19] = numeric::xyzVector< core::Real >(point[0]+1,point[1]+1,point[2]+1);
		neighbor_map[20] = numeric::xyzVector< core::Real >(point[0],point[1]+1,point[2]+1);
		neighbor_map[21] = numeric::xyzVector< core::Real >(point[0]-1,point[1]+1,point[2]+1);

		neighbor_map[22] = numeric::xyzVector< core::Real >(point[0]-1,point[1],point[2]+1);
		neighbor_map[23] = numeric::xyzVector< core::Real >(point[0]-1,point[1],point[2]+1);
		neighbor_map[24] = numeric::xyzVector< core::Real >(point[0]-1,point[1],point[2]+1);

		neighbor_map[25] = numeric::xyzVector< core::Real >(point[0]-1,point[1]-1,point[2]+1);
		neighbor_map[26] = numeric::xyzVector< core::Real >(point[0],point[1]-1,point[2]+1);
		neighbor_map[27] = numeric::xyzVector< core::Real >(point[0]+1,point[1]-1,point[2]+1);
	} else if ( neighborhood_size == 7 ) {
		neighbor_map[1] = numeric::xyzVector< core::Real >(point[0], point[1], point[2]);
		neighbor_map[2] = numeric::xyzVector< core::Real >(point[0]+1, point[1], point[2]);
		neighbor_map[3] = numeric::xyzVector< core::Real >(point[0]-1, point[1], point[2]);
		neighbor_map[4] = numeric::xyzVector< core::Real >(point[0], point[1]+1, point[2]);
		neighbor_map[5] = numeric::xyzVector< core::Real >(point[0], point[1]-1, point[2]);
		neighbor_map[6] = numeric::xyzVector< core::Real >(point[0], point[1], point[2]+1);
		neighbor_map[7] = numeric::xyzVector< core::Real >(point[0], point[1], point[2]-1);
	} else if ( neighborhood_size == 19 ) {
		//share faces
		neighbor_map[1] = numeric::xyzVector< core::Real >(point[0], point[1], point[2]);
		neighbor_map[2] = numeric::xyzVector< core::Real >(point[0]+1, point[1], point[2]);
		neighbor_map[3] = numeric::xyzVector< core::Real >(point[0]-1, point[1], point[2]);
		neighbor_map[4] = numeric::xyzVector< core::Real >(point[0], point[1]+1, point[2]);
		neighbor_map[5] = numeric::xyzVector< core::Real >(point[0], point[1]-1, point[2]);
		neighbor_map[6] = numeric::xyzVector< core::Real >(point[0], point[1], point[2]+1);
		neighbor_map[7] = numeric::xyzVector< core::Real >(point[0], point[1], point[2]-1);

		//share edges
		neighbor_map[8] = numeric::xyzVector< core::Real >(point[0]+1, point[1]+1, point[2]);
		neighbor_map[9] = numeric::xyzVector< core::Real >(point[0]+1, point[1]-1, point[2]);
		neighbor_map[10] = numeric::xyzVector< core::Real >(point[0]-1, point[1]+1, point[2]);
		neighbor_map[11] = numeric::xyzVector< core::Real >(point[0]-1, point[1]-1, point[2]);
		neighbor_map[12] = numeric::xyzVector< core::Real >(point[0]+1, point[1], point[2]+1);
		neighbor_map[13] = numeric::xyzVector< core::Real >(point[0]+1, point[1], point[2]-1);
		neighbor_map[14] = numeric::xyzVector< core::Real >(point[0]-1, point[1], point[2]+1);
		neighbor_map[15] = numeric::xyzVector< core::Real >(point[0]-1, point[1], point[2]-1);
		neighbor_map[16] = numeric::xyzVector< core::Real >(point[0], point[1]+1, point[2]+1);
		neighbor_map[17] = numeric::xyzVector< core::Real >(point[0], point[1]+1, point[2]-1);
		neighbor_map[18] = numeric::xyzVector< core::Real >(point[0], point[1]-1, point[2]+1);
		neighbor_map[19] = numeric::xyzVector< core::Real >(point[0], point[1]-1, point[2]-1);

	}

	return neighbor_map;

}

//checks if a point is a vector of points
bool
LigandAligner::is_point_in_search_group( utility::vector1< std::pair < numeric::xyzVector< core::Real >, core::Real > > selected_points, numeric::xyzVector< core::Real > point_to_compare ) {
	for ( core::Size i = 1; i <= selected_points.size(); i++ ) {
		if ( point_to_compare[0] - selected_points[i].first[0] == 0 && point_to_compare[1] - selected_points[i].first[1] == 0 && point_to_compare[2] - selected_points[i].first[2] == 0 ) {
			return true;
		}
	}
	return false;
}

bool
LigandAligner::check_voxel_distance_to_receptor( numeric::xyzVector< core::Real > voxel, core::pose::Pose const & pose, core::Size const ligid ) {
	bool is_close = false;
	core::Vector resi_atm( 0.0 );
	core::Real x_diff, y_diff, z_diff, distance;
	for ( core::Size ires = 1; ires < ligid; ++ires ) {
		core::conformation::Residue resi = pose.residue( ires );
		for ( core::Size iatm = 1; iatm <= resi.nheavyatoms(); ++iatm ) {
			resi_atm = resi.xyz( iatm );
			x_diff = voxel[0] - resi_atm[0];
			y_diff = voxel[1] - resi_atm[1];
			z_diff = voxel[2] - resi_atm[2];
			distance = sqrt( pow( x_diff, 2 ) + pow( y_diff, 2 ) + pow( z_diff, 2 ) );

			if ( distance < 2.5 ) {
				is_close = true;
			}
		}
	}

	return is_close;
}
// apply
void
LigandAligner::apply(
	LigandConformer & lig
) {
	bool const debug( TR.Debug.visible() );

	// since we are going to minimize a bunch, we work in the context of a minipose
	//    throughout this function rather than converting back and forth

	// store how many structures has been passed through here for biased sampling
	istruct_++;

	core::pose::PoseOP fullpose( utility::pointer::make_shared< core::pose::Pose >() );
	lig.to_pose( fullpose );

	// make a working copy having full pose but redefined moving scs consistent with sf_
	LigandConformer ligwork = LigandConformer( fullpose, lig.ligand_ids(), movable_scs_ );
	ligwork.set_sample_ring_conformers( sample_ring_conformers_ );

	LigandConformer minilig;
	minilig.set_sample_ring_conformers( sample_ring_conformers_ );

	core::pose::PoseOP minipose( utility::pointer::make_shared< core::pose::Pose >() );
	ligwork.to_minipose( minipose, minilig );

	utility::vector1< core::Size > ligids = minilig.ligand_ids();

	// align the ligand CoM to the target
	numeric::xyzVector<core::Real> comTgt(0.0,0.0,0.0);
	for ( core::Size i=1; i<=target_.natoms(); ++i ) {
		comTgt += target_.coord(i);
	}
	comTgt /= target_.natoms();

	// set up grid minimizer
	//   similar to the logic in GridScorer::minimization_loop
	//   but ONLY ligand can move
	core::optimization::MinimizerOptions minopt( "lbfgs_armijo_nonmonotone", 1e-2, true, false, false );

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
	core::Size lig_jump = minilig.get_jumpid();
	mm->set_jump( lig_jump, true );
	if ( !faster_ ) {
		for ( auto ligid : ligids ) {
			mm->set_chi( ligid, true );
			mm->set_bb( ligid, true );
		}
	}

	core::optimization::MinimizerMap min_map;
	min_map.setup( *minipose, *mm );

	core::optimization::Multivec dofs( min_map.nangles() );

	GriddedAtomTreeMultifunc f_i( minilig, *minipose, *sf_, min_map );

	auto start = std::chrono::steady_clock::now();

	// monte carlo
	core::Real score_best=1e6;
	core::pose::Pose minipose_best = *minipose;

	//5*50*3 takes 2sec per model; 5*10*3 1sec ~1sec
	core::Size cycles1 = 20, minsteps1 = 0, repeats1 = 0; // very fast... ~0.001 s per cycle
	if ( lig.has_density_map() ) {
		cycles1 = 20; minsteps1 = 10; repeats1 = 1;
	}
	if ( prealigned_input_ ) {
		cycles1 = 5; minsteps1 = 10; repeats1 = 1;
	}
	core::Size cycles2 = 5, minsteps2 = 50, repeats2 = 3; // adjusted for init-pool gen
	core::Size minsteps_med = 30; // proper minimization b/w stage 1~2

	if ( use_pharmacophore() ) {
		if ( faster_ ) {
			cycles1 = 50; repeats2 = 0; // no stage2
		} else {
			cycles1 = 300; minsteps2 = 10;
		}
	}

	// STAGE 1: randomize w/ increased weight on constraints
	minopt.max_iter( minsteps1 );
	core::optimization::Minimizer min_short( f_i, minopt );

	ConstraintInfo cst_source( *minipose, ligids, use_pharmacophore(), false );
	if ( use_pharmacophore() ) {
		// Caution: should code in updating logic if ligand changes; currently supporting only PHdock mode
		cst_source.map_phores( target_.phores(), false/*update*/ );
		cst_source.select_phore_match( istruct_ );
		// Logic with PHdock:
		// figure out best match within phores at stage1,
		// followed by adding more matches + minimization at stage2
	}

	utility::vector1< std::pair< core::Size, core::Size > > marked_pairs, marked_pairs_best;
	utility::vector1< core::Size > SrcPriorIDs, TgtPriorIDs;
	core::Real const w_prior( 3.0 ); // bonus for phore-matches at stage1

	core::pose::Pose minipose0 = *minipose; //copy of input
	for ( core::Size i=1; i<=cycles1; ++i ) {
		if ( refine_input_ && !prealigned_input_ ) {
			// recover initial
			*minipose = minipose0;
			perturb_lig(*minipose, ligids );
		} else if ( !prealigned_input_ ) {
			randomize_lig(*minipose, ligids, comTgt);
		}

		// superimpose b/w phoreCOMs at receptor & ligand
		if ( use_pharmacophore() ) {
			cst_source.align_to_current_phore_match( *minipose, target_, ligids, marked_pairs,
				SrcPriorIDs, TgtPriorIDs );
		}

		if ( minsteps1 > 0 ) {
			// multiple cstgen+min cycles
			for ( core::Size rpt=1; rpt<=repeats1; ++rpt ) {
				// setup constraints + min; use pharmacophore info if provided for prior match at stage1
				set_constraints(*minipose, ligids, marked_pairs, w_prior, SrcPriorIDs, TgtPriorIDs );

				// the minimizer call
				min_map.copy_dofs_from_pose( *minipose, dofs );
				min_short.run( dofs );
				min_map.reset_jump_rb_deltas( *minipose, dofs );
				min_map.copy_dofs_to_pose( *minipose, dofs );
			}
		} else {
			set_constraints(*minipose, ligids, marked_pairs, w_prior, SrcPriorIDs, TgtPriorIDs );
		}

		core::Real score = sf_->score( *minipose, minilig );
		if ( lig.has_density_map() ) {
			score = sf_->score_init( *minipose, minilig, false ); //hard initial scoring, otherwise overfit into density
		}

		if ( score < score_best ) {
			/*TR << "[iter " << i << "] score = " << score << "/" << score_best
			<< ", npairs " << marked_pairs.size() << std::endl;*/
			score_best = score;
			minipose_best = *minipose;
			marked_pairs_best = marked_pairs;
		}
		*minipose = minipose_best;
	}


	if ( debug ) minipose->dump_pdb("stage1."+std::to_string(istruct_)+".best.pdb" );

	// below: just for reporting
	// use cached info from pose
	*minipose = minipose_best;
	core::Real wcst = (*(sf_->get_sfxn()))[core::scoring::coordinate_constraint];
	core::scoring::ScoreFunctionOP sfxn_cst
		= core::scoring::ScoreFunctionFactory::create_score_function("empty.wts");
	sfxn_cst->set_weight( core::scoring::coordinate_constraint, wcst );
	core::Real cst = sfxn_cst->score( *minipose );

	// STAGE 2: perturb, longer min
	minopt.max_iter( minsteps_med );
	core::optimization::Minimizer min_med( f_i, minopt );

	// minimize w/ hard-cst to the marked_pairs before going to stage2
	TR.Debug << "Stage1, score " << score_best << "(cst=" << cst
		<< ") marked pairs(" << marked_pairs_best.size() << "): ";
	for ( core::Size ipair = 1; ipair <= marked_pairs_best.size(); ++ipair ) {
		TR << " " << marked_pairs_best[ipair].first << "-" << marked_pairs_best[ipair].second;
	}
	TR << std::endl;

	if ( use_pharmacophore() ) {
		set_hard_constraint_on_marked( *minipose, ligids, marked_pairs_best );
		min_map.copy_dofs_from_pose( *minipose, dofs );
		min_med.run( dofs );
		min_map.reset_jump_rb_deltas( *minipose, dofs );
		min_map.copy_dofs_to_pose( *minipose, dofs );
		//if( debug ) minipose->dump_pdb("stage1."+std::to_string(istruct_)+".min.pdb" );

		// substitute back to top-out cst before moving to stage2
		set_constraints(*minipose, ligids, marked_pairs );

		// reporting purpose
		if ( lig.has_density_map() ) {
			score_best = sf_->score_init( *minipose, minilig, false ); //hard initial scoring
		} else {
			score_best = sf_->score( *minipose, minilig );
		}
		cst = sfxn_cst->score( *minipose );
		TR.Debug << "Stage1.5, score " << score_best << "(cst=" << cst << ") " << std::endl;
	}

	minopt.max_iter( minsteps2 );
	core::optimization::Minimizer min_long( f_i, minopt );

	for ( core::Size i=1; i<=cycles2; ++i ) {
		if ( i>1 && !prealigned_input_ ) { perturb_lig(*minipose, ligids); }

		// multiple cstgen+min cycles
		for ( core::Size rpt=1; rpt<=repeats2; ++rpt ) {
			// setup constraints + min
			set_constraints(*minipose, ligids, marked_pairs);

			// the minimizer call
			min_map.copy_dofs_from_pose( *minipose, dofs );
			min_long.run( dofs );
			min_map.reset_jump_rb_deltas( *minipose, dofs );
			min_map.copy_dofs_to_pose( *minipose, dofs );
		}

		core::Real score = 0.0;
		if ( lig.has_density_map() ) {
			score = sf_->score_init( *minipose, minilig, false ); //hard initial scoring
		} else {
			score = sf_->score( *minipose, minilig );
		}

		if ( score < score_best ) {
			//TR << "[iter "<<macro<<".2." << i << "] score = " << score << std::endl;
			score_best = score;
			minipose_best = *minipose;
			marked_pairs_best = marked_pairs;
		}
	}
	// use cached info from pose
	*minipose = minipose_best;
	cst = sfxn_cst->score( *minipose );

	if ( debug ) minipose->dump_pdb("stage2."+std::to_string(istruct_)+".best.pdb" );
	set_constraints(*minipose, ligids, marked_pairs_best );

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end-start;
	//std::chrono::duration<double> diff2 = end-stage1;

	TR.Debug << "Stage2, score "<< score_best << "(cst=" << cst;
	for ( core::Size ipair = 1; ipair <= marked_pairs_best.size(); ++ipair ) {
		TR.Debug << " " << marked_pairs_best[ipair].first << "-" << marked_pairs_best[ipair].second;
	}
	TR.Debug << "), took " << (diff).count() << " seconds " << std::endl;

	// update through ligwork
	ligwork.update_conf_from_minipose( minipose );
	ligwork.to_pose( fullpose );

	lig.update_conf( fullpose );
	lig.score( score_best );
}

}
}
}
