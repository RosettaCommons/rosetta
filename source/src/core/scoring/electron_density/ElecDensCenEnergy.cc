// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ElecDensEnergy.cc
/// @brief  Scoring a structure's fit to electron density
/// @author Frank DiMaio


// Unit headers
#include <core/scoring/electron_density/ElecDensCenEnergy.hh>
#include <core/scoring/electron_density/ElecDensCenEnergyCreator.hh>
#include <basic/options/option.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/OneToAllEnergyContainer.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/statistics/functions.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>


// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/keys/edensity.OptionKeys.gen.hh>


#include <basic/Tracer.hh>

#include <core/scoring/electron_density/ElectronDensity.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif
// C++

namespace core {
namespace scoring {
namespace electron_density {


/// @details This must return a fresh instance of the ElecDensCenEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ElecDensCenEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new ElecDensCenEnergy );
}

ScoreTypes
ElecDensCenEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( elec_dens_whole_structure_ca );
	return sts;
}


using namespace core::scoring::methods;

static thread_local basic::Tracer TR( "core.scoring.electron_density.ElecDensEnergy" );

inline core::Real SQ( core::Real N ) { return N*N; }

methods::LongRangeEnergyType
ElecDensCenEnergy::long_range_type() const { return elec_dens_cen_energy; }

/// c-tor
ElecDensCenEnergy::ElecDensCenEnergy() : parent( methods::EnergyMethodCreatorOP( new ElecDensCenEnergyCreator ) ) {}


/// clone
EnergyMethodOP
ElecDensCenEnergy::clone() const {
	return EnergyMethodOP( new ElecDensCenEnergy( *this ) );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

bool
ElecDensCenEnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const {
	return ( pose.residue( res1 ).aa() == core::chemical::aa_vrt || pose.residue( res2 ).aa() == core::chemical::aa_vrt );
}


void
ElecDensCenEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & /* sf */) const
{
	//if (!pose.is_fullatom()) return;

	core::conformation::symmetry::SymmetryInfoCOP symminfo(0);
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation()).Symmetry_Info();
	}
	structure_score = core::scoring::electron_density::getDensityMap().matchCentroidPose( pose, symminfo , true );

	core::scoring::electron_density::getDensityMap().compute_symm_rotations( pose, symminfo );
}


void
ElecDensCenEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	// Do we have a map?
	if ( ! core::scoring::electron_density::getDensityMap().isMapLoaded() ) {
		utility_exit_with_message("Density scoring function called but no map loaded.");
	}

	// make sure the root of the FoldTree is a virtual atom and is followed by a jump
	kinematics::Edge const &root_edge ( *pose.fold_tree().begin() );
	int virt_res_idx = root_edge.start();
	conformation::Residue const &root_res( pose.residue( virt_res_idx ) );

	pose_is_proper = true;
	if ( root_res.type().name() != "VRT" || root_edge.label() < 0 ) {
		utility_exit_with_message("Fold tree is not set properly for density scoring!");
		//pose_is_proper = false;
	}

	// create LR energy container
	LongRangeEnergyType const & lr_type( long_range_type() );
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		OneToAllEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::OneToAllEnergyContainer > ( lrc ) );
		// make sure size or root did not change
		if ( dec->size() != pose.total_residue() || dec->fixed() != virt_res_idx ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		TR.Debug << "Creating new LRE container (" << pose.total_residue() << ")" << std::endl;
		LREnergyContainerOP new_dec( new OneToAllEnergyContainer( virt_res_idx, pose.total_residue(),  elec_dens_whole_structure_ca ) );
		energies.set_long_range_container( lr_type, new_dec );
	}

	// allocate space for per-AA stats
	int nres = pose.total_residue();
	core::scoring::electron_density::getDensityMap().set_nres( nres );

	// grab symminfo (if defined) from the pose
	// make a copy
	core::conformation::symmetry::SymmetryInfoCOP symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation() ).Symmetry_Info();
	}

	// do the actual matching here; split scores among individual residues
	structure_score = core::scoring::electron_density::getDensityMap().matchCentroidPose( pose, symminfo , false );

	TR.Debug << "ElecDensCenEnergy::setup_for_scoring() returns CC = " << structure_score << std::endl;

	// # of scoring residues
	nreses = 0;
	for ( int i=1; i<=(int)pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() && pose.residue(i).atom_name(2) == " CA " ) {
			nreses++;
		}
	}

	core::scoring::electron_density::getDensityMap().compute_symm_rotations( pose, symminfo );
}


void
ElecDensCenEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const {
	return;
}


void
ElecDensCenEnergy::finalize_total_energy(
	pose::Pose const & ,
	ScoreFunction const &,
	EnergyMap &
) const {
	return;
}


void
ElecDensCenEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & , //pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using namespace numeric::statistics;

	if ( rsd1.aa() != core::chemical::aa_vrt && rsd2.aa() != core::chemical::aa_vrt ) return;
	if ( rsd1.aa() == core::chemical::aa_vrt && rsd2.aa() == core::chemical::aa_vrt ) return;

	core::Real cc = structure_score;
	Real z_CC = cc / 0.1;
	Real p_null = 0.5 * errfc( z_CC/sqrt(2.0) );
	Real edensScore = log ( p_null );

	emap[ elec_dens_whole_structure_ca ] += edensScore;
	return;
}


void
ElecDensCenEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & /*sfxn*/,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	using namespace numeric::statistics;

	if ( !pose_is_proper ) return;
	//if (!pose.is_fullatom()) return;

	// derivative
	int resid = id.rsd();
	int atmid = id.atomno();

	// fn only defined for CA's
	if ( pose.residue(resid).aa() != core::chemical::aa_vrt &&
			( !pose.residue(resid).is_protein() || pose.residue(resid).atom_name(atmid) != " CA ") ) return;

	numeric::xyzVector<core::Real> X = pose.xyz(id);
	numeric::xyzVector< core::Real > dCCdx(0,0,0);
	numeric::xyzMatrix< core::Real > R = numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1);

	// if we're symmetric, but _not_ scoring the symmetric complex,
	//    we need to scale derivatives by the score
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		//////////////////////
		// SYMMETRIC CASE
		//////////////////////
		core::conformation::symmetry::SymmetryInfoCOP symminfo =
			(dynamic_cast<const core::conformation::symmetry::SymmetricConformation &>(pose.conformation()).Symmetry_Info());
		core::Size nsubunits = symminfo->subunits();
		core::Size nres_per = symminfo->num_independent_residues();
		bool remapSymm = basic::options::option[ basic::options::OptionKeys::edensity::score_symm_complex ]();

		if ( pose.residue(resid).aa() == core::chemical::aa_vrt && remapSymm )  {

			// derivative is only defined for the 'ORIG' atom in the virtual
			if ( atmid != 2 ) return;

			// if not a branch node, we're done
			utility::vector1< core::kinematics::Edge > edges_i = pose.fold_tree().get_outgoing_edges(resid), edges_j;
			int nchildren = edges_i.size();
			if ( nchildren < 2 ) return;

			// odd case ... if this vrt is controlled by a cloned jump then we dont have to compute derivs
			//    for now we'll just check the parent, but perhaps we should trace root->here?
			if ( !pose.fold_tree().is_root( resid ) ) {
				core::kinematics::Edge edge_incoming = pose.fold_tree().get_residue_edge(resid);
				if ( ! symminfo->jump_is_independent( edge_incoming.label() ) ) {
					return;
				}
			}

			// if any child jumps are cloned such that the clone jump start is ALSO a branching vrt
			//   then we need to aggregate these as well
			// TO DO: if a jump is _fixed_ we should move downstream to the next movable jump
			// This isn't a problem with symm files generated by the perl script, however,
			//    it may be a problem in hand-coded ones
			utility::vector1< core::Size > vrtclones;
			for ( int j=1; j<=nchildren; ++j ) {
				int basejump = edges_i[j].label();
				if ( ! symminfo->jump_is_independent( basejump ) ) {
					basejump = symminfo->jump_follows( basejump );
				}
				utility::vector1< core::Size > jumpclones = symminfo->jump_clones( basejump );
				for ( int k=0; k<=(int)jumpclones.size(); ++k ) {
					core::Size upstream = pose.fold_tree().jump_edge( k==0 ? basejump : jumpclones[k] ).start();
					edges_j = pose.fold_tree().get_outgoing_edges(upstream);
					if ( edges_j.size() > 1 && std::find( vrtclones.begin(), vrtclones.end(), upstream ) == vrtclones.end() ) {
						vrtclones.push_back( upstream );
					}
				}
			}

			// loop over all clones of this VRT
			for ( int i=1; i<=(int)vrtclones.size(); ++i ) {
				edges_i = pose.fold_tree().get_outgoing_edges(vrtclones[i]);

				// STEP 1: subtract children's contribution
				for ( int j=1; j<=nchildren; ++j ) {
					int downstream = edges_i[j].stop();
					utility::vector1<int> mapping_j;
					numeric::xyzMatrix< core::Real > R_j;
					core::scoring::electron_density::getDensityMap().get_symmMap( downstream , mapping_j, R_j );

					for ( int k=1; k<=(int)nsubunits; ++k ) {
						if ( mapping_j[k] == 0 ) continue;   // subunit k is not under child j

						// loop over all atms in reses in subunit k
						for ( int l=1, l_end=nres_per; l<=l_end; ++l ) {
							// there is a mapping from k->mapping_j[k]
							int source_res = (k-1)*nres_per+l;
							int target_res = (mapping_j[k]-1)*nres_per+l;

							if ( !pose.residue(source_res).is_protein() || pose.residue(source_res).atom_name(2)!=" CA " ) continue;

							numeric::xyzVector<core::Real> X_lm_src = pose.residue(source_res).atom(2).xyz();
							numeric::xyzVector<core::Real> X_lm_tgt = pose.residue(target_res).atom(2).xyz();

							core::scoring::electron_density::getDensityMap().dCCdx_cen( source_res, X_lm_src, pose, dCCdx );
							Real CC = structure_score;
							Real z_CC = CC / 0.1;
							Real p_null = 0.5 * errfc( z_CC/sqrt(2.0) );
							numeric::xyzVector< core::Real > dEdx = 0.5 * ( 1.0 / p_null ) *
								(-2.0/sqrt(M_PI)) *
								exp(-SQ( z_CC/sqrt(2.0) )) *
								1/sqrt(0.02) *
								nreses * R_j * dCCdx / ((core::Real)nsubunits);

							numeric::xyzVector<core::Real> atom_x = X_lm_tgt;
							numeric::xyzVector<core::Real> const f2( dEdx );
							numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
							Vector const f1( atom_x.cross( atom_y ) );

							F1 -= weights[ elec_dens_whole_structure_ca ] * f1;
							F2 -= weights[ elec_dens_whole_structure_ca ] * f2;
						}
					}
				}

				// STEP 2: add my contribution
				utility::vector1<int> mapping_i;
				numeric::xyzMatrix< core::Real > R_i;
				core::scoring::electron_density::getDensityMap().get_symmMap( vrtclones[i] , mapping_i, R_i );
				for ( int k=1; k<=(int)nsubunits; ++k ) {
					if ( mapping_i[k] == 0 ) continue;   // subunit k is not under child j

					// loop over all atms in reses in subunit k
					for ( int l=1, l_end=nres_per; l<=l_end; ++l ) {
						// there is a mapping from k->mapping_i[k]
						int source_res = (k-1)*nres_per+l;
						int target_res = (mapping_i[k]-1)*nres_per+l;

						if ( !pose.residue(source_res).is_protein() || pose.residue(source_res).atom_name(2)!=" CA " ) continue;

						numeric::xyzVector<core::Real> X_lm_src = pose.residue(source_res).atom(2).xyz();
						numeric::xyzVector<core::Real> X_lm_tgt = pose.residue(target_res).atom(2).xyz();

						core::scoring::electron_density::getDensityMap().dCCdx_cen( source_res, X_lm_src, pose, dCCdx );
						Real CC = structure_score;
						Real z_CC = CC / 0.1;
						Real p_null = 0.5 * errfc( z_CC/sqrt(2.0) );
						numeric::xyzVector< core::Real > dEdx = 0.5 * ( 1.0 / p_null ) *
							(-2.0/sqrt(M_PI)) *
							exp(-SQ( z_CC/sqrt(2.0) )) *
							1/sqrt(0.02) *
							nreses * R_i * dCCdx / ((core::Real)nsubunits);

						numeric::xyzVector<core::Real> atom_x = X_lm_tgt;
						numeric::xyzVector<core::Real> const f2( dEdx );
						numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
						Vector const f1( atom_x.cross( atom_y ) );

						F1 += weights[ elec_dens_whole_structure_ca ] * f1;
						F2 += weights[ elec_dens_whole_structure_ca ] * f2;
					}
				}
			}
		} else { // NON-VRT
			if ( ! symminfo->bb_is_independent( resid ) ) return;

			if ( remapSymm ) {
				utility::vector1< Size > myClones = symminfo->bb_clones(resid);
				for ( int i=0; i<=(int)myClones.size(); ++i ) {
					numeric::xyzVector<core::Real> X_i = (i==0) ? X : pose.xyz( id::AtomID( 2, myClones[i] ) );
					core::scoring::electron_density::getDensityMap().dCCdx_cen( (i==0) ? resid : myClones[i], X_i, pose, dCCdx );

					// get R
					core::scoring::electron_density::getDensityMap().get_R( symminfo->subunit_index( (i==0) ? resid : myClones[i] ), R );

					Real CC = structure_score;
					Real z_CC = CC / 0.1;
					Real p_null = 0.5 * errfc( z_CC/sqrt(2.0) );

					// divide by the number of subunits since rosetta will scale up later (??)
					numeric::xyzVector< core::Real > dEdx = 0.5 * ( 1.0 / p_null ) *
						(-2.0/sqrt(M_PI)) *
						exp(-SQ( z_CC/sqrt(2.0) )) *
						1/sqrt(0.02) *
						nreses * R * dCCdx / ((core::Real)nsubunits);

					numeric::xyzVector<core::Real> atom_x = X;
					numeric::xyzVector<core::Real> const f2( dEdx );
					numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
					Vector const f1( atom_x.cross( atom_y ) );

					F1 += weights[ elec_dens_whole_structure_ca ] * f1;
					F2 += weights[ elec_dens_whole_structure_ca ] * f2;
				}
			} else {
				Real CC = structure_score;
				Real z_CC = CC / 0.1;
				Real p_null = 0.5 * errfc( z_CC/sqrt(2.0) );

				core::scoring::electron_density::getDensityMap().dCCdx_cen( resid, X, pose, dCCdx );
				numeric::xyzVector< core::Real > dEdx = ( 1.0 / p_null ) *
					0.5 *
					(-2.0/sqrt(M_PI)) *
					exp(-SQ( z_CC/sqrt(2.0) )) *
					1/sqrt(0.02) *
					nreses *
					dCCdx / ((core::Real)nsubunits);

				numeric::xyzVector<core::Real> atom_x = X;
				numeric::xyzVector<core::Real> const f2( dEdx );
				numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;   // a "fake" atom in the direcion of the gradient
				Vector const f1( atom_x.cross( atom_y ) );

				F1 += weights[ elec_dens_whole_structure_ca ] * f1;
				F2 += weights[ elec_dens_whole_structure_ca ] * f2;
			}
		}
	} else {
		//////////////////////
		// ASYMMETRIC CASE
		//////////////////////

		// the derivative of edens _score_ as a function of dCC/dx
		Real CC = structure_score;
		Real z_CC = CC / 0.1;
		Real p_null = 0.5 * errfc( z_CC/sqrt(2.0) );

		core::scoring::electron_density::getDensityMap().dCCdx_cen( resid, X, pose, dCCdx );

		// divide by the number of subunits since rosetta will scale up later (??)
		numeric::xyzVector< core::Real > dEdx = ( 1.0 / p_null ) *
			0.5 *
			(-2.0/sqrt(M_PI)) *
			exp(-SQ( z_CC/sqrt(2.0) )) *
			1/sqrt(0.02) *
			nreses *
			dCCdx;

		numeric::xyzVector<core::Real> atom_x = X;
		numeric::xyzVector<core::Real> const f2( dEdx );
		numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;   // a "fake" atom in the direcion of the gradient
		Vector const f1( atom_x.cross( atom_y ) );

		F1 += weights[ elec_dens_whole_structure_ca ] * f1;
		F2 += weights[ elec_dens_whole_structure_ca ] * f2;
	}
}
core::Size
ElecDensCenEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
