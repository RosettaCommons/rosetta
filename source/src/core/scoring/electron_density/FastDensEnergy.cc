// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FastDensEnergy.cc
/// @brief  Scoring a structure's fit to electron density
/// @author Frank DiMaio


// Unit headers
#include <core/scoring/electron_density/FastDensEnergy.hh>
#include <core/scoring/electron_density/FastDensEnergyCreator.hh>
#include <basic/options/option.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/OneToAllEnergyContainer.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/statistics/functions.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/PDBInfo.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <core/chemical/AtomType.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <utility/vector1.hh>

#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif

using basic::T;
using basic::Error;
using basic::Warning;

// C++

namespace core {
namespace scoring {
namespace electron_density {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.electron_density.FastDensEnergy" );

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

methods::EnergyMethodOP
FastDensEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & opts
) const {
	return methods::EnergyMethodOP( new FastDensEnergy(opts) );
}

ScoreTypes
FastDensEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( elec_dens_fast );
	return sts;
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

using namespace core::scoring::methods;
inline core::Real SQ( core::Real N ) { return N*N; }


methods::LongRangeEnergyType
FastDensEnergy::long_range_type() const { return elec_dens_fast_energy; }

/// c-tor
FastDensEnergy::FastDensEnergy( methods::EnergyMethodOptions const & opts)
: parent( methods::EnergyMethodCreatorOP( new FastDensEnergyCreator ) ) {
	scoreSymmComplex_ = basic::options::option[ basic::options::OptionKeys::edensity::score_symm_complex ]();

	sc_scale_byres_ = opts.get_density_sc_scale_byres();
	if ( sc_scale_byres_.size() < core::chemical::num_canonical_aas ) {
		sc_scale_byres_.resize(core::chemical::num_canonical_aas, 1.0);
	}
}


/// clone
EnergyMethodOP FastDensEnergy::clone() const {
	return EnergyMethodOP( new FastDensEnergy( *this ) );
}

bool
FastDensEnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const {
	return ( pose.residue( res1 ).aa() == core::chemical::aa_vrt || pose.residue( res2 ).aa() == core::chemical::aa_vrt );
}

bool
FastDensEnergy::pose_is_setup_for_density_scoring( pose::Pose const & pose) const {
	kinematics::Edge const &root_edge ( *pose.fold_tree().begin() );
	int virt_res_idx = root_edge.start();
	conformation::Residue const &root_res( pose.residue( virt_res_idx ) );
	if ( root_res.aa() != core::chemical::aa_vrt || root_edge.label() < 0 ) {
		return false;
	}
	return true;
}


void
FastDensEnergy::setup_for_derivatives( pose::Pose & pose , ScoreFunction const & /* sf */) const {
	// grab symminfo (if defined) from the pose
	core::conformation::symmetry::SymmetryInfoCOP symminfo(0);
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >(
			pose.conformation()).Symmetry_Info();
	}
	core::scoring::electron_density::getDensityMap().compute_symm_rotations( pose, symminfo );
}

void
FastDensEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	using namespace methods;

	if ( ! core::scoring::electron_density::getDensityMap().isMapLoaded() ) {
		utility_exit_with_message("Density scoring function called but no map loaded.");
	}

	if ( ! pose_is_setup_for_density_scoring( pose ) ) {
		TR.Warning << "Warning! Fold tree is not set properly for density scoring!" << std::endl;
		return;
	}

	kinematics::Edge const &root_edge ( *pose.fold_tree().begin() );
	int virt_res_idx = root_edge.start();

	// b factor adjustment
	if ( !pose_has_nonzero_Bs( pose ) ) {
		Real effB = core::scoring::electron_density::getDensityMap().getEffectiveBfactor();

		if ( pose.pdb_info() != NULL ) {
			TR.Debug << "Reset B factors to effective value: " << effB << std::endl;
			for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
				for ( core::Size j=1; j<=pose.residue(i).natoms(); ++j ) {
					pose.pdb_info()->temperature( i, j, effB );
				}
			}
		}
	}

	// create LR energy container (if needed)
	LongRangeEnergyType const & lr_type( long_range_type() );
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		OneToAllEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::OneToAllEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.total_residue() || dec->fixed() != virt_res_idx ) {
			create_new_lre_container = true;  // size or root change; recompute
		}
	}

	if ( create_new_lre_container ) {
		TR.Debug << "Creating new one-to-all energy container (" << pose.total_residue() << ")" << std::endl;
		LREnergyContainerOP new_dec( new OneToAllEnergyContainer( virt_res_idx, pose.total_residue(),  elec_dens_fast ) );
		energies.set_long_range_container( lr_type, new_dec );
	}

	// grab symminfo (if defined) from the pose
	// make a copy
	core::conformation::symmetry::SymmetryInfoCOP symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation() ).Symmetry_Info();
	}

	core::scoring::electron_density::getDensityMap().compute_symm_rotations( pose, symminfo );
}


void
FastDensEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	using namespace numeric::statistics;

	// necessary?
	if ( rsd1.aa() != core::chemical::aa_vrt && rsd2.aa() != core::chemical::aa_vrt ) return;
	if ( rsd1.aa() == core::chemical::aa_vrt && rsd2.aa() == core::chemical::aa_vrt ) return;

	if ( ! pose_is_setup_for_density_scoring( pose ) ) return; // already warned in setup

	conformation::Residue const &rsd (rsd1.aa() == core::chemical::aa_vrt? rsd2 : rsd1 );
	Size r = rsd.seqpos();

	core::Real scalefact = 1.0;
	if ( rsd.aa() <= core::chemical::num_canonical_aas ) {
		scalefact = sc_scale_byres_[(int)rsd.aa()];
	}

	// grab symminfo (if defined) from the pose
	core::conformation::symmetry::SymmetryInfoCOP symminfo(0);
	core::Size nsubunits = 1;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation()).Symmetry_Info();
		nsubunits = symminfo->subunits();
		if ( ! symminfo->bb_is_independent( r ) ) return;
	}

	core::Real cc = core::scoring::electron_density::getDensityMap().matchResFast( r, rsd, pose, symminfo, scalefact );
	Real edensScore = -cc;

	if ( symminfo && scoreSymmComplex_ ) {
		edensScore /= ((core::Real)nsubunits);
		utility::vector1< core::Size > bbclones = symminfo->bb_clones( r );
		for ( int i=1; i<=(int)bbclones.size(); ++i ) {
			cc = core::scoring::electron_density::getDensityMap().matchResFast(
				bbclones[i], pose.residue(bbclones[i]), pose, symminfo, scalefact );
			edensScore -= cc / ((core::Real)nsubunits);
		}
	}

	emap[ elec_dens_fast ] += edensScore;
	return;
}


void
FastDensEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	using namespace numeric::statistics;

	// necessary?
	if ( rsd1.aa() != core::chemical::aa_vrt && rsd2.aa() != core::chemical::aa_vrt ) return;
	if ( !pose_is_setup_for_density_scoring( pose ) ) return; // already warned in setup

	conformation::Residue const & res = (rsd2.aa() == core::chemical::aa_vrt) ? rsd1 : rsd2;
	utility::vector1< DerivVectorPair > &r_atom_derivs = (rsd2.aa() == core::chemical::aa_vrt) ? r1_atom_derivs : r2_atom_derivs;

	// apply sc scaling here
	core::Real weight = weights[ elec_dens_fast ];
	core::Real sc_scale = 1.0;
	if ( res.aa() <= core::chemical::num_canonical_aas ) sc_scale = sc_scale_byres_[(int)res.aa()];

	bool is_symmetric = core::pose::symmetry::is_symmetric(pose);
	core::Size resid = res.seqpos();

	core::conformation::symmetry::SymmetryInfoCOP symminfo;
	core::Size nsubunits = 1;
	core::Size nres_per = 1;

	if ( is_symmetric ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation &>(pose.conformation()).Symmetry_Info();
		nsubunits = symminfo->subunits();
		nres_per = symminfo->num_independent_residues();

		// special logic for symm virtuals with score symm complex
		if ( res.aa() == core::chemical::aa_vrt && scoreSymmComplex_ )  {

			numeric::xyzVector< core::Real > dCCdx(0,0,0);

			// if not a branch node, we're done
			utility::vector1< core::kinematics::Edge > edges_i = pose.fold_tree().get_outgoing_edges(resid), edges_j;
			int nchildren = edges_i.size();
			if ( nchildren < 2 ) return;

			// odd case ... if this vrt is controlled by a cloned jump then we dont have to compute derivs
			//    for now we'll just check the parent, but perhaps we should trace root->here?
			if ( !pose.fold_tree().is_root( resid ) ) {
				core::kinematics::Edge edge_incoming = pose.fold_tree().get_residue_edge(resid);
				if ( ! symminfo->jump_is_independent( edge_incoming.label() ) ) return;
			}

			// if any child jumps are cloned such that the clone jump start is ALSO a branching vrt
			//   then we need to aggregate these as well
			// TO DO: if a jump is _fixed_ we should move downstream to the next movable jump
			// This isn't a problem with symm files generated by the perl script, however, it may be a problem in hand-coded ones
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

							core::Real sc_scale_i = 1.0;
							if ( pose.residue(source_res).aa() <= core::chemical::num_canonical_aas ) {
								sc_scale_i = sc_scale_byres_[(int)pose.residue(source_res).aa()];
							}

							for ( int m=1, m_end=pose.residue(target_res).nheavyatoms(); m<=m_end; ++m ) {
								numeric::xyzVector<core::Real> X_lm_src = pose.residue(source_res).atom(m).xyz();
								numeric::xyzVector<core::Real> X_lm_tgt = pose.residue(target_res).atom(m).xyz();

								core::scoring::electron_density::getDensityMap().dCCdx_fastRes
									( m, source_res, X_lm_src, pose.residue(source_res), pose, dCCdx );
								numeric::xyzVector< core::Real > dEdx = -1*R_j * dCCdx / ((core::Real)nsubunits);

								numeric::xyzVector<core::Real> atom_x = X_lm_tgt;
								numeric::xyzVector<core::Real> const f2( dEdx );
								numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
								Vector const f1( atom_x.cross( atom_y ) );

								// define for "ORIG" atom which is 2 (?)
								if ( m <= (int)pose.residue(target_res).last_backbone_atom() ) {
									r_atom_derivs[ 2 ].f1() -= weight * f1;
									r_atom_derivs[ 2 ].f2() -= weight * f2;
								} else {
									r_atom_derivs[ 2 ].f1() -= sc_scale_i * weight * f1;
									r_atom_derivs[ 2 ].f2() -= sc_scale_i * weight * f2;
								}
							}
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

						core::Real sc_scale_i = 1.0;
						if ( pose.residue(source_res).aa() <= core::chemical::num_canonical_aas ) {
							sc_scale_i = sc_scale_byres_[(int)pose.residue(source_res).aa()];
						}

						for ( int m=1, m_end=pose.residue(target_res).nheavyatoms(); m<=m_end; ++m ) {
							numeric::xyzVector<core::Real> X_lm_src = pose.residue(source_res).atom(m).xyz();
							numeric::xyzVector<core::Real> X_lm_tgt = pose.residue(target_res).atom(m).xyz();

							core::scoring::electron_density::getDensityMap().dCCdx_fastRes
								( m, source_res, X_lm_src, pose.residue(source_res), pose, dCCdx );

							numeric::xyzVector< core::Real > dEdx = -1*R_i * dCCdx / ((core::Real)nsubunits);
							numeric::xyzVector<core::Real> atom_x = X_lm_tgt;
							numeric::xyzVector<core::Real> const f2( dEdx );
							numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
							Vector const f1( atom_x.cross( atom_y ) );

						// define for "ORIG" atom which is 2 (?)
							if ( m <= (int)pose.residue(target_res).last_backbone_atom() ) {
								r_atom_derivs[ 2 ].f1() += weight * f1;
								r_atom_derivs[ 2 ].f2() += weight * f2;
							} else {
								r_atom_derivs[ 2 ].f1() += sc_scale_i * weight * f1;
								r_atom_derivs[ 2 ].f2() += sc_scale_i * weight * f2;
							}
						}
					}
				}
			}
		} // if ( res.aa() == core::chemical::aa_vrt && scoreSymmComplex_)
	} // if (is_symmetric)

	if (res.aa() == core::chemical::aa_vrt) return;

	// we are protein ... loop over all residues
	for ( int i=1; i<=(int)res.nheavyatoms(); ++i ) {
		numeric::xyzVector<core::Real> X = res.xyz(i);
		numeric::xyzVector< core::Real > dCCdx(0,0,0);
		numeric::xyzMatrix< core::Real > R = numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1);

		utility::vector1< conformation::Atom > dummyAtmList;
		if ( is_symmetric ) {
			if ( ! symminfo->bb_is_independent( resid ) ) return;

			if ( scoreSymmComplex_ ) {
				// SYMMETRIC CASE, REMAP
				utility::vector1< Size > myClones = symminfo->bb_clones(resid);
				for ( int j=0; j<=(int)myClones.size(); ++j ) {
					numeric::xyzVector<core::Real> X_j = (j==0) ? X : pose.xyz( id::AtomID( i, myClones[j] ) );
					core::scoring::electron_density::getDensityMap().dCCdx_fastRes
						( i, (j==0)?resid:myClones[j], X_j, pose.residue((j==0)?resid:myClones[j]), pose, dCCdx );

					// get R
					core::scoring::electron_density::getDensityMap().get_R( symminfo->subunit_index( (j==0) ? resid : myClones[j] ), R );

					numeric::xyzVector< core::Real > dEdx = -1*R*dCCdx / ((core::Real)nsubunits);

					numeric::xyzVector<core::Real> atom_x = X;
					numeric::xyzVector<core::Real> const f2( dEdx );
					numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
					Vector const f1( atom_x.cross( atom_y ) );

					if ( i <= (int)pose.residue(resid).last_backbone_atom() ) {
						r_atom_derivs[ i ].f1() += weight * f1;
						r_atom_derivs[ i ].f2() += weight * f2;
					} else {
						r_atom_derivs[ i ].f1() += sc_scale * weight * f1;
						r_atom_derivs[ i ].f2() += sc_scale * weight * f2;
					}
				}
			} else {
				// SYMMETRIC CASE, NO REMAP
				core::scoring::electron_density::getDensityMap().dCCdx_fastRes( i, resid, X, pose.residue(resid), pose, dCCdx );
				numeric::xyzVector< core::Real > dEdx = -1.0 * dCCdx;
				numeric::xyzVector<core::Real> atom_x = X;
				numeric::xyzVector<core::Real> const f2( dEdx );
				numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
				Vector const f1( atom_x.cross( atom_y ) );

				if ( i <= (int)pose.residue(resid).last_backbone_atom() ) {
					r_atom_derivs[ i ].f1() += weight * f1;
					r_atom_derivs[ i ].f2() += weight * f2;
				} else {
					r_atom_derivs[ i ].f1() += sc_scale * weight * f1;
					r_atom_derivs[ i ].f2() += sc_scale * weight * f2;
				}
			}
		} else {
			// ASYMMETRIC CASE
			core::scoring::electron_density::getDensityMap().dCCdx_fastRes( i, resid, X, pose.residue(resid), pose, dCCdx );
			numeric::xyzVector< core::Real > dEdx = -1.0 * dCCdx;
			numeric::xyzVector<core::Real> atom_x = X;
			numeric::xyzVector<core::Real> const f2( dEdx );
			numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
			Vector const f1( atom_x.cross( atom_y ) );

			if ( i <= (int)pose.residue(resid).last_backbone_atom() ) {
				r_atom_derivs[ i ].f1() += weight * f1;
				r_atom_derivs[ i ].f2() += weight * f2;
			} else {
				r_atom_derivs[ i ].f1() += sc_scale * weight * f1;
				r_atom_derivs[ i ].f2() += sc_scale * weight * f2;
			}
		}
	} // for each heavyatom
}


core::Size
FastDensEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}
