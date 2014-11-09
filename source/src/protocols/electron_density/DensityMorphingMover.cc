// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Set up morphing with electron density map
/// @author Yifan Song

#include <protocols/electron_density/DensityMorphingMover.hh>
#include <protocols/electron_density/DensityMorphingMoverCreator.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/kinematics/MoveMap.hh>

// utility
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

static thread_local basic::Tracer TR( "protocols.electron_density.DensityMorphingMover" );

namespace protocols {
namespace electron_density {

DensityMorphingMover::DensityMorphingMover() : Mover(){
	init();
}

DensityMorphingMover::~DensityMorphingMover(){
}

void DensityMorphingMover::init() {
	search_radius_ = -1;
	frag_size_ = 11;
	cst_weight_ = 1.;
	bound_width_ = 0.5;
	coord_dev_factor_ = 1.;
}

core::id::AtomID find_best_anchor(core::pose::Pose & pose) {
	core::Size nres = pose.total_residue();
	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
									dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );

		symm_info = SymmConf.Symmetry_Info();
		nres = symm_info->num_independent_residues();
	}

	// find CoM
	numeric::xyzVector<core::Real> sum_xyz(0.0);
	numeric::xyzVector<core::Real> anchor_xyz(0.0);
	core::Real natom = 0.0;
	for ( core::Size ires = 1; ires <= nres; ++ires ) {
		if ( pose.residue_type(ires).has("CA") ) {
			core::Size iatom = pose.residue_type(ires).atom_index("CA");
			sum_xyz += pose.residue(ires).xyz(iatom);
			natom += 1.;
		}
		if (natom > 1e-3) {
			anchor_xyz = sum_xyz / natom;
		}
	}
	core::Real min_dist2 = -1;
	core::Size best_anchor = 0;
	for ( core::Size ires = 1; ires <= nres; ++ires ) {
		if ( pose.residue_type(ires).has("CA") ) {
			core::Size iatom = pose.residue_type(ires).atom_index("CA");
			core::Real dist2 = pose.residue(ires).xyz(iatom).distance_squared(anchor_xyz);
			if (min_dist2 < 0 || dist2 < min_dist2) {
				min_dist2 = dist2;
				best_anchor = ires;
			}
		}
	}

	core::Size best_anchor_atom = pose.residue_type(best_anchor).atom_index("CA");
	return core::id::AtomID(best_anchor_atom, best_anchor);
}

utility::vector1<core::id::AtomID> DensityMorphingMover::collect_fragment_atom_ids(core::pose::Pose const & pose, core::Size ires, core::Size extend_residues) {
	utility::vector1<core::id::AtomID> atom_ids;

	for (int jres = (int) ires ;
		 jres <= (int) (ires + extend_residues);
		 ++jres) {
		if (jres < 1 || jres > (int) nres_) break;

		if ((Size)jres != ires) {
			if (! pose.residue(jres).is_polymer_bonded((Size) (jres-1))) break;
		}

		for (Size jatom=1; jatom<=pose.residue(jres).nheavyatoms(); ++jatom) {
			core::id::AtomID atomid(jatom, jres);
			atom_ids.push_back(atomid);
		}
	}

	for (int jres = (int) ires - 1;
		 jres >= (int) (ires - extend_residues);
		 --jres) {
		if (jres < 1 || jres > (int) nres_) break;

		if (! pose.residue(jres).is_polymer_bonded((Size) (jres+1))) break;

		for (Size jatom=1; jatom<=pose.residue(jres).nheavyatoms(); ++jatom) {
			core::id::AtomID atomid(jatom, jres);
			atom_ids.push_back(atomid);
		}
	}
	return atom_ids;
}


void DensityMorphingMover::apply(core::pose::Pose & pose) {

	core::scoring::electron_density::ElectronDensity & edensity_map = core::scoring::electron_density::getDensityMap();
	core::Real reso(edensity_map.getResolution());

	nres_ = pose.total_residue();
	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );

		symm_info = SymmConf.Symmetry_Info();
		nres_ = symm_info->num_independent_residues();
	}

	core::id::AtomID best_anchor = find_best_anchor(pose);

	for (core::Size ires = 1; ires<=nres_; ++ires) {
		if ( pose.residue(ires).aa() == core::chemical::aa_vrt ) continue;
		utility::vector1 <core::id::AtomID> atom_ids = collect_fragment_atom_ids(pose, ires, (frag_size_ - 1)/2);

		ObjexxFCL::FArray3D< double > calculated_density(edensity_map.data().u1(), edensity_map.data().u2(), edensity_map.data().u3());
		ObjexxFCL::FArray3D< double > inv_mask(edensity_map.data().u1(), edensity_map.data().u2(), edensity_map.data().u3());
		ObjexxFCL::FArray3D< double > mask(edensity_map.data().u1(), edensity_map.data().u2(), edensity_map.data().u3());

		edensity_map.compute_rho(pose, atom_ids, calculated_density, inv_mask);
		for (core::Size i=0; i<edensity_map.data().size(); ++i) mask[i] = 1.-inv_mask[i];

		core::scoring::electron_density::ElectronDensity temp_edensity_map(edensity_map);
		temp_edensity_map.set_data(mask);
		//temp_edensity_map.writeMRC("mask.mrc");

		temp_edensity_map.set_data(calculated_density);
		//temp_edensity_map.writeMRC("calculated_density.mrc");

		ObjexxFCL::FArray3D< double > density_copy(edensity_map.data().u1(), edensity_map.data().u2(), edensity_map.data().u3());
		for (core::Size i=0; i<edensity_map.data().size(); ++i) density_copy[i] = edensity_map.data()[i];
		numeric::xyzVector < double > best_pos(0,0,0);
		best_pos = edensity_map.match_fragment( calculated_density, mask, density_copy, search_radius_);

		if ( pose.residue_type(ires).has("CA") ) {
			core::Size iatom = pose.residue_type(ires).atom_index("CA");
			numeric::xyzVector < double > cst_pos (pose.residue(ires).xyz(iatom)+best_pos);
			using namespace core::scoring::constraints;
			using core::id::AtomID;
			using namespace ObjexxFCL::format;

			core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, coord_dev_factor_*reso) );
			pose.add_constraint( core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new CoordinateConstraint(
														  AtomID(iatom,ires), best_anchor, pose.residue(ires).xyz(iatom),
														  fx ) ) ) );

			TR.Debug << "Constraint added to residue " << ires << ", atom " << iatom << F(8,3,best_pos[0]) << F(8,3,best_pos[1]) << F(8,3,best_pos[2]) << F(8,3,cst_pos[0]) << F(8,3,cst_pos[1]) << F(8,3,cst_pos[2]) << std::endl;

			/*
			pose.add_constraint( new core::scoring::constraints::CoordinateConstraint(
																			core::id::AtomID(iatom,ires), best_anchor, pose.residue(ires).xyz(iatom)+best_pos,
																			new core::scoring::constraints::BoundFunc( 0, reso/2., 2.*reso, "xyz" )) );
			 */
		}

	}

	/*
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	if (scorefxn->get_weight(core::scoring::coordinate_constraint) < 1e-6) {
		scorefxn->set_weight( core::scoring::coordinate_constraint, 1.0 );
	}
	(*scorefxn)(pose);
	scorefxn->show(TR, pose);
	TR << std::endl;

	core::optimization::MinimizerOptions options_lbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	core::optimization::CartesianMinimizer minimizer;

	core::kinematics::MoveMap mm;
	mm.set_bb  ( true );
	mm.set_chi ( true );
	mm.set_jump( true );
	minimizer.run( pose, mm, *scorefxn, options_lbfgs );
	 */
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
DensityMorphingMover::parse_my_tag(
		TagCOP const tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		moves::Movers_map const &,
		Pose const &)
{
	if ( tag->hasOption("frag_size") ) {
		frag_size_ = tag->getOption<core::Real>("frag_size");
	}
	if ( tag->hasOption("search_radius") ) {
		search_radius_ = tag->getOption<core::Real>("search_radius");
	}
	if ( tag->hasOption("coord_dev_factor") ) {
		coord_dev_factor_ = tag->getOption<core::Real>("coord_dev_factor");
	}
	if ( tag->hasOption("bound_width") ) {
		bound_width_ = tag->getOption<core::Real>("bound_width");
	}
	if ( tag->hasOption("cst_weight") ) {
		cst_weight_ = tag->getOption<core::Real>("cst_weight");
	}

	//parse_task_operations( tag, datamap, filters, movers, pose );
}

protocols::moves::MoverOP
DensityMorphingMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DensityMorphingMover );
}

std::string
DensityMorphingMoverCreator::keyname() const
{
	return DensityMorphingMoverCreator::mover_name();
}

std::string
DensityMorphingMoverCreator::mover_name()
{
	return "DensityMorphing";
}

}
}
