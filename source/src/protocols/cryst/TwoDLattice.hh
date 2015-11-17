// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Zibo Chen


#ifndef INCLUDED_protocols_cryst_TwoDLattice_hh
#define INCLUDED_protocols_cryst_TwoDLattice_hh

#include <protocols/cryst/wallpaper.hh>

#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/scoring/electron_density/util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>


#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>

#include <utility/excn/Exceptions.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>


namespace protocols {
namespace cryst {


class MakeLayerMover : public protocols::moves::Mover {
private:
	WallpaperGroup wg_;
	core::Real contact_dist_;
	bool moving_lattice_;

public:
	MakeLayerMover() : contact_dist_(24.0), moving_lattice_(true) {}

	void
	set_moving_lattice( bool val ) { moving_lattice_ = val; }

	Size
	place_near_origin (
		core::pose::Pose & pose
	);

	void
	add_monomers_to_layer(
		core::pose::Pose const & monomer_pose,
		core::pose::Pose & pose,
		utility::vector1<Size> const & monomer_anchors,
		utility::vector1<Size> & monomer_jumps,
		core::Size rootres
	);

	void
	detect_connecting_subunits(
		core::pose::Pose const & monomer_pose,
		core::pose::Pose const & pose,
		utility::vector1<Size> & monomer_anchors,
		core::Size &basesubunit
	);

	void
	build_layer_of_virtuals(
		core::pose::Pose & pose,
		utility::vector1<Size> &Ajumps,
		utility::vector1<Size> &Bjumps,
		utility::vector1<Size> &subunit_anchors,
		core::Size &basesubunit
	);

	void
	setup_xtal_symminfo(
		core::pose::Pose const & pose,
		core::Size const num_monomers,
		core::Size const num_virtuals,
		core::Size const base_monomer,
		core::Size const nres_monomer,
		utility::vector1<Size> const &Ajumps,
		utility::vector1<Size> const &Bjumps,
		utility::vector1<Size> const &monomer_jumps,
		core::conformation::symmetry::SymmetryInfo & symminfo
	);


	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const {
		return "MakeLayerMover";
	}

	virtual protocols::moves::MoverOP clone() const { return protocols::moves::MoverOP(new MakeLayerMover(*this)); }
	virtual protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP(new MakeLayerMover()); }

	void
	wallpaper_group( WallpaperGroup const & wg_in ) {
		wg_ = wg_in;
	}

	void
	contact_dist( core::Real const setting ) {
		contact_dist_ = setting;
	}

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const & ,
		moves::Movers_map const & ,
		core::pose::Pose const & pose );

private:

	core::conformation::ResidueOP make_vrt( core::Vector O, core::Vector X, core::Vector Y ) {
		core::conformation::ResidueOP vrtrsd
			( core::conformation::ResidueFactory::create_residue(
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( "VRT" ) ) );
		vrtrsd->set_xyz("ORIG",O);
		vrtrsd->set_xyz("X",O+X);
		vrtrsd->set_xyz("Y",O+Y);
		return vrtrsd;
	}
};

}
}

#endif
