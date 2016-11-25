// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_cryst_refinable_lattice_hh
#define INCLUDED_protocols_cryst_refinable_lattice_hh

#include <protocols/cryst/spacegroup.hh>

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
#include <basic/options/keys/cryst.OptionKeys.gen.hh>


namespace protocols {
namespace cryst {

class UpdateCrystInfo : public protocols::moves::Mover {
public:
	void apply( core::pose::Pose & pose ) override;

	// XRW TEMP  std::string get_name() const override {
	// XRW TEMP   return "UpdateCrystInfo";
	// XRW TEMP  }

	protocols::moves::MoverOP clone() const override { return protocols::moves::MoverOP(new UpdateCrystInfo(*this)); }
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP(new UpdateCrystInfo()); }

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const & ,
		moves::Movers_map const & ,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


};


class DockLatticeMover : public protocols::moves::Mover {
private:
	core::scoring::ScoreFunctionOP sf_, sf_vdw_;

	std::map< core::Size, core::conformation::symmetry::SymDof > symdofs_;
	core::Size SUBjump_;  // jump ids, may not correspond to dofs
	core::Real rot_mag_, trans_mag_;
	core::Real temp_;
	core::Size ncycles_;
	core::Real monomer_bump_;
	bool fullatom_, design_;

public:
	DockLatticeMover();
	DockLatticeMover(core::scoring::ScoreFunctionOP sf_cen_in);

	void
	init(core::pose::Pose & pose);

	void
	perturb_trial( core::pose::Pose & pose );

	void
	slide_lattice( core::pose::Pose & pose );

	void
	min_lattice( core::pose::Pose & pose );

	void
	modify_lattice( core::pose::Pose & pose, core::Real mag );

	void apply( core::pose::Pose & pose ) override;

	void set_rot_mag(core::Real rot_mag_in) { rot_mag_=rot_mag_in;}
	void set_trans_mag(core::Real trans_mag_in) { trans_mag_=trans_mag_in;}
	void set_ncycles(core::Size ncycles_in) { ncycles_=ncycles_in;}
	void set_fullatom(bool val) { fullatom_=val; }
	void set_scorefunction(core::scoring::ScoreFunctionOP sf_new) { sf_=sf_new->clone(); }
	void set_design(bool val) { design_=val; }

	void set_temp(core::Real val) { temp_=val; }

	// XRW TEMP  std::string get_name() const override {
	// XRW TEMP   return "DockLatticeMover";
	// XRW TEMP  }

	protocols::moves::MoverOP clone() const override { return protocols::moves::MoverOP(new DockLatticeMover(*this)); }
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP(new DockLatticeMover()); }

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const & ,
		moves::Movers_map const & ,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


};


class MakeLatticeMover : public protocols::moves::Mover {
private:
	Spacegroup sg_;
	core::Real contact_dist_;
	bool refinable_lattice_;

	utility::vector1< numeric::xyzMatrix<core::Real> > allRs_;
	utility::vector1< numeric::xyzVector<core::Real> > allTs_;

public:
	MakeLatticeMover() {
		using namespace basic::options;
		refinable_lattice_ = option[ OptionKeys::cryst::refinable_lattice]();
		contact_dist_ = option[ OptionKeys::cryst::interaction_shell]();
	}

	void
	set_refinable_lattice( bool val ) {
		refinable_lattice_ = val;
	}

	// get R and T from subunit index
	void
	getRT( int idx, numeric::xyzMatrix<core::Real> &R, numeric::xyzVector<core::Real> &T ) {
		R = allRs_[idx];
		T = allTs_[idx];
	}

	// lowerlevel access to sg_
	utility::vector1<core::kinematics::RT> const &symmops() const { return sg_.symmops(); }
	Spacegroup const &sg() const { return sg_; }


	Size
	place_near_origin (
		core::pose::Pose & pose
	);

	void
	add_monomers_to_lattice(
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
	build_lattice_of_virtuals(
		core::pose::Pose & pose,
		numeric::xyzVector<int> EXTEND,
		utility::vector1<Size> &Ajumps,
		utility::vector1<Size> &Bjumps,
		utility::vector1<Size> &Cjumps,
		utility::vector1<Size> &subunit_anchors,
		core::Size &basesubunit
	);

	void
	setup_xtal_symminfo(
		core::pose::Pose & pose,
		core::Size const num_monomers,
		core::Size const num_virtuals,
		core::Size const base_monomer,
		core::Size const nres_monomer,
		utility::vector1<Size> const &Ajumps,
		utility::vector1<Size> const &Bjumps,
		utility::vector1<Size> const &Cjumps,
		utility::vector1<Size> const &monomer_jumps,
		core::conformation::symmetry::SymmetryInfo & symminfo
	);


	void apply( core::pose::Pose & pose ) override;

	// XRW TEMP  std::string get_name() const override {
	// XRW TEMP   return "MakeLatticeMover";
	// XRW TEMP  }

	protocols::moves::MoverOP clone() const override { return protocols::moves::MoverOP(new MakeLatticeMover(*this)); }
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP(new MakeLatticeMover()); }

	void
	spacegroup( Spacegroup const & sg_in ) {
		sg_ = sg_in;
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
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	core::conformation::ResidueOP make_vrt( core::Vector O, core::Vector X, core::Vector Y, bool inv=false ) {
		core::conformation::ResidueOP vrtrsd;
		if ( inv ) {
			vrtrsd = core::conformation::ResidueOP( core::conformation::ResidueFactory::create_residue(
				core::chemical::ChemicalManager::get_instance()->residue_type_set(
				core::chemical::FA_STANDARD )->name_map( "INV_VRT" ) ) );
		} else {
			vrtrsd = core::conformation::ResidueOP( core::conformation::ResidueFactory::create_residue(
				core::chemical::ChemicalManager::get_instance()->residue_type_set(
				core::chemical::FA_STANDARD )->name_map( "VRT" ) ) );
		}
		vrtrsd->set_xyz("ORIG",O);
		vrtrsd->set_xyz("X",O+X);
		vrtrsd->set_xyz("Y",O+Y);
		return vrtrsd;
	}
};

}
}

#endif
