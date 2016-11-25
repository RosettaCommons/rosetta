// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.hh
/// @brief Add SiteConstraints from the Epitope to the Paratope and from the Paratope to the Epitope
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_constraints_ParatopeEpitopeSiteConstraintMover_HH
#define INCLUDED_protocols_antibody_constraints_ParatopeEpitopeSiteConstraintMover_HH

#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.fwd.hh>

#include <core/scoring/constraints/SiteConstraint.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>

#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/moves/Mover.hh>


namespace protocols {
namespace antibody {
namespace constraints {

/// @brief Add SiteConstraints from the Epitope to the Paratope and from the Paratope to the Epitope.
/// Will only add the constraint if not already present.
/// @details
/// If no paratope interface residues are given, detects the epitope at 10 A from antibody chain(s).
/// Optionally constrain to paratope CDRs or a specific set of paratope residues.
/// Uses a Linear Harmonic at 0, 1, 10 by default.  Which means epitope will have penalty at greater than 10 A.
/// Linear Harmonic distance tolerance (last number) is set at the interface distance.
///
class ParatopeEpitopeSiteConstraintMover : public protocols::moves::Mover {


public:

	ParatopeEpitopeSiteConstraintMover();
	ParatopeEpitopeSiteConstraintMover(AntibodyInfoCOP ab_info);
	ParatopeEpitopeSiteConstraintMover(AntibodyInfoCOP ab_info, utility::vector1<CDRNameEnum> paratope_cdrs);
	ParatopeEpitopeSiteConstraintMover(AntibodyInfoCOP ab_info, utility::vector1<CDRNameEnum> paratope_cdrs, utility::vector1<bool> epitope_residues);

	~ParatopeEpitopeSiteConstraintMover();

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		moves::Movers_map const & movers,
		Pose const & pose
	) override;

	void
	apply(core::pose::Pose & pose) override;

	void
	remove(core::pose::Pose & pose);

	//void
	//remove(core::pose::Pose & pose, core::Size resnum);

	void
	constrain_to_paratope_cdrs(utility::vector1<CDRNameEnum> const & paratope_cdrs);

	void
	constrain_to_paratope_cdrs(utility::vector1<bool> const & paratope_cdrs);


	void
	constrain_to_paratope_residues(utility::vector1<bool> const & paratope_residues);

	/// @Brief Manually set the epitope residues via PDB Numbering
	void
	constrain_to_epitope_residues(utility::vector1<design::PDBNumbering> const & epitope_residues, core::pose::Pose const & pose);

	/// @Brief Manually set the epitope residues via pose Numbering
	void
	constrain_to_epitope_residues(utility::vector1<bool> const & epitope_residues);


	void
	set_constraint_func(core::scoring::func::FuncOP constraint_func);


	/// @brief The interface distance for antigen epitope auto-detection and the distance at which the default
	///  at which the default flat-harmonic constraint will give a penalty.  10A default.
	void
	set_interface_distance(core::Real const distance);

	void
	set_defaults();

	// XRW TEMP  std::string
	// XRW TEMP  get_name() const {
	// XRW TEMP   return "ParatopeEpitopeSiteConstraintMover";
	// XRW TEMP  }

	utility::vector1<bool>
	get_epitope_residues() const {
		return epitope_residues_;
	}

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	core::scoring::constraints::SiteConstraintOP
	setup_constraints(core::pose::Pose const & pose, core::Size residue, utility::vector1<bool> const & residues) const;

	utility::vector1<bool>
	paratope_residues_from_cdrs(core::pose::Pose const & pose, utility::vector1<bool> const & paratope_cdrs) const;

private:

	AntibodyInfoCOP ab_info_;
	utility::vector1<bool> paratope_residues_;
	utility::vector1<bool> epitope_residues_;
	utility::vector1<bool> paratope_cdrs_;

	//std::map< core::Size, vector1<core::scoring::constraints::AmbiguousConstraintOP > > constraint_map_;
	core::scoring::func::FuncOP current_func_;

	core::Real interface_distance_;


};


}
}
}

#endif //#ifndef INCLUDED_protocols/antibody_design_ParatopeEpitopeSiteConstraintMover_HH

