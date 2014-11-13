// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/motif/reference_frames.hh
/// @brief  compute motif-related reference frames for residues, chemical groups, etc
/// @author Will Sheffler

#ifndef INCLUDED_core_pose_motif_reference_frames_hh
#define INCLUDED_core_pose_motif_reference_frames_hh

// Utility headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyzTransform.hh>
#include <core/pose/Pose.fwd.hh>

namespace core {
namespace pose {
namespace motif {

numeric::xyzTransform<core::Real>
get_sidechain_reference_frame(
	Pose const & pose,
	core::Size const & ir
);
numeric::xyzTransform<core::Real>
get_backbone_reference_frame(
	Pose const & pose,
	core::Size const & ir
);
numeric::xyzTransform<core::Real>
get_nterminal_peptide_bond_reference_frame(
	Pose const & pose,
	core::Size const & ir
);
numeric::xyzTransform<core::Real>
get_cterminal_peptide_bond_reference_frame(
	Pose const & pose,
	core::Size const & ir
);
inline
numeric::xyzTransform<core::Real>
get_peptide_bond_reference_frame(
	Pose const & pose,
	core::Size const & ir,
	bool const & n_or_c
){
	if(n_or_c) return get_nterminal_peptide_bond_reference_frame(pose,ir);
	else       return get_cterminal_peptide_bond_reference_frame(pose,ir);
}


numeric::xyzTransform<core::Real>
get_backbone_reference_frame(
	numeric::xyzVector<Real> const & N,
	numeric::xyzVector<Real> const & CA,
	numeric::xyzVector<Real> const & C
);
numeric::xyzTransform<core::Real>
get_nterminal_peptide_bond_reference_frame(
	numeric::xyzVector<Real> const & H,
	numeric::xyzVector<Real> const & N,
	numeric::xyzVector<Real> const & CA
);
numeric::xyzTransform<core::Real>
get_cterminal_peptide_bond_reference_frame(
	numeric::xyzVector<Real> const & O,
	numeric::xyzVector<Real> const & C,
	numeric::xyzVector<Real> const & CA
);


numeric::xyzTransform<core::Real> get_frame_ala(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_cys(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_asp(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_glu(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_phe(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_gly(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_his(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_ile(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_lys(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_leu(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_met(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_asn(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_pro(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_gln(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_arg(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_ser(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_thr(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_val(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_trp(Pose const & pose, core::Size const & ir );
numeric::xyzTransform<core::Real> get_frame_tyr(Pose const & pose, core::Size const & ir );

utility::vector1<core::id::AtomID>
get_cterminal_peptide_bond_reference_frame_atomids(
	Pose const & pose,
	core::Size const & ir,
	bool extra=false
);
utility::vector1<core::id::AtomID>
get_nterminal_peptide_bond_reference_frame_atomids(
	Pose const & pose,
	core::Size const & ir,
	bool extra=false
);
inline
utility::vector1<core::id::AtomID>
get_peptide_bond_reference_frame_atomids(
	Pose const & pose,
	core::Size const & ir,
	bool const & n_or_c,
	bool extra=false
){
	if(n_or_c) return get_nterminal_peptide_bond_reference_frame_atomids(pose,ir,extra);
	else       return get_cterminal_peptide_bond_reference_frame_atomids(pose,ir,extra);
}

utility::vector1<core::id::AtomID>
get_backbone_reference_frame_atomids(
	Pose const & pose,
	core::Size const & ir
);

utility::vector1<core::id::AtomID>
get_backbone_reference_frame_atomids_with_downstream(
	Pose const & pose,
	core::Size const & ir
);

utility::vector1<core::id::AtomID>
get_sidechain_reference_frame_atomids(
	Pose const & pose,
	core::Size const & ir
);
utility::vector1<core::id::AtomID>
get_sidechain_reference_frame_atomids_with_downstream(
	Pose const & pose,
	core::Size const & ir
);

utility::vector1<core::id::AtomID> get_atoms_ala(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_cys(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_asp(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_glu(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_phe(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_gly(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_his(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_ile(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_lys(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_leu(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_met(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_asn(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_pro(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_gln(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_arg(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_ser(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_thr(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_val(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_trp(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_tyr(Pose const & pose, core::Size const & ir );

utility::vector1<core::id::AtomID> get_atoms_ala_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_cys_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_asp_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_glu_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_phe_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_gly_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_his_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_ile_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_lys_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_leu_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_met_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_asn_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_pro_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_gln_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_arg_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_ser_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_thr_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_val_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_trp_downstream(Pose const & pose, core::Size const & ir );
utility::vector1<core::id::AtomID> get_atoms_tyr_downstream(Pose const & pose, core::Size const & ir );


} // motif
} // pose
} // core

#endif // INCLUDED_core_pose_util_HH
