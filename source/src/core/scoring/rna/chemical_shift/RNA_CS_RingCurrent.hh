// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file   core/scoring/rna/chemical_shift/RNA_CS_RingCurrent.hh
/// @brief
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)



#ifndef INCLUDED_core_scoring_rna_chemical_shift_RNA_CS_RingCurrent_HH
#define INCLUDED_core_scoring_rna_chemical_shift_RNA_CS_RingCurrent_HH

#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>

#include <core/scoring/rna/chemical_shift/RNA_CS_Parameters.hh>

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {


numeric::xyzVector< core::Real > 
ring_pos( conformation::Residue const & rsd, RNA_CS_residue_parameters const & rna_cs_rsd_params, core::Size const ring_ID ); //NOTE: ONLY INCLUDE HEAVY ATOMS (no hydrogens!)


core::Real
ring_current_effect( numeric::xyzVector< core::Real > const & atom_xyz, conformation::Residue const & source_rsd, RNA_CS_residue_parameters const & rna_cs_rsd_params );

numeric::xyzVector< core::Real >
get_ring_current_deriv( numeric::xyzVector< core::Real > const & CS_data_atom_xyz, 
										  conformation::Residue const & source_rsd, 
											core::Size const source_ring_ID,
											RNA_CS_residue_parameters const & source_rsd_CS_params );


} // chemical_shift
} // rna
} // scoring
} // core

#endif
