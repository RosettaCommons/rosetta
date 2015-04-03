// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file   core/scoring/rna/chemical_shift/RNA_CS_Util.hh
/// @brief
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


#ifndef INCLUDED_core_scoring_rna_chemical_shift_RNA_CS_Util_HH
#define INCLUDED_core_scoring_rna_chemical_shift_RNA_CS_Util_HH

#include <core/types.hh>

#include <core/chemical/AA.hh>


#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>


#include <core/scoring/rna/chemical_shift/RNA_CS_Parameters.hh>


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {


inline 
core::Size 
dround( core::Real var )
{
   return ( core::Size( var + 0.5 ) );
}


numeric::xyzMatrix< core::Real > const
get_rna_base_coordinate_system_from_CS_params( core::conformation::Residue const & rsd, RNA_CS_residue_parameters const & rna_cs_rsd_params );


} //chemical_shift
} //rna
} //scoring
} //core


#endif
