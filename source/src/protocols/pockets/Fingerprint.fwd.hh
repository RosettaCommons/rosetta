// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/Fingerprint.fwd.hh
/// @brief  protocols::pockets::Fingerprint forward declarations header
/// @author Ragul Gowthaman


#ifndef INCLUDED_protocols_pockets_Fingerprint_fwd_hh
#define INCLUDED_protocols_pockets_Fingerprint_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace pockets {

// Forward
class FingerprintBase;
class NonPlaidFingerprint;
class PlaidFingerprint;

// Types
typedef utility::pointer::owning_ptr< NonPlaidFingerprint >  NonPlaidFingerprintOP;
typedef utility::pointer::owning_ptr< NonPlaidFingerprint const >  NonPlaidFingerprintCOP;
typedef utility::pointer::owning_ptr< PlaidFingerprint >  PlaidFingerprintOP;
typedef utility::pointer::owning_ptr< PlaidFingerprint const >  PlaidFingerprintCOP;


} // namespace pockets
} // namespace protocols

#endif
