// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/loops/loop_closure/ccd/RamaCheck.fwd.hh
/// @brief  Forward declaration of RamaCheck classes (RamaCheckBase, RamaCheck1B, RamaCheck2B)
/// @author Brian D. Weitzner


#ifndef INCLUDED_protocols_loops_loop_closure_ccd_RamaCheck_FWD_HH
#define INCLUDED_protocols_loops_loop_closure_ccd_RamaCheck_FWD_HH


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

// Forward declaration and typedefs for RamaCheck base class
class RamaCheckBase;

typedef  utility::pointer::shared_ptr< RamaCheckBase >  RamaCheckBaseOP;
typedef  utility::pointer::shared_ptr< RamaCheckBase const >  RamaCheckBaseCOP;


// Forward declaration and typedefs for one-body (neighbor-independent) Ramachandran scores
class RamaCheck1B;

typedef  utility::pointer::shared_ptr< RamaCheck1B >  RamaCheck1BOP;
typedef  utility::pointer::shared_ptr< RamaCheck1B const >  RamaCheck1BCOP;

// Forward declaration and typedefs for two-body (neighbor-dependent) Ramachandran scores
class RamaCheck2B;

typedef  utility::pointer::shared_ptr< RamaCheck2B >  RamaCheck2BOP;
typedef  utility::pointer::shared_ptr< RamaCheck2B const >  RamaCheck2BCOP;

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loop_closure_ccd_RamaCheck_FWD_HH
