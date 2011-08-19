/// @file   protocols/comparative_modeling/AlignmentClustering.fwd.hh
///
/// @brief
/// @author TJ Brunette

#ifndef INCLUDED_protocols_comparative_modeling_AlignmentClustering_fwd_hh
#define INCLUDED_protocols_comparative_modeling_AlignmentClustering_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace comparative_modeling {

  class AlignmentCluster;
  typedef utility::pointer::owning_ptr< AlignmentCluster >  AlignmentClusterOP;
  typedef utility::pointer::owning_ptr< AlignmentCluster const >  AlignmentClusterCOP;

  class AlignmentClustering;
  typedef utility::pointer::owning_ptr< AlignmentClustering >  AlignmentClusteringOP;
  typedef utility::pointer::owning_ptr< AlignmentClustering const >  AlignmentClusteringCOP;

} // comparative_modeling
} // protocols
#endif
