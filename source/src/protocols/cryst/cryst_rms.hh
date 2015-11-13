/// @file
/// @brief

#ifndef INCLUDED_protocols_cryst_cryst_rms_hh
#define INCLUDED_protocols_cryst_cryst_rms_hh

#include <core/pose/Pose.hh>

namespace protocols {
namespace cryst {

core::Size
get_nres_asu( core::pose::Pose const & pose );

core::Real
crystRMSfast (core::pose::Pose &pose_native, core::pose::Pose &pose_decoy);

core::Real
crystRMS (core::pose::Pose &pose_native, core::pose::Pose &pose_decoy);

}
}

#endif
