
#ifndef INCLUDED_protocols_jd2_LazySilentFileJobInputter_HH
#define INCLUDED_protocols_jd2_LazySilentFileJobInputter_HH

#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentFileData.hh>

namespace protocols{
namespace jd2 {

  class LazySilentFileJobInputter : public protocols::jd2::JobInputter
  {
  public:

    LazySilentFileJobInputter();

    virtual ~LazySilentFileJobInputter();

    virtual void pose_from_job( core::pose::Pose & pose, JobOP job );

    virtual core::io::silent::SilentStruct const& struct_from_job( JobOP job );

    virtual void fill_jobs( Jobs & jobs );

    virtual JobInputterInputSource::Enum input_source() const;
    
    core::io::silent::SilentFileData const& silent_file_data() const { return sfd_; };


  private:
    core::io::silent::SilentFileData sfd_;
  };

} //jd2
} //protocols


#endif //INCLUDED_protocols_jd2_LazySilentFileJobInputter_HH
