// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TracerToFile.hh
/// @brief Class for a tracer that writes all output to a file.
/// @author Matt O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDE_basic_tracer_to_file_hh
#define INCLUDE_basic_tracer_to_file_hh



#include <basic/Tracer.hh>

#include <utility/pointer/owning_ptr.hh>

#include <fstream>

namespace basic {


class TracerToFile : public basic::otstream
{
public:
  TracerToFile( std::string const & file_name );
  virtual ~TracerToFile();


protected:
  void t_flush( std::string const & s );

private:
  std::ofstream file_;
};



} // namespace basic


#endif // INCLUDE_basic_tracer_to_file_hh
