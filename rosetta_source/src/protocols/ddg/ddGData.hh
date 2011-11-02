// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef INCLUDED_protocols_ddg_ddGData_hh
#define INCLUDED_protocols_ddg_ddGData_hh

// Unit header
#include <protocols/ddg/ddGData.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <fstream>

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <map>


namespace protocols{
namespace ddG{

using namespace core;
using namespace scoring;

class ddGData {
public:
	//constructors
	ddGData();

	ddGData(std::string filename);

	//virtual ~ddGData();

	//iterator functions
	bool end();
	void get_next_filenames();
	utility::vector1< pose::Pose > read_mut_data();
	utility::vector1< pose::Pose > read_wt_data();
	core::Real read_exp_data();

private:

	std::string file_to_read_;
	std::string curr_mut_filename;
	std::string curr_wt_filename;
	core::Real experimental;
	bool curr_mut_read;
	bool curr_wt_read;
	std::ifstream inputstream;
};//class ddGData


// the following functions added by flo may '11 to read in ddg prediction output files
// so that this info can be used in other protocols

/// @brief small helper class  that stores the ddGs for mutations
/// at a given position. camel case gets weird when trying to write
/// words containing ddG...
class PositionDdGInfo : public utility::pointer::ReferenceCount {

public:

  PositionDdGInfo(
    core::Size seqpos,
    core::chemical::AA wt_aa
   );

  ~PositionDdGInfo();

  void
  add_mutation_ddG(
    core::chemical::AA aa,
    core::Real ddG
  );

  core::Size
  seqpos() const {
    return seqpos_; }

  core::chemical::AA
  wt_aa() const {
    return wt_aa_; }

  std::map< core::chemical::AA, core::Real > const &
  mutation_ddGs() const {
    return mutation_ddGs_; }


private:
  core::Size seqpos_;
  core::chemical::AA wt_aa_;
  std::map< core::chemical::AA, core::Real > mutation_ddGs_;

};


/// @brief function that reads in a ddg predictions out file
/// and returns the info in it as a map of PositionDdGInfo
const std::map< core::Size, PositionDdGInfoOP >
read_ddg_predictions_file( std::string filename );


}//namespace ddG
} //namespace protocols

#endif
