// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/loops/loop_closure/loophash/LoopHashLoopClosureMover.cc
/// @brief loop closure using (tyka's) loophash library
/// @author Sachko Honda (honda@apl.washington.edu)

// project headers
#include <devel/loophash_loopclosure/LoopHashLoopClosureMover.hh>
#include <devel/loophash_loopclosure/LoopHashLoopClosureMoverCreator.hh>

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <core/types.hh>

#include <devel/init.hh>
#include <basic/options/option.hh>

// Utility Headers

// Unit Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/remodel/RemodelMover.hh>
#include <protocols/forge/methods/util.hh>

#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>

// C++ headers
#include <string>
#include <fstream>

namespace devel{
namespace loophash_loopclosure{

static basic::Tracer TR ("devel.loophash_loopclosure.LoopHashLoopClosureMover" );

LoopHashLoopClosureMoverCreator::LoopHashLoopClosureMoverCreator()
{}
LoopHashLoopClosureMoverCreator::~LoopHashLoopClosureMoverCreator()
{}
protocols::moves::MoverOP
LoopHashLoopClosureMoverCreator::create_mover() const
{
	return new LoopHashLoopClosureMover();
}

std::string
LoopHashLoopClosureMoverCreator::keyname() const
{
	return LoopHashLoopClosureMoverCreator::mover_name();
}
std::string
LoopHashLoopClosureMoverCreator::mover_name() 
{
	return "LoopHashLoopClosureMover";
}

LoopHashLoopClosureMover::LoopHashLoopClosureMover()
{
}
LoopHashLoopClosureMover::~LoopHashLoopClosureMover()
{}
void
LoopHashLoopClosureMover::apply( core::pose::Pose & pose )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if( !option[ remodel::blueprint ].user() ) {
		TR.Error << "remodel:blueprint must be specified" << std::endl;
		exit(1);
	}

	// Use hashloop always because this mover is all about it.
	option[ remodel::RemodelLoopMover::use_loop_hash ]( true );

  if( !option[ lh::db_path ].user() ){
		TR.Error << "loophash_db_path (path to the loophash library) must be specified." << std::endl;
		exit(1);
  }

	if( !option[ remodel::lh_ex_limit ] > 4 ){
		TR.Warning << "lh_ex_limit may be too big and cause segmentation error: " << option[ remodel::lh_ex_limit ] << std::endl;
	}

	remodel_ = new protocols::forge::remodel::RemodelMover();
	remodel_->apply(pose);
}
std::string
LoopHashLoopClosureMover::get_name() const
{
	return "LoopHashLoopClosureMover";
}
protocols::moves::MoverOP
LoopHashLoopClosureMover::clone() const {
	return new LoopHashLoopClosureMover( *this );
}

protocols::moves::MoverOP
LoopHashLoopClosureMover::fresh_instance() const {
	return new LoopHashLoopClosureMover();
}
const std::string
LoopHashLoopClosureMover::make_blueprint( const core::pose::Pose& pose, const std::string& loop_insert_instruction ) const {
	static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );

	std::string bpname = "dummy_blueprint.bp";
	std::ofstream bp( bpname.c_str() );
	runtime_assert( bp.good() );
    int loop_len = 0;
    size_t iChar=0;
    char cur_char = loop_insert_instruction[iChar];
    for ( size_t i=1; i<= pose.total_residue(); ++i ) {
      core::conformation::Residue const & rsd( pose.residue(i) );
      if ( chains[rsd.chain()] == cur_char ) {
			  if( ( loop_insert_instruction.length() > iChar + 1 && loop_insert_instruction[iChar + 1] < 'A' )
						&& ( i+1 < pose.total_residue() && chains[ pose.residue(i+1).chain() ] != cur_char  ) ) {
						bp << i << " A L " << std::endl;
					
        }
        else {
					bp << i << " A . " << std::endl;
				}
      }
      else if( loop_insert_instruction.length() < iChar + 1 ) {
        bp << i << " A . " << std::endl;
      }
      else {
        cur_char = loop_insert_instruction[++iChar];

        if( cur_char >= '1' && cur_char <= '9' ){
          std::string alphanum = "";
          while( cur_char >= '1' && cur_char <= '9' ){
            alphanum += cur_char;
            cur_char = loop_insert_instruction[++iChar];
          }
          loop_len = atoi(alphanum.c_str());
          for( size_t j=0; j<(size_t)loop_len; ++j ) {
            bp << "0 x L" << std::endl;
          }
					bp << i << " A ." << std::endl;
        }
      }
    }
	bp.close();
	return bpname;
}
void LoopHashLoopClosureMover::parse_my_tag(	utility::tag::TagPtr const tag,
																protocols::moves::DataMap & ,
	    	                        protocols::filters::Filters_map const & ,
  	    	                      protocols::moves::Movers_map const &,
    	    	                    core::pose::Pose const & pose ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// Instruction string: e.g.  A6BC6D for closing between A-B and C-D with 6 aa's.
  // Currently, fixed loop length is supported.
  std::string loop_insert_instruction = tag->getOption<std::string>("loop_insert");
	std::string bpname = make_blueprint(pose, loop_insert_instruction);
	option[ OptionKeys::remodel::blueprint ]( bpname );

	bool is_quick_and_dirty = tag->getOption<bool>("quick_and_dirty", true);
	option[ remodel::quick_and_dirty ]( is_quick_and_dirty );

	// Use hashloop always because this mover is all about it.
	option[ remodel::RemodelLoopMover::use_loop_hash ]( true );

	Size num_trajectory = tag->getOption<Size>("num_trajectory", 1);
	option[ remodel::num_trajectory ]( num_trajectory );

	Size num_save_top = tag->getOption<Size>("save_top", 1);
	option[ remodel::save_top ](num_save_top);

	bool no_optH = tag->getOption<bool>( "no_optH", false );
	option[ packing::no_optH ]( no_optH );

	if( !tag->hasOption( "loophash_db_path" ) ){
		TR.Error << "loophash_db_path (path to the loophash library) must be specified." << std::endl;
		exit(1);
  }
	std::string loophash_db_path = tag->getOption<std::string>( "loophash_db_path" );
	option[ lh::db_path ]( loophash_db_path );

	Size loophash_ex_limit = tag->getOption<Size>( "loophash_ex_limit", 4 );
	option[ remodel::lh_ex_limit ]( loophash_ex_limit );
}

} // loophash_loopclosure
} // devel

