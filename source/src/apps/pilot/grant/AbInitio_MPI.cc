// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author grant murphy ( g.s.murphy@gmail.com )


// Headers
#include <protocols/init/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

//#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/fragment/ConstantLengthFragSet.hh>

#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/relax/ClassicRelax.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/moves/OutputMovers.hh>

#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>
#include <numeric/random/random.hh>

#include <utility/excn/Exceptions.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "apps.pilot.grant.AbInitio_MPI" );

namespace core{ namespace options{ namespace OptionKeys{

}}} // core::options::OptionKeys



    ///@brief AbInitio_MPI mover
class AbInitio_MPI : public protocols::moves::Mover {
public:
	AbInitio_MPI() : init_for_input_yet_(false)
	{
	}

	virtual std::string get_name() const { return "AbInitio_MPI"; }


	///@brief init_on_new_input system allows for initializing these details the first time apply() is called.
	virtual
	void
	init_on_new_input( core::pose::Pose & pose) {
		init_for_input_yet_ = true;

		std::string sequence = core::sequence::read_fasta_file( basic::options::option[ basic::options::OptionKeys::in::file::fasta ]()[1] )[1]->sequence();

		core::pose::make_pose_from_sequence( pose, sequence,
			*( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ))
		);

		extended_pose_ = pose;

		fragset_3mer_ = new core::fragment::ConstantLengthFragSet( 3 );
		fragset_3mer_->read_fragment_file( basic::options::option[ basic::options::OptionKeys::in::file::frag3 ].value() );

		fragset_9mer_ = new core::fragment::ConstantLengthFragSet( 9 );
		fragset_9mer_->read_fragment_file( basic::options::option[ basic::options::OptionKeys::in::file::frag9 ].value() );

		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
    movemap->set_bb( true );

		core::scoring::ScoreFunctionOP fullfxn(
                                           core::scoring::get_score_function()
		);

		protocols::abinitio::ClassicAbinitioOP abinitio( new protocols::abinitio::ClassicAbinitio( fragset_3mer_, fragset_9mer_, movemap ));
		abinitio->init( pose );

		abrelax_ = new protocols::moves::SequenceMover();
		abrelax_->add_mover( abinitio );

		//		protocols::moves::PDBDumpMoverOP dump1 = new protocols::moves::PDBDumpMover( "after_abinitio.pdb" );
		//		abrelax_->add_mover( dump1 );

		protocols::simple_moves::SwitchResidueTypeSetMoverOP switch_to_full_atom = new protocols::simple_moves::SwitchResidueTypeSetMover( "fa_standard" );
		abrelax_->add_mover( switch_to_full_atom );

		//		protocols::moves::PDBDumpMoverOP dump2 = new protocols::moves::PDBDumpMover( "after_switch.pdb" );
		//		abrelax_->add_mover( dump2 );


		protocols::relax::RelaxProtocolBaseOP relax( new protocols::relax::ClassicRelax( fullfxn ));
		abrelax_->add_mover( relax );


  }

  virtual ~AbInitio_MPI(){};

  virtual
  void
  apply( core::pose::Pose & pose ){

		if( !init_for_input_yet_ ){ // we will only go thru this loop on the first apply of this protocol
			init_on_new_input( pose );
		}

		pose = extended_pose_;
		abrelax_->apply( pose );

	}

private:
	bool init_for_input_yet_;
	core::fragment::ConstantLengthFragSetOP fragset_3mer_;
	core::fragment::ConstantLengthFragSetOP fragset_9mer_;
	protocols::moves::SequenceMoverOP abrelax_;
	core::pose::Pose extended_pose_;

};

typedef utility::pointer::owning_ptr< AbInitio_MPI > AbInitio_MPIOP;

int main( int argc, char* argv[] )
{
	try {

	using basic::options::option;
	using namespace basic::options::OptionKeys;

	protocols::init::init( argc, argv);

	protocols::jd2::JobDistributor::get_instance()->go(new AbInitio_MPI);
	TR << "************************d**o**n**e**************************************" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
