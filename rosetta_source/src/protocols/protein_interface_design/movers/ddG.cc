// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/ddG.cc
/// @brief implementation of the ddG class for computing interface delta dGs
/// @author Sarel Fleishman (sarelf@u.washington.edu)


#include <core/types.hh>
//#include <protocols/moves/ResidueMover.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pose/Pose.hh>
// Auto-header: duplicate removed #include <core/scoring/Energies.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <numeric/xyzVector.hh>
#include <protocols/protein_interface_design/movers/DesignRepackMover.hh>
#include <protocols/protein_interface_design/movers/ddG.hh>

#include <protocols/moves/RigidBodyMover.hh>

#include <ObjexxFCL/format.hh>

// C++ headers
#include <map>

#include <basic/Tracer.hh>

//Auto Headers
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/tag/Tag.hh> // REQUIRED FOR WINDOWS


namespace protocols {
namespace protein_interface_design {
namespace movers {

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.protein_interface_design.movers.ddG" );

using namespace core;
using namespace protocols::protein_interface_design;
using namespace core::scoring;

ddG::ddG( core::scoring::ScoreFunctionCOP scorefxn_in, core::Size const jump/*=1*/, bool const symmetry /*=false*/ )
{
	scorefxn_ = new core::scoring::ScoreFunction( *scorefxn_in );
	rb_jump_ = jump;
	symmetry_ = symmetry;
}

ddG::~ddG() {}

/// @details a private function for storing energy values
void
ddG::fill_energy_vector( pose::Pose const & pose, std::map< ScoreType, core::Real > & energy_map )
{
	using namespace core::scoring;
	for ( int i=1; i<= n_score_types; ++i ) {
		if ( (*scorefxn_)[ ScoreType(i) ] != 0.0 && ScoreType(i) != pro_close ) {
			energy_map.insert( std::make_pair( ScoreType( i ), (*scorefxn_)[ ScoreType(i)] * pose.energies().total_energies()[ ScoreType( i ) ] ) );
		}
	}
}

/// @details output the ddG values
void
ddG::report_ddG( std::ostream & out ) const
{

	using ObjexxFCL::fmt::F;
	using ObjexxFCL::fmt::LJ;
	out << "-----------------------------------------\n";
	out << " Scores                       Wghtd.Score\n";
	out << "-----------------------------------------\n";
	std::map< ScoreType, Real >::const_iterator unbound_it=unbound_energies_.begin();
	for( std::map< ScoreType, Real >::const_iterator bound_it=bound_energies_.begin();
			 bound_it!=bound_energies_.end();
			 ++bound_it ) {
			 if( std::abs( unbound_it->second ) > 0.001 || std::abs( bound_it->second ) > 0.001 )
				out << ' ' << LJ( 24, bound_it->first ) << ' ' << F( 9,3, bound_it->second - unbound_it->second )<<'\n';
		++unbound_it;
	}
	out << "-----------------------------------------\n";
	out << "Sum ddg: "<< sum_ddG()<<std::endl;
}

/// @details returns the total ddG
Real
ddG::sum_ddG() const
{
	Real sum_energy(0.0);

	std::map< ScoreType, Real >::const_iterator unbound_it=unbound_energies_.begin();
	for( std::map< ScoreType, Real >::const_iterator bound_it=bound_energies_.begin();
			 bound_it!=bound_energies_.end();
			 ++bound_it ) {
		sum_energy += bound_it->second - unbound_it->second;
		++unbound_it;
	}
	return sum_energy;
}

/// @details compute the energy of the repacked complex in the bound and unbound states
void
ddG::calculate( pose::Pose const & pose_in )
{
	using namespace pack;
	using namespace protocols::moves;

	if ( symmetry_ ) {
		symm_ddG( pose_in );
		return;
	}

	repack_partner1_ = true;
	repack_partner2_ = true;
	design_partner1_ = false;
	design_partner2_ = false;
	pose::Pose pose = pose_in;

	//setup_packer_and_movemap( pose );
	task_ = core::pack::task::TaskFactory::create_packer_task( pose );
	task_->initialize_from_command_line().or_include_current( true );
	core::pack::task::operation::RestrictToRepacking rpk;
	rpk.apply( pose, *task_ );
	core::pack::task::operation::NoRepackDisulfides nodisulf;
	nodisulf.apply( pose, *task_ );
	protocols::toolbox::task_operations::RestrictToInterface rti( rb_jump_, 8.0 /*interface_distance_cutoff_*/ );
	rti.apply( pose, *task_ );

	pack::pack_rotamers( pose, *scorefxn_, task_ );
	(*scorefxn_)( pose );
	fill_energy_vector( pose, bound_energies_ );

	RigidBodyTransMoverOP translate( new RigidBodyTransMover( pose, rb_jump_ ) );
	translate->step_size( 1000.0 );
	translate->apply( pose );
	pack::pack_rotamers( pose, *scorefxn_, task_ );
	(*scorefxn_)( pose );
	fill_energy_vector( pose, unbound_energies_ );
}

// @details compute the energy of a repacked symmetrical complex in bound and unbound states
void
ddG::symm_ddG( pose::Pose const & pose_in )
{
	using namespace pack;
	using namespace protocols::moves;


	pose::Pose pose = pose_in;

	assert( core::pose::symmetry::is_symmetric( pose ));
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

  std::map< Size, core::conformation::symmetry::SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	// convert to symetric scorefunction
	scorefxn_ = new scoring::symmetry::SymmetricScoreFunction( scorefxn_ );

	//setup_packer_and_movemap( pose );
	task_ = core::pack::task::TaskFactory::create_packer_task( pose );
	task_->initialize_from_command_line().or_include_current( true );
	core::pack::task::operation::RestrictToRepacking rpk;
	rpk.apply( pose, *task_ );
	core::pack::task::operation::NoRepackDisulfides nodisulf;
	nodisulf.apply( pose, *task_ );
	protocols::toolbox::task_operations::RestrictToInterface rti( rb_jump_, 8.0 /*interface_distance_cutoff_*/ );
	rti.apply( pose, *task_ );

	pack::symmetric_pack_rotamers( pose, *scorefxn_, task_ );
	(*scorefxn_)( pose );
	fill_energy_vector( pose, bound_energies_ );

	RigidBodyDofSeqTransMoverOP translate( new RigidBodyDofSeqTransMover( dofs ) );
	translate->step_size( 1000.0 );
	translate->apply( pose );
	pack::symmetric_pack_rotamers( pose, *scorefxn_, task_ );
	(*scorefxn_)( pose );
	fill_energy_vector( pose, unbound_energies_ );

}

std::string
ddG::get_name() const {
	return "ddG";
}

protocols::moves::MoverOP
ddG::clone() const {
	return (protocols::moves::MoverOP) new ddG( *this );
}

} //movers
} //protein_interface_design
} //protocols
