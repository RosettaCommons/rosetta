// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/ralford/make_canonical_helix.cc
/// @brief Creates an ideal a-helix from sequence
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/pose/annotated_sequence.hh>

#include <core/types.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>

// Utility Headers
#include <cstdlib>

static basic::Tracer TR( "apps.pilot.ralford.make_canonical_helix" );

class MakeCanonicalHelix : public protocols::moves::Mover {

public:

	MakeCanonicalHelix() :
		protocols::moves::Mover(),
		phi_( -57.0 ),
		psi_( -47.0 ),
		omega_( 175.0 )
	{
		//register_options();
		//init_options();
	}

	MakeCanonicalHelix( MakeCanonicalHelix const & src ) :
		protocols::moves::Mover( src ),
		phi_( src.phi_ ),
		psi_( src.psi_ ),
		omega_( src.omega_ )
	{}

	virtual ~MakeCanonicalHelix() {}

	virtual std::string get_name() const { return "MakeCanonicalHelix"; }

	virtual
	void
	apply( core::pose::Pose & pose ) {

		TR << "Trasforming current pose coordinates into an ideal alpha helix" << std::endl;

		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			pose.set_phi( i, phi_ );
			pose.set_psi( i, psi_ );
			pose.set_omega( i, omega_ );
		}
	}

private:

	core::Real phi_;
	core::Real psi_;
	core::Real omega_;

};

////////////////////////////////////////////////////////////////////////////////////

typedef utility::pointer::shared_ptr< MakeCanonicalHelix > MakeCanonicalHelixOP;
typedef utility::pointer::shared_ptr< MakeCanonicalHelix const  > MakeCanonicalHelixCOP;

////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{

		devel::init(argc, argv);
		MakeCanonicalHelixOP make_helix( new MakeCanonicalHelix );
		protocols::jd2::JobDistributor::get_instance()->go(make_helix);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "Caught Exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

