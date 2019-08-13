// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/plot_scoreterm.cc
/// @brief Plot a given score term as a function to two degrees of freedom for an arbitrary residue type.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>

// core headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreTypeManager.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

// C++ headers
#include <sstream>

static basic::Tracer TR("plot_scoreterm");



OPT_KEY( String, residue_type )
OPT_KEY( String, dim1 )
OPT_KEY( String, dim2 )
OPT_KEY( Real, phi )
OPT_KEY( Real, psi )
OPT_KEY( Real, omega )
OPT_KEY( Real, omega_minusone )
OPT_KEY( Real, theta )
OPT_KEY( Real, mu )
OPT_KEY( Real, chi1 )
OPT_KEY( Real, chi2 )
OPT_KEY( Real, chi3 )
OPT_KEY( Real, chi4 )
OPT_KEY( Real, stepsize )
OPT_KEY( String, scoreterm )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

}

// Enums:

enum ALLOWED_DOFS {
	PHI=1,
	PSI,
	THETA,
	MU,
	OMEGA,
	OMEGA_MINUSONE,
	CHI1,
	CHI2,
	CHI3,
	CHI4
};

// Prototypes:

ALLOWED_DOFS determine_allowed_dofs( int dofno );

void initialize_pose();

void set_dof( core::pose::Pose &pose, ALLOWED_DOFS const dof, core::Real const &setting);

void plot_scoreterm();


//Functions:

ALLOWED_DOFS determine_allowed_dofs( int dofno ) {
	using namespace basic::options;

	runtime_assert( dofno==1 || dofno==2 );
	std::string const dof(
		dofno == 1
		?
		option[ basic::options::OptionKeys::dim1 ]()
		:
		option[ basic::options::OptionKeys::dim2 ]()
	);

	if ( !dof.compare( "PHI" ) ) { return PHI; }
	else if ( !dof.compare("PSI") ) { return PSI; }
	else if ( !dof.compare("THETA") ) { return THETA; }
	else if ( !dof.compare("MU") ) { return MU; }
	else if ( !dof.compare("OMEGA") ) { return OMEGA; }
	else if ( !dof.compare("OMEGA_MINUSONE") ) { return OMEGA_MINUSONE; }
	else if ( !dof.compare("CHI1") ) { return CHI1; }
	else if ( !dof.compare("CHI2") ) { return CHI2; }
	else if ( !dof.compare("CHI3") ) { return CHI3; }
	else if ( !dof.compare("CHI4") ) { return CHI4; }

	utility_exit_with_message( "Invalid dof \"" + dof + "\" specified for dimension " + std::to_string(dofno) + "." );

	return PHI; //'Cause I'm tryin' to keep the customer satisfied, satisfied!
}

void initialize_pose( core::pose::Pose & pose ) {
	using namespace basic::options;

	std::string const resname( option[ OptionKeys::residue_type ]() );

	core::pose::make_pose_from_sequence( pose, "GX[" + resname + "]G", "fa_standard", true );
	pose.set_phi( 2, option[OptionKeys::phi]() );
	pose.set_psi( 2, option[OptionKeys::psi]() );
	if ( pose.residue_type(2).is_beta_aa() || pose.residue_type(2).is_oligourea() ) {
		pose.set_theta( 2, option[OptionKeys::theta]() );
		if ( pose.residue_type(2).is_oligourea() ) {
			pose.set_mu( 2, option[OptionKeys::mu]() );
		}
	}
	pose.set_omega( 2, option[OptionKeys::omega]() );
	pose.set_omega( 1, option[OptionKeys::omega_minusone]() );
	if ( pose.residue_type(2).nchi() >= 1 ) { pose.set_chi( 1, 2, option[OptionKeys::chi1]() ); }
	if ( pose.residue_type(2).nchi() >= 2 ) { pose.set_chi( 2, 2, option[OptionKeys::chi2]() ); }
	if ( pose.residue_type(2).nchi() >= 3 ) { pose.set_chi( 3, 2, option[OptionKeys::chi3]() ); }
	if ( pose.residue_type(2).nchi() >= 4 ) { pose.set_chi( 4, 2, option[OptionKeys::chi4]() ); }
}

void set_dof( core::pose::Pose &pose, ALLOWED_DOFS const dof, core::Real const &setting) {
	switch( dof ) {
	case PHI :
		pose.set_phi(2, setting);
		break;
	case PSI :
		pose.set_psi(2, setting);
		break;
	case THETA :
		pose.set_theta(2, setting);
		break;
	case MU :
		pose.set_mu(2, setting);
		break;
	case OMEGA :
		pose.set_omega(2, setting);
		break;
	case OMEGA_MINUSONE :
		pose.set_omega(1, setting);
		break;
	case CHI1 :
		pose.set_chi(1, 2, setting);
		break;
	case CHI2 :
		pose.set_chi(2, 2, setting);
		break;
	case CHI3 :
		pose.set_chi(3, 2, setting);
		break;
	case CHI4 :
		pose.set_chi(4, 2, setting);
		break;
	}
}

void plot_scoreterm() {
	using namespace basic::options;

	ALLOWED_DOFS const dof1( determine_allowed_dofs(1) );
	ALLOWED_DOFS const dof2( determine_allowed_dofs(2) );

	core::scoring::ScoreType const score_type( core::scoring::ScoreTypeManager::score_type_from_name( option[OptionKeys::scoreterm]() ) ); //Errors if string isn't parseable.

	core::pose::Pose pose;

	initialize_pose(pose);

	core::Real const stepsize( option[OptionKeys::stepsize]() );
	core::Real curscore(0.0);
	std::stringstream outstream;
	outstream << std::endl;
	core::scoring::ScoreFunction sfxn;
	sfxn.set_weight( score_type, 1.0 );
	for ( core::Real i(-180.0); i<180.0; i+=stepsize ) {
		for ( core::Real j(-180.0); j<180.0; j+=stepsize ) {
			set_dof( pose, dof2, i );
			set_dof( pose, dof1, j );
			curscore = sfxn( pose );
			outstream << curscore << "\t";
		}
		outstream << std::endl;
	}

	TR << outstream.str() << std::endl;
}

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();

		NEW_OPT( residue_type, "The full name of the residue type", "TYR" );
		NEW_OPT( dim1, "The DoF on the first axis (one of PHI, PSI, THETA, MU, OMEGA, OMEGA_MINUSONE, CHI1, CHI2, CHI3, CHI4).", "CHI1" );
		NEW_OPT( dim2, "The DoF on the second axis (one of PHI, PSI, THETA, MU, OMEGA, OMEGA_MINUSONE, CHI1, CHI2, CHI3, CHI4).", "CHI2" );
		NEW_OPT( phi, "Value for phi.", -60.0 );
		NEW_OPT( psi, "Value for psi.", -40.0 );
		NEW_OPT( omega, "Value for omega.", 180.0 );
		NEW_OPT( omega_minusone, "Value for previous omega.", 180.0 );
		NEW_OPT( theta, "Value for theta.", 180.0 );
		NEW_OPT( mu, "Value for mu.", 180.0 );
		NEW_OPT( chi1, "Value for chi1.", 180.0 );
		NEW_OPT( chi2, "Value for chi2.", 180.0 );
		NEW_OPT( chi3, "Value for chi3.", 180.0 );
		NEW_OPT( chi4, "Value for chi4.", 180.0 );
		NEW_OPT( stepsize, "Spacing of sampled grid, in degrees.", 5.0 );
		NEW_OPT( scoreterm, "The score term to plot.  Defaults to fa_dun.", "fa_dun" );

		devel::init( argc, argv );
		plot_scoreterm();

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
