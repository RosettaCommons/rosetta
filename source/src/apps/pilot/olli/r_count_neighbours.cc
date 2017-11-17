// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

#include <protocols/abinitio/Templates.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <devel/init.hh>

#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/BoundConstraint.hh>


#include <protocols/jd2/JobDistributor.hh>


// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/exit.hh>


// Unit headers

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.hh>
//#include <core/scoring/CachedData.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>


// Project headers
// Auto-header: duplicate removed #include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>


// Utility headers


//Auto Headers
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/mpistream.hh>


static basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;
using namespace abinitio;
using namespace jumping;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace scoring;
using namespace scoring::constraints;


OPT_KEY( Integer, level ) //viol_level already in AbrelaxApplication definier
OPT_1GRP_KEY( File, in, top )
OPT_KEY( Real, threshold )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	Templates::register_options();
	OPT( in::file::s );
	OPT( constraints::cst_file );
	OPT( out::prefix );
	NEW_OPT( threshold, "what makes a violation", 1 );
	NEW_OPT( level, "how much detail for violation output", 71 );
	NEW_OPT( in::top, "read topology from this file for checking", "");
	//NEW_OPT( viol_type, "work only on these types of constraints", "");
}


Distance const start_sig = 9.8;
Distance const end_sig   = 10.2;

DistanceSquared const start_sig2 = start_sig*start_sig;
DistanceSquared const end_sig2   = end_sig*end_sig;

// Forward
class EnvToolMover;

// Types
typedef  utility::pointer::owning_ptr< EnvToolMover >  EnvToolMoverOP;
typedef  utility::pointer::owning_ptr< EnvToolMover const >  EnvToolMoverCOP;

class EnvToolMover : public moves::Mover {
public:
	virtual void apply( core::pose::Pose& );
	Real sigmoidish_neighbor( DistanceSquared const sqdist ) const;

	std::string const &
	representative_atom_name( chemical::AA const aa ) const;

private:
};

inline Real sqr ( Real x ) {
	return x*x;
}

Real
EnvToolMover::sigmoidish_neighbor( DistanceSquared const sqdist ) const
{
	if ( sqdist > end_sig2 ) {
		return 0.0;
	} else if ( sqdist < start_sig2 ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0  - sqr( (dist - start_sig) / (end_sig - start_sig) ) );
	}
}


/// @details returns const & to static data members to avoid expense
/// of string allocation and destruction.  Do not call this function
/// on non-canonical aas
std::string const &
EnvToolMover::representative_atom_name( chemical::AA const aa ) const
{
	assert( aa >= 1 && aa <= chemical::num_canonical_aas );

	static std::string const cbeta_string(  "CB"  );
	static std::string const sgamma_string( "SG"  );
	static std::string const cgamma_string( "CG"  );
	static std::string const cdelta_string( "CD"  );
	static std::string const czeta_string(  "CZ"  );
	static std::string const calpha_string( "CA"  );
	static std::string const ceps_1_string( "CE1" );
	static std::string const cdel_1_string( "CD1" );
	static std::string const ceps_2_string( "CE2" );
	static std::string const sdelta_string( "SD"  );

	switch ( aa ) {
	case ( chemical::aa_ala ) : return cbeta_string;  break;
	case ( chemical::aa_cys ) : return sgamma_string; break;
	case ( chemical::aa_asp ) : return cgamma_string; break;
	case ( chemical::aa_glu ) : return cdelta_string; break;
	case ( chemical::aa_phe ) : return czeta_string;  break;
	case ( chemical::aa_gly ) : return calpha_string; break;
	case ( chemical::aa_his ) : return ceps_1_string; break;
	case ( chemical::aa_ile ) : return cdel_1_string; break;
	case ( chemical::aa_lys ) : return cdelta_string; break;
	case ( chemical::aa_leu ) : return cgamma_string; break;
	case ( chemical::aa_met ) : return sdelta_string; break;
	case ( chemical::aa_asn ) : return cgamma_string; break;
	case ( chemical::aa_pro ) : return cgamma_string; break;
	case ( chemical::aa_gln ) : return cdelta_string; break;
	case ( chemical::aa_arg ) : return czeta_string;  break;
	case ( chemical::aa_ser ) : return cbeta_string;  break;
	case ( chemical::aa_thr ) : return cbeta_string;  break;
	case ( chemical::aa_val ) : return cbeta_string;  break;
	case ( chemical::aa_trp ) : return ceps_2_string; break;
	case ( chemical::aa_tyr ) : return czeta_string;  break;
	default :
		utility_exit_with_message( "ERROR: Failed to find amino acid " + chemical::name_from_aa( aa ) + " in EnvSmooth::representative_atom_name" );
		break;
	}

	// unreachable
	return calpha_string;
}

void EnvToolMover::apply( core::pose::Pose &pose ) {
	ScoreFunction scfxn;
	scfxn.set_weight( scoring::envsmooth, 1.0 );
	scfxn( pose );

	pose.update_residue_neighbors();
	Size const nres( pose.size() );
	utility::vector1< Real > residue_N;

	// iterate over all the residues in the protein and count their neighbours
	// and save values of E, N, and dEdN
	for ( Size i = 1; i <= nres; ++i ) {

		// get the appropriate residue from the pose.
		conformation::Residue const & rsd( pose.residue(i) );
		// currently this is only for protein residues
		if ( ! rsd.is_protein() ) continue; //return;

		Size const atomindex_i = rsd.atom_index( representative_atom_name( rsd.aa() ));

		core::conformation::Atom const & atom_i = rsd.atom(atomindex_i);

		const Energies & energies( pose.energies() );
		const TwelveANeighborGraph & graph ( energies.twelveA_neighbor_graph() );

		Real countN    =  0.0;

		// iterate across neighbors within 12 angstroms
		for ( utility::graph::Graph::EdgeListConstIter
				ir  = graph.get_node(i)->const_edge_list_begin(),
				ire = graph.get_node(i)->const_edge_list_end();
				ir != ire; ++ir ) {
			Size const j( (*ir)->get_other_ind( i ) );
			conformation::Residue const & rsd_j( pose.residue(j) );
			Size atomindex_j( rsd_j.type().nbr_atom() );

			core::conformation::Atom const & atom_j = rsd_j.atom(atomindex_j);

			Real sqdist = atom_i.xyz().distance_squared(atom_j.xyz());
			countN += sigmoidish_neighbor( sqdist );
		}

		//Real score = 0;
		//Real dscoredN = 0;

		residue_N.push_back( countN );
		std::cout << "ENV:  " << i << "  " << countN << std::endl;
	}


}


void run() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !basic::options::option[ basic::options::OptionKeys::in::path::database ].user() ) {
		basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "/work/olange/minirosetta_database");
	}

	if ( option[ OptionKeys::constraints::no_linearize_bounded ] ) {
		tr.Info << "use fully harmonic potential for BOUNDED " << std::endl;
		ConstraintIO::get_func_factory().add_type("BOUNDED", new scoring::constraints::BoundFunc(0,0,0,1000,"dummy") );
	}
	if ( option[ OptionKeys::constraints::named ] ) {
		tr.Info << "use named constraints in AtomPairConstraint to avoid problems with cutpoint-variants " << std::endl;
		ConstraintIO::get_cst_factory().add_type( new scoring::constraints::NamedAtomPairConstraint( id::NamedAtomID(), id::NamedAtomID(), NULL) );
	}

	EnvToolMoverOP env_tool =  new EnvToolMover;
	protocols::jd2::JobDistributor::get_instance()->go( env_tool );

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
		register_options();
		devel::init( argc, argv );

		try{
			run();
		} catch (utility::excn::Exception& excn ) {
			excn.show( std::cerr );
		}
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


