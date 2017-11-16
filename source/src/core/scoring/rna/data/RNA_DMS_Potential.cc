// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/rna/data/RNA_DMS_Potential.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu

#include <core/scoring/rna/data/RNA_DMS_Potential.hh>
#include <core/pose/rna/RNA_DataInfo.hh>
#include <core/scoring/rna/data/util.hh>
#include <core/scoring/carbon_hbonds/CarbonHBondEnergy.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/PDBInfo.hh>

#include <numeric/conversions.hh>
#include <numeric/constants.hh>

#include <numeric/interpolation/InterpolatedPotential.tmpl.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <utility/tools/make_vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

static basic::Tracer TR( "core.scoring.rna.data.RNA_DMS_Potential" );

using utility::tools::make_vector1;
using utility::vector1;
using namespace basic::options;

//////////////////////////////////////////////////////
//
//              O
//              |
// DMS is CH3-O-S-O-CH3   (dimethyl sulfate)
//              |
//              O
//
// puts a methyl group at N1 of adenosine.
//
// For DMS (dimethyl sulfate) attack at adenosine, three features were discovered
//  to be useful (see rna_features app and analysis at git@github.com:DasLab/DeepChemicalProfiling.git )
//
//   *  whether the N1 is H-bonded/CH-bonded or not.
//   *  whether the space 1.5 A away from the N1 is occupied or not.
//   *  how good of a binding pocket is the space 3.5-5.5 A away from the N1 for
//        a methyl and an oxygen that would come from DMS.
//
// the flag -DMS_separate_features calculates 3 separate scores, as if the features were independent,
//   and adds them up.
//
// TO DO:
//  The following code is insanely compute-intensive -- used below to compute a potential
//   for binding the chemical probe. Need to replace with calls to atom-level
//   scoring terms for atr/rep, geom_sol, etc. Will definitely need to do that
//   to get derivatives!
//
//  Put in separate statistics for DMS of syn-A [?]. Currently giving back score = 0.0,
//   i.e. assuming that DMS distribution for syn-A matches that of total. But I think
//   syn-A typically has very low DMS.
//
//  Work out cytidine potential?
//
//  Replace lookup of potential with an interpolation. There are some great functions to do this
//   in numeric/interpolation/spline. In particular, TricubicSpline looks perfect.
//
//  Double-check whether we should include bin size when normalizing integrals below (e.g., DMS_stats_total); I think
//   these divide out, but I should check.
//
//  -- rhiju
//
//  AMW: Some stuff I've done here to make things a little more stable. The prior
//  implementation as nested vector1s doesn't guarantee e.g. that every 2nd-level
//  vector is the same length - missing data points can cause mis-alignments.
//  I think instead they should be zeroes (and interpolated if necessary). You're
//  also correct ot complain of the lack of a general grid object, and using a
//  pre-dimensioned array helps with that. It'll also help plug into eventual spline
//  work... so I'm going to use a MathNTensor (hey, I made that a couple years ago!).
//
//  I also don't love separately storing axis values and indices, but I understand
//  why it's important (and why a map would probably be a poor performance solution)
//  but I am going to store these values as a std::set (after all, the axis values
//  are, as you're checking, ascending and unique!)
//
//  --AMW
//
//////////////////////////////////////////////////////


namespace core {
namespace scoring {
namespace rna {
namespace data {

//Constructor
RNA_DMS_Potential::RNA_DMS_Potential():
	separate_scores_( option[ OptionKeys::score::rna::DMS_separate_features ]() ),
	occ_dist_( 1.5 ),
	methyl_probe_dist_( 3.5 ),
	oxygen_probe_dist_( 5.5 ),
	occ_shells_( std::make_pair( 2.0, 4.0 ) ),
	DMS_stats_(),
	DMS_potential_()
{
	probe_scorefxn_ = get_probe_scorefxn( true /*soft_rep*/, false /* just_atr_rep */);
}

//Destructor
RNA_DMS_Potential::~RNA_DMS_Potential()
{}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_DMS_Potential::initialize_DMS_potential() {

	numeric::MathNTensor< Real, 3 > not_bonded = read_DMS_stats_file( "scoring/rna/chem_map/dms/ade_N1_not_bonded_logstats.txt" );
	//Read in data file, and fill in private data.
	numeric::MathNTensor< Real, 3 > bonded = read_DMS_stats_file( "scoring/rna/chem_map/dms/ade_N1_bonded_logstats.txt" );

	utility::fixedsizearray1< Size, 4 > n_dimensions;
	n_dimensions[ 1 ] = 2;
	n_dimensions[ 2 ] = not_bonded.n_bins(1);
	n_dimensions[ 3 ] = not_bonded.n_bins(2);
	n_dimensions[ 4 ] = not_bonded.n_bins(3);

	DMS_stats_ = numeric::MathNTensor< Real, 4 >( n_dimensions, 0.0 );
	DMS_stats_.replace_layer( 0, not_bonded );
	DMS_stats_.replace_layer( 1, bonded );

	is_bonded_values_.insert( 0.0 );
	is_bonded_values_.insert( 1.0 );

	figure_out_potential(); // updates DMS_stats_, DMS_potential_, etc.
}

//////////////////////////////////////////////////////////////////////////////////
numeric::MathNTensor< Real, 3 >
RNA_DMS_Potential::read_DMS_stats_file( std::string const & potential_file ) {

	utility::io::izstream stream;
	basic::database::open( stream, potential_file );
	if ( !stream.good() ) utility_exit_with_message( "Unable to open "+potential_file );

	std::string line;

	// AMW:
	getline( stream, line );
	std::istringstream l1( line );
	utility::fixedsizearray1< Size, 3 > dimensions;
	l1 >> dimensions[1] >> dimensions[2] >> dimensions[3];

	numeric::MathNTensor< Real, 3 > DMS_stats( dimensions, 0.0 );

	// check labels
	getline( stream, line );
	std::istringstream l2( line );
	utility::vector1< std::string > labels( 4, "" );
	l2 >> labels[1] >> labels[2] >> labels[3] >> labels[4];
	runtime_assert( labels[1] == "occ" );
	runtime_assert( labels[2] == "Ebind" );
	runtime_assert( labels[3] == "DMS" );
	runtime_assert( labels[4] == "log-stats" );

	// Indices for the MathMatrix aren't obtained by lookup in vector
	// anymore. That formalism both assumes they're in increasing order (look
	// at lookup_idx!) and is only valuable if they might not be.
	// Instead, loop from 0 to dimensions[1], [2], [3]
	Size occ_idx = 0, binding_energy_idx = 0, DMS_idx = 0;
	Real occ, binding_energy, DMS, stats_value, log_stats_value;
	while ( getline( stream, line ) ) {
		if ( occ_idx > dimensions[1] ) {
			utility_exit_with_message( "The number of lines in a DMS file and the label indicating how many to expect are out of sync." );
		}
		std::istringstream l( line );
		l  >> occ >> binding_energy >> DMS >> log_stats_value;
		stats_value = exp( log_stats_value );

		// populate std::set axis values
		occ_values_.insert( occ );
		if ( occ_idx == 0 ) {
			binding_energy_values_.insert( binding_energy );
			if ( binding_energy_idx == 0 ) {
				DMS_values_.insert( DMS );
			}
		}

		DMS_stats( occ_idx, binding_energy_idx, DMS_idx ) = stats_value;

		// Increment and wrap as needed
		++DMS_idx;
		if ( DMS_idx == dimensions[3] ) {
			DMS_idx = 0;
			++binding_energy_idx;
			if ( binding_energy_idx == dimensions[2] ) {
				binding_energy_idx = 0;
				++occ_idx;
			}
		}
	}

	return DMS_stats;
}

//////////////////////////////////////////////////////////////////////////////////
// I could have done this externally, but I like recording
// the idea behind the log-odds score within Rosetta, for future
// reference. -- rhiju
void
RNA_DMS_Potential::figure_out_potential(){
	// indices are: is_bonded, DMS, occ, binding_energy
	//
	// log-odds score is: -kT log P( is_bonded, DMS, occ, binding_energy ) / [ P( DMS ) P( is_bonded, occ, binding_energy ) ]

	// first of all, need to normalize everything to total.
	Real const DMS_stats_total = DMS_stats_.sum();
	DMS_stats_ /= DMS_stats_total;

	p_DMS_ = numeric::MathVector< Real >( DMS_values_.size(), 0.0 );
	p_model_ = numeric::MathTensor< Real >(
		is_bonded_values_.size(),
		occ_values_.size(),
		binding_energy_values_.size(),
		0.0 ); // 'model' = (is_bonded, occ, binding_energy).

	// for separate feature-scores [ not on by default ]
	numeric::MathVector< Real > p_is_bonded( is_bonded_values_.size(), 0.0 );
	numeric::MathVector< Real > p_occ( occ_values_.size(), 0.0 );
	numeric::MathVector< Real > p_binding_energy( binding_energy_values_.size(), 0.0 );

	numeric::MathMatrix< Real > p_is_bonded_DMS( is_bonded_values_.size(), DMS_values_.size(), 0.0 );
	numeric::MathMatrix< Real > p_occ_DMS( occ_values_.size(), DMS_values_.size(), 0.0 );
	numeric::MathMatrix< Real > p_binding_energy_DMS( binding_energy_values_.size(), DMS_values_.size(), 0.0 );

	// AMW todo: implement "slice sums" for cross sections along the axes.
	// fill projections, which give denominator of log-odds score.
	for ( Size h = 0; h < is_bonded_values_.size(); h++ ) {
		for ( Size i = 0; i < occ_values_.size(); i++ ) {
			for ( Size j = 0; j < binding_energy_values_.size(); j++ ) {
				for ( Size k = 0; k < DMS_values_.size(); k++ ) {

					Real const val = DMS_stats_( h, i, j, k );
					p_DMS_( k )                    += val;
					p_model_( h, i, j )        += val;

					// for separate feature scores [not on by default]
					p_is_bonded( h )               += val;
					p_is_bonded_DMS( h, k )      += val;
					p_occ( i )                     += val;
					p_occ_DMS( i, k )            += val;
					p_binding_energy( j )          += val;
					p_binding_energy_DMS( j, k ) += val;
				}
			}
		}
	}

	//DMS_potential_ = DMS_stats_; // values will be replaced
	DMS_potential_ = DMS_stats_; // values will be replaced
	for ( Size h = 0; h < is_bonded_values_.size(); h++ ) {
		for ( Size i = 0; i < occ_values_.size(); i++ ) {
			for ( Size j = 0; j < binding_energy_values_.size(); j++ ) {
				for ( Size k = 0; k < DMS_values_.size(); k++ ) {
					DMS_potential_( h, i, j, k ) =
						-1.0 * log( DMS_stats_( h, i, j, k )/ ( p_model_( h, i, j ) * p_DMS_( k ) ) );
					//std::cout << "potential " << h << " " << i << " " <<  j << " " << k << " is " <<  DMS_potential_( h, i, j, k ) << std::endl;
					//std::cout << "stats is " <<  DMS_stats_( h, i, j, k ) << std::endl;
					//std::cout << "p_model_ is " <<  p_model_( h, i, j ) << std::endl;
					//std::cout << "p_DMS_ is " <<  p_DMS_( k ) << std::endl;
				}
			}
		}
	}

	// separate feature scores [not on by default]
	DMS_potential_is_bonded_ = p_is_bonded_DMS; // values will be replaced
	for ( Size h = 0; h < is_bonded_values_.size(); h++ ) {
		for ( Size k = 0; k < DMS_values_.size(); k++ ) {
			DMS_potential_is_bonded_( h, k ) =
				-1.0 * log( p_is_bonded_DMS( h, k ) /( p_is_bonded( h ) * p_DMS_( k )) );
		}
	}
	DMS_potential_occ_ = p_occ_DMS; // values will be replaced
	for ( Size i = 0; i < occ_values_.size(); i++ ) {
		for ( Size k = 0; k < DMS_values_.size(); k++ ) {
			DMS_potential_occ_( i, k ) =
				-1.0 * log( p_occ_DMS( i, k ) / ( p_occ( i ) * p_DMS_( k ) ) );
		}
	}
	DMS_potential_binding_energy_ = p_binding_energy_DMS; // values will be replaced
	for ( Size j = 0; j < binding_energy_values_.size(); j++ ) {
		for ( Size k = 0; k < DMS_values_.size(); k++ ) {
			DMS_potential_binding_energy_( j, k ) =
				-1.0 * log( p_binding_energy_DMS( j, k ) /( p_binding_energy( j ) * p_DMS_( k ) ) );
		}
	}

	using namespace numeric::interpolation::spline;
	// Train polycubic spline on the potential.
	utility::fixedsizearray1< BorderFlag, 4 > const BORDER( e_Natural );

	utility::fixedsizearray1< double, 4 > START;
	START[1] = *is_bonded_values_.begin();
	START[2] = *occ_values_.begin();
	START[3] = *binding_energy_values_.begin();
	START[4] = *DMS_values_.begin();

	utility::fixedsizearray1< double, 4 >  DELTA;
	DELTA[1] = *std::next(is_bonded_values_.begin()) - *is_bonded_values_.begin();
	DELTA[2] = *std::next(occ_values_.begin()) - *occ_values_.begin();
	DELTA[3] = *std::next(binding_energy_values_.begin()) - *binding_energy_values_.begin();
	DELTA[4] = *std::next(DMS_values_.begin()) - *DMS_values_.begin();

	utility::fixedsizearray1< bool, 4 > const LINCONT( true );
	utility::fixedsizearray1< std::pair< Real, Real >, 4 > const FIRSTBE( std::make_pair( 0.0, 0.0 ) );

	numeric::interpolation::spline::PolycubicSpline< 4 > pcs;
	pcs.train( BORDER, START, DELTA, DMS_potential_, LINCONT, FIRSTBE );

	utility::fixedsizearray1< Size, 4 > dims;
	dims[1] = is_bonded_values_.size();
	dims[2] = occ_values_.size();
	dims[3] = binding_energy_values_.size();
	dims[4] = DMS_values_.size();
	interpolated_potential_.dimension( dims );
	interpolated_potential_.set_bin_width( DELTA );
	interpolated_potential_.set_periodic( false );

	for ( Size h = 0; h < is_bonded_values_.size(); h++ ) {
		for ( Size i = 0; i < occ_values_.size(); i++ ) {
			for ( Size j = 0; j < binding_energy_values_.size(); j++ ) {
				for ( Size k = 0; k < DMS_values_.size(); k++ ) {
					utility::fixedsizearray1< Size, 4 > indices;
					indices[1] = h; indices[2] = i; indices[3] = j; indices[4] = k;
					debug_assert( interpolated_potential_( h, i, j, k ).size() == pcs.get_all_derivs( indices ).size() );
					//for ( Size mm = 1; mm <= interpolated_potential_( h, i, j, k ).size(); ++mm ) {
					interpolated_potential_( h, i, j, k ) = pcs.get_all_derivs( indices );
					//}
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_DMS_Potential::initialize( core::pose::Pose const & pose ) {
#ifdef MULTI_THREADED
	utility_exit_with_message("The RNA_DMS_Potential is fundamentally non-threadsafe, and cannot be used in a multi-threaded context.");
#endif

	// calculate H-bonds -- this is currently very inefficient, as we
	// really only need to look at Hbond status of adenosine N1's.
	//
	hbonds::HBondOptionsOP hbond_options( new hbonds::HBondOptions() );
	hbond_options->use_hb_env_dep( false );
	hbond_set_ = hbonds::HBondSetOP( new hbonds::HBondSet( *hbond_options ) );
	hbonds::fill_hbond_set( pose, false /*calc deriv*/, *hbond_set_ );

	// setup poses in which one residue base can be virtualized, and
	//  with a 'probe' atom that will be used to calculate occupancy/binding energies
	//  for a mock DMS molecule. This is also super-inefficient.
	working_pose_            = pose.clone(); // for add/delete probe residues.
	working_pose_with_probe_ = pose.clone(); // for add/delete probe residues.
	add_probe_to_pose( *working_pose_with_probe_ );

	// clone() above actually does not reset scoring_ flag in Energies object, leading to early exit.
	working_pose_->set_new_energies_object( scoring::EnergiesOP( new scoring::Energies ) );
	working_pose_with_probe_->set_new_energies_object( scoring::EnergiesOP( new scoring::Energies ) );

	// I think this needs to be initialized ...
	if ( DMS_potential_.size() == 0 ) initialize_DMS_potential();

}

////////////////////////////////////////////////////////////////////////////////////
bool
RNA_DMS_Potential::get_features( pose::Pose const & pose,
	Size const i,
	bool & ade_n1_bonded,
	Real & binding_energy,
	Real & occupancy_density )
{
	runtime_assert ( pose.residue( i ).aa() == core::chemical::na_rad );

	// no syn-adenosines (some of these are exposed but not DMS-reactive -- chemical understanding
	// is currently incomplete.
	if ( pose.chi( i ) < 0.0 ) return false;

	// is the N1 sequestered in a hydrogen bond?
	bool ade_n1_hbonded  = check_hbonded(  pose, i, " N1 ", true /*acceptor*/ );
	bool ade_n1_chbonded = check_chbonded( pose, i, " N1 " );
	ade_n1_bonded = ade_n1_hbonded || ade_n1_chbonded;

	Vector const occ_xyz = get_probe_xyz( pose.residue( i ), occ_dist_ /* 1.5 A */);
	occupancy_density = get_occupancy_density( pose, i /*for exclusion*/, occ_xyz, occ_shells_ );

	Vector const methyl_probe_xyz = get_probe_xyz( pose.residue( i ), methyl_probe_dist_ /*3.5 A */ );
	Real const methyl_binding_energy = get_binding_energy( i, methyl_probe_xyz, *probe_scorefxn_ );

	Vector const oxygen_probe_xyz = get_probe_xyz( pose.residue( i ), oxygen_probe_dist_ /*5.5 A */ );
	Real const oxygen_binding_energy = get_binding_energy( i, oxygen_probe_xyz, *probe_scorefxn_ );

	binding_energy = methyl_binding_energy + oxygen_binding_energy;
	return true;
}

//////////////////////////////////////////////////////////////////////////////////
Real
RNA_DMS_Potential::evaluate( core::pose::Pose const & pose,
	pose::rna::RNA_Reactivity const & rna_reactivity ) {

#ifdef MULTI_THREADED
	utility_exit_with_message("The RNA_DMS_Potential is fundamentally non-threadsafe, and cannot be used in a multi-threaded context.");
#endif

	using namespace core::pose::full_model_info;

	if ( DMS_potential_.size() == 0 ) initialize_DMS_potential();

	runtime_assert( rna_reactivity.type() == pose::rna::DMS );
	Size const & pos = rna_reactivity.position();
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	if ( pos < 1 || pos >= full_model_info.full_sequence().size() || full_model_info.full_sequence()[ pos - 1 ] != 'a' ) return 0.0;
	if ( !full_model_info.working_res().has_value( pos ) )   return 0.0;

	// initialize feature values to what would be appropriate for a bulged A
	if ( !full_model_info.res_list().has_value( pos ) ) return 0.0;

	// AMW: note that our framework is robust to the possibility that "ade_n1_bonded"
	// becomes a Real to reflect a continuous confidence in whether the N1 is
	// H-bonded.

	Size const i = full_model_info.full_to_sub( pos );
	if ( pose.residue( i ).has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) return 0.0;
	if ( pose.residue( i ).has_variant_type( chemical::REPLONLY ) ) return 0.0;

	bool ade_n1_bonded( false );
	Real binding_energy( 0.0 ), occupancy_density( 0.0 );
	bool const success = get_features( pose, i, ade_n1_bonded, binding_energy, occupancy_density );
	if ( !success ) return 0.0;

	Real score( 0.0 );
	if ( true ) {
		Size const n1_bond_idx = ade_n1_bonded ? 1 : 0;
		Size const DMS_idx = get_idx( rna_reactivity.value(), DMS_values_ );
		Size const occ_idx = get_idx( occupancy_density, occ_values_ );
		Size const binding_energy_idx = get_idx( binding_energy, binding_energy_values_ );

		if ( separate_scores_ ) {
			Real score_n1_bond        = DMS_potential_is_bonded_     ( n1_bond_idx, DMS_idx );
			Real score_occ            = DMS_potential_occ_           ( occ_idx, DMS_idx );
			Real score_binding_energy = DMS_potential_binding_energy_( binding_energy_idx, DMS_idx );
			score = score_n1_bond + score_occ + score_binding_energy;
		} else {
			score = DMS_potential_( n1_bond_idx, occ_idx,binding_energy_idx, DMS_idx );
			//std::cout << "score " << score << " " << n1_bond_idx << " " <<  occ_idx << " " << binding_energy_idx << " " << DMS_idx << std::endl;
			//std::cout << "near vals " << occupancy_density << " " << binding_energy << " " << rna_reactivity.value() << std::endl;
		}
	} else {
		utility::fixedsizearray1< Real, 4 > dscoredfeat;
		utility::fixedsizearray1< Real, 4 > values;
		values[1] = ade_n1_bonded;
		values[2] = rna_reactivity.value();
		values[3] = occupancy_density;
		values[4] = binding_energy;
		numeric::interpolation::polycubic_interpolation( interpolated_potential_, values, score, dscoredfeat );
	}

	//  TR <<  pose.pdb_info()->number(i) << " ade_n1_bonded " << ade_n1_bonded << " (" << n1_bond_idx << ")" << "   occupancy " << occupancy_density << " (" << occ_idx << ")" << "   binding_energy " << binding_energy << " (" << binding_energy_idx << ") " << "  value " << rna_reactivity.value() << " (" << DMS_idx << ")" << " SCORE " << score << std::endl;

	return score;
}


///////////////////////////////////////////////////////////////////////
bool
RNA_DMS_Potential::check_hbonded( pose::Pose const & pose,
	Size const & i /*residue number*/,
	std::string const & atom_name,
	bool is_acceptor ) const {

	using namespace core::scoring::hbonds;
	// need to know what atoms could be in this nt.
	// go through list of donors and list of acceptors

	// check in HBond List
	for ( Size n = 1; n <= hbond_set_->nhbonds(); n++ ) {
		HBond const & hbond( hbond_set_->hbond( n ) );
		if ( is_acceptor && hbond.acc_res() == i && pose.residue_type(i).atom_name( hbond.acc_atm() ) == atom_name )   return true;
		if ( !is_acceptor && hbond.don_res() == i && pose.residue_type(i).atom_name( hbond.don_hatm() ) == atom_name )   return true;
	}
	return false;

}

///////////////////////////////////////////////////////////////////////////////
Real
RNA_DMS_Potential::get_N1_lonepair_donor_angle(  core::conformation::Residue const & acc_rsd,
	core::conformation::Residue const & don_rsd,
	Size const don_h_atm ) const {
	runtime_assert( acc_rsd.has( " N1 " ) );
	runtime_assert( acc_rsd.aa() == core::chemical::na_rad ); // for now
	core::Vector u = ( acc_rsd.xyz( " N1 " ) - 0.5 * ( acc_rsd.xyz( " C6 " ) + acc_rsd.xyz( " C2 " ) ) ).normalized();
	core::Vector v = ( don_rsd.xyz( don_h_atm  ) - acc_rsd.xyz( " N1 " ) ).normalized();
	return numeric::conversions::degrees( std::acos( dot( u, v ) ) );
}


///////////////////////////////////////////////////////////////////////////////
bool
RNA_DMS_Potential::check_chbonded( pose::Pose const & pose,
	Size const & i /*residue number*/,
	std::string const & atom_name ) const {

	static core::scoring::carbon_hbonds::CarbonHBondEnergy const carbon_hbond_energy;

	core::conformation::Residue const & acc_rsd = pose.residue( i );
	Size const acc_atm( acc_rsd.atom_index( atom_name ) );
	static Real const ENERGY_CUTOFF = -0.4;
	bool found_chbond( false );
	for ( Size j = 1; j <= pose.size(); j++ ) {

		if ( i == j ) continue;

		Real energy( 0.0 ), angle( 0.0 ), res_res_energy( 0.0 );
		core::conformation::Residue const & don_rsd = pose.residue( j );
		for ( Size const don_h_atm : don_rsd.Hpos_apolar() ) {

			/// square-distance check outside the call to get_atom_atom_carbon_hbond_energy
			if ( don_rsd.xyz( don_h_atm ).distance_squared( acc_rsd.xyz( acc_atm ) ) > carbon_hbond_energy.max_dis2() ) continue;
			if ( carbon_hbond_energy.get_atom_atom_carbon_hbond_energy( don_h_atm, don_rsd,
					acc_atm, acc_rsd, energy ) ) {

				// wait, we should probably set an angle based on lone-pair direction.
				// This should probably be *inside* the CH-bond potential.
				angle = get_N1_lonepair_donor_angle( acc_rsd, don_rsd, don_h_atm );
				if ( angle <= 90.0 ) {
					res_res_energy += energy;
					if ( energy < ENERGY_CUTOFF ) found_chbond = true;
				}
			}
		}
	}
	return found_chbond;
}

///////////////////////////////////////////////////////////////////////////////
core::Vector
RNA_DMS_Potential::get_probe_xyz( core::conformation::Residue const & rsd, Distance const probe_dist ) const {
	runtime_assert( rsd.has( " N1 " ) );
	runtime_assert( rsd.name1() == 'a' ); // for now
	// define unit vector.
	core::Vector u = ( rsd.xyz( " N1 " ) - 0.5 * ( rsd.xyz( " C6 " ) + rsd.xyz( " C2 " ) ) ).normalized();
	return rsd.xyz( " N1 " ) + u * probe_dist;
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_DMS_Potential::get_occupancy_densities( utility::vector1< Real > & occupancy_densities,
	pose::Pose const & pose,
	Size const i /*for exclusion*/,
	core::Vector const & probe_xyz,
	utility::vector1< Distance > const & shells ) const {

	utility::vector1< Size > num_atoms_in_shells( shells.size() - 1, 0 );
	Real const max_distance = shells[ shells.size() ];
	for ( Size j = 1; j <= pose.size(); j++ ) {
		if ( i == j ) continue;
		core::conformation::Residue const & rsd = pose.residue( j );
		for ( Size jj = 1; jj <= rsd.natoms(); jj++ ) {
			if ( rsd.is_virtual( jj ) ) continue;
			Distance dist = ( rsd.xyz( jj ) - probe_xyz ).length();
			if ( dist < max_distance ) {
				Size k( 1 );
				for ( k = 1; k < shells.size(); k++ ) {
					if ( dist < shells[k+1] ) break;
				}
				num_atoms_in_shells[ k ]++;
			}
		}
	}

	occupancy_densities.clear();
	for ( Size k = 1; k < shells.size(); k++ ) {
		Real const shell_volume = (4 / 3) * numeric::constants::d::pi * ( pow( shells[k+1], 3 ) - pow( shells[k], 3 ) );
		occupancy_densities.push_back( static_cast<Real>( num_atoms_in_shells[ k ] )/ shell_volume );
	}
}


///////////////////////////////////////////////////////////////////////////////
Real
RNA_DMS_Potential::get_occupancy_density( pose::Pose const & pose,
	Size const i /*for exclusion*/,
	core::Vector const & probe_xyz,
	std::pair< Distance, Distance > const & shells_pair ) const {
	vector1< Distance > shells = make_vector1( shells_pair.first, shells_pair.second );
	vector1< Distance > occupancy_densities;
	get_occupancy_densities( occupancy_densities, pose, i, probe_xyz, shells );
	return occupancy_densities[ 1 ];
}

///////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionOP
RNA_DMS_Potential::get_probe_scorefxn( bool const soft_rep, bool const just_atr_rep ) const {

	using namespace core::scoring;
	using namespace core::scoring::methods;

	ScoreFunctionOP scorefxn( new ScoreFunction );
	EnergyMethodOptions energy_method_options = scorefxn->energy_method_options();
	energy_method_options.hbond_options().use_hb_env_dep( false );
	if ( soft_rep ) energy_method_options.etable_type( FA_STANDARD_SOFT );
	scorefxn->set_energy_method_options( energy_method_options );
	scorefxn->set_weight( fa_atr, 0.21 );
	scorefxn->set_weight( fa_rep, 0.20 );
	if ( !just_atr_rep ) {
		scorefxn->set_weight( fa_stack, 0.13 ); // ?
		scorefxn->set_weight( hbond_sc, 0.17 ); // needed for geom_sol in special cases incl. RNA/protein. a little ridiculous.
		scorefxn->set_weight( geom_sol_fast, 0.17 );
		scorefxn->set_weight( lk_nonpolar, 0.25 );
	}
	return scorefxn;
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_DMS_Potential::update_virtual_base_if_necessary( pose::Pose & pose, Size const i ){
	// virtualize adenosine where probe would stick, or will get weird clashes. Note also that
	// this should be constant contribution to all A's.
	add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_BASE, i );
	for ( Size j = 1; j <= pose.size(); j++ ) {
		if ( i == j ) continue;
		remove_variant_type_from_pose_residue( pose, chemical::VIRTUAL_BASE, j );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_DMS_Potential::add_probe_to_pose( pose::Pose & pose ){
	Size const i = pose.size();
	core::chemical::ResidueTypeSetCOP rsd_set = pose.residue_type_set_for_pose( pose.residue_type( i ).mode() );
	core::chemical::ResidueTypeCOP rsd_type ( rsd_set->get_representative_type_name3 ( " CZ" ) ); // just a carbon atom.
	core::conformation::ResidueOP probe_res = ( core::conformation::ResidueFactory::create_residue ( *rsd_type ) );
	pose.append_residue_by_jump( *probe_res, i );
}

///////////////////////////////////////////////////////////////////////////////
Real
RNA_DMS_Potential::get_binding_energy( Size const i,
	core::Vector const & probe_xyz,
	core::scoring::ScoreFunction const & scorefxn ){

	core::pose::Pose & pose            = *working_pose_;
	core::pose::Pose & pose_with_probe = *working_pose_with_probe_;

	update_virtual_base_if_necessary( pose, i );
	update_virtual_base_if_necessary( pose_with_probe, i );

	Real const score_start = ( scorefxn )( pose );

	Size const probe_res = pose_with_probe.size();
	runtime_assert( pose_with_probe.residue( probe_res ).name3() == " CZ" ); // probe.
	core::conformation::Residue const & probe_rsd = pose_with_probe.residue( pose_with_probe.size() );
	core::Vector const start_xyz = probe_rsd.xyz( " C1 " );
	core::Vector const translation = probe_xyz - start_xyz;
	for ( Size k = 1; k <= probe_rsd.natoms(); k++ ) pose_with_probe.set_xyz( id::AtomID( k, probe_res ), probe_rsd.xyz( k ) + translation );
	Real const score_probe = (scorefxn)( pose_with_probe );

	Real const binding_energy = score_probe - score_start;
	// std::cout << "checking binding energy at " << seqpos << " : " << binding_energy << std::endl;

	return binding_energy;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// TO DO: should interpolate DMS_stats.
utility::vector1< Real >
RNA_DMS_Potential::get_logL_values( pose::Pose const & pose, Size const i /*, utility::vector1< Real > const & DMS_values */  ){

	runtime_assert( working_pose_ != 0 );
	runtime_assert( working_pose_with_probe_ != 0 );

	if ( pose.aa( i ) != chemical::na_rad ) {
		return vector1< Real >( DMS_values_.size(), 0.0 );
	}

	vector1< Real > logL_values;
	bool ade_n1_bonded( false );
	Real binding_energy( 0.0 ), occupancy_density( 0.0 );
	bool const success = get_features( pose, i, ade_n1_bonded, binding_energy, occupancy_density );

	if ( success ) {
		Size const n1_bond_idx = ade_n1_bonded ? 1 : 0;
		//get_bool_idx( ade_n1_bonded, is_bonded_values_ );
		Size const occ_idx = get_idx( occupancy_density, occ_values_ );
		Size const binding_energy_idx = get_idx( binding_energy, binding_energy_values_ );
		for ( Size k = 1; k <= DMS_values_.size(); k++ ) {
			logL_values.push_back( log( DMS_stats_( n1_bond_idx, occ_idx, binding_energy_idx, k ) /
				p_model_( n1_bond_idx, occ_idx, binding_energy_idx ) ) /*for normalization*/ );
		}
	} else {
		// return generic DMS distribution.
		for ( Size k = 1; k <= DMS_values_.size(); k++ ) {
			logL_values.push_back( log( p_DMS_( k ) ) );
		}
	}

	return logL_values;
}


} //data
} //rna
} //scoring
} //core
