// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.cc
/// @brief  ScoreFunction class definition.
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Modified by Sergey Lyskov


// Unit headers
#include <core/scoring/ScoreFunction.hh>

// C/C++ headers
#include <algorithm>
#include <string>
#include <vector>

// External headers
#include <boost/algorithm/string/join.hpp>

// Package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/methods/FreeDOF_Options.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParam.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <core/scoring/Energies.hh>

// Project headers
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/id/TorsionID.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/Pose.hh>

/// Utility headers
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#ifdef PYROSETTA
#include <basic/init.hh>
#endif

#include <numeric/random/DistributionSampler.hh>
#include <ObjexxFCL/format.hh>
#include <utility/io/izstream.hh>

#include <core/id/DOF_ID.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

static THREAD_LOCAL basic::Tracer tr( "core.scoring" );

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

#ifdef APL_TEMP_DEBUG
Size n_minimization_sfxn_evals( 0 );
#endif

///////////////////////////////////////////////////////////////////////////////
ScoreFunction::ScoreFunction()
{
	reset();
}


ScoreFunction::~ScoreFunction() {}

///////////////////////////////////////////////////////////////////////////////
ScoreFunctionOP
ScoreFunction::clone() const
{
	ScoreFunctionOP new_score_function( new ScoreFunction );
	new_score_function->assign(*this);
	return new_score_function;
}

///////////////////////////////////////////////////////////////////////////////
ScoreFunctionOP
ScoreFunction::clone_as_base_class() const
{
	ScoreFunctionOP new_score_function( new ScoreFunction );
	new_score_function->assign(*this);
	return new_score_function;
}

///////////////////////////////////////////////////////////////////////////////
/// @details reset function that can be called either independently or from default
/// ScoreFunction creation
void
ScoreFunction::reset()
{
#ifdef PYROSETTA
		// Sanity check: check if core::init was called already and abort otherwise with helpful message...
		if( !basic::was_init_called() ) utility_exit_with_message("Attempt to initialize ScoreFunction object before core::init was called detected… Have you forgot to call core::init?");
#endif

	score_function_info_current_ = true;
	score_function_info_ = ScoreFunctionInfoOP( new ScoreFunctionInfo );
	any_intrares_energies_ = false;
	energy_method_options_ = methods::EnergyMethodOptionsOP( new methods::EnergyMethodOptions );
	initialize_methods_arrays();
	weights_.clear();
}

///////////////////////////////////////////////////////////////////////////////
/// read info from file


void
ScoreFunction::initialize_from_file( std::string const & filename )
{
	reset();

	add_weights_from_file( filename );
}

/// Format: { term : weight, ... }
std::string ScoreFunction::serialize_weights() const {
	using std::string;
	using std::vector;

	vector<string> tokens;
	ScoreTypes terms = get_nonzero_weighted_scoretypes();
	for ( ScoreTypes::const_iterator i = terms.begin(); i != terms.end(); ++i ) {
		const ScoreType& term = *i;
		if ( has_nonzero_weight(term) ) {
			tokens.push_back(str(boost::format("%1% : %2%") % term % get_weight(term)));
		}
	}
	return "{" + boost::algorithm::join(tokens, ", ") + "}";
}

/// @detail Perturbs each non-zero weight independently by adding a Gaussian
/// noise with u=0 and sd=weight / 8. Weights are not allowed to become negative.
///   - 68% of the time, percent change [0, 12.5)
///   - 26% of the time, percent change [12.5, 25)
///   -  4% of the time, percent change [25...
void ScoreFunction::perturb_weights() {
	using boost::math::normal;
	using core::Real;
	using numeric::random::DistributionSampler;

	// Perturb the weights of the non-zero score terms
	ScoreTypes terms = get_nonzero_weighted_scoretypes();
	for ( ScoreTypes::const_iterator i = terms.begin(); i != terms.end(); ++i ) {
		const ScoreType& term = *i;
		if ( has_nonzero_weight(term) ) {
			Real weight = get_weight(term);

			normal dist(0, weight / 8);
			DistributionSampler<normal> sampler(dist);

			Real perturbed_weight = weight + sampler.sample();
			set_weight(term, std::max(perturbed_weight, 0.0));
			tr.Debug << name_from_score_type(term) << ": "
				<< weight << " => " << perturbed_weight << std::endl;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// read weights/etc from file, which may be in database directory.
/// Does not clear weights beforehand.

void
ScoreFunction::add_weights_from_file( std::string const & filename )
{
	_add_weights_from_file( find_weights_file(filename, ".wts") );
}

///////////////////////////////////////////////////////////////////////////////
/// read weights/etc from file. Does not clear weights beforehand.
/// no lookup in database directory

void
ScoreFunction::_add_weights_from_file( std::string const & filename, bool patch/*=false*/ )
{
	utility::io::izstream data( filename );
	if ( !data.good() ) {
		if ( ! patch ) { utility_exit_with_message( "Unable to open weights file: "+filename ); }
		else { utility_exit_with_message( "Unable to open weights-patch file: "+filename ); }
	}

	std::string line,tag,operation;
	ScoreType score_type;
	Real wt;
	Real real_value;
	Size size_value;
	std::string string_value;
	utility::vector1< std::pair< ScoreType, Real > >  score_types_from_file;
	utility::vector1< std::pair< std::pair< ScoreType, Real >, std::string > > patches_from_file;

	while ( getline( data, line ) ) {
		std::istringstream l( line );
		l >> tag;
		if ( l.fail() || tag[0] == '#' ) continue;

		// //////////// Property Setting Tags ///////////////////////
		if ( tag == "ETABLE" ) {
			l >> tag;
			set_etable( tag );
		} else if ( tag == "METHOD_WEIGHTS" ) {
			l >> score_type;
			if ( l.fail() ) {
				utility_exit_with_message( "unrecognized score_type: "+line );
			}
			utility::vector1< Real > wts;
			while ( !l.fail() ) {
				l >> wt;
				if ( l.fail() ) break;
				wts.push_back( wt );
			}
			energy_method_options_->set_method_weights( score_type, wts );
		} else if ( tag == "STRAND_STRAND_WEIGHTS" ) {
			utility::vector1< int > values;
			int value;
			while ( !l.fail() ) {
				l >> value;
				if ( l.fail() ) break;
				values.push_back( value );
			}
			if ( values.size() != 2 ) {
				utility_exit_with_message( "incorrect number of arguments to STRAND_STRAND_WEIGHTS: " + line );
			}
			energy_method_options_->set_strand_strand_weights( values[1], values[2] );
		} else if ( tag == "NO_PROTEIN_PROTEIN_FA_ELEC" ) {
			energy_method_options_->exclude_protein_protein_fa_elec( true );
		} else if ( tag == "NO_MONOMER_FA_ELEC" ) {
			energy_method_options_->exclude_monomer_fa_elec( true );
		} else if ( tag == "INCLUDE_DNA_DNA" ) {
			energy_method_options_->exclude_DNA_DNA( false );
		} else if ( tag == "NO_HB_ENV_DEP" ) {
			energy_method_options_->hbond_options().use_hb_env_dep( false );
		} else if ( tag == "NO_HB_ENV_DEP_DNA" ) {
			energy_method_options_->hbond_options().use_hb_env_dep_DNA( false );
		} else if ( tag == "NO_SMOOTH_HB_ENV_DEP" ) {
			energy_method_options_->hbond_options().smooth_hb_env_dep( false );
		} else if ( tag == "NO_BB_DONOR_ACCEPTOR_CHECK" ) {
			energy_method_options_->hbond_options().bb_donor_acceptor_check( false );
		} else if ( tag == "COARSE_RNA" ) {
			std::cout << "ATOM_VDW set to COARSE_RNA" << std::endl;
			energy_method_options_->atom_vdw_atom_type_set_name( "coarse_rna" );
		} else if ( tag == "CENTROID_ROT" ) {
			tr << "ATOM_VDW set to CENTROID_ROT" << std::endl;
			energy_method_options_->atom_vdw_atom_type_set_name( "centroid_rot" );
		} else if ( tag == "CENTROID_ROT_MIN" ) {
			tr << "ATOM_VDW set to CENTROID_ROT_MIN" << std::endl;
			energy_method_options_->atom_vdw_atom_type_set_name( "centroid_rot:min" );
		} else if ( tag == "INCLUDE_INTRA_RES_RNA_HB" ) {
			utility_exit_with_message( "INCLUDE_INTRA_RES_RNA_HB no longer an option. Just use hbond_intra or hbond score term..." );
		} else if ( tag == "EXCLUDE_INTRA_RES_RNA_HB" ) {
			energy_method_options_->hbond_options().exclude_intra_res_RNA( true );
		} else if ( tag == "INCLUDE_INTRA_RES_PROTEIN" ) {
			energy_method_options_->exclude_intra_res_protein( false );
		} else if ( tag == "PUT_INTRA_INTO_TOTAL" ) {
			energy_method_options_->put_intra_into_total( true ); // also updates hbond_options
		} else if ( tag == "INCLUDE_HB_DNA_DNA" ) {
			energy_method_options_->hbond_options().exclude_DNA_DNA( false );
		} else if ( tag == "FA_MAX_DIS" ) {
			l >> real_value;
			energy_method_options_->etable_options().max_dis = real_value;
			reset_energy_methods(); // ensure that etable is recomputed.
			score_function_info_current_ = false;
		} else if ( tag == "ENLARGE_H_LJ_WDEPTH" ) {
			energy_method_options_->etable_options().enlarge_h_lj_wdepth = true;
			reset_energy_methods(); // ensure that etable is recomputed.
			score_function_info_current_ = false;
		} else if ( tag == "BOND_ANGLE_CENTRAL_ATOMS_TO_SCORE" ) {
			utility::vector1<std::string> central_atoms;
			std::string central_atom;
			while ( !l.fail() ) {
				l >> central_atom;
				if ( l.fail() ) break;
				central_atoms.push_back( central_atom );
			}
			if ( ! energy_method_options_->bond_angle_residue_type_param_set() ) {
				energy_method_options_->bond_angle_residue_type_param_set( scoring::mm::MMBondAngleResidueTypeParamSetOP( new core::scoring::mm::MMBondAngleResidueTypeParamSet() ) );
			}
			energy_method_options_->bond_angle_residue_type_param_set()->central_atoms_to_score( central_atoms );
		} else if ( tag == "BOND_ANGLE_USE_RESIDUE_TYPE_THETA0" ) {
			if ( ! energy_method_options_->bond_angle_residue_type_param_set() ) {
				energy_method_options_->bond_angle_residue_type_param_set( scoring::mm::MMBondAngleResidueTypeParamSetOP( new core::scoring::mm::MMBondAngleResidueTypeParamSet() ) );
			}
			energy_method_options_->bond_angle_residue_type_param_set()->use_residue_type_theta0( true );
		} else if ( tag == "UNFOLDED_ENERGIES_TYPE" ) {
			std::string type;
			l >> type;
			energy_method_options_->unfolded_energies_type( type );
		} else if ( tag == "SPLIT_UNFOLDED_LABEL_TYPE" ) {
			std::string label_type;
			l >> label_type;
			energy_method_options_->split_unfolded_label_type( label_type );
		} else if ( tag == "SPLIT_UNFOLDED_VALUE_TYPE" ) {
			std::string value_type;
			l >> value_type;
			energy_method_options_->split_unfolded_value_type( value_type );
		} else if ( tag == "FA_ELEC_MIN_DIS" ) {
			l >> real_value;
			energy_method_options_->elec_min_dis( real_value );
		} else if ( tag == "FA_ELEC_MAX_DIS" ) {
			l >> real_value;
			energy_method_options_->elec_max_dis( real_value );
		} else if ( tag == "FA_ELEC_NO_DIS_DEP_DIE" ) {
			energy_method_options_->elec_no_dis_dep_die( true );
		} else if ( tag == "NO_LK_POLAR_DESOLVATION" ) {
			energy_method_options_->etable_options().no_lk_polar_desolvation = true;
			reset_energy_methods(); // ensure that etable is recomputed.
			score_function_info_current_ = false;
		} else if ( tag == "GEOM_SOL_INTERRES_PATH_DISTANCE_CUTOFF" ) {
			l >> size_value;
			energy_method_options_->geom_sol_interres_path_distance_cutoff( size_value );
		} else if ( tag == "GEOM_SOL_INTRARES_PATH_DISTANCE_CUTOFF" ) {
			l >> size_value;
			energy_method_options_->geom_sol_intrares_path_distance_cutoff( size_value );
		} else if ( tag == "RNA_SYN_G_POTENTIAL_BONUS" ) {
			l >> real_value;
			energy_method_options_->rna_options().syn_G_potential_bonus( real_value );
		} else if ( tag == "RNA_TORSION_POTENTIAL" ) {
			l >> string_value;
			energy_method_options_->rna_options().torsion_potential( string_value );
		} else if ( tag == "RNA_SUITENESS_BONUS" ) {
			l >> string_value;
			energy_method_options_->rna_options().suiteness_bonus( string_value );
		} else if ( tag == "FREE_SUITE_BONUS" ) {
			l >> real_value;
			energy_method_options_->free_dof_options().free_suite_bonus( real_value );
		} else if ( tag == "FREE_SUGAR_BONUS" ) {
			l >> real_value;
			energy_method_options_->free_dof_options().free_sugar_bonus( real_value );
		} else if ( tag == "FREE_2HOPRIME_BONUS" ) {
			l >> real_value;
			energy_method_options_->free_dof_options().free_2HOprime_bonus( real_value );
		} else if ( tag == "PACK_PHOSPHATE_PENALTY" ) {
			l >> real_value;
			energy_method_options_->free_dof_options().pack_phosphate_penalty( real_value );
		} else if ( tag == "FREE_SIDE_CHAIN_BONUS" ) {
			l >> real_value;
			energy_method_options_->free_dof_options().free_side_chain_bonus( real_value );
		} else if ( tag == "SCALE_SIDECHAIN_DENSITY_WEIGHTS" ) {
			utility::vector1< core::Real > scale_sc_density;
			core::Real sc_i;
			while ( !l.fail() ) {
				l >> sc_i;
				if ( l.fail() ) break;
				scale_sc_density.push_back( sc_i );
			}
			if ( scale_sc_density.size() == 1 ) scale_sc_density.resize( core::chemical::num_canonical_aas, scale_sc_density[1] );
			if ( scale_sc_density.size() != core::chemical::num_canonical_aas ) {
				utility_exit_with_message( "incorrect number of arguments to SCALE_SIDECHAIN_DENSITY_WEIGHTS: " + line );
			}
			for ( int i=1; i<=(int) core::chemical::num_canonical_aas; ++i ) {
				energy_method_options_->set_density_sc_scale_byres( (core::chemical::AA)i, scale_sc_density[i] );
			}
		} else {

			// //////////// Regular Weights ///////////////////////
			l.str( line );

			if ( ! patch ) {

				// Weights file parsing
				l >> score_type >> wt;
				if ( l.fail() ) {
					utility_exit_with_message( "bad line in file "+filename+":"+line );
				}
				score_types_from_file.push_back( std::make_pair(  score_type, wt ) ); // will apply at end.

			} else {

				// Patch file parsing
				l >> score_type >> operation >> wt;
				if ( l.fail() ) {
					tr.Error << "could not parse line in patch-file: " << line << std::endl;
					continue;
				}
				// will apply at end.
				patches_from_file.push_back( std::make_pair( std::make_pair( score_type, wt ), operation ) );

			} // if ( ! patch )
		} // if ... else if ... else if ...
	} // while ( getline( data, line ) )

	for ( Size n = 1; n <= score_types_from_file.size(); n++ ) {
		ScoreType const & score_type = score_types_from_file[n].first;
		Real const        & weight     = score_types_from_file[n].second;
		set_weight( score_type, weight );
	}
	for ( Size n = 1; n <= patches_from_file.size(); n++ ) {
		ScoreType const & score_type = patches_from_file[n].first.first;
		Real const & weight     = patches_from_file[n].first.second;
		std::string const & operation = patches_from_file[n].second;
		if ( operation == "*=" ) {
			set_weight( score_type, weights_[ score_type ]*weight );
		} else if ( operation == "=" ) {
			set_weight( score_type, weight );
		} else {
			utility_exit_with_message(
				"unrecognized scorefunction patch operation "+operation+" in file: "+filename
			);
		}
	}


}


///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::set_energy_method_options(
	methods::EnergyMethodOptions const &
	energy_method_options_in )
{
	energy_method_options_
		= methods::EnergyMethodOptionsOP( new methods::EnergyMethodOptions( energy_method_options_in ) );
	//*energy_method_options_ = energy_method_options_in;

	// Some of the energy methods only know about these options when
	// they are constructed. So the safest thing is to destroy them and
	// create them again.
	reset_energy_methods();

	score_function_info_current_ = false;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::reset_energy_methods()
{
	//Note from rhiju: This is super hacky, but it was Andrew Leaver-Fay's idea.
	// Zero out the weights ... which
	// destroys the methods. Then return the weights to their
	// non-zero values to construct them again from scratch.
	EnergyMap weights_old( weights_ );
	for ( Size i = 1; i <= n_score_types; ++i ) {
		set_weight( ScoreType( i ), 0.0 );
	}
	for ( Size i = 1; i <= n_score_types; ++i ) {
		set_weight( ScoreType( i ), weights_old.get( ScoreType( i ) ) );
	}

}

void
ScoreFunction::apply_patch_from_file( std::string const & patch_tag )
{
	_add_weights_from_file( find_weights_file(patch_tag, ".wts_patch"), /*patch=*/ true );
}


/// @brief Given a filename (represented by a std::string), set the e_table for this ScoreFunction.
void
ScoreFunction::set_etable(
	std::string const & etable_name
)
{
	energy_method_options_->etable_type( etable_name );
	reset_energy_methods();
	score_function_info_current_ = false;
}


void
ScoreFunction::set_method_weights(
	ScoreType const & t,
	utility::vector1< Real > const & wts
)
{
	//debug_assert( weights_[ t ] == 0.0 );
	energy_method_options_->set_method_weights( t, wts );
	reset_energy_methods();
	score_function_info_current_ = false;
}

////////////////////////////////////////////////
// Don't use this function unless you're mucking with ScoreFunction implementations
// It may discard subtype information.
void
ScoreFunction::assign( ScoreFunction const & src )
{
	if ( this == &src ) return;

	// copy the weights
	weights_ = src.weights_;

	// deep copy of the energy method options
	energy_method_options_ = methods::EnergyMethodOptionsOP( new methods::EnergyMethodOptions( * src.energy_method_options_ ) );

	// copy the methods:
	initialize_methods_arrays(); // clears & sizes the arrays

	for ( AllMethods::const_iterator it = src.all_methods_.begin(),
			ite = src.all_methods_.end(); it != ite; ++it ) {
		add_method( (*it)->clone() );
	}
	update_intrares_energy_status();
	debug_assert( check_methods() );

	/// SL: Fri Jun 29 12:51:25 EDT 2007 @744 /Internet Time/
	score_function_info_current_ = src.score_function_info_current_;
	score_function_info_ = ScoreFunctionInfoOP( new ScoreFunctionInfo( *src.score_function_info_ ) );

	any_intrares_energies_ = src.any_intrares_energies_;
}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::show( std::ostream & out ) const
{
	out << "ScoreFunction::show():\nweights:";
	for ( int i=1; i<= n_score_types; ++i ) {
		if ( weights_[ ScoreType(i) ] != 0.0 ) {
			out << " (" << ScoreType(i) << ' '<< weights_[ScoreType(i)] << ')';
		}
	}
	out << '\n';

	out << "energy_method_options: " << *energy_method_options_ << '\n';
}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::merge( const ScoreFunction & scorefxn_to_be_merged )
{
	for ( int i=1; i<= n_score_types; ++i ) {
		core::Real the_weight;
		the_weight = scorefxn_to_be_merged.get_weight( ScoreType(i) );
		if ( the_weight != 0.0 ) {
			set_weight(  ScoreType(i), the_weight );
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::show( std::ostream & out,  pose::Pose & pose ) const
{
	(*this)(pose); //make sure scores are set
	out << "------------------------------------------------------------\n";
	out << " Scores                       Weight   Raw Score Wghtd.Score\n";
	out << "------------------------------------------------------------\n";
	float sum_weighted=0.0;
	for ( int i=1; i<= n_score_types; ++i ) {
		if ( weights_[ ScoreType(i) ] != 0.0 ) {
			out << ' ' << LJ(24,ScoreType(i)) << ' '<< F(9,3,weights_[ ScoreType(i) ]) << "   "
				<<    F(9,3,pose.energies().total_energies()[ ScoreType(i) ] ) << "   "
				<<    F(9,3, weights_[ ScoreType(i) ] * pose.energies().total_energies()[ ScoreType(i) ] )
				<< '\n';
			sum_weighted += weights_[ ScoreType(i) ] * pose.energies().total_energies()[ ScoreType(i) ];
		}
	}
	out << "---------------------------------------------------\n";
	out << " Total weighted score:                    " << F(9,3,sum_weighted) << '\n';
}

///////////////////////// output as show( os, pose ) but without the pose //////////////////////////
void
ScoreFunction::show_pretty( std::ostream & out ) const {
	out << "---------------------------------------------\n";
	out << " Scores                       Weight   \n";
	out << "---------------------------------------------\n";
	//float sum_weighted=0.0;
	for ( int i=1; i<= n_score_types; ++i ) {
		if ( weights_[ ScoreType(i) ] != 0.0 ) {
			out << ' ' << LJ(24,ScoreType(i)) << ' '<< F(9,3,weights_[ ScoreType(i) ]) << "   "
				<< '\n';
		}
	}
	out << "---------------------------------------------------\n";
}

///////////////////////////////////////////////////////////////////////////////
void
show_detail( std::ostream & out, EnergyMap & energies,  EnergyMap weights )
{
	out << "------------------------------------------------------------\n";
	out << " Scores                       Weight   Raw Score Wghtd.Score\n";
	out << "------------------------------------------------------------\n";
	float sum_weighted=0.0;
	for ( int i=1; i<= n_score_types; ++i ) {
		if ( weights[ ScoreType(i) ] != 0.0 ) {
			out << ' ' << LJ(24,ScoreType(i)) << ' '<< F(9,3,weights[ ScoreType(i) ]) << "   "
				<<    F(9,3,energies[ ScoreType(i) ] ) << "   "
				<<    F(9,3, weights[ ScoreType(i) ] * energies[ ScoreType(i) ] )
				<< '\n';
			sum_weighted += weights[ ScoreType(i) ] * energies[ ScoreType(i) ];
		}
	}
	out << "---------------------------------------------------\n";
	out << " Total weighted score:                    " << F(9,3,sum_weighted) << '\n';
}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::show_line_headers( std::ostream & out ) const
{

	for ( int i=1; i<= n_score_types; ++i ) {
		if ( weights_[ ScoreType(i) ] != 0.0 ) {
			out << ' ' << LJ(16,ScoreType(i)) << ' ';
		}
	}

}
///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::show_line( std::ostream & out,  pose::Pose const & pose ) const
{

	float sum_weighted=0.0;
	for ( int i=1; i<= n_score_types; ++i ) {
		if ( weights_[ ScoreType(i) ] != 0.0 ) {
			sum_weighted += weights_[ ScoreType(i) ] * pose.energies().total_energies()[ ScoreType(i) ];
		}
	}
	out << F(9,3,sum_weighted);

	for ( int i=1; i<= n_score_types; ++i ) {
		if ( weights_[ ScoreType(i) ] != 0.0 ) {
			// OL: added the " " to the output. Otherwise column separation can be missing in silent-files
			// OL: if large numbers come into play
			// OL: reduced F(9,3) to F(8,3) to account for the extra space.
			out << " " << F(8,3, weights_[ ScoreType(i) ] * pose.energies().total_energies()[ ScoreType(i) ] );
		}
	}

}


///////////////////////////////////////////////////////////////////////////////

// to start out, just thinking fullatom energies
//
// NOTE: no freakin rotamer trials inside scoring!

Real
ScoreFunction::operator()( pose::Pose & pose ) const
{
#ifdef APL_TEMP_DEBUG
	if ( n_minimization_sfxn_evals != 0 && ! pose.energies().use_nblist() ) {
		std::cout << "n_minimization_sfxn_evals: " << n_minimization_sfxn_evals << std::endl;
		n_minimization_sfxn_evals = 0;
	}
	if ( pose.energies().use_nblist() ) {
		++n_minimization_sfxn_evals;
	}
#endif

	// completely unnecessary temporary hack to force refold if nec. for profiling
	if ( pose.total_residue() > 0 ) {
		pose.residue( pose.total_residue() );
	}

	//fpd fail if this is called on a symmetric pose
	runtime_assert( !core::pose::symmetry::is_symmetric( pose ) );

	PROF_START( basic::SCORE );
	//std::cout << "ScoreFunction::operator()\n";

	// notify the pose that we are starting a score evaluation.
	// also tells the cached energies object about the scoring
	// parameters, for the purposes of invalidating cached data
	// if necessary
	//
	// at this point the dof/xyz-moved information will be converted
	// to a domain map. Energy/neighbor links between pair-moved residues
	// will be deleted, and 1d energies of res-moved residues will be
	// cleared.
	//
	// further structure modification will be prevented until scoring is
	// completed
	//
	PROF_START( basic::SCORE_BEGIN_NOTIFY );
	pose.scoring_begin( *this );

	//std::cout << "ScoreFunction::operator() 1\n";
	if ( pose.energies().total_energy() != 0.0 ) {
		std::cout << "STARTING SCORE NON-ZERO!" << std::endl;
	}
	PROF_STOP( basic::SCORE_BEGIN_NOTIFY );
	// ensure that the total_energies are zeroed out -- this happens in Energies.scoring_begin()
	//std::cout << "ScoreFunction::operator() 2\n";
	PROF_START( basic::SCORE_SETUP );
	// do any setup necessary
	setup_for_scoring( pose );
	//std::cout << "ScoreFunction::operator() 3\n";
	PROF_STOP( basic::SCORE_SETUP );

	// evaluate the residue-residue energies that only exist between
	// neighboring residues
	PROF_START( basic::SCORE_NEIGHBOR_ENERGIES );

	eval_twobody_neighbor_energies( pose );

	PROF_STOP ( basic::SCORE_NEIGHBOR_ENERGIES );

	// evaluate the residue pair energies that exist between possibly-distant residues
	PROF_START( basic::SCORE_LONG_RANGE_ENERGIES );

	eval_long_range_twobody_energies( pose );

	PROF_STOP ( basic::SCORE_LONG_RANGE_ENERGIES );

	PROF_START( basic::SCORE_ONEBODY_ENERGIES );

	// evaluate the onebody energies -- rama, dunbrack, ...
	eval_onebody_energies( pose );

	PROF_STOP( basic::SCORE_ONEBODY_ENERGIES );

	PROF_START( basic::SCORE_FINALIZE );
	// give energyfunctions a chance update/finalize energies
	// etable nblist calculation is performed here
	for ( AllMethods::const_iterator it=all_methods_.begin(),
			it_end = all_methods_.end(); it != it_end; ++it ) {
		(*it)->finalize_total_energy( pose, *this, pose.energies().finalized_energies() );
	}

	pose.energies().total_energies() += pose.energies().finalized_energies();

	if ( pose.energies().use_nblist() ) {
		pose.energies().total_energies() += pose.energies().minimization_graph()->fixed_energies();
	}

	PROF_STOP( basic::SCORE_FINALIZE );

	//std::cout << "Total energies: ";
	//pose.energies().total_energies().show_if_nonzero_weight( std::cout, weights_ );
	//std::cout << std::endl;
	PROF_START( basic::SCORE_DOT );

	// dot the weights with the scores
	pose.energies().total_energy() = pose.energies().total_energies().dot( weights_ );
	pose.energies().total_energies()[ total_score ] = pose.energies().total_energy();

	PROF_STOP( basic::SCORE_DOT );

	PROF_START( basic::SCORE_END_NOTIFY );

	// notify that scoring is over
	pose.scoring_end( *this );

	PROF_STOP( basic::SCORE_END_NOTIFY );

	PROF_STOP ( basic::SCORE );

	return pose.energies().total_energy();
}

Real
ScoreFunction::score( pose::Pose & pose ) const {
	return (*this)(pose);
}


///////////////////////////////////////////////////////////////////////////////
Real
ScoreFunction::score_by_scoretype(
	pose::Pose & pose,
	ScoreType const t,
	bool const weighted
) const {
	(*this)(pose); //make sure scores are set
	Real score = pose.energies().total_energies()[t];
	if ( weighted ) score *= weights_[t];
	return score;
}
///////////////////////////////////////////////////////////////////////////////
core::Real
ScoreFunction::get_sub_score(
	pose::Pose const & pose,
	utility::vector1< bool > const & residue_mask
) const {
	debug_assert(residue_mask.size() == pose.total_residue());

	// retrieve cached energies object
	Energies const & energies( pose.energies() );
	debug_assert(energies.energies_updated());
	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// should be zeroed at the beginning of scoring
	Real total( 0.0 ); //( energies.totals() );
	for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
		if ( !residue_mask[i] ) continue;

		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const & edge( static_cast< EnergyEdge const &> (**iru) );
			Size const j( edge.get_second_node_ind() );
			if ( !residue_mask[j] ) continue;

			total += edge.dot( weights_ );
		} // nbrs of i
	} // i=1,nres

	///////////////////////////////////////////////////////////////
	// context independent onebody energies

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !residue_mask[i] ) continue;
		EnergyMap const & emap( energies.onebody_energies( i ) );

		total += weights_.dot( emap, ci_1b_types() );
		total += weights_.dot( emap, cd_1b_types() );

		// 2body energy methods are allowed to define 1body intxns //////
		if ( any_intrares_energies_ ) {
			// context independent:
			total += weights_.dot( emap, ci_2b_types() );
			total += weights_.dot( emap, ci_lr_2b_types() );
			// contex dependent:
			total += weights_.dot( emap, cd_2b_types() );
			total += weights_.dot( emap, cd_lr_2b_types() );
		}
	}

	// add in hbonding
	if ( pose.energies().data().has( EnergiesCacheableDataType::HBOND_SET ) ) {

		hbonds::HBondSet const & hbond_set
			( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));

		hbonds::HBondSet hbond_set_excl(hbond_set, residue_mask);
		EnergyMap hbond_emap;
		hbonds::get_hbond_energies(hbond_set_excl, hbond_emap);
		total += weights_[hbond_sr_bb] * hbond_emap[hbond_sr_bb];
		total += weights_[hbond_lr_bb] * hbond_emap[hbond_lr_bb];
	}

	//////////////////////////////////////////////////
	///  Context Independent Long Range 2Body methods

	for ( CI_LR_2B_Methods::const_iterator iter = ci_lr_2b_methods_.begin(),
			iter_end = ci_lr_2b_methods_.end(); iter != iter_end; ++iter ) {
		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( ! residue_mask[ii] ) continue;
			if ( ! lrec->any_upper_neighbors_for_residue( ii ) ) continue;
			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				if ( !residue_mask[rni->upper_neighbor_id()] ) continue;

				EnergyMap emap;
				rni->retrieve_energy( emap ); // pbmod
				total += weights_.dot( emap, ci_lr_2b_types() );
			}
		}
	}
	/////////////////////////////////////////////////////
	///  Context Independent Long Range twobody methods

	for ( CD_LR_2B_Methods::const_iterator iter = cd_lr_2b_methods_.begin(),
			iter_end = cd_lr_2b_methods_.end(); iter != iter_end; ++iter ) {
		LREnergyContainerCOP lrec =
			pose.energies().long_range_container( (*iter)->long_range_type() );
		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( !residue_mask[ii] ) continue;
			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				if ( !residue_mask[rni->upper_neighbor_id()] ) continue;

				EnergyMap emap;
				rni->retrieve_energy( emap ); // pbmod
				total += weights_.dot( emap, cd_lr_2b_types() );
			}
		}
	}

	return total;
}

core::Real
ScoreFunction::get_sub_score(
	pose::Pose & pose,
	utility::vector1< bool > const & residue_mask
) const {
	//If the energies are not up-to-date score the pose
	if ( !pose.energies().energies_updated() ) ( *this)(pose);

	return get_sub_score(const_cast<pose::Pose const &>(pose), residue_mask);
}


///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::get_sub_score(
	pose::Pose const & pose,
	utility::vector1< bool > const & residue_mask,
	EnergyMap & emap
) const {
	debug_assert(residue_mask.size() == pose.total_residue());

	// retrieve cached energies object
	Energies const & energies( pose.energies() );
	debug_assert(energies.energies_updated());

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// should be zeroed at the beginning of scoring
	for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
		if ( !residue_mask[i] ) continue;

		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const & edge( static_cast< EnergyEdge const &> (**iru) );
			Size const j( edge.get_second_node_ind() );
			if ( !residue_mask[j] ) continue;

			EnergyMap edge_emap(edge.fill_energy_map());
			edge_emap *= weights_;
			emap += edge_emap;
		} // nbrs of i
	} // i=1,nres


	///////////////////////////////////////////////////////////////
	// context independent onebody energies

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !residue_mask[i] ) continue;
		EnergyMap onebody_emap(energies.onebody_energies(i));
		onebody_emap *= weights_;

		emap.accumulate( onebody_emap, ci_1b_types() );
		emap.accumulate( onebody_emap, cd_1b_types() );

		// 2body energy methods are allowed to define 1body intxns //////
		if ( any_intrares_energies_ ) {
			emap.accumulate( onebody_emap, ci_2b_types() );
			emap.accumulate( onebody_emap, ci_lr_2b_types() );
			emap.accumulate( onebody_emap, cd_2b_types() );
			emap.accumulate( onebody_emap, cd_lr_2b_types() );
		}
	}

	// add in backbone-backbone hbonding
	if ( pose.energies().data().has( EnergiesCacheableDataType::HBOND_SET ) ) {

		hbonds::HBondSet const & hbond_set
			( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));

		hbonds::HBondSet hbond_set_excl(hbond_set, residue_mask);
		EnergyMap hbond_emap;
		hbonds::get_hbond_energies(hbond_set_excl, hbond_emap);
		emap[hbond_sr_bb] += weights_[hbond_sr_bb] * hbond_emap[hbond_sr_bb];
		emap[hbond_lr_bb] += weights_[hbond_lr_bb] * hbond_emap[hbond_lr_bb];
	}

	//////////////////////////////////////////////////
	///  Context Independent Long Range 2Body methods

	for ( CI_LR_2B_Methods::const_iterator iter = ci_lr_2b_methods_.begin(),
			iter_end = ci_lr_2b_methods_.end(); iter != iter_end; ++iter ) {
		LREnergyContainerCOP lrec =
			pose.energies().long_range_container((*iter)->long_range_type());
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( !residue_mask[ii] ) continue;

			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				if ( !residue_mask[rni->upper_neighbor_id()] ) continue;

				EnergyMap lr_emap;
				rni->retrieve_energy( lr_emap ); // pbmod
				emap += lr_emap;
			}
		}

		/////////////////////////////////////////////////////
		///  Context Independent Long Range twobody methods

		for ( CD_LR_2B_Methods::const_iterator iter = cd_lr_2b_methods_.begin(),
				iter_end = cd_lr_2b_methods_.end(); iter != iter_end; ++iter ) {
			LREnergyContainerCOP lrec =
				pose.energies().long_range_container((*iter)->long_range_type());
			// Potentially O(N^2) operation...
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				if ( !residue_mask[ii] ) continue;
				for ( ResidueNeighborConstIteratorOP
						rni = lrec->const_upper_neighbor_iterator_begin( ii ),
						rniend = lrec->const_upper_neighbor_iterator_end( ii );
						(*rni) != (*rniend); ++(*rni) ) {
					if ( !residue_mask[rni->upper_neighbor_id()] ) continue;

					EnergyMap lr_emap;
					rni->retrieve_energy( lr_emap ); // pbmod
					emap += lr_emap;

				}
			}
		}
	}
}

void
ScoreFunction::get_sub_score(
	pose::Pose & pose,
	utility::vector1< bool > const & residue_mask,
	EnergyMap & emap
) const {
	//If the energies are not up-to-date score the pose
	if ( !pose.energies().energies_updated() ) ( *this)(pose);
	return get_sub_score(const_cast<pose::Pose const &>(pose), residue_mask, emap);
}


///////////////////////////////////////////////////////////////////////////////
Real
ScoreFunction::get_sub_score_exclude_res(
	pose::Pose const & pose,
	utility::vector1< Size > const & exclude_res
) const {

	utility::vector1< bool > residue_mask(pose.total_residue(), true);
	for (
			utility::vector1< Size >::const_iterator
			ii = exclude_res.begin(), ie = exclude_res.end();
			ii != ie; ++ii ) {
		residue_mask[*ii] = false;
	}
	return get_sub_score(pose, residue_mask);
}

///////////////////////////////////////////////////////////////////////////////
Real
ScoreFunction::get_sub_score_exclude_res(
	pose::Pose & pose,
	utility::vector1< Size > const & exclude_res
) const {

	//If the energies are not up-to-date score the pose
	if ( !pose.energies().energies_updated() ) ( *this)(pose);

	utility::vector1< bool > residue_mask(pose.total_residue(), true);
	for (
			utility::vector1< Size >::const_iterator
			ii = exclude_res.begin(), ie = exclude_res.end();
			ii != ie; ++ii ) {
		residue_mask[*ii] = false;
	}
	return get_sub_score(const_cast<pose::Pose const &>(pose), residue_mask);
}


///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::get_sub_score_exclude_res(
	pose::Pose const & pose,
	utility::vector1< core::Size > const & exclude_res,
	EnergyMap & emap
) const {
	utility::vector1< bool > residue_mask(pose.total_residue(), true);
	for (
			utility::vector1< Size >::const_iterator
			ii = exclude_res.begin(), ie = exclude_res.end();
			ii != ie; ++ii ) {
		residue_mask[*ii] = false;
	}
	get_sub_score(pose, residue_mask, emap);
}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::get_sub_score_exclude_res(
	pose::Pose & pose,
	utility::vector1< core::Size > const & exclude_res,
	EnergyMap & emap
) const {

	//If the energies are not up-to-date score the pose
	if ( !pose.energies().energies_updated() ) ( *this)(pose);

	utility::vector1< bool > residue_mask(pose.total_residue(), true);
	for (
			utility::vector1< Size >::const_iterator
			ii = exclude_res.begin(), ie = exclude_res.end();
			ii != ie; ++ii ) {
		residue_mask[*ii] = false;
	}
	get_sub_score(const_cast<pose::Pose const &>(pose), residue_mask, emap);
}


///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::eval_twobody_neighbor_energies(
	pose::Pose & pose
) const
{

	/// Use the EnergyGraph or the MinimizationGraph to direct scoring

	// cached energies object
	Energies & energies( pose.energies() );

	EnergyMap & total_energies( const_cast< EnergyMap & > ( energies.total_energies() ) );

	// the neighbor/energy links
	EnergyGraph & energy_graph( energies.energy_graph() );

	// are we using the atom-atom nblist?
	// renaming for true purpose: are we minimizing -- if so,
	// zero the energies stored on edges in the Energy graph, but do
	// not mark the edges as having had their energies computed
	bool const minimizing( energies.use_nblist() );

	if ( minimizing ) {
		/// When minimizing, do not touch the EnergyGraph -- leave it fixed
		MinimizationGraphCOP g = energies.minimization_graph();
		for ( Size ii = 1; ii < pose.total_residue(); ++ii ) {
			conformation::Residue const & ii_rsd( pose.residue( ii ) );
			for ( core::graph::Graph::EdgeListConstIter
					edge_iter = g->get_node( ii )->const_upper_edge_list_begin(),
					edge_iter_end = g->get_node( ii )->const_upper_edge_list_end();
					edge_iter != edge_iter_end; ++edge_iter ) {
				Size const jj = (*edge_iter)->get_second_node_ind();
				conformation::Residue const & jj_rsd( pose.residue( jj ));
				MinimizationEdge const & minedge( static_cast< MinimizationEdge const & > (**edge_iter) );

				eval_res_pair_energy_for_minedge( minedge, ii_rsd, jj_rsd, pose, *this, total_energies );
			}
		}

	} else {
		EnergyMap tbemap;

		for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
			conformation::Residue const & resl( pose.residue( i ) );
			for ( graph::Graph::EdgeListIter
					iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->upper_edge_list_end();
					iru != irue; ++iru ) {
				EnergyEdge & edge( static_cast< EnergyEdge & > (**iru) );

				Size const j( edge.get_second_node_ind() );
				conformation::Residue const & resu( pose.residue( j ) );
				// the pair energies cached in the link
				///EnergyMap & emap( energy_edge->energy_map() );
				tbemap.zero( cd_2b_types() );
				tbemap.zero( ci_2b_types() );

				// the context-dependent guys can't be cached, so they are always reevaluated
				eval_cd_2b( resl, resu, pose, tbemap );

				if ( edge.energies_not_yet_computed() ) {
					// energies not yet computed? <-> (during minimization w/ nblist) moving rsd pair
					/// TEMP zero portions;

					if ( minimizing ) { // APL -- 5/18/2010 -- DEPRICATED
						// confirm that this rsd-rsd interaction will be included
						// in the atom pairs on the nblist:
						//debug_assert( ( pose.energies().domain_map(i) !=
						//     pose.energies().domain_map(j) ) ||
						//    ( pose.energies().res_moved(i) ) );
						// ensure that these are zeroed, since we will hit them at the
						// end inside the nblist calculation
						// they almost certainly should be, since the energies have not
						// yet been computed...
						eval_ci_2b( resl, resu, pose, tbemap );
						edge.store_active_energies( tbemap );
						// do not mark energies as computed!!!!!!!!!!!!!
					} else {
						eval_ci_2b( resl, resu, pose, tbemap );
						edge.store_active_energies( tbemap );
						edge.mark_energies_computed();
					}
				} else {
					/// Read the CI energies from the edge, as they are still valid;
					for ( Size ii = 1; ii <= ci_2b_types().size(); ++ii ) {
						tbemap[ ci_2b_types()[ ii ]] = edge[ ci_2b_types()[ ii ] ];
					}

					/// Save the freshly computed CD energies on the edge
					edge.store_active_energies( tbemap, cd_2b_types() );
				}

				total_energies.accumulate( tbemap, ci_2b_types() );
				total_energies.accumulate( tbemap, cd_2b_types() );
			} // nbrs of i
		} // i=1,nres
	} // not minimizing
}

void
ScoreFunction::eval_long_range_twobody_energies( pose::Pose & pose ) const
{
	EnergyMap & total_energies(const_cast< EnergyMap & > (pose.energies().total_energies()));

	bool const minimizing( pose.energies().use_nblist() );
	if ( minimizing ) return; // long range energies are handled as part of the 2-body energies in the minimization graph

	for ( CI_LR_2B_Methods::const_iterator iter = ci_lr_2b_methods_.begin(),
			iter_end = ci_lr_2b_methods_.end(); iter != iter_end; ++iter ) {

		LREnergyContainerOP lrec = pose.energies().nonconst_long_range_container( (*iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( ResidueNeighborIteratorOP
					rni = lrec->upper_neighbor_iterator_begin( ii ),
					rniend = lrec->upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				EnergyMap emap;
				if ( ! rni->energy_computed() ) {
					(*iter)->residue_pair_energy(
						pose.residue(ii),
						pose.residue( rni->upper_neighbor_id() ),
						pose, *this, emap );
					rni->save_energy( emap );

					/// DANGER DANGER DANGER.  use_nblist() now means "In the process of a minimization". There is
					/// no such thing as a non-neighborlist minimization.  This will confuse people at some point.
					/// We ought to have some "are we minimizing currently" flag who's meaning can be decoupled
					/// from the neighborlist idea.
					if ( ! pose.energies().use_nblist() ) {
						rni->mark_energy_computed();
					}
				} else {
					rni->retrieve_energy( emap ); // pbmod
				}
				total_energies += emap;
			}
		}

	}

	for ( CD_LR_2B_Methods::const_iterator iter = cd_lr_2b_methods_.begin(),
			iter_end = cd_lr_2b_methods_.end(); iter != iter_end; ++iter ) {

		LREnergyContainerOP lrec
			= pose.energies().nonconst_long_range_container( (*iter)->long_range_type() );

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( ResidueNeighborIteratorOP
					rni = lrec->upper_neighbor_iterator_begin( ii ),
					rniend = lrec->upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {

				EnergyMap emap;
				(*iter)->residue_pair_energy(
					pose.residue(ii),
					pose.residue( rni->upper_neighbor_id() ),
					pose, *this, emap );

				rni->save_energy( emap );
				//rni->mark_energy_computed();
				total_energies += emap;

			}
		}

	}
}


///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::eval_onebody_energies( pose::Pose & pose ) const
{
	// context independent onebody energies
	Energies & energies( pose.energies() );
	EnergyMap & totals( energies.total_energies() );

	bool const minimizing( energies.use_nblist() );

	if ( minimizing ) {
		MinimizationGraphCOP mingraph = energies.minimization_graph();
		debug_assert( mingraph );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			conformation::Residue const & rsd = pose.residue( ii );

			MinimizationNode const & iiminnode =  * mingraph->get_minimization_node( ii );
			eval_res_onebody_energies_for_minnode( iiminnode, rsd, pose, *this, totals );

		}

	} else {
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			EnergyMap & emap( energies.onebody_energies( i ) );

			// 1body intxns ///////////////////////////
			if ( energies.res_moved( i ) ) {
				// have to recalculate
				emap.clear(); // should already have been done when the domain map was internalized
				eval_ci_1b( pose.residue(i), pose, emap );
			}

			emap.zero( cd_1b_types() ); // cant be cached
			eval_cd_1b( pose.residue(i), pose, emap );


			// 2body energy methods are allowed to define 1body intxns ///////////////////////////
			if ( any_intrares_energies_ ) {
				// context independent:
				if ( energies.res_moved( i ) ) {
					eval_ci_intrares_energy( pose.residue(i), pose, emap );
				}
				// context dependent
				emap.zero( cd_2b_types() ); // cant be cached (here only intrares are relevant)
				emap.zero( cd_lr_2b_types() ); // cant be cached (here only intrares are relevant)
				eval_cd_intrares_energy( pose.residue(i), pose, emap );
			}

			totals += emap;

			energies.reset_res_moved( i ); // mark one body energies as having been calculated
			//std::cout << "totals: "<<  i  << totals;
		}
	}
}

void
ScoreFunction::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	if ( ! any_intrares_energies_ ) return;
	for ( TWO_B_Methods::const_iterator iter = ci_2b_intrares_.begin(),
			iter_end = ci_2b_intrares_.end(); iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy( rsd, pose, *this, emap );
	}
	for ( TWO_B_Methods::const_iterator iter = cd_2b_intrares_.begin(),
			iter_end = cd_2b_intrares_.end(); iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy( rsd, pose, *this, emap );
	}
}

void
ScoreFunction::evaluate_rotamer_intrares_energies(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	utility::vector1< core::PackerEnergy > & energies
) const
{
	if ( ! any_intrares_energies_ ) return;
	for ( TWO_B_Methods::const_iterator iter = ci_2b_intrares_.begin(),
			iter_end = ci_2b_intrares_.end(); iter != iter_end; ++iter ) {
		(*iter)->evaluate_rotamer_intrares_energies( set, pose, *this, energies );
	}
	for ( TWO_B_Methods::const_iterator iter = cd_2b_intrares_.begin(),
			iter_end = cd_2b_intrares_.end(); iter != iter_end; ++iter ) {
		(*iter)->evaluate_rotamer_intrares_energies( set, pose, *this, energies );
	}
}

void
ScoreFunction::evaluate_rotamer_intrares_energy_maps(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	utility::vector1< EnergyMap > & emaps
) const
{
	if ( ! any_intrares_energies_ ) return;
	for ( TWO_B_Methods::const_iterator iter = ci_2b_intrares_.begin(),
			iter_end = ci_2b_intrares_.end(); iter != iter_end; ++iter ) {
		(*iter)->evaluate_rotamer_intrares_energy_maps( set, pose, *this, emaps );
	}
	for ( TWO_B_Methods::const_iterator iter = cd_2b_intrares_.begin(),
			iter_end = cd_2b_intrares_.end(); iter != iter_end; ++iter ) {
		(*iter)->evaluate_rotamer_intrares_energy_maps( set, pose, *this, emaps );
	}
}

void
ScoreFunction::eval_ci_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	if ( ! any_intrares_energies_ ) return;
	for ( TWO_B_Methods::const_iterator iter = ci_2b_intrares_.begin(),
			iter_end = ci_2b_intrares_.end(); iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy( rsd, pose, *this, emap );
	}
}

void
ScoreFunction::eval_cd_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	if ( ! any_intrares_energies_ ) return;
	for ( TWO_B_Methods::const_iterator iter = cd_2b_intrares_.begin(),
			iter_end = cd_2b_intrares_.end(); iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy( rsd, pose, *this, emap );
	}
}


///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::eval_ci_1b(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CI_1B_Methods::const_iterator iter = ci_1b_methods_.begin(),
			iter_end = ci_1b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->residue_energy( rsd, pose, emap );
	}
}

void
ScoreFunction::eval_cd_1b(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CD_1B_Methods::const_iterator iter = cd_1b_methods_.begin(),
			iter_end = cd_1b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->residue_energy( rsd, pose, emap );
	}
}


///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::eval_ci_2b(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy( rsd1, rsd2, pose, *this, emap );
	}

}


void
ScoreFunction::eval_ci_2b_bb_bb(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->backbone_backbone_energy( rsd1, rsd2, pose, *this, emap );
	}

}


void
ScoreFunction::eval_ci_2b_bb_sc(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->backbone_sidechain_energy( rsd1, rsd2, pose, *this, emap );
	}

}

void
ScoreFunction::eval_ci_2b_sc_sc(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->sidechain_sidechain_energy( rsd1, rsd2, pose, *this, emap );
	}

}


///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::eval_cd_2b(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CD_2B_Methods::const_iterator
			iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy( rsd1, rsd2, pose, *this, emap );
	}
}

void
ScoreFunction::eval_cd_2b_bb_bb(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->backbone_backbone_energy( rsd1, rsd2, pose, *this, emap );
	}

}


void
ScoreFunction::eval_cd_2b_bb_sc(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->backbone_sidechain_energy( rsd1, rsd2, pose, *this, emap );
	}

}

void
ScoreFunction::eval_cd_2b_sc_sc(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->sidechain_sidechain_energy( rsd1, rsd2, pose, *this, emap );
	}

}


/// @brief score the sidechain from rsd1 against the entirety of rsd2
///
/// low fidelity description of what it would be like to replace
/// the sidechain of rsd1 at a particular position -- only those scoring functions
/// that specify they should be included in the bump check will be evaluated.
/// The inaccuracy in this function is based on the change in LK solvation types
/// and/or charges for backbone atoms when residues are mutated from one amino
/// acid type to another (or from any other residue with backbone atoms)
void
ScoreFunction::bump_check_full(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->bump_energy_full( rsd1, rsd2, pose, *this, emap );
	}
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->bump_energy_full( rsd1, rsd2, pose, *this, emap );
	}

}


/// @details low fidelity description of what it would be like to replace
/// the sidechain of rsd1 at a particular position knowing that
/// rsd2's sidechain might also change sychronously.
/// only those scoring functions that specify they should be included in the
/// bump check will be evaluated.
/// The inaccuracy in this function is based on the change in LK solvation types
/// and/or charges for backbone atoms when residues are mutated from one amino
/// acid type to another (or from any other residue with backbone atoms)
void
ScoreFunction::bump_check_backbone(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->bump_energy_backbone( rsd1, rsd2, pose, *this, emap );
	}
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->bump_energy_backbone( rsd1, rsd2, pose, *this, emap );
	}

}

void
ScoreFunction::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	FArray2D< core::PackerEnergy > & energy_table
) const
{
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->evaluate_rotamer_pair_energies(
			set1, set2, pose, *this, weights(), energy_table );
	}
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->evaluate_rotamer_pair_energies(
			set1, set2, pose, *this, weights(), energy_table );
	}

}

/// @details Determines whether two positions could have non-zero interaction energy
bool
ScoreFunction::are_they_neighbors(
	pose::Pose const & pose,
	Size const pos1,
	Size const pos2
) const
{
	conformation::Residue const & rsd1( pose.residue( pos1 ) );
	conformation::Residue const & rsd2( pose.residue( pos2 ) );

	Real const intxn_radius( rsd1.nbr_radius() + rsd2.nbr_radius() + max_atomic_interaction_cutoff() );

	return ( ( intxn_radius > 0 &&
		rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) < intxn_radius * intxn_radius ) ||
		( any_lr_residue_pair_energy( pose, pos1, pos2 ) ) );
}


bool
ScoreFunction::any_lr_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const
{
	for ( LR_2B_Methods::const_iterator iter = lr_2b_methods_.begin(),
			iter_end = lr_2b_methods_.end(); iter != iter_end; ++iter ) {
		if ( (*iter)->defines_residue_pair_energy( pose, res1, res2 ) ) {
			return true;
		}
	}

	return false;
}

ScoreFunction::AllMethodsIterator
ScoreFunction::all_energies_begin() const
{
	return all_methods_.begin();
}

ScoreFunction::AllMethodsIterator
ScoreFunction::all_energies_end() const
{
	return all_methods_.end();
}

ScoreFunction::LR_2B_MethodIterator
ScoreFunction::long_range_energies_begin() const
{
	return lr_2b_methods_.begin();
}

ScoreFunction::LR_2B_MethodIterator
ScoreFunction::long_range_energies_end() const
{
	return lr_2b_methods_.end();
}

ScoreFunction::TWO_B_MethodIterator
ScoreFunction::ci_2b_intrares_begin() const
{
	return ci_2b_intrares_.begin();
}

ScoreFunction::TWO_B_MethodIterator
ScoreFunction::ci_2b_intrares_end() const
{
	return ci_2b_intrares_.end();
}

ScoreFunction::TWO_B_MethodIterator
ScoreFunction::cd_2b_intrares_begin() const
{
	return cd_2b_intrares_.begin();
}

ScoreFunction::TWO_B_MethodIterator
ScoreFunction::cd_2b_intrares_end() const
{
	return cd_2b_intrares_.end();
}

ScoreFunction::CI_2B_Methods::const_iterator
ScoreFunction::ci_2b_begin() const
{
	return ci_2b_methods_.begin();
}

ScoreFunction::CI_2B_Methods::const_iterator
ScoreFunction::ci_2b_end() const
{
	return ci_2b_methods_.end();
}

ScoreFunction::CD_2B_Methods::const_iterator
ScoreFunction::cd_2b_begin() const
{
	return cd_2b_methods_.begin();
}

ScoreFunction::CD_2B_Methods::const_iterator
ScoreFunction::cd_2b_end() const
{
	return cd_2b_methods_.end();
}


ScoreFunction::CI_LR_2B_MethodIterator
ScoreFunction::ci_lr_2b_methods_begin() const
{
	return ci_lr_2b_methods_.begin();
}

ScoreFunction::CI_LR_2B_MethodIterator
ScoreFunction::ci_lr_2b_methods_end() const
{
	return ci_lr_2b_methods_.end();
}

ScoreFunction::CD_LR_2B_MethodIterator
ScoreFunction::cd_lr_2b_methods_begin() const
{
	return cd_lr_2b_methods_.begin();
}

ScoreFunction::CD_LR_2B_MethodIterator
ScoreFunction::cd_lr_2b_methods_end() const
{
	return cd_lr_2b_methods_.end();
}

ScoreFunction::WS_MethodIterator
ScoreFunction::ws_methods_begin() const
{
	return ws_methods_.begin();
}

ScoreFunction::WS_MethodIterator
ScoreFunction::ws_methods_end() const
{
	return ws_methods_.end();
}

/// @details Evaluate the energies between a set of rotamers and a single other (background) residue
void
ScoreFunction::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set1,
	conformation::Residue const & residue2,
	pose::Pose const & pose,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->evaluate_rotamer_background_energies(
			set1, residue2, pose, *this, weights(), energy_vector );
	}
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->evaluate_rotamer_background_energies(
			set1, residue2, pose, *this, weights(), energy_vector );
	}
}


///////////////////////////////////////////////////////////////////////////////
/// @details set the weight
void
ScoreFunction::set_weight( ScoreType const & t, Real const & setting )
{
	using namespace methods;

	Real const old_weight( weights_[t] );
	weights_[t] = setting;

	if ( setting == Real(0.0) ) {
		if ( old_weight != Real(0.0) ) {

			EnergyMethodCOP method( methods_by_score_type_[t] );
			debug_assert( method );

			// check to see if there's a non-zero weight for this method
			bool has_nonzero_weight( false );
			for ( ScoreTypes::const_iterator
					iter     = method->score_types().begin(),
					iter_end = method->score_types().end(); iter != iter_end; ++iter ) {
				if ( weights_[ *iter ] != 0.0 ) {
					has_nonzero_weight = true;
					break;
				}
			}
			if ( ! has_nonzero_weight ) {
				score_function_info_current_ = false; // cached score_function_info_ object no longer up-to-date
				remove_method( utility::pointer::const_pointer_cast< EnergyMethod > (method) ); // cast to non-const, though, all this does is delete the object
			}
		} else {
			// Nothing needs doing.  Early return.  Be careful of intrares status and check_methods() for early returns!!!
			return;
		}
	} else {
		if ( old_weight == Real(0.0) ) {
			// do we already have a method that evaluates this score??
			if ( ! methods_by_score_type_[ t ] != 0  &&
					( t != python ) )  { // sheffler

				score_function_info_current_ = false; // cached score_function_info_ object no longer up-to-date
				debug_assert( std::abs( old_weight ) < 1e-2 ); // sanity

				// get pointer to the energyfunction that evaluates this score
				EnergyMethodOP method ( ScoringManager::get_instance()->energy_method( t, *energy_method_options_));
				add_method( method );
			} else {
				if ( t == python ) return; // escape for python defined scores -- intentially avoid update_intrares_status
			}
		} else {
			// only adjusted a weight that was already non-zero.  Nothing else needs doing.
			// Early return.  Be careful of intrares status and check_methods() for early returns!!!
			debug_assert( old_weight != Real( 0.0 ) || weights_[ t ] != Real( 0.0 ) );
			return;
		}
	}

	// arriving here, we know that a weight went either to or from zero.
	// perform some neccessary bookkeeping.
	debug_assert( old_weight == Real( 0.0 ) || weights_[ t ] == Real( 0.0 ) );

	update_intrares_energy_status();
	debug_assert( check_methods() );
}

///////////////////////////////////////////////////////////////////////////////
/// @details set the weight if zero
void
ScoreFunction::set_weight_if_zero(const ScoreType& t, const Real& setting) {
	if ( get_weight(t) == 0 ) {
		set_weight(t, setting);
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details get the weight
Real
ScoreFunction::get_weight( ScoreType const & t) const
{
	return  weights_[t];
}


///////////////////////////////////////////////////////////////////////////////
/// @details private
bool
ScoreFunction::check_methods() const
{
	using namespace methods;

	EnergyMap counter;
	// counter indexes scoretypes to the number of energy methods. Each scoretype should only
	// map to one energy method.
	for ( AllMethods::const_iterator it=all_methods_.begin(),
			it_end = all_methods_.end(); it != it_end; ++it ) {
		for ( ScoreTypes::const_iterator t=(*it)->score_types().begin(),
				t_end=(*it)->score_types().end(); t != t_end; ++t ) {
			//++counter[ *t ]; //++ not defined for floats, right?
			counter[ *t ] += 1;
		}
	}
	for ( EnergyMap::const_iterator it=counter.begin(), it_end=counter.end();
			it != it_end; ++it ) {
		// total hack!!!
		ScoreType const t( static_cast< ScoreType >( it - counter.begin() + 1));
		Real const count( *it );
		if ( count > 1.5 ) {
			utility_exit_with_message("multiple methods for same score type!" );
		} else if ( count < 0.5 && std::abs( weights_[ t ] ) > 1e-2 ) {
			utility_exit_with_message(
				"no method for scoretype " + name_from_score_type( t ) + " w/nonzero weight!"
			);
		}
	}
	return true;
}

void
ScoreFunction::update_intrares_energy_status()
{
	any_intrares_energies_ = false;
	ci_2b_intrares_.clear();
	cd_2b_intrares_.clear();
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		if ( (*iter)->defines_intrares_energy(weights()) ) {
			any_intrares_energies_ = true;
			cd_2b_intrares_.push_back( (*iter) );
		}
	}
	for ( CD_LR_2B_Methods::const_iterator iter = cd_lr_2b_methods_.begin(),
			iter_end = cd_lr_2b_methods_.end(); iter != iter_end; ++iter ) {
		if ( (*iter)->defines_intrares_energy(weights()) ) {
			any_intrares_energies_ = true;
			cd_2b_intrares_.push_back( (*iter) );
		}
	}

	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		if ( (*iter)->defines_intrares_energy( weights() ) ) {
			any_intrares_energies_ = true;
			ci_2b_intrares_.push_back( (*iter) );
		}
	}
	for ( CI_LR_2B_Methods::const_iterator iter = ci_lr_2b_methods_.begin(),
			iter_end = ci_lr_2b_methods_.end(); iter != iter_end; ++iter ) {
		if ( (*iter)->defines_intrares_energy( weights() ) ) {
			any_intrares_energies_ = true;
			ci_2b_intrares_.push_back( (*iter) );
		}
	}
}

void
ScoreFunction::setup_for_scoring(
	pose::Pose & pose
) const
{

	if ( ! pose.energies().use_nblist() ) {
		for ( AllMethods::const_iterator it=all_methods_.begin(),
				it_end = all_methods_.end(); it != it_end; ++it ) {
			(*it)->setup_for_scoring( pose, *this );
		}
	} else {
		debug_assert( pose.energies().minimization_graph() );
		MinimizationGraphOP mingraph = pose.energies().minimization_graph();

		/// 1. setup_for_scoring for residue singles
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			mingraph->get_minimization_node( ii )->setup_for_scoring( pose.residue( ii ), pose, *this );
		}

		/// 2. inter-residue setup for derivatives
		for ( graph::Graph::EdgeListIter
				edgeit = mingraph->edge_list_begin(), edgeit_end = mingraph->edge_list_end();
				edgeit != edgeit_end; ++edgeit ) {
			MinimizationEdge & minedge = static_cast< MinimizationEdge & > ( (**edgeit) );
			minedge.setup_for_scoring(
				pose.residue(minedge.get_first_node_ind()),
				pose.residue(minedge.get_second_node_ind()),
				pose, *this );
		}

		/// 3. Whole-pose-context energies should be allowed to setup for scoring now.
		for ( MinimizationGraph::Energies::const_iterator
				iter     = mingraph->whole_pose_context_enmeths_begin(),
				iter_end = mingraph->whole_pose_context_enmeths_end();
				iter != iter_end; ++iter ) {
			(*iter)->setup_for_scoring( pose, *this );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const & residues_repacking,
	utility::vector1< bool > const & residues_designing
) const
{
	for ( AllMethods::const_iterator iter=all_methods_.begin(),
			iter_end= all_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_packing( pose, residues_repacking, residues_designing );
	}
}


void
ScoreFunction::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & set
) const
{
	for ( AllMethods::const_iterator iter=all_methods_.begin(),
			iter_end= all_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->prepare_rotamers_for_packing( pose, set );
	}
}


void
ScoreFunction::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	for ( AllMethods::const_iterator iter=all_methods_.begin(),
			iter_end= all_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->update_residue_for_packing( pose, resid );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::setup_for_derivatives(
	pose::Pose & pose
) const
{
	//std::cout << "ScoreFunction::setup_for_derivatives" << std::endl;
	//for ( AllMethods::const_iterator iter=all_methods_.begin(),
	//  iter_end= all_methods_.end(); iter != iter_end; ++iter ) {
	// (*iter)->setup_for_derivatives( pose, *this );
	//}

	debug_assert( pose.energies().minimization_graph() );
	MinimizationGraphOP mingraph = pose.energies().minimization_graph();

	/// 1. setup_for_derivatives for residue singles
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		MinimizationNode & iiminnode( * mingraph->get_minimization_node( ii ));
		iiminnode.setup_for_derivatives( pose.residue(ii), pose, *this );
	}

	/// 2. inter-residue setup for derivatives
	for ( graph::Graph::EdgeListIter
			edgeit = mingraph->edge_list_begin(), edgeit_end = mingraph->edge_list_end();
			edgeit != edgeit_end; ++edgeit ) {
		MinimizationEdge & minedge = static_cast< MinimizationEdge & > ( (**edgeit) );
		minedge.setup_for_derivatives(
			pose.residue(minedge.get_first_node_ind()),
			pose.residue(minedge.get_second_node_ind()),
			pose, *this );
	}

	/// 3. Whole-pose-context energies should be allowed to setup for derivatives now.
	for ( MinimizationGraph::Energies::const_iterator
			iter     = mingraph->whole_pose_context_enmeths_begin(),
			iter_end = mingraph->whole_pose_context_enmeths_end();
			iter != iter_end; ++iter ) {
		(*iter)->setup_for_derivatives( pose, *this );
	}

}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::finalize_after_derivatives(
	pose::Pose & pose
) const
{
	//std::cout << "ScoreFunction::setup_for_derivatives" << std::endl;
	for ( AllMethods::const_iterator iter=all_methods_.begin(),
			iter_end= all_methods_.end(); iter != iter_end; ++iter ) {
		(*iter)->finalize_after_derivatives( pose, *this );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::setup_for_minimizing(
	pose::Pose & pose,
	kinematics::MinimizerMapBase const & min_map
) const
{

	/// 1. Initialize the nodes in the energy graph with one-body
	/// and the intra-residue two-body energy methods

	/// 2. Initialize the short-ranged edges in the energy graph with
	/// the short-ranged two-body energy methods.

	/// 3. Initialize any additional edges with the long-range
	/// two-body energies

	/// 4. Drop any edges from the minimization graph that produce no
	/// energies, and call setup-for-minimization on all edges that remain

	/// 5. Setup whole-structure energies and energy methods that opt-out
	/// of the minimization-graph control over derivative evaluation.


	MinimizationGraphOP g( new MinimizationGraph( pose.total_residue() ) );
	std::list< methods::EnergyMethodCOP > eval_derivs_with_pose_enmeths;
	for ( AllMethods::const_iterator iter=all_methods_.begin(),
			iter_end= all_methods_.end(); iter != iter_end; ++iter ) {
		if ( (*iter)->defines_high_order_terms( pose ) || (*iter)->minimize_in_whole_structure_context( pose ) ) {
			eval_derivs_with_pose_enmeths.push_back( *iter );
		}
	}

	EnergyMap fixed_energies; // portions of the score function that will not change over the course of minimization.

	kinematics::DomainMap const & domain_map = min_map.domain_map();
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		setup_for_minimizing_for_node( * g->get_minimization_node( ii ), pose.residue( ii ),
			min_map, pose, true, fixed_energies );
	}

	//EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	g->copy_connectivity( pose.energies().energy_graph() );
	// 2. Now initialize the edges of the minimization graph using the edges of the EnergyGraph;
	// The energy graph should be up-to-date before this occurs
	for ( core::graph::Graph::EdgeListIter
			edge_iter = g->edge_list_begin(),
			edge_iter_end = g->edge_list_end(),
			ee_edge_iter = pose.energies().energy_graph().edge_list_begin();
			edge_iter != edge_iter_end; ++edge_iter, ++ee_edge_iter ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();
		debug_assert( node1 == (*ee_edge_iter)->get_first_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		debug_assert( node2 == (*ee_edge_iter)->get_second_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering

		MinimizationEdge & minedge( static_cast< MinimizationEdge & > (**edge_iter) );
		// domain map check here?
		bool const res_moving_wrt_eachother(
			domain_map( node1 ) == 0 ||
			domain_map( node2 ) == 0 ||
			domain_map( node1 ) != domain_map( node2 ) );

		setup_for_minimizing_sr2b_enmeths_for_minedge(
			pose.residue( node1 ), pose.residue( node2 ),
			minedge, min_map, pose, res_moving_wrt_eachother, true,
			static_cast< EnergyEdge const * > (*ee_edge_iter), fixed_energies );
	}

	/// 3. Long range energies need time to get included into the graph, which may require the addition of new edges
	/// 3a: CILR2B energies
	///    i.   Iterate across all ci2b long range energy methods
	///    ii.  Iterate across all residue pairs indicated in the long-range energy containers
	///    iii. If two residues have the same non-zero domain-map coloring, then accumulate their interaction energy and move on
	///    iv.  otherwise, find the corresponding minimization-graph edge for this residue pair
	///    v.   adding a new edge if necessary,
	///    vi.  and prepare the minimization data for this edge

	for ( LR_2B_MethodIterator
			iter = long_range_energies_begin(),
			iter_end = long_range_energies_end();
			iter != iter_end; ++iter ) {

		// NO! go ahead an include these terms in the minimization graph for the sake of scoring
		/// if ( (*iter)->minimize_in_whole_structure_context( pose ) ) continue;

		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue;

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				Size const jj = rni->upper_neighbor_id();
				bool const res_moving_wrt_eachother(
					domain_map( ii ) == 0 ||
					domain_map( jj ) == 0 ||
					domain_map( ii ) != domain_map( jj ) );
				setup_for_lr2benmeth_minimization_for_respair(
					pose.residue( ii ), pose.residue( jj ), *iter, *g, min_map, pose,
					res_moving_wrt_eachother, true, rni, fixed_energies );
			}
		}
	}

	/// 4. Call setup_for_minimizing on each edge that has active twobody energies, and drop
	/// all other edges.
	for ( core::graph::Graph::EdgeListIter edge_iter = g->edge_list_begin(),
			edge_iter_end = g->edge_list_end(); edge_iter != edge_iter_end; /* no increment */ ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();

		MinimizationEdge & minedge( static_cast< MinimizationEdge & > (**edge_iter) );

		core::graph::Graph::EdgeListIter edge_iter_next( edge_iter );
		++edge_iter_next;

		if ( minedge.any_active_enmeths() ) {
			minedge.setup_for_minimizing( pose.residue(node1), pose.residue(node2), pose, *this, min_map );
		} else {
			/// The edge will not contribute anything to scoring during minimization,
			/// so delete it from the graph, so we don't have to pay the expense of traversing
			/// through it.
			g->delete_edge( *edge_iter );
		}
		edge_iter = edge_iter_next;
	}

	/// 5.  Whole structure energies and energies that are opting out of the MinimizationGraph
	/// routines get a chance to setup for minimizing (using the entire pose as context) and
	for ( std::list< methods::EnergyMethodCOP >::const_iterator
			iter     = eval_derivs_with_pose_enmeths.begin(),
			iter_end = eval_derivs_with_pose_enmeths.end();
			iter != iter_end; ++iter ) {
		(*iter)->setup_for_minimizing( pose, *this, min_map );
		g->add_whole_pose_context_enmeth( *iter );
	}


	//std::cout << "Fixed energies: ";
	//fixed_energies.show_if_nonzero_weight( std::cout, weights_ );
	//std::cout << std::endl;

	g->set_fixed_energies( fixed_energies );
	pose.energies().set_minimization_graph( g );
}

/// @details Note that energy methods should be added to the MinimizationGraph
/// reguardless of whether or not they have been modernized to use the
/// eval_residue_derivatives/eval_residue_pair_derivatives machinery.
/// (Grandfathered EnergyMethods return "false" in their version of
/// minimize_in_whole_structure_context() and have their derivatives still
/// evaluated in the one-atom-at-a-time scheme.)  The reason these
/// grandfathered energy methods should be added to the MinimizationGraph is
/// to ensure that their energies are still computed during minimization-score-
/// function evaluation.
void
ScoreFunction::setup_for_minimizing_for_node(
	MinimizationNode & min_node,
	conformation::Residue const & rsd,
	kinematics::MinimizerMapBase const & min_map,
	pose::Pose & pose, // context
	bool accumulate_fixed_energies,
	EnergyMap & fixed_energies
) const
{
	kinematics::DomainMap const & domain_map( min_map.domain_map() );
	Size seqpos( rsd.seqpos() );
	/// 1a: context-independent 1body energies
	for ( CI_1B_Methods::const_iterator iter = ci_1b_methods_.begin(),
			iter_end = ci_1b_methods_.end(); iter != iter_end; ++iter ) {

		min_node.add_onebody_enmeth( *iter, rsd, pose, domain_map( seqpos ) );
		if ( domain_map( seqpos ) != 0 &&  accumulate_fixed_energies ) {
			fixed_energies.accumulate( pose.energies().onebody_energies( seqpos ), (*iter)->score_types() );
		}
	}
	/// 1b: context-dependent 1body energies
	for ( CD_1B_Methods::const_iterator iter = cd_1b_methods_.begin(),
			iter_end = cd_1b_methods_.end(); iter != iter_end; ++iter ) {
		// domain map check here?
		// continue; -- BUG!
		min_node.add_onebody_enmeth( *iter, rsd, pose, domain_map( seqpos ) );
	}
	/// 1c: context-independent 2body energies
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		// domain map check here?
		//std::cout << "seqpos " << seqpos << " " << domain_map( seqpos ) <<  std::endl;
		min_node.add_twobody_enmeth( *iter, rsd, pose, weights_, domain_map( seqpos ) );
		if ( domain_map( seqpos ) != 0 && accumulate_fixed_energies ) {
			//std::cout << "accumulating one body energies for " << seqpos << std::endl;
			fixed_energies.accumulate( pose.energies().onebody_energies( seqpos ), (*iter)->score_types() );
		}
	}

	/// 1d: context-dependent 2body energies
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		min_node.add_twobody_enmeth( *iter, rsd, pose, weights_, domain_map( seqpos ) );
	}
	/// 1e: all long-range 2body energies -- maybe this should be divided into ci and cd loops?
	for ( LR_2B_Methods::const_iterator iter = lr_2b_methods_.begin(),
			iter_end = lr_2b_methods_.end(); iter != iter_end; ++iter ) {
		min_node.add_twobody_enmeth( *iter, rsd, pose, weights_, domain_map( seqpos ) );
		if ( domain_map( seqpos ) == 0 && accumulate_fixed_energies && (*iter)->method_type() == methods::ci_lr_2b ) {
			fixed_energies.accumulate( pose.energies().onebody_energies( seqpos ), (*iter)->score_types() );
		}
	}
	min_node.setup_for_minimizing( rsd, pose, *this, min_map );
}

void
ScoreFunction::reinitialize_minnode_for_residue(
	MinimizationNode & min_node,
	conformation::Residue const & rsd,
	kinematics::MinimizerMapBase const & min_map,
	pose::Pose & pose // context
) const
{
	min_node.update_active_enmeths_for_residue( rsd, pose, weights_, min_map.domain_map()( rsd.seqpos() ) );
	min_node.setup_for_minimizing( rsd, pose, *this, min_map );
	/*
	kinematics::DomainMap const & domain_map( min_map.domain_map() );
	Size seqpos( rsd.seqpos() );
	/// 1a: context-independent 1body energies
	for ( CI_1B_Methods::const_iterator iter = ci_1b_methods_.begin(),
	iter_end = ci_1b_methods_.end(); iter != iter_end; ++iter ) {
	if ( domain_map( seqpos ) == 0 && ! (*iter)->minimize_in_whole_structure_context( pose ) ) {
	if ( (*iter)->defines_score_for_residue( rsd )) {
	(*iter)->setup_for_minimizing_for_residue(
	rsd, pose, *this, min_map,
	min_node.res_min_data() );
	}
	}
	}
	/// 1b: context-dependent 1body energies
	for ( CD_1B_Methods::const_iterator iter = cd_1b_methods_.begin(),
	iter_end = cd_1b_methods_.end(); iter != iter_end; ++iter ) {
	// domain map check here?
	if ( ! (*iter)->minimize_in_whole_structure_context( pose ) && (*iter)->defines_score_for_residue( rsd )) {
	(*iter)->setup_for_minimizing_for_residue(
	rsd, pose, *this, min_map,
	min_node.res_min_data() );
	}
	}
	/// 1c: context-independent 2body energies
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
	iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
	if ( ! (*iter)->minimize_in_whole_structure_context( pose ) ) {
	/// Only provide the setup opportunity to enmeths that are not relying on whole-structure
	/// minimization logic
	(*iter)->setup_for_minimizing_for_residue(
	rsd, pose, *this, min_map,
	min_node.res_min_data() );
	}
	}

	/// 1d: context-dependent 2body energies
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
	iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
	// domain map check here?
	if ( ! (*iter)->minimize_in_whole_structure_context( pose ) ) {
	(*iter)->setup_for_minimizing_for_residue(
	rsd, pose, *this, min_map,
	min_node.res_min_data() );
	}
	}
	/// 1e: all long-range 2body energies -- maybe this should be divided into ci and cd loops?
	for ( LR_2B_Methods::const_iterator iter = lr_2b_methods_.begin(),
	iter_end = lr_2b_methods_.end(); iter != iter_end; ++iter ) {
	// domain map check here? -- should separate CI and CD LR2B energies
	if ( ! (*iter)->minimize_in_whole_structure_context( pose ) ) {
	(*iter)->setup_for_minimizing_for_residue(
	rsd, pose, *this, min_map,
	min_node.res_min_data() );
	}
	}*/

}


void
ScoreFunction::setup_for_minimizing_sr2b_enmeths_for_minedge(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	MinimizationEdge & min_edge,
	kinematics::MinimizerMapBase const &,
	pose::Pose & pose,
	bool const res_moving_wrt_eachother,
	bool accumulate_fixed_energies,
	EnergyEdge const * energy_edge,
	EnergyMap & fixed_energies,
	Real const // apl -- remove this parameter edge_weight
) const
{
	debug_assert( ! accumulate_fixed_energies || energy_edge );

	/// First, the context-independent 2body energies
	if ( res_moving_wrt_eachother ) {
		for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
				iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
			min_edge.add_twobody_enmeth( *iter, res1, res2, pose, res_moving_wrt_eachother );
		}
	} else if ( accumulate_fixed_energies ) {
		energy_edge->add_to_energy_map( fixed_energies, ci_2b_types() );
	}

	/// Second, the context-dependent 2body energies
	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		if ( ! min_edge.add_twobody_enmeth( *iter, res1, res2, pose, res_moving_wrt_eachother )
				&& accumulate_fixed_energies )  {
			energy_edge->add_to_energy_map( fixed_energies, (*iter)->score_types() );
		}
	}
}

void
ScoreFunction::setup_for_lr2benmeth_minimization_for_respair(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	methods::LongRangeTwoBodyEnergyCOP lr2benergy,
	MinimizationGraph & g,
	kinematics::MinimizerMapBase const &,
	pose::Pose & pose,
	bool const res_moving_wrt_eachother,
	bool accumulate_fixed_energies,
	ResidueNeighborConstIteratorOP rni,
	EnergyMap & fixed_energies,
	Real const edge_weight,
	Real const edge_dweight
) const
{
	Size const seqpos1( res1.seqpos() );
	Size const seqpos2( res2.seqpos() );

	if ( res_moving_wrt_eachother &&
			lr2benergy->defines_score_for_residue_pair( res1, res2, res_moving_wrt_eachother ) ) {

		/// 3. iv Find the edge in the graph
		MinimizationEdge * minedge( static_cast< MinimizationEdge * > (g.find_edge( seqpos1, seqpos2)) );
		/// 3. v  Create a new edge if necessary
		if ( ! minedge ) minedge = static_cast< MinimizationEdge * > (g.add_edge( seqpos1, seqpos2));
		/// 3. vi. Now initialize this edge

		minedge->add_twobody_enmeth( lr2benergy, res1, res2, pose, res_moving_wrt_eachother );
		minedge->weight( edge_weight );
		minedge->dweight( edge_dweight );
	} else if ( accumulate_fixed_energies ) {

		/// Even if edge_weight != 1.0, don't scale the energies stored in the residue-neighbor-iterator,
		/// because, in the case of symmetric scoring, these energies have already been scaled.
		rni->accumulate_energy( fixed_energies );
	}

}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::eval_npd_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	Vector & F1,
	Vector & F2
) const
{
	//for ( AllMethods::const_iterator iter=all_methods_.begin(),
	//  iter_end= all_methods_.end(); iter != iter_end; ++iter ) {
	// (*iter)->eval_atom_derivative( atom_id, pose, domain_map, *this, weights_, F1, F2 );
	//}

	debug_assert( pose.energies().minimization_graph() );
	MinimizationGraphCOP mingraph = pose.energies().minimization_graph();

	/// Whole-pose-context energies should have their contribution calculated here.
	for ( MinimizationGraph::Energies::const_iterator
			iter     = mingraph->whole_pose_context_enmeths_begin(),
			iter_end = mingraph->whole_pose_context_enmeths_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_atom_derivative( atom_id, pose, domain_map, *this, weights_, F1, F2 );
	}

}

///////////////////////////////////////////////////////////////////////////////
Real
ScoreFunction::eval_dof_derivative(
	id::DOF_ID const & dof_id,
	id::TorsionID const & torsion_id,
	pose::Pose const & pose
) const
{

	/*for ( AllMethods::const_iterator iter=all_methods_.begin(),
	iter_end= all_methods_.end(); iter != iter_end; ++iter ) {
	//  for ( CI_1B_Methods::const_iterator iter = ci_1b_methods_.begin(),
	//    iter_end = ci_1b_methods_.end(); iter != iter_end; ++iter ) {
	deriv += (*iter)->eval_dof_derivative( dof_id, torsion_id, pose, *this, weights_ );
	}*/

	debug_assert( pose.energies().minimization_graph() );
	MinimizationGraphCOP mingraph = pose.energies().minimization_graph();

	Size const rsdno = torsion_id.valid() ? torsion_id.rsd() : dof_id.atom_id().rsd();
	conformation::Residue const & rsd = pose.residue( rsdno );

	MinimizationNode const & minnode = * mingraph->get_minimization_node( rsdno );

	return eval_dof_deriv_for_minnode( minnode, rsd, pose, dof_id, torsion_id, *this, weights_ );
}



///////////////////////////////////////////////////////////////////////////////
ScoreFunctionInfoOP
ScoreFunction::info() const
{
	if ( ! score_function_info_current_ ) {
		score_function_info_->initialize_from( *this );
		score_function_info_current_ = true;
	}
	return ScoreFunctionInfoOP( new ScoreFunctionInfo( *score_function_info_ ) );
}

ScoreFunction::AllMethods const & ScoreFunction::all_methods() const {
	return all_methods_;
}


///////////////////////////////////////////////////////////////////////////////
/// @details private -- handles setting the derived data

void
ScoreFunction::add_method( methods::EnergyMethodOP method )
{
	using namespace methods;

	all_methods_.push_back( method );

	// mapping from ScoreType to EnergyMethod
	for ( ScoreTypes::const_iterator
			iter     = method->score_types().begin(),
			iter_end = method->score_types().end(); iter != iter_end; ++iter ) {
		methods_by_score_type_[ *iter ] = method;
		score_types_by_method_type_[ method->method_type() ].push_back( *iter );
	}

	// PHIL replace these with utility::down_cast when it's working
	switch ( method->method_type() ) {
	case ci_2b :
		// there must be an easier way to cast these owning pointers!
		ci_2b_methods_.push_back
			( utility::pointer::static_pointer_cast< core::scoring::methods::ContextIndependentTwoBodyEnergy > ( method ) );
		break;

	case cd_2b :
		cd_2b_methods_.push_back
			( utility::pointer::static_pointer_cast< core::scoring::methods::ContextDependentTwoBodyEnergy > ( method ) );
		break;

	case ci_1b :
		ci_1b_methods_.push_back
			( utility::pointer::static_pointer_cast< core::scoring::methods::ContextIndependentOneBodyEnergy > ( method ) );
		break;

	case cd_1b :
		cd_1b_methods_.push_back
			( utility::pointer::static_pointer_cast< core::scoring::methods::ContextDependentOneBodyEnergy > ( method ) );
		break;

	case ci_lr_2b :
		ci_lr_2b_methods_.push_back(
			utility::pointer::static_pointer_cast< core::scoring::methods::ContextIndependentLRTwoBodyEnergy > ( method ) );
		lr_2b_methods_.push_back(
			utility::pointer::static_pointer_cast< core::scoring::methods::LongRangeTwoBodyEnergy > ( method ) );
		break;

	case cd_lr_2b :
		cd_lr_2b_methods_.push_back(
			utility::pointer::static_pointer_cast< core::scoring::methods::ContextDependentLRTwoBodyEnergy > ( method ) );
		lr_2b_methods_.push_back(
			utility::pointer::static_pointer_cast< core::scoring::methods::LongRangeTwoBodyEnergy > ( method ) );
		break;

	case ws :
		ws_methods_.push_back(
			utility::pointer::static_pointer_cast< core::scoring::methods::WholeStructureEnergy > ( method ) );
		break;

	default :
		utility_exit_with_message( "unrecognized method type " );
	} // switch
}


/// @details  Add a scoring method that's not necessarily included in the core library
/// @note  Currently only works for methods associated with a single specific score_type
void
ScoreFunction::add_extra_method(
	ScoreType const & new_type,
	Real const new_weight,
	methods::EnergyMethod const & new_method
)
{
	if ( weights_[ new_type]  != Real( 0.0 ) ||
			new_weight == Real(0.0) ||
			new_method.score_types().size() != 1 ||
			new_method.score_types().front() != new_type ) {
		utility_exit_with_message( "bad call to ScoreFunction::add_extra_method" );
	}
	weights_[ new_type ] = new_weight;
	add_method( new_method.clone() );
}


/// @details  Add a scoring method that's not necessarily included in the core library
///
void
ScoreFunction::add_extra_method(
	std::map< ScoreType, Real > const & new_weights,
	methods::EnergyMethod const & new_method
)
{
	for ( std::map< ScoreType, Real >::const_iterator it= new_weights.begin(); it != new_weights.end(); ++it ) {
		ScoreType const new_type( it->first );
		Real const new_weight( it->second );
		if ( ( weights_[ new_type ] != Real(0.0) ) ||
				( std::find( new_method.score_types().begin(), new_method.score_types().end(), new_type ) ==
				new_method.score_types().end() ) ) {
			utility_exit_with_message( "bad call to ScoreFunction::add_extra_method" );
		}
		weights_[ new_type ] = new_weight;
	}
	add_method( new_method.clone() );
}


///////////////////////////////////////////////////////////////////////////////
/// private -- handles setting the derived data

// This can/should be moved to utility/vector1.hh
template< class T >
inline
void
vector1_remove( utility::vector1< T > & v, T const & t )
{
	debug_assert( std::find( v.begin(), v.end(), t ) != v.end() );
	v.erase( std::find( v.begin(), v.end(), t ) );
}

void
ScoreFunction::remove_method( methods::EnergyMethodOP method )
{
	using namespace methods;
	vector1_remove( all_methods_, methods::EnergyMethodCOP(method) );

	// mapping from ScoreType to EnergyMethod
	for ( ScoreTypes::const_iterator
			iter     = method->score_types().begin(),
			iter_end = method->score_types().end(); iter != iter_end; ++iter ) {
		methods_by_score_type_[ *iter ] = 0;
		vector1_remove( score_types_by_method_type_[ method->method_type() ], *iter );
	}

	bool rebuild_lr_methods( false );
	// PHIL replace these with utility::down_cast when it's working
	switch ( method->method_type() ) {
	case ci_2b :
		// there must be an easier way to cast these owning pointers!
		vector1_remove( ci_2b_methods_,
			ContextIndependentTwoBodyEnergyOP(( utility::pointer::static_pointer_cast< core::scoring::methods::ContextIndependentTwoBodyEnergy > ( method ) )));
		break;

	case cd_2b :
		vector1_remove( cd_2b_methods_,
			ContextDependentTwoBodyEnergyOP(( utility::pointer::static_pointer_cast< core::scoring::methods::ContextDependentTwoBodyEnergy > ( method ) ) ));
		break;

	case ci_1b :
		vector1_remove( ci_1b_methods_,
			ContextIndependentOneBodyEnergyOP(( utility::pointer::static_pointer_cast< core::scoring::methods::ContextIndependentOneBodyEnergy > ( method ) )));
		break;

	case cd_1b :
		vector1_remove( cd_1b_methods_,
			ContextDependentOneBodyEnergyOP(( utility::pointer::static_pointer_cast< core::scoring::methods::ContextDependentOneBodyEnergy > ( method ) )));
		break;

	case ci_lr_2b :
		rebuild_lr_methods = true;
		vector1_remove( ci_lr_2b_methods_,
			ContextIndependentLRTwoBodyEnergyOP(utility::pointer::static_pointer_cast< core::scoring::methods::ContextIndependentLRTwoBodyEnergy > ( method )));
		break;

	case cd_lr_2b :
		rebuild_lr_methods = true;
		vector1_remove( cd_lr_2b_methods_,
			ContextDependentLRTwoBodyEnergyOP( utility::pointer::static_pointer_cast< core::scoring::methods::ContextDependentLRTwoBodyEnergy > ( method )));
		break;

	case ws :
		vector1_remove( ws_methods_,
			WholeStructureEnergyOP( utility::pointer::static_pointer_cast< core::scoring::methods::WholeStructureEnergy > ( method )));
		break;

	default :
		utility_exit_with_message( "unrecognized method type " );
	} // switch

	if ( rebuild_lr_methods )  {
		lr_2b_methods_.clear();
		for ( CD_LR_2B_Methods::const_iterator iter = cd_lr_2b_methods_.begin(),
				iter_end = cd_lr_2b_methods_.end(); iter != iter_end; ++iter ) {
			lr_2b_methods_.push_back( *iter );
		}
		for ( CI_LR_2B_Methods::const_iterator iter = ci_lr_2b_methods_.begin(),
				iter_end = ci_lr_2b_methods_.end(); iter != iter_end; ++iter ) {
			lr_2b_methods_.push_back( *iter );
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
ScoreFunction::initialize_methods_arrays()
{
	// clear (may not be necessary)
	all_methods_.clear();
	methods_by_score_type_.clear();
	score_types_by_method_type_.clear();

	ci_2b_methods_.clear();
	cd_2b_methods_.clear();
	ci_1b_methods_.clear();
	cd_1b_methods_.clear();
	ci_lr_2b_methods_.clear();
	cd_lr_2b_methods_.clear();
	lr_2b_methods_.clear();
	ws_methods_.clear();

	// some need resizing
	methods_by_score_type_.resize( n_score_types, 0 );
	//std::fill( methods_by_score_type_.begin(), methods_by_score_type_.end(),  ); // wipe old pointers
	for ( Size ii = 1; ii <= n_score_types; ++ii ) { methods_by_score_type_[ ii ] = 0; }
	score_types_by_method_type_.resize( methods::n_energy_method_types );
	for ( Size ii = 1; ii <= methods::n_energy_method_types; ++ii ) {
		score_types_by_method_type_[ ii ].clear();
	}
}

/// @details determines the furthest-reach of the short-range two body energies
/// as well as the whole-structure energies, which are allowed to require
/// that the EnergyGraph have edges of a minimum length.
Distance
ScoreFunction::max_atomic_interaction_cutoff() const {
	Distance max_cutoff = 0;

	for ( CD_2B_Methods::const_iterator iter = cd_2b_methods_.begin(),
			iter_end = cd_2b_methods_.end(); iter != iter_end; ++iter ) {
		if ( (*iter)->atomic_interaction_cutoff() > max_cutoff ) {
			max_cutoff = (*iter)->atomic_interaction_cutoff();
		}
	}
	for ( CI_2B_Methods::const_iterator iter = ci_2b_methods_.begin(),
			iter_end = ci_2b_methods_.end(); iter != iter_end; ++iter ) {
		if ( (*iter)->atomic_interaction_cutoff() > max_cutoff ) {
			max_cutoff = (*iter)->atomic_interaction_cutoff();
		}
	}
	for ( WS_Methods::const_iterator iter = ws_methods_.begin(),
			iter_end = ws_methods_.end(); iter != iter_end; ++iter ) {
		if ( (*iter)->atomic_interaction_cutoff() > max_cutoff ) {
			max_cutoff = (*iter)->atomic_interaction_cutoff();
		}
	}

	return max_cutoff;
}

/// @brief scoring function fills in the context graphs that its energy methods require
///
/// input vector should be false and have num_context_graph_types slots.  Each method
/// ors its required context graphs
void
ScoreFunction::indicate_required_context_graphs(
	utility::vector1< bool > & context_graphs_required
) const
{
	for ( AllMethods::const_iterator it=all_methods_.begin(),
			it_end = all_methods_.end(); it != it_end; ++it ) {
		(*it)->indicate_required_context_graphs(context_graphs_required);
	}
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// @brief Utility function to locate a weights or patch file, either with a fully qualified path,
/// in the local directory, or in the database. Names may be passed either with or without the
/// optional extension.
std::string
find_weights_file(std::string name, std::string extension/*=".wts"*/) {
	utility::io::izstream data1( name );
	if ( data1.good() ) {
		return name;
	} else {
		utility::io::izstream data2( name + extension );
		if ( data2.good() ) {
			return name + extension;
		} else {
			utility::io::izstream data3(  basic::database::full_name( "scoring/weights/"+name+extension, /*warn=*/false )  );
			if ( data3.good() ) {
				return basic::database::full_name( "scoring/weights/"+name+extension );
			} else {
				utility::io::izstream data4(  basic::database::full_name( "scoring/weights/"+name, false )  );
				if ( data4.good() ) {
					return basic::database::full_name( "scoring/weights/"+name );
				} else {
					utility_exit_with_message( "Unable to open weights/patch file. None of (./)" + name + " or " +
						"(./)" + name + extension + " or " +
						basic::database::full_name( "scoring/weights/"+name, false )  + " or " +
						basic::database::full_name( "scoring/weights/"+name+extension, false )  + " exist"  );
					return "invalid"; // To make the compiler happy - should never reach here.
				}
			}
		}
	}
}

/// @brief check order of methods
bool
ScoreFunction::check_methods_in_right_order( ScoreType const & score_type_in_first_method,
	ScoreType const & score_type_in_second_method ) const
{
	bool first_method_found( false );
	for ( AllMethods::const_iterator it=all_methods_.begin(),
			it_end = all_methods_.end(); it != it_end; ++it ) {

		if ( ( (*it)->score_types() ).has_value( score_type_in_first_method ) ) first_method_found = true;

		if ( ( (*it)->score_types() ).has_value( score_type_in_second_method ) &&  !first_method_found ) return false;

	}
	return true;
}

bool
ScoreFunction::ready_for_nonideal_scoring() const
{
	// if any of cart_bonded terms are on and pro_close is off, return true
	if ( (has_nonzero_weight( cart_bonded ) ||
			has_nonzero_weight( cart_bonded_angle ) ||
			has_nonzero_weight( cart_bonded_length ) ||
			has_nonzero_weight( cart_bonded_torsion )) &&
			has_zero_weight( pro_close ) ) {
		return true;
	}

	// if any of the mm_ terms are on, return true
	if ( has_nonzero_weight( mm_lj_intra_rep ) ||
			has_nonzero_weight( mm_lj_intra_atr ) ||
			has_nonzero_weight( mm_lj_inter_rep ) ||
			has_nonzero_weight( mm_lj_inter_atr ) ||
			has_nonzero_weight( mm_twist ) ||
			has_nonzero_weight( mm_bend ) ||
			has_nonzero_weight( mm_stretch ) ) {
		return true;
	}

	// if any of the mm_ terms are on, return true
	if ( has_nonzero_weight( bond_geometry ) ||
			has_nonzero_weight( rna_bond_geometry ) ) {
		return true;
	}

	return false;
}

} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::ScoreFunction::save( Archive & arc ) const {
	arc( CEREAL_NVP( name_ ) ); // std::string
	arc( CEREAL_NVP( weights_ ) ); // EnergyMap
	arc( CEREAL_NVP( energy_method_options_ ) ); // methods::EnergyMethodOptionsOP

	// The following data members are basically derived from the ones above, and
	// therefore, should not be serialized
	// EXEMPT ci_2b_methods_ cd_2b_methods_ ci_1b_methods_ cd_1b_methods_
	// EXEMPT ci_lr_2b_methods_ cd_lr_2b_methods_
	// EXEMPT lr_2b_methods_ ws_methods_ all_methods_ methods_by_score_type_
	// EXEMPT score_types_by_method_type_ score_function_info_current_
	// EXEMPT score_function_info_ any_intrares_energies_
	// EXEMPT ci_2b_intrares_ cd_2b_intrares_
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::ScoreFunction::load( Archive & arc ) {
	arc( name_ ); // std::string
	arc( weights_ ); // EnergyMap
	arc( energy_method_options_ ); // methods::EnergyMethodOptionsOP

	// The following data members are basically derived from the ones above, and
	// therefore, should not be serialized
	// EXEMPT ci_2b_methods_ cd_2b_methods_ ci_1b_methods_ cd_1b_methods_
	// EXEMPT ci_lr_2b_methods_ cd_lr_2b_methods_
	// EXEMPT lr_2b_methods_ ws_methods_ all_methods_ methods_by_score_type_
	// EXEMPT score_types_by_method_type_ score_function_info_current_
	// EXEMPT score_function_info_ any_intrares_energies_
	// EXEMPT ci_2b_intrares_ cd_2b_intrares_

	for ( core::Size ii = 1; ii <= n_score_types; ++ii ) {
		set_weight( ScoreType( ii ), weights_[ ScoreType( ii ) ] );
	}

}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::ScoreFunction );
CEREAL_REGISTER_TYPE( core::scoring::ScoreFunction )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_ScoreFunction )
#endif // SERIALIZATION
