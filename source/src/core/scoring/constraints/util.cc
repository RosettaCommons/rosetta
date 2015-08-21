// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/util.cc
/// @brief utility functions for defining and using constraints.
/// @author James Thompson

#include <core/scoring/constraints/util.hh>
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/options/option.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AmbiguousNMRConstraint.hh>

#include <core/conformation/Residue.hh>

#include <numeric/random/random_permutation.hh>
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray2D.hh>

// option key includes
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace constraints {


static thread_local basic::Tracer tr( "core.scoring.constraints.util" );

/// @brief Returns the weighted value of a normal distribution evaluated
///  with the given mean, sd, and x values. Returns zero if the weight is
///  less than 1e-10.
Real logdgaussian( Real x, Real mean, Real sd, Real weight )
{
	using namespace std;
	Real r = x - mean;
	static Real sqrt_2pi = 2.50662721600161;

	Real answer = log(weight / (sd * sqrt_2pi)) - ( r * r / (2*sd*sd) );
	return answer;
}

/// @brief Returns the weighted value of a normal distribution evaluated
///  with the given mean, sd, and x values. Returns zero if the weight is
///  less than 1e-10.
Real logdgaussian_deriv( Real x, Real mean, Real sd, Real )
{
	using namespace std;
	Real r = x - mean;
	Real answer =  - (r / (sd*sd) );
	return answer;
}

/// @brief Returns the weighted value of a normal distribution evaluated
///  with the given mean, sd, and x values. Returns zero if the weight is
///  less than 1e-10.
Real dgaussian( Real x, Real mean, Real sd, Real weight ) {
	if ( weight < 1e-10 ) {
		return 0;
	}

	using namespace std;
	Real r = x - mean;
	static Real sqrt_2pi = 2.50662721600161;
	Real answer = weight * (1 / (sd * sqrt_2pi)) * exp( -1 * r * r / (2*sd*sd) );
	return answer;
}

/// @brief Returns the weighted derivative of a normal distribution evaluated
/// with the given mean, sd, and x values. Returns zero if the weight is less
/// than 1e-10.
Real gaussian_deriv( Real x, Real mean, Real sd, Real weight ) {
	debug_assert( weight >= 0.0 && weight <= 1.0 );

	if ( weight < 1e-10 ) {
		return 0;
	}

	using namespace std;
	Real r = std::abs( x - mean );
	Real answer = dgaussian( x, mean, sd , weight ) * ( 1 / ( sd * sd ) ) * r;
	return answer;
}

/// @brief Returns the weighted value of an exponential distribution
/// evaluated with the given anchor, rate, and x values. Returns zero if the
/// weight is less than 1e-10.
Real dexponential( Real x, Real anchor, Real rate, Real weight ) {
	debug_assert( weight >= 0.0 && weight <= 1.0 );
	if ( weight < 1e-10 ) {
		return 0;
	}

	using namespace std;
	Real r = std::abs( x - anchor );
	Real answer = weight * rate * exp( -1 * rate * r );
	return answer;
}

/// @brief Returns the weighted derivative of a log-exponential distribution
/// evaluated with the given anchor, rate, and x values. Returns zero if the
/// weight is less than 1e-10.
Real exponential_deriv( Real x, Real anchor, Real rate, Real weight ) {
	if ( weight < 1e-10 ) {
		return 0;
	}

	using namespace std;
	Real r = std::abs( x - anchor );
	return weight  * rate * rate * exp( -1 * rate * r );
}

Real linear_interpolate(
	Real const x_val,
	Real const x1,
	Real const x2,
	Real const y1,
	Real const y2
) {

	if ( x_val == x1 ) return y1;
	if ( x_val == x2 ) return y2;

	// calculate slope
	Real slope = ( y2 - y1 ) / ( x2 - x1 );
	// walk along line
	return (x_val - x1) * slope;
}

void
cull_violators(
	ConstraintCOPs const& target_list,
	ConstraintCOPs &culled_list,
	pose::Pose const& filter_pose,
	core::Real threshold
) {
	culled_list.clear();
	for ( ConstraintCOPs::const_iterator it = target_list.begin(),
			eit = target_list.end(); it != eit; ++it ) {
		if ( (*it)->show_violations( tr.Debug, filter_pose, 1, threshold ) == 0 ) {
			culled_list.push_back( *it );
		}
	}
}

/////////////////////////////////////////////////////////////////////////
//////////////This block is for add-constraints-from-command-line utilities
//////////////They REPLACE whatever constraint values may exist
/////////////////////////////////////////////////////////////////////////

////////// Centroid constraints
std::string get_cst_file_option(){
	using namespace basic::options;
	utility::vector1< std::string>  cst_files = option[ OptionKeys::constraints::cst_file ]();
	core::Size choice=1;
	if ( cst_files.size() > 1 ) choice=core::Size( numeric::random::rg().random_range( 1,cst_files.size() ) );
	tr.Info << "Constraint choice: " << cst_files[choice] << std::endl;
	return cst_files[choice];
}

/// @details add constraints from command line to POSE ONLY, if cst file is supplied by user.  Overwrites any constraints already in the Pose.  Assumed to be "centroid constraints" using the cst_file flag (will work fine with fa constraints, but uses the not-fa command line option).
void add_constraints_from_cmdline_to_pose( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace core::scoring::constraints;
	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		ConstraintSetOP cstset_ = ConstraintIO::get_instance()->read_constraints( get_cst_file_option() ,ConstraintSetOP( new ConstraintSet ), pose );
		pose.constraint_set( cstset_ );
	}
}

/// @details add constraints from command line to SCOREFUNCTION ONLY, if cst file is supplied by user.  Overwrites any constraint weights already in the scorefunction.  Assumed to be "centroid constraints" using the cst_weight flag (will work fine with fa constraints, but uses the not-fa command line option).
void add_constraints_from_cmdline_to_scorefxn( core::scoring::ScoreFunction &scorefxn_  ) {
	using namespace basic::options;
	if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		scorefxn_.set_weight( atom_pair_constraint,  option[ OptionKeys::constraints::cst_weight ]() );
		scorefxn_.set_weight( angle_constraint,      option[ OptionKeys::constraints::cst_weight ]() );
		scorefxn_.set_weight( dihedral_constraint,   option[ OptionKeys::constraints::cst_weight ]() );
		scorefxn_.set_weight( coordinate_constraint, option[ OptionKeys::constraints::cst_weight ]() );
	}
}

/// @details add constraints from command line to SCOREFUNCTION AND POSE BOTH, if cst file is supplied by user.  Overwrites any constraints present in either.  Assumed to be "centroid constraints" using cst_file & cst_weight (will work fine with fa constraints, but uses the not-fa command line options).
void add_constraints_from_cmdline( core::pose::Pose & pose, core::scoring::ScoreFunction &scorefxn_  ) {
	add_constraints_from_cmdline_to_pose( pose );
	add_constraints_from_cmdline_to_scorefxn( scorefxn_ );
}

////////// FA constraints
std::string get_cst_fa_file_option() {
	using namespace basic::options;
	utility::vector1< std::string> cst_files
		= option[ OptionKeys::constraints::cst_fa_file ]();
	core::Size choice=1;
	if ( cst_files.size() > 1 ) choice=core::Size( numeric::random::rg().random_range( 2,cst_files.size() ) );
	tr.Info << "Constraint choice: " << cst_files[choice] << std::endl;
	return cst_files[choice];
}

/// @details add constraints from command line to POSE ONLY, if cst file is supplied by user.  Overwrites any constraints already in the Pose.  Assumed to be "fullatom constraints" using the cst_fa_file flag (will work fine with centroid constraints otherwise)
void add_fa_constraints_from_cmdline_to_pose( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace core::scoring::constraints;
	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		ConstraintSetOP cstset_ = ConstraintIO::get_instance()->read_constraints(
			get_cst_fa_file_option(), ConstraintSetOP( new ConstraintSet ), pose
		);
		pose.constraint_set( cstset_ );
	}
}

/// @details add constraints from command line to SCOREFUNCTION ONLY, if cst file is supplied by user.  Overwrites any constraint weights already in the scorefunction.  Assumed to be "fullatom constraints" because it uses the cst_fa_weight flag, but will work with centroids otherwise.
void add_fa_constraints_from_cmdline_to_scorefxn( core::scoring::ScoreFunction & scorefxn_  ) {
	using namespace basic::options;
	if ( option[ OptionKeys::constraints::cst_fa_weight ].user() ) {
		scorefxn_.set_weight( atom_pair_constraint,  option[ OptionKeys::constraints::cst_fa_weight ]() );
		scorefxn_.set_weight( angle_constraint,      option[ OptionKeys::constraints::cst_fa_weight ]() );
		scorefxn_.set_weight( dihedral_constraint,   option[ OptionKeys::constraints::cst_fa_weight ]() );
		scorefxn_.set_weight( coordinate_constraint, option[ OptionKeys::constraints::cst_fa_weight ]() );
	}
}


/// @details add constraints from command line to SCOREFUNCTION AND POSE BOTH, if cst file is supplied by user.  Overwrites any constraints present in either.  Assumed to be "fullatom constraints" because it uses the cst_fa_file & cst_fa_weight flags, but will work with centroids otherwise.
void add_fa_constraints_from_cmdline(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction & scorefxn_
) {
	add_fa_constraints_from_cmdline_to_pose( pose );
	add_fa_constraints_from_cmdline_to_scorefxn( scorefxn_ );
}

/////////////////////////////////////////////////////////////////////////
//////////////This block is for add-constraints-from-command-line utilities
//////////////They MERGE constraints from cmdline with whatever preexists
/////////////////////////////////////////////////////////////////////////

////////// Centroid constraints

/// @details merge centroid constraints from commandline (cst_file) to pose: read the user-specified constraints file, and merge those to the Pose's ConstraintSet.  Creates a ConstraintSet normally, then merges it to the Pose's set.
void merge_constraints_from_cmdline_to_pose( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace core::scoring::constraints;
	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		ConstraintSetCOP const new_cstset = ConstraintIO::get_instance()->read_constraints( get_cst_file_option(), ConstraintSetOP( new ConstraintSet ), pose );
		pose.add_constraints( new_cstset->get_all_constraints() );
	}
}

/// @details weight merging helper function - only sets previously-nonzero weights!
void merge_csts_to_scorefunction( core::Real const weight, core::scoring::ScoreFunction &scorefxn_ ){

	if ( scorefxn_.has_zero_weight( atom_pair_constraint  ) ) scorefxn_.set_weight( atom_pair_constraint , weight );
	else {
		tr.Warning
			<< "WARNING: merge_csts_to_scorefunction not overwriting atom_pair_constraint weight with file weight "
			<< weight
			<< ", instead leaving untouched pre-existing scorefunction weight "
			<< scorefxn_.get_weight(atom_pair_constraint) << std::endl;
	}

	if ( scorefxn_.has_zero_weight( angle_constraint      ) ) scorefxn_.set_weight( angle_constraint     , weight );
	else {
		tr.Warning
			<< "WARNING: merge_csts_to_scorefunction not overwriting angle_constraint weight with file weight "
			<< weight
			<< ", instead leaving untouched pre-existing scorefunction weight "
			<< scorefxn_.get_weight(angle_constraint) << std::endl;
	}

	if ( scorefxn_.has_zero_weight( dihedral_constraint   ) ) scorefxn_.set_weight( dihedral_constraint  , weight );
	else {
		tr.Warning
			<< "WARNING: merge_csts_to_scorefunction not overwriting dihedral_constraint weight with file weight "
			<< weight
			<< ", instead leaving untouched pre-existing scorefunction weight "
			<< scorefxn_.get_weight(dihedral_constraint) << std::endl;
	}

	if ( scorefxn_.has_zero_weight( coordinate_constraint ) ) scorefxn_.set_weight( coordinate_constraint, weight );
	else {
		tr.Warning
			<< "WARNING: merge_csts_to_scorefunction not overwriting coordinate_constraint weight with file weight "
			<< weight
			<< ", instead leaving untouched pre-existing scorefunction weight "
			<< scorefxn_.get_weight(coordinate_constraint) << std::endl;
	}


	return;
}


/// @details "merge" constraint weights for scorefunction.  I don't know what "merge" means here, so it will only modify the zero weights; nonzero weights are untouched.  Reads -cst_weight from command line, and sets nonzero constraint weights to the new value
void merge_constraints_from_cmdline_to_scorefxn( core::scoring::ScoreFunction &scorefxn_  ) {
	using namespace basic::options;
	if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		merge_csts_to_scorefunction( option[ OptionKeys::constraints::cst_weight ].value(), scorefxn_);
	}

	return;
}

/// @details perform both pose and scorefunction updates in merge mode for centroid flags -cst_file and -cst_weight
void merge_constraints_from_cmdline( core::pose::Pose & pose, core::scoring::ScoreFunction &scorefxn_  ) {
	merge_constraints_from_cmdline_to_pose( pose );
	merge_constraints_from_cmdline_to_scorefxn( scorefxn_ );
}

////////// FA constraints

/// @details merge fullatom constraints from commandline (cst_fa_file) to pose: read the user-specified constraints file, and merge those to the Pose's ConstraintSet.  Creates a ConstraintSet normally, then merges it to the Pose's set.
void merge_fa_constraints_from_cmdline_to_pose( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace core::scoring::constraints;
	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		ConstraintSetCOP const new_cstset = ConstraintIO::get_instance()->read_constraints( get_cst_fa_file_option(), ConstraintSetOP( new ConstraintSet ), pose );
		pose.add_constraints( new_cstset->get_all_constraints() );
	}
}

/// @details "merge" constraint weights for scorefunction.  I don't know what "merge" means here, so it will only modify the zero weights; nonzero weights are untouched.  Reads -cst_fa_weight from command line, and sets nonzero constraint weights to the new value
void merge_fa_constraints_from_cmdline_to_scorefxn( core::scoring::ScoreFunction &scorefxn_  ) {
	using namespace basic::options;
	if ( option[ OptionKeys::constraints::cst_fa_weight ].user() ) {
		merge_csts_to_scorefunction( option[ OptionKeys::constraints::cst_fa_weight ].value(), scorefxn_);
	}

	return;
}

/// @details perform both pose and scorefunction updates in merge mode for centroid flags -cst_fa_file and -cst_fa_weight
void merge_fa_constraints_from_cmdline(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction & scorefxn_
) {
	merge_fa_constraints_from_cmdline_to_pose( pose );
	merge_fa_constraints_from_cmdline_to_scorefxn( scorefxn_ );
}

///////////////////////////////////////////////////////////////////////////////////////////////
void
add_coordinate_constraints( pose::Pose & pose, Real const coord_sdev /* = 10.0 */, bool include_sc /* true */) {

	using namespace core::id;
	using namespace core::conformation;
	using namespace core::scoring::constraints;

	Size const my_anchor( pose.fold_tree().root() ); //Change to use the root of the current foldtree as done by Rocco in AtomCoordinateCstMover - JAB.

	ConstraintSetOP cst_set = pose.constraint_set()->clone();

	Size const nres( pose.total_residue() );
	for ( Size i=1; i<= nres;  ++i ) {

		Residue const & i_rsd( pose.residue(i) );

		core::Size last_atom = i_rsd.last_backbone_atom();
		if ( include_sc ) {
			last_atom = i_rsd.nheavyatoms();
		}
		for ( Size ii = 1; ii <= last_atom; ++ii ) {

			func::FuncOP f( new func::HarmonicFunc( 0.0, coord_sdev ) );
			cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( AtomID(ii,i), AtomID(1,my_anchor), i_rsd.xyz(ii), f ) ) ) );
		}
	}

	pose.constraint_set( cst_set );

}

void
add_coordinate_constraints( pose::Pose & pose, core::Size const start_res, core::Size const end_res, Real const coord_sdev /* = 10.0 */, bool include_sc /* true*/) {
	using namespace core::id;
	using namespace core::conformation;
	using namespace core::scoring::constraints;

	Size const my_anchor( pose.fold_tree().root() ); //Change to use the root of the current foldtree as done by Rocco in AtomCoordinateCstMover - JAB.

	ConstraintSetOP cst_set = pose.constraint_set()->clone();

	for ( Size i=start_res; i<= end_res;  ++i ) {

		Residue const & i_rsd( pose.residue(i) );

		core::Size last_atom = i_rsd.last_backbone_atom();
		if ( include_sc ) {
			last_atom = i_rsd.nheavyatoms();
		}
		for ( Size ii = 1; ii<= last_atom; ++ii ) {

			func::FuncOP f( new func::HarmonicFunc( 0.0, coord_sdev ) );
			cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( AtomID(ii,i), AtomID(1,my_anchor), i_rsd.xyz(ii), f ) ) ) );
		}
	}

	pose.constraint_set( cst_set );
}


void
remove_constraints_of_type(core::pose::Pose & pose, std::string const & type) {
	utility::vector1< ConstraintCOP > all_csts = pose.constraint_set()->get_all_constraints();
	for ( core::Size i=1; i<=all_csts.size(); ++i ) {
		if ( all_csts[i]->type() == type ) {
			pose.remove_constraint(all_csts[i], true);
		}
	}
}

void
remove_constraints_of_type(core::pose::Pose & pose, std::string const & type, core::Size const start_res, core::Size const end_res){
	utility::vector1< ConstraintCOP > all_csts = pose.constraint_set()->get_all_constraints();
	for ( core::Size i=1; i<=all_csts.size(); ++i ) {
		if ( all_csts[i]->type() == type ) {
			utility::vector1< core::Size > residues = all_csts[i]->residues();

			for ( core::Size x = 1; x <= residues.size(); ++x ) {
				if ( (start_res <= residues[x]) && (residues[x] <= end_res) ) {
					pose.remove_constraint(all_csts[i], true);
					break;
				}
			}
		}
	}
}


void
remove_nonbb_constraints( pose::Pose & pose) {

	using namespace core::id;
	using namespace core::conformation;
	using namespace core::scoring::constraints;

	utility::vector1<ConstraintCOP> cst_set = pose.constraint_set()->get_all_constraints();

	for ( ConstraintCOPs::const_iterator it = cst_set.begin(), eit = cst_set.end(); it!=eit; ++it ) {
		for ( core::Size i=1; i<=(*it)->natoms(); i++ ) {
			if ( !pose.residue((*it)->atom(i).rsd()).atom_is_backbone((*it)->atom(i).atomno()) ) {
				pose.remove_constraint(*it, false);
				break;
			} //remove non-backbone constraint
		}//loop through residues in constraints
	}//loop through constraint
}
bool combinable( Constraint const& cst, utility::vector1< Size > exclude_res ) {
	if ( exclude_res.size() == 0 ) return true;
	utility::vector1< Size > pos_list( cst.residues() );
	for ( core::Size i(1); i <= pos_list.size(); ++i ) {
		Size const seqpos( pos_list[i] );
		runtime_assert( seqpos <= exclude_res.size() );
		if ( !exclude_res[ seqpos ] ) {
			return true;
		}
	}
	return false;
}

/// @brief combine constraints randomly into AmbiguousConstraints N -> 1 this greatly decreases the odds to have a wrong constraint
void choose_effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp, ConstraintCOPs& in ) {
	for ( utility::vector1< ConstraintCOP >::const_iterator it = in.begin(); it != in.end(); ++it ) {
		Constraint& cst = const_cast< Constraint& >( (**it) );
		cst.choose_effective_sequence_separation( sp, numeric::random::rg() );
	}
}


inline core::Size bin_by_seq_separation( core::Size sep ) {
	if ( sep < 5 ) return 1;
	if ( sep < 20 ) return 2;
	if ( sep < 50 ) return 3;
	return 4;
}

/// @brief combine constraints randomly into AmbiguousConstraints N -> 1 this greatly decreases the odds to have a wrong constraint
void combine_constraints(
	ConstraintCOPs& in,
	core::Size combine_ratio,
	utility::vector1< bool > exclude_res,
	core::kinematics::ShortestPathInFoldTree const& sp_in
) {

	tr.Info << " combine constraints " << combine_ratio << " -> 1 " << std::endl;
	using namespace scoring::constraints;
	if ( combine_ratio <= 1 ) return;

	tr.Trace << "figure out sequence-separation bins" << std::endl;
	//first bin constraints by effective sequence separation  --- combine within bins
	typedef std::map< core::Size, ConstraintCOPs > SeqSepMap;
	SeqSepMap seq_sep_map;
	for ( utility::vector1< ConstraintCOP >::const_iterator it = in.begin(); it != in.end(); ++it ) {
		Size seq_bin = bin_by_seq_separation( (*it)->effective_sequence_separation( sp_in ) );
		seq_sep_map[ seq_bin ].push_back( *it );
	}

	tr.Trace << "combine within bins..."<< std::endl;
	// combine constraints within bins
	utility::vector1< ConstraintCOP > out;
	for ( SeqSepMap::iterator bin_it=seq_sep_map.begin(); bin_it != seq_sep_map.end(); ++bin_it ) {
		// random permutation within bin... ensures random combination
		random_permutation( bin_it->second, numeric::random::rg() ); //I don't think a single pass of pairwise exchanges is enough to randomize the vector.
		random_permutation( bin_it->second, numeric::random::rg() );
		random_permutation( bin_it->second, numeric::random::rg() );

		// combine bin
		for ( utility::vector1< ConstraintCOP >::const_iterator it = bin_it->second.begin(); it != bin_it->second.end();
				//DO NOT INCREMENT, already incremented in next loop -- gives segfaults otherwise
				) {
			Size ct( combine_ratio );
			MultiConstraintOP combined_cst( new AmbiguousConstraint );
			for ( ; ct > 0 && it != bin_it->second.end(); ++it ) {
				tr.Trace << " add constraint " << ct << std::endl;
				//check if constraint is combinable:
				if ( combinable( **it, exclude_res ) ) {
					combined_cst->add_individual_constraint( *it );
					--ct;
				} else {
					out.push_back( *it ); //keep uncombined constraints around
				}
			}
			//fill up with more constraints if ct is not 0 yet.
			if ( ct > 0 ) {
				tr.Trace << " fill last Ambiguous constraint " << ct << std::endl;
				for ( ConstraintCOPs::const_iterator it2 = bin_it->second.begin(); ct > 0 && it2 != bin_it->second.end(); ++it2 ) {
					if ( combinable( **it2, exclude_res ) ) {
						--ct;
						combined_cst->add_individual_constraint( *it2 );
					}
				}
			} // ct > 0
			combined_cst->choose_effective_sequence_separation( sp_in, numeric::random::rg() );
			out.push_back( combined_cst );
		}
	}  // combination within bin
	in = out;
}

//helper function for skip_redundant_cosntraints
void count_constraint( ConstraintCOP cst, bool redundant, ObjexxFCL::FArray2D_int& count_matrix, Size influence_width, Size total_residue ){

	// figure out if it's inter-res, residue_pair, or 3+body
	utility::vector1< int > pos_list( cst->residues() );

	if ( pos_list.size() != 2 ) {
		tr.Error << "problems understanding constraint in skip_redundant_constraints ... ignore and keep this one" << std::endl;
		for ( utility::vector1< int >::const_iterator it = pos_list.begin(); it != pos_list.end(); ++it ) {
			tr.Debug << " resid: " << *it << std::endl;
		}
		return;
	}
	//  count_matrix( pos_list[1],pos_list[2] ) = std::max(  count_matrix( pos_list[1],pos_list[2] ), redundant ? 1 : 2 );
	//  count_matrix( pos_list[2],pos_list[1] ) = count_matrix( pos_list[1],pos_list[2] );
	for ( int i=pos_list[1]-influence_width; i<= pos_list[1]+ (int) influence_width ; ++i ) {
		for ( int j=pos_list[2]-influence_width; j<= pos_list[2]+ (int) influence_width ; ++j ) {
			if ( i < 1 || j < 1 || i > (int) total_residue || j > (int) total_residue ) continue;
			count_matrix( i, j ) = std::max( count_matrix( i,j ), redundant ? 1 : 2 );
			count_matrix( j, i ) = count_matrix( i,j );
		}
	}
}

//helper function for skip_redundant_cosntraints
bool keep_constraint( ConstraintCOP cst, bool redundant, ObjexxFCL::FArray2D_int& count_matrix, Size influence_width, Size total_residue ) {

	// figure out if it's inter-res, residue_pair, or 3+body
	utility::vector1< int > pos_list( cst->residues() );

	if ( pos_list.size() != 2 ) {
		tr.Error << "problems understanding constraint in skip_redundant_constraints ... ignore and keep this one" << std::endl;
		return true;
	}
	bool keep( false );

	if ( ( !redundant && count_matrix( pos_list[1],pos_list[2] ) == 3 ) || count_matrix( pos_list[1],pos_list[2] ) == ( redundant ? 1 : 2 ) ) {
		keep = true;
		for ( int i=pos_list[1]-influence_width; i<=  pos_list[1]+ (int) influence_width ; ++i ) {
			for ( int j=pos_list[2]-influence_width; j<= (int) pos_list[2]+ (int)influence_width ; ++j ) {
				if ( i < 1 || j < 1 || i > (int) total_residue || j > (int) total_residue ) continue;
				count_matrix( i, j ) += 2;
				count_matrix( j, i ) = count_matrix( i,j );
			}
		}
	}
	// count_matrix( pos_list[1],pos_list[2] ) += 2;
	//   count_matrix( pos_list[2],pos_list[1] ) = count_matrix( pos_list[1],pos_list[2] );
	return keep;
}

void skip_redundant_constraints( ConstraintCOPs& in, Size total_residue, Size influence_width ) {
	if ( influence_width == 0 ) return;
	--influence_width; //1 means now same residue
	tr.Info << " skip redundant constraints... starting with " << in.size() << " constraints" << std::endl;
	using namespace scoring::constraints;
	utility::vector1< ConstraintCOP > out;
	ObjexxFCL::FArray2D_int count_matrix( total_residue, total_residue, 0 );
	using namespace numeric::random;
	random_permutation( in, numeric::random::rg() ); //I don't think a single pass of pairwise exchanges is enough to randomize the vector.
	random_permutation( in, numeric::random::rg() );
	random_permutation( in, numeric::random::rg() );

	for ( utility::vector1< ConstraintCOP >::const_iterator it = in.begin(); it != in.end(); ++it ) {
		AmbiguousNMRConstraintCOP cst_in_casted;
		cst_in_casted = utility::pointer::dynamic_pointer_cast< AmbiguousNMRConstraint const >( *it );
		if ( cst_in_casted ) {
			tr.Debug << "casted to AmbiguousNMRConstraint: " << std::endl;
			for ( utility::vector1< ConstraintCOP >::const_iterator multi_it = cst_in_casted->member_constraints().begin(); multi_it != cst_in_casted->member_constraints().end(); ++multi_it ) {
				count_constraint( *multi_it, cst_in_casted->member_constraints().size() > 1 , count_matrix, influence_width, total_residue );
			}
		} else {
			count_constraint( *it, false, count_matrix, influence_width, total_residue );
		}
	}

	if ( tr.Trace.visible() ) {
		for ( Size i=1; i<=total_residue; ++i ) {
			for ( Size j=i+1; j<=total_residue; ++j ) {
				tr.Trace << i << " " << j << "   " << count_matrix( i, j ) << std::endl;
			}
		}
	}

	for ( utility::vector1< ConstraintCOP >::const_iterator it = in.begin(); it != in.end(); ++it ) {

		AmbiguousNMRConstraintCOP cst_in_casted;
		cst_in_casted = utility::pointer::dynamic_pointer_cast< AmbiguousNMRConstraint const >( *it );
		bool keep( false );
		if ( cst_in_casted ) {
			for ( utility::vector1< ConstraintCOP >::const_iterator multi_it = cst_in_casted->member_constraints().begin(); multi_it != cst_in_casted->member_constraints().end(); ++multi_it ) {
				keep |= keep_constraint( *multi_it,  cst_in_casted->member_constraints().size() > 1, count_matrix, influence_width, total_residue );
			}
		} else {
			keep = keep_constraint( *it, false, count_matrix, influence_width, total_residue );
		}

		if ( keep ) {
			out.push_back( *it );
		}
	}

	in = out;

	if ( tr.Trace.visible() ) {
		for ( Size i=1; i<=total_residue; ++i ) {
			for ( Size j=i+1; j<=total_residue; ++j ) {
				tr.Trace << i << " " << j << "   " << count_matrix( i, j ) << std::endl;
			}
		}
	}
	tr.Info << "remaining non-redundant constraints " << in.size() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////
void drop_constraints( ConstraintCOPs& in, core::Real drop_rate ) {
	utility::vector1< ConstraintCOP > out;
	for ( utility::vector1< ConstraintCOP >::const_iterator it = in.begin(); it != in.end(); ++it ) {
		if ( numeric::random::rg().uniform() >= drop_rate ) {
			out.push_back( *it );
		}
	}
	in = out;
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// @details example of how to go through a pose constraint set and print out stuff.
void
print_atom_pair_constraints( pose::Pose const & pose, std::ostream & out /* = std::cout */ ){
	ConstraintSetCOP cst_set = pose.constraint_set();
	typedef ResidueConstraints::const_iterator ResiduePairConstraintsIterator;
	// should probably do intra-residue too. Oh well.
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		for ( ResiduePairConstraintsIterator
				iter = cst_set->residue_pair_constraints_begin( n ),
				iter_end = cst_set->residue_pair_constraints_end( n );
				iter  != iter_end; ++iter ) {
			ConstraintsOP csts = iter->second;
			for ( ConstraintCOPs::const_iterator it=csts->begin(), ite = csts->end();
					it != ite;
					++it ) {
				Constraint const & cst( **it );
				if ( cst.type() == "AtomPair" ) cst.show_def( out, pose );
			}
		}
	}
}

} // namespace constraints
} // namespace scoring
} // namespace core
