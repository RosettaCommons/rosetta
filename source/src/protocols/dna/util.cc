// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.cc
/// @author ashworth

#include <protocols/dna/util.hh>
#include <protocols/dna/DnaChains.hh>
#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/DNAParameters.hh>

#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/graph/Graph.hh>
#include <basic/options/option.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh> // file_exists, create_directory
#include <utility/vector1.hh>
using utility::vector1;
#include <utility/vector0.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
using utility::string_split;

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
typedef numeric::xyzVector< core::Real > xyzVec;


#include <algorithm> // std::min
#include <iostream>
#include <sstream>

// option key includes

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>


//#include <fstream> no, use zstreams

namespace protocols {
namespace dna {

using namespace core;
using namespace conformation;
using namespace chemical;
using namespace basic::options;
using namespace pack;
using namespace rotamer_set;
using namespace scoring;
using namespace ObjexxFCL::format;

static THREAD_LOCAL basic::Tracer TR( "protocols.dna.util", basic::t_info );

typedef utility::vector1< std::string > Strings;

/// @details checks c-beta (except glycine) to base atom-atom distances, not including ribose or phosphate backbone.
/// @author ashworth
bool
close_to_dna(
	Residue const & pres,
	Residue const & dres,
	Real threshold,
	bool base_only /* = false */
)
{
	// TR << "pres " << pres << " dres " << dres << std::endl;
	// iterate over dna base ('sidechain') atoms, check for distance to protein sidechain takeoff point
	Atoms::const_iterator baseatom = ( base_only ? dres.sidechainAtoms_begin() : dres.atom_begin() );
	for ( Atoms::const_iterator end( dres.heavyAtoms_end() ); baseatom != end; ++baseatom ) {
		if ( baseatom->xyz().distance_squared( pres.nbr_atom_xyz() ) < threshold ) return true;
	}
	return false;
}

/// @details arginine rotamer sweep at a protein residue to see if it should be considered a (potentially) 'dna-contacting' residue
/// @author ashworth
Real
argrot_dna_dis2(
	pose::Pose const & pose,
	Size presid,
	Residue const & pres,
	Residue const & dres,
	Real threshold,
	bool base_only /* = false */
)
{
	using namespace pack;
	using namespace scoring;
	using namespace task;

	// TR << "Arg rot screen for " << pres << " " << presid << " vs " << dres << std::endl;

	PackerTaskOP ptask( TaskFactory::create_packer_task( pose ) );
	ptask->set_bump_check( false );
	ptask->temporarily_set_pack_residue( presid, true );

	// use ex1 rotamers
	ResidueLevelTask & restask( ptask->nonconst_residue_task( presid ) );
	restask.or_ex1( true );

	// restrict to arginine (not currently necessary if calling build...concrete directly
	vector1< bool > keep_aas( num_canonical_aas, false );
	keep_aas[ aa_arg ] = true;
	restask.restrict_absent_canonical_aas( keep_aas );

	// mostly if not completely irrelevant here, but required as an argument for building rotamers
	std::string weights_tag("dna");
	if ( option[ OptionKeys::score::weights ].user() ) {
		weights_tag = option[ OptionKeys::score::weights ]();
	}
	ScoreFunctionOP scrfxn( ScoreFunctionFactory::create_score_function( weights_tag ) );
	// unnecessary here, yet also required
	graph::GraphOP dummygraph( new graph::Graph( pose.total_residue() ) );

	RotamerSetFactory rsf;
	RotamerSetOP rotset( rsf.create_rotamer_set( pres ) );
	rotset->set_resid( presid );
	rotset->build_rotamers( pose, *scrfxn, *ptask, dummygraph, false );

	// TR(basic::t_debug) << "arg screen w/ " << rotset->num_rotamers() << " rots" << std::endl;

	// add a bump check here first

	Real shortest_dis2(10000), dis2;

	for ( Rotamers::const_iterator rotamer( rotset->begin() ); rotamer != rotset->end(); ++rotamer ) {
		if ( (*rotamer)->aa() != aa_arg ) {
			// for packer safety, RotamerSet will add in a native rotamer if it didn't actually build any rotamers
			if ( rotset->num_rotamers() == 1 ) continue;
			TR << "warning non-arg rotamer " << (*rotamer)->aa() << std::endl;
			runtime_assert( false );
		}

		Atoms::const_iterator prot_begin( (*rotamer)->sidechainAtoms_begin() ),
			prot_end( (*rotamer)->heavyAtoms_end() ),
			dna_end( dres.heavyAtoms_end() );
		Atoms::const_iterator dna_begin =
			( base_only ? dres.sidechainAtoms_begin() : dres.atom_begin() );

		dis2 = contact_distance2( prot_begin, prot_end, dna_begin, dna_end, threshold );
		if ( dis2 < shortest_dis2 ) shortest_dis2 = dis2;
		if ( shortest_dis2 < threshold ) return shortest_dis2;
	}
	return shortest_dis2;
}

/// @details distance check for contact between two sets of atoms
/// @author ashworth
Real
contact_distance2(
	Atoms::const_iterator a_begin,
	Atoms::const_iterator a_end,
	Atoms::const_iterator b_begin,
	Atoms::const_iterator b_end,
	Real threshold // default is 0.0
)
{
	Real shortest_dis2(10000), dis2;

	for ( Atoms::const_iterator atm_a( a_begin ); atm_a != a_end; ++atm_a ) {
		for ( Atoms::const_iterator atm_b( b_begin ); atm_b != b_end; ++atm_b ) {

			dis2 = atm_a->xyz().distance_squared( atm_b->xyz() );
			if ( dis2 < shortest_dis2 ) shortest_dis2 = dis2;
			if ( shortest_dis2 < threshold ) return shortest_dis2; // early exit mode
		}
	}
	return shortest_dis2;
}

/// @details A sanity check for the arginine rotamer screen. Can prevent the design of positions that are best left alone because they are too far away along the helical axis ('laterally').
/// @author ashworth
Real z_axis_dist(
	Residue const & pres,
	Residue const & dres
)
{
	using namespace scoring::dna;

	xyzVec const Z( get_z_axis( dres, get_y_axis(dres,1) ) ),
		prot( pres.nbr_atom_xyz() ),
		dna( dres.xyz( dres.first_sidechain_atom() ) );

	// vector from first protein sidechain atom to DNA N1/N9
	xyzVec vec( prot - dna );
	// return scalar projection onto DNA helical axis
	return std::abs( dot(vec,Z) );
}

/// @brief also consider using the dna_base_partner function below
/// @author ashworth
std::string dna_comp_name_str( std::string const & dna ) {
	if ( dna == "ADE" ) return "THY";
	if ( dna == "CYT" ) return "GUA";
	if ( dna == "GUA" ) return "CYT";
	if ( dna == "THY" ) return "ADE";
	if ( dna == " DA" ) return " DT";
	if ( dna == " DC" ) return " DG";
	if ( dna == " DG" ) return " DC";
	if ( dna == " DT" ) return " DA";
	utility_exit_with_message( "Bad DNA name " + dna );
	return "NONE";
}

/// @brief intended to convert any DNA "threeletter code" into the full three-letter code. Note that this does not (necessarily) return the same thing as residue_type::name3 (which returns "  N" format as of Dec 2008)
std::string dna_full_name3( std::string const & name3 )
{
	if ( name3 == " DA" || name3 == "ADE" ) return "ADE";
	if ( name3 == " DC" || name3 == "CYT" ) return "CYT";
	if ( name3 == " DG" || name3 == "GUA" ) return "GUA";
	if ( name3 == " DT" || name3 == "THY" ) return "THY";
	if ( name3 == "  A" ) return "RAD";
	if ( name3 == "  C" ) return "RCY";
	if ( name3 == "  G" ) return "RGU";
	if ( name3 == "  U" ) return "URA";
	return name3;
}

/// helper function
chemical::AA
dna_base_partner( chemical::AA const & na )
{
	using namespace chemical;

	switch( na ) {
	case na_ade :
		return na_thy;
	case na_thy :
		return na_ade;
	case na_gua :
		return na_cyt;
	case na_cyt :
		return na_gua;
	default :
		utility_exit_with_message( "Bad DNA aa "+chemical::name_from_aa(na) );
	}
	// return na_ade;
	return aa_unk;
}

/// @details DnaChains version, adapted from pbradley's code.  More paranoid geometry checks, in order to allow highly distorted basepairs without making mistakes
/// @author ashworth
void
find_basepairs(
	pose::Pose const & pose,
	DnaChains & dna_chains,
	bool include_unpaired // defaults to true
)
{
	using namespace scoring::dna;

	TR << "\nFinding basepairs:\n";

	Real const max_d( 4.0 );
	Size const nres( pose.total_residue() );

	dna_chains.clear();
	runtime_assert( dna_chains.empty() );

	std::map< AA, AA > base_partner;
	base_partner[ na_ade ] = na_thy;
	base_partner[ na_thy ] = na_ade;
	base_partner[ na_gua ] = na_cyt;
	base_partner[ na_cyt ] = na_gua;

	std::map< AA, std::string > hbond_atom;
	hbond_atom[ na_ade ] = "N1";
	hbond_atom[ na_thy ] = "N3";
	hbond_atom[ na_gua ] = "N1";
	hbond_atom[ na_cyt ] = "N3";

	//ja RNA support
	base_partner[ na_rad ] = na_ura;
	base_partner[ na_ura ] = na_rad;
	base_partner[ na_rgu ] = na_rcy;
	base_partner[ na_rcy ] = na_rgu;

	hbond_atom[ na_rad ] = "N1";
	hbond_atom[ na_ura ] = "N3";
	hbond_atom[ na_rgu ] = "N1";
	hbond_atom[ na_rcy ] = "N3";

	//std::map< Size, Size > partner; // temporary

	for ( Size i(1); i <= nres; ++i ) {
		Residue const & i_rsd( pose.residue(i) );
		AA const & i_aa( i_rsd.aa() );
		if ( !i_rsd.is_DNA() ) continue;
		// the following is false for already-added bottom strand residue indices
		if ( dna_chains.contains(i) ) continue;

		// hbond atom, base y-axis, base z-axis
		xyzVec const hbatm_xyz_i( i_rsd.xyz( hbond_atom[ i_aa ] ) ),
			base_yaxis_i( get_y_axis( i_rsd, 1 /*strand*/ ) );
		xyzVec const base_zaxis_i( get_z_axis( i_rsd, base_yaxis_i ) );

		bool paired( false );
		// check for a basepairing partner
		Real bestdotsum(0.);
		Size best_j(0);
		for ( Size j(i+1); j <= nres; ++j ) {

			Residue const & j_rsd( pose.residue(j) );
			AA const & j_aa( j_rsd.aa() );
			if ( !j_rsd.is_DNA() ) continue;
			// no Watson-Crick check here on purpose

			// -distance check-
			xyzVec const hbatm_xyz_j( j_rsd.xyz( hbond_atom[ j_aa ] ) );

			Real d( hbatm_xyz_i.distance( hbatm_xyz_j ) );
			if ( d >= max_d ) continue;

			// -geometry check-
			xyzVec const base_yaxis_j( get_y_axis( j_rsd, 2 ) );
			xyzVec const base_zaxis_j( get_z_axis( j_rsd, base_yaxis_j ) ),
				hb_vec( ( hbatm_xyz_i - hbatm_xyz_j ).normalized() );

			Real const
				// base y-axes parallel?
				ydot( std::abs( dot( base_yaxis_i, base_yaxis_j ))),
				// hbond vector parallel to base y axis?
				dothbyi( std::abs( dot( base_yaxis_i, hb_vec ))),
				dothbyj( std::abs( dot( base_yaxis_j, hb_vec ))),
				// hbond vector perpendicular to base z axis?
				dothbzi( std::abs( dot( base_zaxis_i, hb_vec ))),
				dothbzj( std::abs( dot( base_zaxis_j, hb_vec )));

			Real const dotsum( 2*ydot + dothbyi + dothbyj - dothbzi - dothbzj );
			int pdbi(i), pdbj(j);
			if ( pose.pdb_info() ) {
				pdbi = pose.pdb_info()->number(i);
				pdbj = pose.pdb_info()->number(j);
			}
			int verbosity(0); // to do: learn to use Tracer properly
			if ( verbosity >= 2 ) {
				TR << "basepair geom "
					<< pdbi << " vs. " << pdbj << " dis " << d
					<< " ydot " << ydot
					<< " hbydots " << dothbyi << " " << dothbyj
					<< " hbzdots " << dothbzi << " " << dothbzj;
			}

			if ( dotsum < bestdotsum || ydot < 0.8 ||
					dothbyi < 0.6 || dothbyj < 0.6 ||
					dothbzi > 0.5 || dothbzj > 0.5 ) {
				if ( verbosity >= 2 ) TR << '\n';
				continue;
			}
			if ( verbosity >= 2 ) TR << " acceptable" << '\n';

			// -complementarity check-
			if ( j_aa != base_partner[ i_aa ] ) {
				std::cerr << "Warning: nucleic acids " << i_rsd.name3() << " " <<
					pdbi << " and " << j_rsd.name3() << " " <<
					pdbj << " have basepaired geometry, but are not " <<
					"complementary types" << '\n';
				continue; // skip non-canonical basepairs for now
			}
			// passed: save j as optimal basepair partner
			bestdotsum = dotsum;
			best_j = j;
		} // end j loop
		if ( bestdotsum != 0. ) {
			dna_chains[ i ] = DnaPosition( i, best_j );
			paired = true;
		}
		if ( paired || !include_unpaired ) continue;
		// include unpaired dna position
		dna_chains[ i ] = DnaPosition( i );
	} // end i loop
	dna_chains.print( pose, TR );
}

/// @details populates a set of all possible sequence combinations over a given range of positions. recursive.
/// @author ashworth
void
make_sequence_combinations(
	utility::vector1< Size >::const_iterator seqset_iter,
	utility::vector1< Size > const & seq_indices,
	task::PackerTaskCOP ptask,
	ResTypeSequence & sequence,
	ResTypeSequences & sequences
)
{
	using namespace task;

	Size resid( *seqset_iter );
	ResidueLevelTask const & restask( ptask->residue_task( resid ) );

	for ( ResidueLevelTask::ResidueTypeCOPListConstIter type( restask.allowed_residue_types_begin() );
			type != restask.allowed_residue_types_end(); ++type ) {
		// ignore adduct variant types for now (probably hydrated)
		if ( (*type)->has_variant_type( chemical::ADDUCT_VARIANT ) ) continue;
		sequence[ resid ] = *type;

		if ( seqset_iter == seq_indices.end() - 1 ) sequences.push_back( sequence );
		else make_sequence_combinations( seqset_iter+1, seq_indices, ptask, sequence, sequences );
	}
}

/// @brief make a list of all single mutants from a base sequence
/// @author ashworth
void
make_single_mutants(
	ResTypeSequence const & sequence,
	task::PackerTaskCOP ptask,
	ResTypeSequences & sequences
)
{
	using namespace task;
	for ( ResTypeSequence::const_iterator it( sequence.begin() ); it != sequence.end(); ++it ) {
		Size index( it->first );
		ResidueLevelTask const & rtask( ptask->residue_task(index) );
		for ( ResidueLevelTask::ResidueTypeCOPListConstIter type( rtask.allowed_residue_types_begin() );
				type != rtask.allowed_residue_types_end(); ++type ) {
			// ignore adduct variant types for now (probably hydrated)
			if ( (*type)->has_variant_type( chemical::ADDUCT_VARIANT ) ) continue;
			if ( (*type)->aa() == it->second->aa() ) continue; // avoid duplicating input sequence
			ResTypeSequence mutant( sequence );
			mutant[ index ] = *type;
			sequences.push_back( mutant );
		}
	}
}

void
design_residues_list(
	std::list< PositionType > & design_residues,
	pose::Pose const & pose,
	task::PackerTask const & ptask
)
{
	Size nres( pose.total_residue() );
	for ( Size index(1); index <= nres; ++index ) {
		if ( pose.residue_type(index).is_DNA() ) {
			if ( !ptask.residue_task(index).has_behavior("TARGET") &&
					!ptask.residue_task(index).has_behavior("SCAN") &&
					!ptask.residue_task(index).being_designed() ) continue;
		} else if ( !ptask.pack_residue( index ) ) continue;
		design_residues.push_back(
			PositionType( index, pose.residue_type( index ).get_self_ptr(), ptask.design_residue( index ) ) );
	}
}

// (relevant typdefs are in fwd.hh)
std::ostream & operator << ( std::ostream & os, ResTypeSequence const & seq )
{
	for ( ResTypeSequence::const_iterator pos( seq.begin() ); pos != seq.end(); ++pos ) {
		if ( pos != seq.begin() ) os << ", ";
		os << pos->first << "-" << pos->second->name1();
	}
	return os;
}

std::string seq_to_str( ResTypeSequence const & seq ) {
	std::string str;
	for ( ResTypeSequence::const_iterator pos( seq.begin() ); pos != seq.end(); ++pos ) {
		str += pos->second->name1();
	}
	return str;
}

std::ostream & operator << ( std::ostream & os, ResTypeSequences const & seqs )
{
	for ( ResTypeSequences::const_iterator seq( seqs.begin() ); seq != seqs.end(); ++seq ) {
		os << *seq << '\n';
	}
	return os;
}

std::string seq_pdb_str(
	ResTypeSequence const & seq,
	pose::Pose const & pose
)
{
	std::ostringstream os;
	for ( ResTypeSequence::const_iterator pos( seq.begin() ); pos != seq.end(); ++pos ) {
		Size const index( pos->first );
		if ( index < 1 || index > pose.total_residue() ) {
			assert(false);
			continue;
		}
		if ( pos != seq.begin() ) os << ",";
		if ( pose.pdb_info() ) {
			os << pose.pdb_info()->chain( index ) << "." << pose.pdb_info()->number( index );
		} else {
			os << pose.chain( index ) << "." << index;
		}
		os << "." << dna_full_name3( pos->second->name3() );
	}
	return os.str();
}

void print_sequence_pdb_nums(
	ResTypeSequence const & seq,
	pose::Pose const & pose,
	std::ostream & os
)
{
	os << seq_pdb_str( seq, pose ) << '\n';
}

void print_sequences_pdb_nums(
	ResTypeSequences const & seqs,
	pose::Pose const & pose,
	std::ostream & os
)
{
	for ( ResTypeSequences::const_iterator seq( seqs.begin() ); seq != seqs.end(); ++seq ) {
		print_sequence_pdb_nums( *seq, pose, os );
	}
}

/// @details for packing a single DNA sequence out of a multi-DNA-sequence RotamerSet
/// @author ashworth
void
restrict_dna_rotamers(
	RotamerSetsCOP rotamer_sets,
	ResTypeSequence const & seq,
	utility::vector0<int> & rot_to_pack
)
{
	rot_to_pack.clear();
	Size const nrot( rotamer_sets->nrotamers() );
	for ( Size roti(1); roti <= nrot; ++roti ) {

		Size const rotpos( rotamer_sets->res_for_rotamer(roti) );
		ResidueType const & rot_type( rotamer_sets->rotamer(roti)->type() );

		ResTypeSequence::const_iterator seqindex( seq.find( rotpos ) );
		if ( seqindex != seq.end() ) {
			// compare only the name3's on order to allow variants
			std::string seq_typename( (seqindex->second)->name3() ),
				rot_typename( rot_type.name3() );
			if ( seq_typename != rot_typename ) continue;
		}
		rot_to_pack.push_back( roti );
	}
	Size const rots_off( nrot - rot_to_pack.size() );
	TR << "Fixing DNA rotamers: " << rots_off
		<< " out of " << nrot << " rotamers disabled." << std::endl;
}

/// @details for packing a single sequence out of a RotamerSets that (potentially) represents sequence variability
/// @author ashworth
void
restrict_to_single_sequence(
	rotamer_set::RotamerSetsCOP rotamer_sets,
	vector1< ResidueTypeCOP > const & single_sequence,
	utility::vector0< int > & rot_to_pack
)
{
	rot_to_pack.clear();
	Size const nrot( rotamer_sets->nrotamers() );
	for ( Size roti(1); roti <= nrot; ++roti ) {
		Size const rotpos( rotamer_sets->res_for_rotamer(roti) );
		ResidueType const & rot_type( rotamer_sets->rotamer(roti)->type() );
		// a comparison operator is not defined for the ResidueType class
		// compare names here
		// name3 comparison should allow variants
		std::string seq_typename( ( single_sequence[ rotpos ] )->name3() ),
			rot_typename( rot_type.name3() );
		if ( seq_typename != rot_typename ) continue;
		rot_to_pack.push_back( roti );
	}
	Size const rots_off( nrot - rot_to_pack.size() );
	TR << "Fixing rotamers for a single sequence: " << rots_off
		<< " out of " << nrot << " rotamers disabled." << std::endl;
}

/// @details
/// @author ashworth
void
substitute_residue(
	pose::Pose & pose,
	Size index,
	ResidueType const & new_type
)
{
	Residue const & existing( pose.residue( index ) );
	ResidueOP new_res( ResidueFactory::create_residue( new_type, existing, pose.conformation() ) );
	new_res->set_chi( 1, existing.chi(1) );
	pose.replace_residue( index, *new_res, false );
}

// @brief
// @author ashworth
void
write_checkpoint( pose::Pose & pose, Size iter )
{
	if ( ! option[ OptionKeys::dna::design::checkpoint ].user() ) return;
	std::string fileroot( option[ OptionKeys::dna::design::checkpoint ]() );

	TR << "writing dna mode checkpoint files..." << '\n';

	// write out current Pose
	std::string pdbname( fileroot + ".pdb.checkpoint" );
	utility::io::ozstream pdbout( pdbname.c_str() );
	pose.dump_pdb( pdbout );
	pdbout.close();

	// write checkpoint file
	std::string checkpointname( fileroot + ".checkpoint" );
	utility::io::ozstream out( checkpointname.c_str() );
	// std::ofstream out( checkpointname.c_str() );

	if ( !out ) {
		std::cerr << "trouble opening file " << checkpointname
			<< " for writing... skipping checkpoint" << std::endl;
		runtime_assert( false ); // die here in debug mode
		return;
	}

	// here iter should refer to the last complete iteration
	out << "Iteration " << iter << '\n' << pdbname << '\n';
	out.close();

	TR << "wrote " << pdbname << ", " << checkpointname << std::endl;
}

// @brief
// @author ashworth
void
load_checkpoint( pose::Pose & pose, Size & iter )
{
	if ( ! option[ OptionKeys::dna::design::checkpoint ].user() ) return;
	std::string fileroot( option[ OptionKeys::dna::design::checkpoint ]() );

	utility::io::izstream file;
	std::string filename( fileroot + ".checkpoint" );
	file.open( filename.c_str() );
	if ( !file ) return;

	TR << "Reading DNA design checkpoint info from " << filename << '\n';

	std::string /*line, */word, pdbfile;
	// get iteration
	Size last_iter;
	file >> word >> last_iter >> skip; // first line
	if ( ( word != "Iteration" ) ) return;
	file >> pdbfile >> skip;
	file.close();

	if ( option[ OptionKeys::out::pdb_gz ]() ) pdbfile += ".gz";
	pose::Pose temp_pose;
	core::import_pose::pose_from_pdb( temp_pose, filename );

	pose = temp_pose;
	// here iter should refer to the last complete iteration
	iter = last_iter + 1;

	TR << "loaded " << pdbfile << " for iteration " << iter << std::endl;
}

// @brief make sure that old checkpoint files will not be accidentally reused
// @author ashworth
void
checkpoint_cleanup()
{
	if ( ! option[ OptionKeys::dna::design::checkpoint ].user() ) return;
	std::string fileroot( option[ OptionKeys::dna::design::checkpoint ]() );

	std::list< std::string > filenames;
	filenames.push_back( fileroot + ".checkpoint" );
	filenames.push_back( fileroot + ".pdb.checkpoint" );

	for ( std::list< std::string >::const_iterator filename( filenames.begin() );
			filename != filenames.end(); ++filename ) {
		if ( utility::file::file_exists( *filename ) ) {
			std::string nameold( *filename + ".old" );
			std::rename( (*filename).c_str(), nameold.c_str() );
		}
	}
}

/// @brief loads command-line dna design definitions (shorthand alternative to using resfile)
/// option value is string vector
///   i.e. -dna_defs C.-6 C.-5
///   or   -dna_defs C.-6.GUA C.-5.CYT
/// @author ashworth
void
load_dna_design_defs_from_strings(
	DnaDesignDefOPs & defs,
	Strings const & str_defs
)
{
	for ( Strings::const_iterator str_def( str_defs.begin() ), end( str_defs.end() );
			str_def != end; ++str_def ) {
		defs.push_back( DnaDesignDefOP( new DnaDesignDef( *str_def ) ) );
	}
}

void
load_dna_design_defs_from_file(
	DnaDesignDefOPs & defs,
	std::string const & filename,
	std::string const & pdb_prefix
)
{
	std::string stripped_prefix( pdb_prefix );

	TR << "Getting dna_defs from file " << filename;
	if ( ! stripped_prefix.empty() ) {
		stripped_prefix = string_split( stripped_prefix, '/' ).back();
		TR << " for " << stripped_prefix;
	}
	TR << '\n';

	utility::io::izstream defs_file( filename.c_str() );
	std::string line;
	while ( getline( defs_file, line ) ) {
		utility::vector1< std::string > words( string_split( line, ' ' ) );
		// multiple pdbs may be specified in this file:
		// only match lines beginning with stripped_prefix, if specified
		if ( ! stripped_prefix.empty() && words.front() != stripped_prefix ) continue;
		Strings str_defs;
		str_defs.insert( str_defs.begin(), words.begin()+1, words.end() );
		load_dna_design_defs_from_strings( defs, str_defs );
	}
}

void
load_dna_design_defs_from_options(
	DnaDesignDefOPs & defs,
	std::string pdb_prefix /* = std::string() */
)
{
	if ( option[ OptionKeys::dna::design::dna_defs ].user() ) {
		// list of defs for a single pdb
		Strings str_defs( option[ OptionKeys::dna::design::dna_defs ]().vector() );
		load_dna_design_defs_from_strings( defs, str_defs );
	} else if ( option[ OptionKeys::dna::design::dna_defs_file ].user() ) {
		// file containing lists of defs, with format 'pdbcode def def def'
		load_dna_design_defs_from_file(
			defs,
			option[ OptionKeys::dna::design::dna_defs_file ](),
			pdb_prefix
		);
	}
	TR.flush();
}

void
add_constraints_from_file(
	pose::Pose & pose
)
{
	using namespace scoring::constraints;

	std::string cst_file;
	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		cst_file = option[ OptionKeys::constraints::cst_file ]().front();
	} else return;

	ConstraintSetOP cst_set =
		ConstraintIO::get_instance()->read_constraints_new( cst_file, ConstraintSetOP( new ConstraintSet ), pose );

	pose.constraint_set( cst_set );
}

kinematics::FoldTree
make_base_pair_aware_fold_tree ( pose::Pose const & pose )
{

	Size const nres( pose.total_residue() );

	pose::PDBInfoCOP pdb_data( pose.pdb_info() );

	// Identify DNA duplexed regions
	protocols::dna::DNAParameters dna_info( pose );
	// dna_info.calculate( pose );

	Size num_chains( 1 );
	utility::vector1< Size > chain_start;
	utility::vector1< Size > chain_end;
	utility::vector1< Size > chain_type;

	chain_start.push_back( 1 );
	for ( Size resid = 1 ; resid < nres ; ++resid ) {
		if ( pdb_data->chain( resid ) != pdb_data->chain( resid + 1 ) ) {
			chain_end.push_back( resid );
			chain_start.push_back( resid + 1 );
			num_chains++;
		}
	}
	chain_end.push_back( nres );

	// Allocate the FArrays to call FoldTree::tree_from_cuts_and_jumps
	Size num_cuts( num_chains - 1 );
	ObjexxFCL::FArray1D_int cut_positions( num_cuts, 0 );
	ObjexxFCL::FArray2D_int jump_pairs( 2, num_cuts, 0 );
	Size jump_pair_count( 1 );

	// We can fill the cut info now
	for ( Size cut_num = 1 ; cut_num < chain_end.size() ; ++cut_num ) {
		cut_positions( cut_num ) = chain_end[ cut_num ];
	}

	Size const amino( 1 );
	Size const bped_dna( 2 );
	Size const non_bped_dna( 3 );

	utility::vector1< Size > protein_root( num_chains, 0 );
	utility::vector1< Size > closest_base( num_chains, 0 );

	// Analyze each chain
	for ( Size this_chain = 1 ; this_chain <= num_chains ; ++this_chain ) {
		TR << "Working on chain " << this_chain << std::endl;

		// Check for amino acid chain
		if ( pose.residue( chain_start[ this_chain ] ).is_protein() ) {
			chain_type.push_back( amino );
			TR << "Found 1 initial segments for chain " << this_chain << std::endl;
			TR << "Chain " << this_chain << " segment 1 start res " << chain_start[ this_chain ] <<
				" end res " << chain_end[ this_chain ] << " of type 1" << std::endl;


			// Find the DNA base with the closest C1' atom to some Calpha in this protein
			Real best_dist( 9999.0 );
			for ( Size prot_res = chain_start[ this_chain ] ; prot_res <= chain_end[ this_chain ] ; ++prot_res ) {
				for ( Size dna_res = 1 ; dna_res <= nres ; ++dna_res ) {
					if ( !pose.residue( dna_res ).is_DNA() ) continue;
					Real check_dist = pose.residue( prot_res ).xyz( "CA" ).distance_squared( pose.residue( dna_res ).xyz( "C1'" ) );
					if ( check_dist < best_dist ) {
						best_dist = check_dist;
						protein_root[ this_chain ] = prot_res;
						closest_base[ this_chain ] = dna_res;
					}
				}
			}

			TR << "Protein closest approach is res " << protein_root[ this_chain ]<< " with base " << closest_base[ this_chain ] << " with distance " << std::sqrt( best_dist ) << std::endl;

			// Note this pair as a jump
			if ( protein_root[this_chain] < closest_base[this_chain] ) {
				jump_pairs( 1, jump_pair_count ) = protein_root[this_chain];
				jump_pairs( 2, jump_pair_count ) = closest_base[this_chain];
			} else {
				jump_pairs( 2, jump_pair_count ) = protein_root[this_chain];
				jump_pairs( 1, jump_pair_count ) = closest_base[this_chain];
			}
			jump_pair_count++;

			continue;
		}

		// Bail if it's something other than protein or DNA
		if ( !pose.residue( chain_start[ this_chain ] ).is_DNA() ) {
			std::cerr << "Bad call to make_basepair_aware_fold_tree() with non-protein, non-DNA type" << std::endl;
			utility_exit_with_message( "make_base_aware_fold_tree() takes only protein, DNA!" );
		}

		chain_type.push_back( bped_dna );

		// Break up this DNA into segments

		Size num_segments( 1 );
		utility::vector1< Size > segment_start;
		utility::vector1< Size > segment_end;
		utility::vector1< Size > segment_type;

		segment_start.push_back( chain_start[ this_chain ] );
		if ( dna_info.find_partner( chain_start[ this_chain ] ) != 0 ) {
			segment_type.push_back( bped_dna );
		} else {
			segment_type.push_back( non_bped_dna );
		}
		for ( Size resid = chain_start[this_chain]  ; resid < chain_end[ this_chain ] ; ++resid ) {
			// Check for difference
			bool this_bped( dna_info.find_partner( resid ) != 0 );
			bool next_bped( dna_info.find_partner( resid + 1 ) != 0 );
			// Switch from base paired to not base paired
			if ( this_bped && !next_bped ) { // Switch from base paired to not base paired
				segment_end.push_back( resid );
				segment_start.push_back( resid + 1 );
				segment_type.push_back( non_bped_dna );
			} else if ( !this_bped && next_bped ) { // Switch from not base paired to base paired
				segment_end.push_back( resid );
				segment_start.push_back( resid + 1 );
				segment_type.push_back( bped_dna );
			} else if ( !this_bped && !next_bped ) { // No switch - do nothing
				continue;
			} else if ( pdb_data->chain( dna_info.find_partner( resid ) ) !=
					pdb_data->chain( dna_info.find_partner( resid + 1 ) ) ) { // Both base-paired, but to different strands
				segment_end.push_back( resid );
				segment_start.push_back( resid + 1 );
				segment_type.push_back( bped_dna );
			}
		}
		segment_end.push_back( chain_end[ this_chain ] );

		num_segments = segment_start.size();

		// Let's see what we have
		TR << "Found " << num_segments << " initial segments for chain " << this_chain << std::endl;
		for ( Size this_segment = 1 ; this_segment <= num_segments ; ++this_segment ) {
			TR << "Chain " << this_chain << " segment " << this_segment << " start res " << segment_start[ this_segment ] <<
				" end res " << segment_end[ this_segment ] << " of type " << segment_type[ this_segment ] << std::endl;
		}


		// Record mid-points of base paired segments
		utility::vector1< Size > bp_middle( num_segments, 0 );

		for ( Size this_segment = 1 ; this_segment <= num_segments ; ++this_segment ) {
			if ( segment_type[ this_segment ] == bped_dna ) {
				bp_middle[ this_segment ] = ( segment_start[ this_segment ] + segment_end[ this_segment ] ) / 2;
			}
		}

		// Merge non-base paired segments into base paired segments if possible

		Size num_processed( num_segments );

		if ( num_segments  > 1 ) {
			for ( Size this_segment = 1 ; this_segment <= num_segments ; ++this_segment ) {
				if ( segment_type[ this_segment ] == non_bped_dna  && num_segments > 1 ) {
					num_processed--;
					if ( this_segment == 1 ) { // Just merge with the next segment
						segment_start[ this_segment + 1 ] = segment_start[ this_segment ];
					} else if ( this_segment == num_segments ) { // Just merge with the previous segment
						segment_end[ this_segment - 1 ] = segment_end[ this_segment ];
					} else {  // Must be between two base paired segments - split evenly
						// Handle single residue segment
						if ( segment_start[ this_segment ] == segment_end[ this_segment ] ) {
							// Just give it to the previous
							segment_end[ this_segment - 1 ] = segment_end[ this_segment ];
						} else {
							// Divide in half
							Size split_pos(  ( segment_start[ this_segment ] + segment_end[ this_segment ] ) / 2 );
							segment_end[ this_segment - 1 ] = split_pos;
							segment_start[ this_segment + 1 ] = split_pos + 1;
						}
					}
				}
			}
		}

		// Need to handle the case of a single segment chain of non-base paired DNA - it needs a jump somewhere
		if ( num_segments == 1 && chain_type[ this_chain ] == non_bped_dna ) {
			// Find the amino acid with the closest Calpha atom to some C1' atom this chain
			Real best_dist( 9999.0 );
			for ( Size dna_res = chain_start[ this_chain ] ; dna_res <= chain_end[ this_chain ] ; ++dna_res ) {
				for ( Size prot_res = 1 ; prot_res <= nres ; ++prot_res ) {
					if ( !pose.residue( prot_res ).is_protein() ) continue;
					Real check_dist = pose.residue( prot_res ).xyz( "CA" ).distance_squared( pose.residue( dna_res ).xyz( "C1'" ) );
					if ( check_dist < best_dist ) {
						best_dist = check_dist;
						protein_root[ this_chain ] = prot_res;
						closest_base[ this_chain ] = dna_res;
					}
				}
			}

			TR << "Unpaired DNA closest approach is res " << closest_base[ this_chain ]<< " with amino acid " << protein_root[ this_chain ] << " with distance " << std::sqrt( best_dist ) << std::endl;

			// Note this pair as a jump
			if ( protein_root[this_chain] < closest_base[this_chain] ) {
				jump_pairs( 1, jump_pair_count ) = protein_root[this_chain];
				jump_pairs( 2, jump_pair_count ) = closest_base[this_chain];
			} else {
				jump_pairs( 2, jump_pair_count ) = protein_root[this_chain];
				jump_pairs( 1, jump_pair_count ) = closest_base[this_chain];
			}
			jump_pair_count++;
		}


		// Let's see what we have
		TR << "Found " << num_processed << " final segments for chain " << this_chain << std::endl;
		Size accum_count( 0 );
		for ( Size this_segment = 1 ; this_segment <= num_segments ; ++this_segment ) {
			if ( segment_type[ this_segment ] == bped_dna ) {
				accum_count++;
				TR << "Chain " << this_chain << " segment " << accum_count << " start res " << segment_start[ this_segment ] <<
					" end res " << segment_end[ this_segment ] << " of type " << segment_type[ this_segment ] << std::endl;

				// retrieve the mid-point base pair and its partner
				Size mid_partner = dna_info.find_partner( bp_middle[ this_segment ] );

				// Store the jump info if this chain is the lower number (to avoid adding twice)
				// Also check to make sure these chains haven't already been connected.  This
				// can happen in two strands are base-paired at the ends but have a non-bp-ed
				// bubble in between.
				if ( ( bp_middle[ this_segment ] < mid_partner ) &&
						not_already_connected( pose, jump_pair_count - 1, pose.pdb_info()->chain( bp_middle[ this_segment] ), pose.pdb_info()->chain( mid_partner ), jump_pairs ) ) {
					TR << "Making jump between " << bp_middle[ this_segment ] << " and " << mid_partner << std::endl;
					jump_pairs( 1, jump_pair_count ) = bp_middle[ this_segment ];
					jump_pairs( 2, jump_pair_count ) = mid_partner;
					jump_pair_count++;
				}
			}
		}
	}

	kinematics::FoldTree ft( nres );

	ft.tree_from_jumps_and_cuts( nres, num_cuts, jump_pairs, cut_positions, 1 );

	return ft;
}

bool
not_already_connected(
	pose::Pose const & pose,
	Size const num_jumps,
	char const this_chain,
	char const other_chain,
	ObjexxFCL::FArray2D_int & jump_pairs
)
{

	for ( Size i = 1 ; i <= num_jumps ; ++i ) {

		// Get chain ids for residues involved in jumps
		char const jump_chain1( pose.pdb_info()->chain( jump_pairs( 1, i ) ) );
		char const jump_chain2( pose.pdb_info()->chain( jump_pairs( 2, i ) ) );

		// Check versus pre-existing jump ( both ways )

		if ( ( jump_chain1 == this_chain && jump_chain2 == other_chain ) ||
				( jump_chain2 == this_chain && jump_chain1 == other_chain ) ) {
			return false;
		}
	}

	return true;
}


void
set_base_segment_chainbreak_constraints(
	pose::Pose & pose,
	core::Size const start_base,
	core::Size const end_base
)
{
	using namespace scoring::constraints;
	using namespace id;
	using numeric::conversions::radians;

	pose::PDBInfoCOP pdb_data( pose.pdb_info() );

	// Size const nres( pose.total_residue() );

	// From Phil
	Real const O3_P_distance( 1.608 );
	Real const O3_angle( 119.8 );
	Real const  P_angle( 103.4 );
	Real const  OP2_angle( 108.23 );

	Real const distance_stddev( 0.3 ); // amber is 0.0659
	Real const angle_stddev_degrees( 35 ); // amber is 8.54 (P angle), 5.73 (O3 angle)

	core::scoring::func::FuncOP const distance_func( new core::scoring::func::HarmonicFunc( O3_P_distance, distance_stddev ) );
	core::scoring::func::FuncOP const O3_angle_func( new core::scoring::func::HarmonicFunc( radians( O3_angle ), radians( angle_stddev_degrees ) ) );
	core::scoring::func::FuncOP const  P_angle_func( new core::scoring::func::HarmonicFunc( radians(  P_angle ), radians( angle_stddev_degrees ) ) );
	core::scoring::func::FuncOP const OP2_angle_func( new core::scoring::func::HarmonicFunc( radians(  OP2_angle ), radians( angle_stddev_degrees ) ) );

	assert( start_base <= end_base );

	// First the start base
	if ( !pose.residue_type( start_base ).is_lower_terminus() ) {
		conformation::Residue const & rsd1( pose.residue( start_base-1 ) );
		conformation::Residue const & rsd2( pose.residue( start_base   ) );

		// Setup constraints to close bb

		AtomID const C3_id( rsd1.atom_index( "C3'" ), start_base - 1 );
		AtomID const O3_id( rsd1.atom_index( "O3'" ), start_base - 1 );
		AtomID const  P_id( rsd2.atom_index( "P"   ), start_base );
		AtomID const O5_id( rsd2.atom_index( "O5'" ), start_base );
		AtomID const OP2_id( rsd2.atom_index( "OP2" ), start_base );

		// distance from O3' to P
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint( O3_id, P_id, distance_func ) ) ) );
		// angle at O3'
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func ) ) ) );
		// angle at P
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) ) ) );
		// another angle at P - try not to get goofy geometries
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) ) ) );
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AngleConstraint( O3_id, P_id, OP2_id,  OP2_angle_func ) ) ) );
	}

	// Next the end base
	if ( !pose.residue_type( end_base ).is_upper_terminus() ) {
		conformation::Residue const & rsd1( pose.residue( end_base   ) );
		conformation::Residue const & rsd2( pose.residue( end_base+1 ) );

		// Setup constraints to close bb

		AtomID const C3_id( rsd1.atom_index( "C3'" ), end_base );
		AtomID const O3_id( rsd1.atom_index( "O3'" ), end_base );
		AtomID const  P_id( rsd2.atom_index( "P"   ), end_base + 1 );
		AtomID const O5_id( rsd2.atom_index( "O5'" ), end_base + 1 );
		AtomID const OP2_id( rsd2.atom_index( "OP2" ), end_base + 1 );

		// distance from O3' to P
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint( O3_id, P_id, distance_func ) ) ) );
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AngleConstraint( O3_id, P_id, OP2_id,  OP2_angle_func ) ) ) );
		// angle at O3'
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func ) ) ) );
		// angle at P
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) ) ) );
	}

}


} // namespace dna
} // namespace protocols
