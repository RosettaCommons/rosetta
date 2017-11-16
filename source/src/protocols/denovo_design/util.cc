// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/util.cc
/// @brief util functions for denovo design of structures
/// @details
/// @author Tom Linsky


//Unit Headers
#include <protocols/denovo_design/util.hh>

//Project Headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>

//Protocol Headers
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

//Core Headers
#include <core/chemical/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/LinearChainbreakEnergy.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

//Boost/ObjexxFCL Headers
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

//C++ Headers

static basic::Tracer TR("protocols.denovo_design.components.util");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Tells whether the two given poses are identical based on # resides and dihedrals
bool same_pose( core::pose::Pose const & pose1, core::pose::Pose const & pose2 )
{
	if ( pose1.size() != pose2.size() ) {
		return false;
	}

	for ( core::Size i = 1; i <= pose1.size(); ++i ) {
		if ( pose1.residue(i).name() != pose2.residue(i).name() ) {
			return false;
		}
		if ( pose1.residue(i).is_protein() && !pose2.residue(i).is_protein() ) {
			return false;
		}
		if ( !pose1.residue(i).is_protein() && pose2.residue(i).is_protein() ) {
			return false;
		}
		if ( !pose1.residue(i).is_protein() && !pose2.residue(i).is_protein() ) {
			continue;
		}

		if ( std::abs(pose1.phi(i) - pose2.phi(i)) > 0.0001 ) {
			return false;
		}
		if ( std::abs(pose1.psi(i) - pose2.psi(i)) > 0.0001 ) {
			return false;
		}
		if ( std::abs(pose1.omega(i) - pose2.omega(i)) > 0.0001 ) {
			return false;
		}
	}
	return true;
}

core::kinematics::FoldTree
remove_all_jump_atoms( core::kinematics::FoldTree const & orig )
{
	core::kinematics::FoldTree ft;
	for ( auto const & e : orig ) {
		if ( e.is_jump() && e.has_atom_info() ) {
			core::kinematics::Edge newedge = e;
			newedge.start_atom() = "";
			newedge.stop_atom() = "";
			ft.add_edge( newedge );
		} else {
			ft.add_edge( e );
		}
	}
	TR.Debug << "Removed all jump atoms, fold tree=" << ft << std::endl;
	return ft;
}


/// @brief removes atoms missing from the current pose
core::kinematics::FoldTree
remove_missing_jump_atoms( core::pose::Pose const & pose, core::kinematics::FoldTree const & orig )
{
	core::kinematics::FoldTree ft;
	for ( auto const & e : orig ) {
		if ( e.is_jump() && e.has_atom_info() ) {
			core::kinematics::Edge newedge = e;
			if ( ( !pose.residue(e.start()).has(e.start_atom()) ) || ( !pose.residue(e.stop()).has(e.stop_atom()) ) ) {
				newedge.start_atom() = "";
				newedge.stop_atom() = "";
			}
			ft.add_edge( newedge );
		} else {
			ft.add_edge( e );
		}
	}
	TR.Debug << "Removed jump atoms, fold tree=" << ft << std::endl;
	return ft;
}

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
/// @details If keep_chirality is true, the D-amino acids are mutated to D-alanine.
void construct_poly_ala_pose(
	core::pose::Pose & pose,
	bool const keep_disulf,
	core::select::residue_selector::ResidueSubset const & set1,
	core::select::residue_selector::ResidueSubset const & set2,
	bool const keep_chirality
) {
	runtime_assert_string_msg( set1.size() == set2.size(), "Error in protocols::denovo_design::util::construct_poly_ala_pose(): The two residue subsets given must be the same size, and they are not." );
	std::set< core::Size > res_set;
	for ( core::Size resid=1; resid<=set1.size(); ++resid ) {
		if ( set1[ resid ] ) res_set.insert( resid );
	}
	for ( core::Size resid=1; resid<=set2.size(); ++resid ) {
		if ( set2[ resid ] ) res_set.insert( resid );
	}
	construct_poly_ala_pose( pose, keep_disulf, res_set, keep_chirality );
}

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
/// @details If keep_chirality is true, the D-amino acids are mutated to D-alanine.
/// @note Updated for beta-3-amino acids by VKM on 26 November 2016.
void construct_poly_ala_pose(
	core::pose::Pose & pose,
	bool const keep_disulf,
	std::set< core::Size > const & res_set,
	bool const keep_chirality
) {
	utility::vector1< core::Size > positions;
	utility::vector1< core::Size > d_positions;
	utility::vector1< core::Size > beta_positions;
	utility::vector1< core::Size > d_beta_positions;

	// remove jump atoms, which may cause problems in mutating residues
	pose.fold_tree( remove_all_jump_atoms( pose.fold_tree() ) );

	for ( core::Size resid : res_set ) {
		if ( !pose.residue( resid ).is_protein() ) continue;
		if ( !pose.residue_type( resid ).is_alpha_aa() && !pose.residue_type( resid ).is_beta_aa() ) continue;

		if ( pose.residue_type(resid).is_alpha_aa() ) {
			if ( !keep_chirality || !pose.residue( resid ).type().is_d_aa() ) positions.push_back( resid );
			else if ( keep_chirality && pose.residue( resid ).type().is_d_aa() ) d_positions.push_back( resid );
		} else if ( pose.residue_type(resid).is_beta_aa() ) {
			if ( !keep_chirality || !pose.residue( resid ).type().is_d_aa() ) beta_positions.push_back( resid );
			else if ( keep_chirality && pose.residue( resid ).type().is_d_aa() ) d_beta_positions.push_back( resid );
		}

		if ( !keep_disulf && ( pose.residue( resid ).type().is_disulfide_bonded() ) ) {
			core::Size const bonded_partner( core::conformation::get_disulf_partner( pose.conformation(), resid ) );
			core::conformation::break_disulfide( pose.conformation(), resid, bonded_partner );
			if ( pose.residue_type(resid).is_alpha_aa() ) {
				if ( !pose.residue( resid ).type().is_d_aa() || !keep_chirality ) {
					protocols::simple_moves::MutateResidue mut( resid, "ALA" );
					mut.apply( pose );
				} else {
					protocols::simple_moves::MutateResidue mut( resid, "DALA" );
					mut.apply( pose );
				}
			} else if ( pose.residue_type(resid).is_beta_aa() ) {
				if ( !pose.residue( resid ).type().is_d_aa() || !keep_chirality ) {
					protocols::simple_moves::MutateResidue mut( resid, "B3A" );
					mut.apply( pose );
				} else {
					protocols::simple_moves::MutateResidue mut( resid, "DB3A" );
					mut.apply( pose );
				}
			}

			if ( pose.residue_type(bonded_partner).is_alpha_aa() ) {
				if ( !pose.residue(bonded_partner).type().is_d_aa() || !keep_chirality ) {
					protocols::simple_moves::MutateResidue mut( bonded_partner, "ALA" );
					mut.apply( pose );
				} else {
					protocols::simple_moves::MutateResidue mut( bonded_partner, "DALA" );
					mut.apply( pose );
				}
			} else if ( pose.residue_type(bonded_partner).is_beta_aa() ) {
				if ( !pose.residue(bonded_partner).type().is_d_aa() || !keep_chirality ) {
					protocols::simple_moves::MutateResidue mut( bonded_partner, "B3A" );
					mut.apply( pose );
				} else {
					protocols::simple_moves::MutateResidue mut( bonded_partner, "DB3A" );
					mut.apply( pose );
				}
			}
		}
	}

	protocols::toolbox::pose_manipulation::construct_poly_ala_pose(
		pose, positions,
		false, // bool keep_pro,
		true, // bool keep_gly,
		keep_disulf ); // bool keep_disulfide_cys
	if ( keep_chirality && d_positions.size() > 0 ) {
		protocols::toolbox::pose_manipulation::construct_poly_d_ala_pose(
			pose, d_positions,
			false, // bool keep_pro,
			true, // bool keep_gly,
			keep_disulf ); // bool keep_disulfide_cys
	}
	protocols::toolbox::pose_manipulation::construct_poly_beta_ala_pose(
		pose, beta_positions,
		false, // bool keep_pro,
		true, // bool keep_gly,
		keep_disulf ); // bool keep_disulfide_cys
	if ( keep_chirality && d_positions.size() > 0 ) {
		protocols::toolbox::pose_manipulation::construct_poly_d_beta_ala_pose(
			pose, d_beta_positions,
			false, // bool keep_pro,
			true, // bool keep_gly,
			keep_disulf ); // bool keep_disulfide_cys
	}
}

core::pose::PoseOP
construct_dummy_pose( std::string const & restype_name )
{
	return construct_dummy_pose( restype_name, 2 );
}

core::pose::PoseOP
construct_dummy_pose( std::string const & restype_name, core::Size const length )
{
	core::chemical::ResidueTypeSetCOP typeset = core::chemical::rsd_set_from_cmd_line().lock();
	core::chemical::ResidueType const & typ = typeset->name_map( restype_name );
	return construct_dummy_pose( typ, length );
}

core::pose::PoseOP
construct_dummy_pose( core::chemical::ResidueType const & restype )
{
	return construct_dummy_pose( restype, 2 );
}

core::pose::PoseOP
construct_dummy_pose( core::chemical::ResidueType const & restype, core::Size length )
{
	debug_assert( length );
	core::pose::PoseOP newp = core::pose::PoseOP( new core::pose::Pose() );
	core::conformation::ResidueOP newres = core::conformation::ResidueFactory::create_residue( restype );
	newp->append_residue_by_jump( *newres, 1 );
	for ( core::Size i=1; i<=length-1; ++i ) {
		newp->append_polymer_residue_after_seqpos( *newres, newp->size(), true );
		newp->set_omega( newp->size()-1, 180.0 );
	}
	core::pose::add_lower_terminus_type_to_pose_residue( *newp, 1 );
	core::pose::add_upper_terminus_type_to_pose_residue( *newp, newp->size() );
	return newp;
}

core::select::residue_selector::ResidueSelectorCOP
get_residue_selector( basic::datacache::DataMap const & data, std::string const & name )
{
	core::select::residue_selector::ResidueSelectorCOP selector;
	try {
		selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", name );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream error_msg;
		error_msg << "Failed to find ResidueSelector named '" << name << "' from the Datamap.\n";
		error_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
	}
	debug_assert( selector );
	return selector;
}

std::string
abego_str( utility::vector1< std::string > const & abego )
{
	std::string ab = "";
	for ( auto a=abego.begin(), enda=abego.end(); a!=enda; ++a ) {
		debug_assert( a->size() == 1 );
		ab += *a;
	}
	debug_assert( ab.size() == abego.size() );
	return ab;
}

utility::vector1< std::string >
abego_vector( std::string const & ab )
{
	utility::vector1< std::string > abego;
	for ( char c : ab ) {
		std::string res_abego = "";
		res_abego += c;
		abego.push_back( res_abego );
	}
	debug_assert( abego.size() == ab.size() );
	return abego;
}

// gets a remark line, pasting multiple lines together if necessary
std::string
get_remark_line( core::io::Remarks::const_iterator & it_rem, core::io::Remarks::const_iterator const & end )
{
	std::string line = it_rem->value;
	int num = it_rem->num;
	while ( ( it_rem != end ) && ( line.substr( line.size()-1, 1 ) == "#" ) ) {
		++it_rem;
		debug_assert( num == it_rem->num );
		line = line.substr(0, line.size()-1) + it_rem->value;
	}
	boost::algorithm::trim( line );
	return line;
}

/// @brief adds a remark to a Remarks object, splitting it into multiple remarks if it is too long
void add_remark( core::io::Remarks & remarks, core::Size const num, std::string const & str_val )
{
	core::Size i = 0;
	while ( i < str_val.size() ) {
		core::io::RemarkInfo remark;
		remark.num = num;
		int chunksize = str_val.size() - i;
		bool add_pound = false;
		if ( chunksize >= 59 ) {
			chunksize = 58;
			add_pound = true;
		}
		remark.value = str_val.substr( i, chunksize );
		if ( add_pound ) {
			remark.value += "#";
		}
		remarks.push_back( remark );
		i += chunksize;
	}
}

// helper function to calculate stop of loop without overlap
core::Size
loop_stop_without_overlap( core::pose::Pose const & pose, core::Size stopres, core::Size const overlap )
{
	// calculate start component without overlap
	for ( core::Size i = 1; i <= overlap; ++i ) {
		// see if this loop is lower-terminal
		if ( ( stopres == 1 ) || //loop starts at first residue
				( ! pose.residue( stopres-1 ).is_protein() ) || //residue before start is not protein
				( pose.chain( stopres-1 ) != pose.chain( stopres ) ) || // residues before start are on another chain
				( pose.residue( stopres ).is_lower_terminus() ) ) { // start of residue is lower terminus
			break;
		}
		--stopres;
	}
	return stopres;
}

// helper function to calculate stop residue of loop without overlap
core::Size
loop_start_without_overlap( core::pose::Pose const & pose, core::Size startres, core::Size const overlap )
{
	// calculate start component without overlap
	for ( core::Size i = 1; i <= overlap; ++i ) {
		// see if this loop is upper-terminal
		if ( ( startres == pose.size() ) || // loop end at last residue
				( !pose.residue( startres+1 ).is_protein() ) || // residue after end is not protein
				( pose.chain( startres+1 ) != pose.chain( startres ) ) || // residues before start is other chain
				( pose.residue( startres ).is_upper_terminus() ) ) { // explicit terminus variant @ end of loop
			break;
		}
		++startres;
	}
	return startres;
}

/// @brief given a residue, rebuilds all missing atoms
void rebuild_missing_atoms( core::pose::Pose & pose, core::Size const resi )
{
	core::conformation::Residue const & res = pose.residue( resi );
	bool any_missing = false;
	utility::vector1< bool > missing;
	for ( core::Size i=1; i<=res.type().natoms(); ++i ) {
		if ( res.has( res.type().atom_name(i) ) ) {
			missing.push_back( false );
		} else {
			missing.push_back( true );
			any_missing = true;
		}
	}

	if ( any_missing ) {
		TR << "Rebuilding residue " << resi << std::endl;
		core::conformation::Residue newres =  pose.residue(resi);
		newres.fill_missing_atoms( missing, pose.conformation() );
		pose.conformation().replace_residue( resi, newres, false );
	} else {
		TR << "All atoms present" << std::endl;
	}
}

/// @brief helper function that looks for the given residue in a fold tree and returns the jump that controls its 6D-DoFs
int find_jump_rec(
	core::kinematics::FoldTree const & ft,
	core::Size const residue )
{
	debug_assert( residue > 0 );
	debug_assert( residue <= ft.nres() );

	// if this residue is root, jump is 0
	if ( ft.root() == residue ) {
		return 0;
	}

	// search for jump edges that contains this residue
	for ( auto const & e : ft ) {
		if ( ( e.label() > 0 ) && ( e.stop() == residue ) ) {
			return e.label();
		}
	}

	// search upstream for peptide edges containing this residue
	for ( auto e = ft.begin(); e != ft.end(); ++e ) {
		if ( e->label() < 0 ) {
			if ( ( ( e->start() < residue ) && ( residue <= e->stop() ) ) ||
					( ( e->stop() <= residue ) && ( residue < e->start() ) ) ) {
				return find_jump_rec( ft, e->start() );
			}
		}
	}
	return -1;
}

/// @brief inserts the peptide edges to accomodate the new jump edge given
void insert_peptide_edges( core::kinematics::FoldTree & ft, core::kinematics::Edge const & jedge )
{
	debug_assert( jedge.label() > 0 );
	core::Size pos1 = jedge.start();
	core::Size pos2 = jedge.stop();
	utility::vector1< core::kinematics::Edge > new_edges;
	utility::vector1< core::kinematics::Edge > remove_edges;
	core::kinematics::FoldTree const & ft_const = ft; // ugly hack needed to prevent gcc from trying to use the protected non-const FoldTree::begin() method...
	for ( auto const & it : ft_const ) {
		if ( it.label() != core::kinematics::Edge::PEPTIDE ) continue;
		if ( it.start() <= pos1 && it.stop() >= pos1 ) {
			//disallow edges to self
			if ( it.start() != pos1 ) {
				new_edges.push_back( core::kinematics::Edge( it.start(), pos1, core::kinematics::Edge::PEPTIDE ) );
			}
			if ( it.stop() != pos1 ) {
				new_edges.push_back( core::kinematics::Edge( pos1, it.stop(), core::kinematics::Edge::PEPTIDE ) );
			}
			remove_edges.push_back( it );
		} else if ( it.stop() <= pos1 && it.start() >= pos1 ) { // edges not always in sequential order (eg - jump in middle of chain)
			//disallow edges to self
			if ( it.start() != pos1 ) {
				new_edges.push_back( core::kinematics::Edge( pos1, it.start(), core::kinematics::Edge::PEPTIDE ) );
			}
			if ( it.stop() != pos1 ) {
				new_edges.push_back( core::kinematics::Edge( it.stop(), pos1, core::kinematics::Edge::PEPTIDE ) );
			}
			remove_edges.push_back( it );
		}
		if ( it.start() <= pos2 && it.stop() >= pos2 ) {
			if ( it.start() != pos2 ) {
				new_edges.push_back( core::kinematics::Edge( it.start(), pos2, core::kinematics::Edge::PEPTIDE ) );
			}
			if ( it.stop() != pos2 ) {
				new_edges.push_back( core::kinematics::Edge( pos2, it.stop(), core::kinematics::Edge::PEPTIDE ) );
			}
			remove_edges.push_back( it );
		} else if ( it.stop() <= pos2 && it.start() >= pos2 ) { // the backwards version of the above
			//TR.Debug << "start-pos2-stop " << start <<" " << pos2 << " " << stop << std::endl;
			if ( it.start() != pos2 ) {
				new_edges.push_back( core::kinematics::Edge( pos2, it.start(), core::kinematics::Edge::PEPTIDE ) );
			}
			if ( it.stop() != pos2 ) {
				new_edges.push_back( core::kinematics::Edge( it.stop(), pos2, core::kinematics::Edge::PEPTIDE ) );
			}
			remove_edges.push_back( it );
		}
	}
	for ( auto & remove_edge : remove_edges ) {
		TR.Debug << "Removing edge " << remove_edge << std::endl;
		ft.delete_edge( remove_edge );
	}
	for ( auto & new_edge : new_edges ) {
		TR.Debug << "Adding edge " << new_edge << std::endl;
		ft.add_edge( new_edge );
	}
	TR.Debug << "FT after peptide edge redo: " << ft << std::endl;
}

// @brief parses a string containing single integers and ranges. Returns a vector of all possible values.
utility::vector1< core::Size >
parse_length_string( std::string const & len_str )
{
	utility::vector1< core::Size > retval;
	utility::vector1< std::string > const str_residues( utility::string_split( len_str , ',' ) );
	for ( core::Size i=1; i<=str_residues.size(); ++i ) {
		if ( str_residues[i] == "" ) continue;
		utility::vector1< std::string > const ranges( utility::string_split( str_residues[i], ':' ) );
		if ( ranges.size() == 1 ) {
			retval.push_back( boost::lexical_cast< core::Size >( ranges[1] ) );
		} else if ( ranges.size() == 2 ) {
			core::Size const start( boost::lexical_cast< core::Size >( ranges[1] ) );
			core::Size const end( boost::lexical_cast< core::Size >( ranges[2] ) );
			for ( core::Size i=start; i<=end; ++i ) {
				retval.push_back( i );
			}
		} else {
			throw utility::excn::EXCN_Msg_Exception( "Invalid length input: " + len_str );
		}
	}
	return retval;
}

/// @brief given a number 0 <= x < 1, calculate an integer M <= x <= N
/// NOTE THAT THIS FUNCTION MODIFIES THE PARAMETER
core::Size
extract_int( core::Real & num, core::Size const m, core::Size const n )
{
	debug_assert( num < 1 );
	debug_assert( num >= 0 );
	if ( n == m ) {
		return m;
	}
	debug_assert( n > m );
	TR.Debug << "num: " << num << " m: " << m << " n: " << n << std::endl;
	core::Size const len( n-m+1 );
	num *= len;
	int const val( static_cast< int >(num) );
	num -= val;
	TR.Debug << "num: " << num << " len: " << len << " val: " << val << std::endl;
	debug_assert( val >= 0 );
	debug_assert( val < static_cast< int >(n) );
	return val + m;
}

std::pair< std::string, std::string >
parse_strand_pair( std::string const & strand_pair_str )
{
	utility::vector1< std::string > const strandvec = utility::string_split( strand_pair_str, ',' );
	if ( strandvec.size() != 2 ) {
		std::stringstream err;
		err << "The given strand pair (" << strand_pair_str << ") has more or less than 2 values. It must have exactly two.  Edge strands can be represented by blanks (e.g. \"\",sheet.s2 )." << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	return std::make_pair( strandvec[ 1 ], strandvec[ 2 ] );
}

core::Size
count_bulges( components::StructureData const & perm, std::string const & segment )
{
	std::string const & secstruct = perm.segment( segment ).ss();
	std::string const & abego = perm.segment( segment ).abego();
	if ( secstruct.size() != abego.size() ) {
		std::stringstream msg;
		msg << "ERROR: secstruct " << secstruct << ", size of " << secstruct.size()
			<< " not the same as the size of the abego (" << abego << ") of the segment: "
			<< abego.size() << std::endl;
		utility_exit_with_message( msg.str() );
	}

	core::Size bulges = 0;
	// iterate through vector and string simultaneously
	std::string::const_iterator ab = abego.begin();
	std::string::const_iterator ss = secstruct.begin();
	for ( ; ss!=secstruct.end() && ab!=abego.end(); ++ss, ++ab ) {
		if ( ( *ss == 'E' ) && ( *ab == 'A' ) ) {
			++bulges;
		}
	}
	return bulges;
}

/// @brief dumps a pose into another pose as a new chain
void
add_chain_from_pose( core::pose::PoseCOP to_add, core::pose::PoseOP combined )
{
	runtime_assert( to_add );
	runtime_assert( combined );

	if ( ! to_add->size() ) return;

	if ( combined->size() ) {
		// here we want an anchor equal to the root of the fold tree
		core::Size const anchor_res = 1;
		combined->conformation().buffer_signals();
		combined->conformation().insert_conformation_by_jump( to_add->conformation(), combined->size()+1, combined->num_jump()+2, anchor_res, combined->num_jump()+1 );
		combined->conformation().unblock_signals();
	} else {
		*combined = *to_add;
	}
	// copy remarks
	if ( to_add->pdb_info() ) {
		for ( auto const & r : to_add->pdb_info()->remarks() ) {
			TR.Debug << "Copying remark to new pose: " << r.value << std::endl;
			if ( !combined->pdb_info() ) {
				combined->pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( *combined, true ) ) );
			}
			debug_assert( combined->pdb_info() );
			combined->pdb_info()->remarks().push_back( r );
		}
	}
	TR << "Added segment to pose of length " << to_add->size() << std::endl;
}

/// @brief adds residues from template_pose to pose.  If new_chain == true, creates covalent bond
void
add_residues_to_pose(
	core::pose::Pose & pose,
	core::pose::Pose const & template_pose,
	bool const new_chain )
{
	if ( !template_pose.size() ) return;

	if ( pose.empty() ) {
		pose = template_pose;
		return;
	}

	pose.conformation().buffer_signals();
	if ( new_chain ) {
		core::Size const anchor_res = 1;
		pose.conformation().insert_conformation_by_jump(
			template_pose.conformation(),
			pose.size() + 1,
			pose.num_jump() + 2,
			anchor_res,
			pose.num_jump() + 1 );
	} else { // new_chain == false
		for ( core::Size resid=1; resid<=template_pose.size(); ++resid ) {
			pose.append_polymer_residue_after_seqpos( template_pose.residue( resid ), pose.size(), false );
		}
	}
	pose.conformation().unblock_signals();
	TR.Debug << "Added pose segment of length " << template_pose.size() << std::endl;
}

core::Size
get_resid(
	components::StructureData const & perm,
	std::string const & resid_str )
{
	std::stringstream target_stream( resid_str );
	core::Size const target_resid = boost::lexical_cast< core::Size >( perm.substitute_variables( target_stream ) );
	if ( ( target_resid == 0 ) || ( target_resid > perm.pose_length() ) ) {
		std::stringstream msg;
		msg << "get_resid(): Invalid residue string (" << resid_str
			<< ") was given. A bad target resid (" << target_resid
			<< ") was produced." << std::endl;
		msg << perm << std::endl;
		utility_exit_with_message( msg.str() );
	}
	debug_assert( target_resid );
	debug_assert( target_resid <= perm.pose_length() );
	return target_resid;
}

/// @brief evaluate linear chainbreak at a position
/// @remarks If necessary, will evaluate using a copy of the Pose with a cut
///  fold tree.  If cutpoint variants are present at chainbreak, will use
///  existing variants and not modify them.  If cutpoint variants are not
///  found will add them and then remove them once calculation is finished.
///  Eventually this should be merged into protocols/forge/methods/linear_chainbreak.cc
///  However, the changes I needed to make to that file break certain parts
///  of remodel
core::Real
linear_chainbreak(
	core::pose::Pose const & pose,
	core::Size const pos
)
{
	using core::kinematics::FoldTree;
	using core::pose::Pose;
	using core::scoring::EnergyMap;
	using core::scoring::ScoreFunction;
	using core::scoring::methods::LinearChainbreakEnergy;

	if ( pose.fold_tree().num_cutpoint() == 0 ) {
		return core::Real( 0.0 );
	}
	debug_assert( pos > 0 );
	debug_assert( pos < pose.size() );

	Pose scratch = pose;

	FoldTree ft;
	ft.add_edge( 1, scratch.size(), core::kinematics::Edge::PEPTIDE );
	ft.new_jump( pos, pos + 1, pos );

	protocols::forge::methods::add_cutpoint_variants( scratch, pos );

	EnergyMap emap;
	ScoreFunction const fx; // dummy, needed for function call
	LinearChainbreakEnergy energy;
	energy.finalize_total_energy( scratch, fx, emap );

	return emap[ core::scoring::linear_chainbreak ];
}

core::kinematics::FoldTree
slide_jump(
	core::kinematics::FoldTree const & ft_orig,
	core::Size const jump_idx,
	core::Size const new_start,
	core::Size const new_stop )
{
	core::kinematics::FoldTree ft;
	TR << "Sliding jump " << jump_idx << " to " << new_start << " --> " << new_stop << " in " << ft_orig << std::endl;
	debug_assert( jump_idx > 0 );
	debug_assert( jump_idx <= ft_orig.num_jump() );
	core::kinematics::Edge jedge = ft_orig.jump_edge( jump_idx );
	auto begin_edge = ft_orig.begin();
	debug_assert( begin_edge != ft_orig.end() );

	if ( ft_orig.root() == jedge.start() ) {
		for ( auto e=ft_orig.begin(); e!=ft_orig.end(); ++e ) {
			if ( ( e->start() == ft_orig.root() ) && ( *e != jedge ) &&
					( static_cast< core::Size >( e->stop() ) != new_start ) &&
					( static_cast< core::Size >( e->stop() ) != new_stop ) ) {
				begin_edge = e;
				break;
			}
		}
	}
	TR << "Selected new begin edge: " << *begin_edge << std::endl;

	ft.add_edge( *begin_edge );
	for ( auto e=ft_orig.begin(); e!=ft_orig.end(); ++e ) {
		if ( (e != begin_edge) && (*e != jedge) ) {
			ft.add_edge( *e );
		}
	}
	ft.add_edge( jedge );
	TR.Debug << "Sliding jump " << jump_idx << " to " << new_start << " --> " << new_stop << " in " << ft << std::endl;
	ft.slide_jump( jump_idx, new_start, new_stop );
	debug_assert( ft.check_fold_tree() );
	return ft;
}

void
add_cutpoints( core::pose::Pose & pose, components::StructureData const & sd )
{
	for ( auto s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
		core::Size const cut = sd.segment( *s ).cutpoint();
		if ( cut ) {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, cut );
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, cut+1 );
			pose.conformation().declare_chemical_bond( cut, pose.residue( cut ).atom_name( pose.residue( cut ).upper_connect_atom() ),
				cut + 1, pose.residue( cut + 1 ).atom_name( pose.residue( cut + 1 ).lower_connect_atom() ) );

		}
	}
}

/// @brief Given a symmetric pose and a secstruct for the asymmetric unit, constructs and
///        returns a secondary structure string compatible with the symmetric pose based on
///        the secondary structure of the asymmetric unit
std::string
symmetric_secstruct( core::pose::Pose const & pose, std::string const & asymm_secstruct )
{
	using core::conformation::Conformation;
	using core::conformation::symmetry::SymmetricConformation;
	using core::conformation::symmetry::SymmetryInfo;
	using core::kinematics::FoldTree;

	Conformation const & conf( pose.conformation() );
	SymmetricConformation const & symm_conf( dynamic_cast< SymmetricConformation const & >( conf ) );
	SymmetryInfo const & symm_info( *symm_conf.Symmetry_Info() );
	core::Size const nres_subunit( symm_info.num_independent_residues() );
	core::Size const nsubunits( symm_info.subunits() );
	core::Size const num_nonvrt( symm_info.num_total_residues_without_pseudo() );

	if ( nres_subunit != asymm_secstruct.size() ) {
		std::stringstream msg;
		msg << "protocols::denovo_design::symmetric_secstruct(): The secondary structure for the asymmetric unit ("
			<< asymm_secstruct << ") has length (" << asymm_secstruct.size()
			<< ") that differs from the length of each symmetric subunit (" << nres_subunit << ")" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	std::stringstream symm_secstruct;
	// add asymm unit secstruct for each subunit
	for ( core::Size subunit=1; subunit<=nsubunits; ++subunit ) {
		symm_secstruct << asymm_secstruct;
	}

	// add secstruct for virtuals from the pose
	for ( core::Size resid=num_nonvrt+1; resid<=pose.size(); ++resid ) {
		symm_secstruct << pose.secstruct( resid );
	}

	if ( symm_secstruct.str().size() != pose.size() ) {
		std::stringstream msg;
		msg << "protocols::denovo_design::symmetric_secstruct(): The generated secondary structure for the symmetric pose ("
			<< symm_secstruct.str() << ") has length (" << symm_secstruct.str().size()
			<< ") that differs from the length of the symmetric pose (" << pose.size() << ")" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	return symm_secstruct.str();
}


/// @brief Given a symmetric pose, and a fold tree for the asymmetric unit, constructs and
///        returns a symmetric fold tree while preserving the topology of the aysmmetric
///        unit's fold tree
core::kinematics::FoldTree
symmetric_fold_tree( core::pose::Pose const & pose, core::kinematics::FoldTree const & asymm_ft )
{
	using core::conformation::Conformation;
	using core::conformation::symmetry::SymmetricConformation;
	using core::conformation::symmetry::SymmetryInfo;
	using core::kinematics::FoldTree;

	Conformation const & conf( pose.conformation() );
	SymmetricConformation const & symm_conf( dynamic_cast< SymmetricConformation const & >( conf ) );
	SymmetryInfo const & symm_info( *symm_conf.Symmetry_Info() );
	core::Size const nres_subunit( symm_info.num_independent_residues() );
	core::Size const nsubunits( symm_info.subunits() );
	core::Size const num_nonvrt( symm_info.num_total_residues_without_pseudo() );
	core::Size const root = pose.fold_tree().root();

	if ( num_nonvrt + 1 != root ) {
		std::stringstream msg;
		msg << "FoldTreeFromFoldGraphMover::symmetric_fold_tree(): The residue after num_nonvrt ("
			<< num_nonvrt + 1 << ") should be the root, but the root is " << root << std::endl;
		msg << "Pose fold tree: " << pose.fold_tree() << std::endl;
		msg << "Asymm fold tree: " << asymm_ft << std::endl;
		utility_exit_with_message( msg.str() );
	}

	if ( nsubunits != ( pose.size() - num_nonvrt ) ) {
		std::stringstream msg;
		msg << "FoldTreeFromFoldGraphMover::symmetric_fold_tree(): number of subunits ("
			<< nsubunits << ") does not match the number of virtuals ("
			<< pose.size() - num_nonvrt << ")!" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	// Build symmetric fold tree
	FoldTree symm_ft;

	// Count jumps and resids as we add them
	utility::vector1< int > added_jumps;
	int cur_jump = asymm_ft.num_jump() + 1;
	TR.Debug << "Cur_jump=" << cur_jump << std::endl;
	core::Size cur_pose_resid = 1;

	// add root so there is a vertex in there
	symm_ft.add_edge( cur_pose_resid, cur_pose_resid, core::kinematics::Edge::PEPTIDE );

	// add subunit jumps from virtuals
	for ( core::Size resid=num_nonvrt+1; resid<=pose.size(); ++resid ) {
		debug_assert( cur_pose_resid <= num_nonvrt );
		TR.Debug << "inserting asymm fold tree, resid=" << cur_pose_resid
			<< " jump=" << cur_jump << " cur_ft=" << symm_ft << std::endl;
		symm_ft.insert_fold_tree_by_jump( asymm_ft, cur_pose_resid, cur_jump, cur_pose_resid );
		if ( resid != root ) added_jumps.push_back( cur_jump );
		++cur_jump;
		cur_pose_resid += nres_subunit;
	}

	// at this point, the next residue should be the root
	if ( cur_pose_resid != root ) {
		std::stringstream msg;
		msg << "FoldTreeFromFoldGraphMover::symmetric_fold_tree(): root of new symmetric fold tree ("
			<< cur_pose_resid << ") does not match root of pose fold tree (" << root << std::endl;
		msg << "Pose fold tree: " << pose.fold_tree() << std::endl;
		msg << "Asymm fold tree: " << asymm_ft << std::endl;
		msg << "Symm fold tree: " << symm_ft << std::endl;
		utility_exit_with_message( msg.str() );
	}

	utility::vector1< int >::const_iterator jump = added_jumps.begin();
	// add root jumps from virtuals
	for ( core::Size resid=root; resid<=pose.size(); ++resid ) {
		debug_assert( pose.residue( resid ).aa() == core::chemical::aa_vrt );
		if ( resid == root ) continue;
		if ( jump == added_jumps.end() ) {
			std::stringstream msg;
			msg << "FoldTreeFromFoldGraphMover::symmetric_fold_tree(): list of jumps added ("
				<< added_jumps << ") does not match virtual residues, which range from "
				<< num_nonvrt + 1 << " to " << pose.size() << std::endl;
			utility_exit_with_message( msg.str() );
		}
		symm_ft.jump_edge( *jump ).start() = resid;
		symm_ft.add_edge( root, resid, cur_jump );
		++jump;
		++cur_jump;
	}

	symm_ft.delete_extra_vertices();

	TR << "Orig FT " << asymm_ft << std::endl;
	TR << "Created symmetric FT " << symm_ft << std::endl;
	return symm_ft;
}

/// @brief Given a symmetric pose and a ResidueSubset for the asymmetric unit, constructs
///        and returns a residue subset compatible with the symmetric pose based on the
///        given asymmetric unit residue subset
core::select::residue_selector::ResidueSubset
symmetric_residue_subset( core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & asymm_subset )
{
	using core::conformation::Conformation;
	using core::conformation::symmetry::SymmetricConformation;
	using core::conformation::symmetry::SymmetryInfo;
	using core::kinematics::FoldTree;

	Conformation const & conf( pose.conformation() );
	SymmetricConformation const & symm_conf( dynamic_cast< SymmetricConformation const & >( conf ) );
	SymmetryInfo const & symm_info( *symm_conf.Symmetry_Info() );
	core::Size const nres_subunit( symm_info.num_independent_residues() );
	core::Size const nsubunits( symm_info.subunits() );
	//core::Size const num_nonvrt( symm_info.num_total_residues_without_pseudo() );

	if ( nres_subunit != asymm_subset.size() ) {
		std::stringstream msg;
		msg << "protocols::denovo_design::symmetric_residue_subset(): The residue subset given for the asymmetric unit has length ("
			<< asymm_subset.size() << ") that differs from the length of each symmetric subunit ("
			<< nres_subunit << ")" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	core::select::residue_selector::ResidueSubset subset( pose.size(), false );
	core::Size cur_resid = 1;
	for ( core::Size subunit=1; subunit<=nsubunits; ++subunit ) {
		for ( core::Size sub_resid=1; sub_resid<=nres_subunit; ++sub_resid, ++cur_resid ) {
			subset[ cur_resid ] = asymm_subset[ sub_resid ];
		}
	}
	return subset;
}

/// @brief Computes secondary structure string from the given motifs
/// @param[in]  motif_str  Motif string to be parsed (e.g. "5EB-2LG-5EB")
/// @param[out] secstruct  Secondary structure string to be cleared and filled
/// @param[out] abego      ABEGO string to be cleared and filled
void
parse_motif_string( std::string const & motif_str, std::string & secstruct, std::string & abego )
{
	secstruct.clear();
	abego.clear();

	utility::vector1< std::string > const motifs = utility::string_split( motif_str, '-' );
	for ( auto const & motif : motifs ) {
		// here, we can accept "3LX" or "3:LX"
		std::string motif_seg = "";
		for ( char c : motif ) {
			if ( c == ' ' ) continue;
			if ( c == '\t' ) continue;
			if ( c == '\n' ) continue;
			if ( c == ':' ) continue;
			motif_seg += c;
		}

		if ( motif_seg.empty() ) continue;

		char const ss_type( motif_seg[motif_seg.size()-2] );
		if ( (ss_type != 'H') && (ss_type != 'L') && (ss_type != 'E') ) {
			TR.Error << "Segment::parse_motif(): Invalid SS type in motif " << motif_seg << std::endl;
			utility_exit();
		}

		char const abego_type( motif_seg[ motif_seg.size() - 1 ] );
		if ( abego_type > 'Z' || abego_type < 'A' ) {
			TR.Error << "Segment::parse_motif(): Invalid abego type in motif " << motif_seg << std::endl;
			utility_exit();
		}

		int const len( utility::string2int( motif_seg.substr( 0, motif_seg.size()-2 ) ) );

		std::string const secstruct_m( len, ss_type );
		std::string const abego_m( len, abego_type );
		secstruct += secstruct_m;
		abego += abego_m;
	}
}


} // protocols
} // denovo_design
