/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/util.cc
/// @brief util functions for denovo design of structures
/// @details
/// @author Tom Linsky


//Unit Headers
#include <protocols/denovo_design/util.hh>

//Project Headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/Segment.hh>

//Protocol Headers
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

//Core Headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

//Boost/ObjexxFCL Headers
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

//C++ Headers

static THREAD_LOCAL basic::Tracer TR("protocols.denovo_design.components.util");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Tells whether the two given poses are identical based on # resides and dihedrals
bool same_pose( core::pose::Pose const & pose1, core::pose::Pose const & pose2 )
{
	if ( pose1.total_residue() != pose2.total_residue() ) {
		return false;
	}

	for ( core::Size i = 1; i <= pose1.total_residue(); ++i ) {
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
remove_jump_atoms( core::kinematics::FoldTree const & orig )
{
	core::kinematics::FoldTree ft;
	for ( core::kinematics::FoldTree::EdgeList::const_iterator e=orig.begin(); e!=orig.end(); ++e ) {
		if ( e->is_jump() && e->has_atom_info() ) {
			core::kinematics::Edge newedge = *e;
			newedge.start_atom() = "";
			newedge.stop_atom() = "";
			ft.add_edge( newedge );
		} else {
			ft.add_edge( *e );
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
	std::set< core::Size > const & set1,
	std::set< core::Size > const & set2,
	bool const keep_chirality
) {
	std::set< core::Size > setunion = set1;
	for ( std::set< core::Size >::const_iterator r=set2.begin(), endr=set2.end(); r!=endr; ++r ) {
		setunion.insert( *r );
	}
	construct_poly_ala_pose( pose, keep_disulf, setunion, keep_chirality );
}

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
/// @details If keep_chirality is true, the D-amino acids are mutated to D-alanine.
void construct_poly_ala_pose(
	core::pose::Pose & pose,
	bool const keep_disulf,
	std::set< core::Size > const & res_set,
	bool const keep_chirality
) {
	utility::vector1< core::Size > positions;
	utility::vector1< core::Size > d_positions;

	// remove jump atoms, which may cause problems in mutating residues
	pose.fold_tree( remove_jump_atoms( pose.fold_tree() ) );

	for ( core::Size i=1, endi=pose.total_residue(); i<=endi; ++i ) {
		if ( pose.residue(i).is_protein() && ( res_set.find(i) != res_set.end() ) ) {
			if ( !keep_chirality || !pose.residue(i).type().is_d_aa() ) positions.push_back(i);
			else if ( keep_chirality && pose.residue(i).type().is_d_aa() ) d_positions.push_back(i);
		}
		if ( !keep_disulf && ( pose.residue(i).type().is_disulfide_bonded() ) ) {
			if ( pose.residue(i).type().is_l_aa() || !keep_chirality ) {
				protocols::simple_moves::MutateResidue mut( i, "ALA" );
				mut.apply( pose );
			} else {
				protocols::simple_moves::MutateResidue mut( i, "DALA" );
				mut.apply( pose );
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
	assert( selector );
	return selector;
}

std::string
abego_str( utility::vector1< std::string > const & abego )
{
	std::string ab = "";
	for ( utility::vector1< std::string >::const_iterator a=abego.begin(), enda=abego.end(); a!=enda; ++a ) {
		assert( a->size() == 1 );
		ab += *a;
	}
	assert( ab.size() == abego.size() );
	return ab;
}

utility::vector1< std::string >
abego_vector( std::string const & ab )
{
	utility::vector1< std::string > abego;
	for ( std::string::const_iterator c=ab.begin(), endc=ab.end(); c!=endc; ++c ) {
		std::string res_abego = "";
		res_abego += *c;
		abego.push_back( res_abego );
	}
	assert( abego.size() == ab.size() );
	return abego;
}

// gets a remark line, pasting multiple lines together if necessary
std::string get_remark_line( core::io::Remarks::const_iterator & it_rem, core::io::Remarks::const_iterator const & end )
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
		if ( ( startres == pose.total_residue() ) || // loop end at last residue
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
	int const residue )
{
	debug_assert( residue > 0 );
	debug_assert( residue <= static_cast<int>(ft.nres()) );

	// if this residue is root, jump is 0
	if ( ft.root() == residue ) {
		return 0;
	}

	// search for jump edges that contains this residue
	for ( core::kinematics::FoldTree::const_iterator e = ft.begin(); e != ft.end(); ++e ) {
		if ( ( e->label() > 0 ) && ( e->stop() == residue ) ) {
			return e->label();
		}
	}

	// search upstream for peptide edges containing this residue
	for ( core::kinematics::FoldTree::const_iterator e = ft.begin(); e != ft.end(); ++e ) {
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
	assert( jedge.label() > 0 );
	int pos1 = jedge.start();
	int pos2 = jedge.stop();
	utility::vector1< core::kinematics::Edge > new_edges;
	utility::vector1< core::kinematics::Edge > remove_edges;
	core::kinematics::FoldTree const & ft_const = ft; // ugly hack needed to prevent gcc from trying to use the protected non-const FoldTree::begin() method...
	for ( core::kinematics::FoldTree::const_iterator it=ft_const.begin(); it != ft_const.end(); ++it ) {
		if ( it->label() != core::kinematics::Edge::PEPTIDE ) continue;
		if ( it->start() <= pos1 && it->stop() >= pos1 ) {
			//disallow edges to self
			if ( it->start() != pos1 ) {
				new_edges.push_back( core::kinematics::Edge( it->start(), pos1, core::kinematics::Edge::PEPTIDE ) );
			}
			if ( it->stop() != pos1 ) {
				new_edges.push_back( core::kinematics::Edge( pos1, it->stop(), core::kinematics::Edge::PEPTIDE ) );
			}
			remove_edges.push_back( *it );
		} else if ( it->stop() <= pos1 && it->start() >= pos1 ) { // edges not always in sequential order (eg - jump in middle of chain)
			//disallow edges to self
			if ( it->start() != pos1 ) {
				new_edges.push_back( core::kinematics::Edge( pos1, it->start(), core::kinematics::Edge::PEPTIDE ) );
			}
			if ( it->stop() != pos1 ) {
				new_edges.push_back( core::kinematics::Edge( it->stop(), pos1, core::kinematics::Edge::PEPTIDE ) );
			}
			remove_edges.push_back( *it );
		}
		if ( it->start() <= pos2 && it->stop() >= pos2 ) {
			if ( it->start() != pos2 ) {
				new_edges.push_back( core::kinematics::Edge( it->start(), pos2, core::kinematics::Edge::PEPTIDE ) );
			}
			if ( it->stop() != pos2 ) {
				new_edges.push_back( core::kinematics::Edge( pos2, it->stop(), core::kinematics::Edge::PEPTIDE ) );
			}
			remove_edges.push_back( *it );
		} else if ( it->stop() <= pos2 && it->start() >= pos2 ) { // the backwards version of the above
			//TR.Debug << "start-pos2-stop " << start <<" " << pos2 << " " << stop << std::endl;
			if ( it->start() != pos2 ) {
				new_edges.push_back( core::kinematics::Edge( pos2, it->start(), core::kinematics::Edge::PEPTIDE ) );
			}
			if ( it->stop() != pos2 ) {
				new_edges.push_back( core::kinematics::Edge( it->stop(), pos2, core::kinematics::Edge::PEPTIDE ) );
			}
			remove_edges.push_back( *it );
		}
	}
	for ( utility::vector1< core::kinematics::Edge >::iterator it=remove_edges.begin(); it!=remove_edges.end(); ++it ) {
		TR.Debug << "Removing edge " << *it << std::endl;
		ft.delete_edge( *it );
	}
	for ( utility::vector1< core::kinematics::Edge >::iterator it=new_edges.begin(); it!=new_edges.end(); ++it ) {
		TR.Debug << "Adding edge " << *it << std::endl;
		ft.add_edge( *it );
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
	assert( n > m );
	TR.Debug << "num: " << num << " m: " << m << " n: " << n << std::endl;
	core::Size const len( n-m+1 );
	num *= len;
	int const val( static_cast< int >(num) );
	num -= val;
	TR.Debug << "num: " << num << " len: " << len << " val: " << val << std::endl;
	assert( val >= 0 );
	assert( val < static_cast< int >(n) );
	return val + m;
}

/// @brief copies rotamers from the pose "src" into the permutation "dest"
/// no backbone changes are made
/// if detect_disulf flag is on, disulfides will be re-detected
void
copy_rotamers( components::StructureData & dest, core::pose::Pose const & src )
{
	debug_assert( dest.pose_length() == src.total_residue() );

	for ( core::Size r=1; r<=src.total_residue(); ++r ) {
		dest.replace_residue( r, src.residue(r), true );
	}
	// re-detect disulfides
	core::scoring::ScoreFunctionOP sfx = core::scoring::get_score_function();
	dest.detect_disulfides( sfx );
}

std::pair< std::string, std::string >
parse_strand_pair( std::string const & strand_pair_str )
{
	StringVec const strandvec = utility::string_split( strand_pair_str, ',' );
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
	std::string const & ss = perm.segment( segment ).ss();
	utility::vector1< std::string > const & abego = perm.segment( segment ).abego();
	core::Size bulges = 0;
	for ( core::Size resid=perm.segment( segment ).start(); resid<=perm.segment( segment ).stop(); ++resid ) {
		if ( ( ss[ resid - 1 ] == 'E' ) && ( abego[ resid ] == "A" ) ) {
			++bulges;
		}
	}
	return bulges;
}

std::string
strand_pair_str(
	components::StructureData const & perm,
	std::string const & strand1,
	std::string const & strand2,
	std::map< std::string, core::Size > const & name_to_strandnum,
	bool use_register_shift )
{
	if ( strand1.empty() || strand2.empty() ) {
		return "";
	}

	std::stringstream out;
	bool in_order = true;
	core::Size const snum1 = name_to_strandnum.find( strand1 )->second;
	core::Size const snum2 = name_to_strandnum.find( strand2 )->second;
	if ( snum1 < snum2 ) {
		out << snum1 << '-' << snum2;
	} else if ( snum1 > snum2 ) {
		out << snum2 << '-' << snum1;
		in_order = false;
	} else if ( snum1 == snum2 ) {
		std::stringstream ss;
		ss << "Two strands have the same number " << snum1 << " -- check your input! name_to_strandnum=" << name_to_strandnum << " strandnames=" << strand1 << "," << strand2 << std::endl;
		throw utility::excn::EXCN_BadInput( ss.str() );
	}

	int orientation = 1;
	if ( perm.has_data_int( strand2, "orientation" ) ) orientation = perm.get_data_int( strand2, "orientation" );

	int prevorientation = 1;
	if ( perm.has_data_int( strand1, "orientation" ) ) prevorientation = perm.get_data_int( strand1, "orientation" );

	core::Size const bulges1 = count_bulges( perm, strand1 );
	core::Size const bulges2 = count_bulges( perm, strand2 );

	bool parallel = true;
	out << '.';
	if ( orientation == prevorientation ) {
		out << 'P';
	} else {
		parallel = false;
		out << 'A';
	}
	out << '.';

	int shift = 0;
	if ( perm.has_data_int( strand2, "shift" ) ) shift = perm.get_data_int( strand2, "shift" );

	if ( use_register_shift ) {
		// determine "Nobu-style" register shift
		if ( parallel && in_order ) {
			// i + 1's shift is the register shift
			out << shift;
		} else if ( parallel && ! in_order ) {
			// i + 1's shift is simply made negative
			out << -shift;
		} else if ( ! parallel && in_order ) {
			// I'm pretty sure we just use i+1's shift here
			out << shift;
		} else {
			// most complicated
			int const len1 = static_cast< int >( perm.segment( strand2 ).length() ) - bulges2;
			int const len2 = static_cast< int >( perm.segment( strand1 ).length() ) - bulges1;
			out << len2 - len1 - shift;
			TR.Debug << "len1=" << len1 << " bulges1=" << bulges1 << " len2=" << len2 << " bulges2=" << bulges2 << " shift=" << shift << " value=" << len2 - len1 - shift << std::endl;
		}
	} else {
		out << 99;
	}
	return out.str();
}

std::string
get_strandpairings(
	components::StructureData const & perm,
	bool const use_register_shift )
{
	// initiate strand pairing check
	std::map< std::string, core::Size > name_to_strandnum;
	StringVec strandnames;
	core::Size strandcount = 0;
	for ( StringList::const_iterator s=perm.segments_begin(); s!=perm.segments_end(); ++s ) {
		// scan to see if this is a strand
		bool strand = false;
		for ( std::string::const_iterator ch = perm.segment( *s ).ss().begin(); ch != perm.segment( *s ).ss().end(); ++ch ) {
			if ( *ch == 'E' ) {
				strand = true;
				break;
			}
		}
		if ( !strand ) {
			continue;
		}

		++strandcount;
		name_to_strandnum[ *s ] = strandcount;
		strandnames.push_back( *s );
	}
	TR.Debug << "name_to_strandnum=" << name_to_strandnum << " name_to_compidx=" << strandnames << std::endl;

	// now build sheet topology string
	std::stringstream sheet_str;
	std::set< std::pair< std::string, std::string > > visited;
	for ( StringVec::const_iterator s=strandnames.begin(); s!=strandnames.end(); ++s ) {
		if ( !perm.has_data_str( *s, "paired_strands" ) ) {
			continue;
		}
		std::pair< std::string, std::string > const spair =
			parse_strand_pair( perm.get_data_str( *s, "paired_strands" ) );
		TR.Debug << "Found strand pair " << spair.first << " : " << spair.second << std::endl;

		if ( visited.find( std::make_pair( spair.first, *s ) ) == visited.end() ) {
			visited.insert( std::make_pair( spair.first, *s ) );
			visited.insert( std::make_pair( *s, spair.first ) );
			std::string const pair_str = strand_pair_str( perm, spair.first, *s, name_to_strandnum, use_register_shift );
			if ( !pair_str.empty() ) {
				if ( sheet_str.str().size() ) {
					sheet_str << ';';
				}
				sheet_str << pair_str;
			}
		}

		if ( visited.find( std::make_pair( *s, spair.second ) ) == visited.end() ) {
			visited.insert( std::make_pair( spair.second, *s ) );
			visited.insert( std::make_pair( *s, spair.second ) );
			std::string const pair_str = strand_pair_str( perm, *s, spair.second, name_to_strandnum, use_register_shift );
			if ( !pair_str.empty() ) {
				if ( sheet_str.str().size() ) {
					sheet_str << ';';
				}
				sheet_str << strand_pair_str( perm, *s, spair.second, name_to_strandnum, use_register_shift );
			}
		}
	}
	return sheet_str.str();
}

/// @brief dumps a pose into another pose as a new chain
void
add_chain_from_pose( core::pose::PoseCOP to_add, core::pose::PoseOP combined )
{
	runtime_assert( to_add );
	runtime_assert( combined );

	if ( ! to_add->total_residue() ) return;

	if ( combined->total_residue() ) {
		// here we want an anchor equal to the root of the fold tree
		core::Size const anchor_res = 1;
		combined->conformation().buffer_signals();
		combined->conformation().insert_conformation_by_jump( to_add->conformation(), combined->total_residue()+1, combined->num_jump()+2, anchor_res, combined->num_jump()+1 );
		combined->conformation().unblock_signals();
	} else {
		*combined = *to_add;
	}
	// copy remarks
	if ( to_add->pdb_info() ) {
		for ( core::io::Remarks::const_iterator r=to_add->pdb_info()->remarks().begin(); r!=to_add->pdb_info()->remarks().end(); ++r ) {
			TR.Debug << "Copying remark to new pose: " << r->value << std::endl;
			if ( !combined->pdb_info() ) {
				combined->pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( *combined, true ) ) );
			}
			debug_assert( combined->pdb_info() );
			combined->pdb_info()->remarks().push_back( *r );
		}
	}
	TR << "Added segment to pose of length " << to_add->total_residue() << std::endl;
}

} // protocols
} // denovo_design

//////////////////////////////////////////////////////////////////////////
/// Output operators for std classes                                   ///
//////////////////////////////////////////////////////////////////////////

namespace std {

/// @brief outputs a matrix
std::ostream &
operator<<( std::ostream & os, numeric::xyzMatrix< core::Real > const & mat ) {
	os << "[ [" << mat.xx() << ", " << mat.xy() << ", " << mat.xz() << "]" << std::endl;
	os << "  [" << mat.yx() << ", " << mat.yy() << ", " << mat.yz() << "]" << std::endl;
	os << "  [" << mat.zx() << ", " << mat.zy() << ", " << mat.zz() << "] ]";
	return os;
}

/// @brief outputs a set
std::ostream &
operator<<( std::ostream & os, std::set< int > const & set ) {
	os << "[ ";
	for ( std::set< int >::const_iterator it=set.begin(); it != set.end(); ++it ) {
		os << *it << " ";
	}
	os << "]";
	return os;
}

/// @brief outputs a set
std::ostream &
operator<<( std::ostream & os, std::set< core::Size > const & set ) {
	os << "[ ";
	for ( std::set< core::Size >::const_iterator it=set.begin(); it != set.end(); ++it ) {
		os << *it << " ";
	}
	os << "]";
	return os;
}

/// @brief outputs a list of sizes
std::ostream & operator<<( std::ostream & os, std::list< core::Size > const & list )
{
	os << "[ ";
	for ( std::list< core::Size >::const_iterator c=list.begin(), end=list.end(); c != end; ++c ) {
		os << *c << " ";
	}
	os << "]";
	return os;
}

/// @brief outputs a list of strings
std::ostream &
operator<<( std::ostream & os, std::list< std::string > const & list ) {
	os << "[ ";
	for ( std::list< std::string >::const_iterator c=list.begin(), end=list.end(); c != end; ++c ) {
		os << *c << " ";
	}
	os << "]";
	return os;
}

/// @brief outputs a set
std::ostream &
operator<<( std::ostream & os, std::set< std::string > const & set ) {
	os << "[ ";
	for ( std::set< std::string >::const_iterator it=set.begin(); it != set.end(); ++it ) {
		os << *it << " ";
	}
	os << "]";
	return os;
}

/// @brief outputs a map
std::ostream &
operator<<( std::ostream & os, std::map< core::Size, core::Size > const & map ) {
	os << "{";
	std::map< core::Size, core::Size >::const_iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) {
		if ( it != map.begin() ) {
			os << ", ";
		}
		os << " " << it->first << ":" << it->second;
	}
	os << " }";
	return os;
}

/// @brief outputs a map
std::ostream &
operator<<( std::ostream & os, std::map< std::string, core::Size > const & map ) {
	os << "{";
	std::map< std::string, core::Size >::const_iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) {
		if ( it != map.begin() ) {
			os << ", ";
		}
		os << " " << it->first << ":" << it->second;
	}
	os << " }";
	return os;
}

/// @brief outputs a map
std::ostream &
operator<<( std::ostream & os, std::map< std::pair< std::string, std::string >, core::Size > const & map ) {
	os << "{";
	std::map< std::pair< std::string, std::string >, core::Size >::const_iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) {
		if ( it != map.begin() ) {
			os << ", ";
		}
		os << " " << it->first.first << "__" << it->first.second << ":" << it->second;
	}
	os << " }";
	return os;
}

/// @brief outputs a map
std::ostream &
operator<<( std::ostream & os, std::map< std::string, core::Real > const & map ) {
	os << "{";
	std::map< std::string, core::Real >::const_iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) {
		if ( it != map.begin() ) {
			os << ", ";
		}
		os << " " << it->first << ":" << it->second;
	}
	os << " }";
	return os;
}

/// @brief outputs a vector
std::ostream &
operator<<( std::ostream & os, numeric::xyzVector< core::Real > const & vec )
{
	os << "{ " << vec.x() << ", " << vec.y() << ", " << vec.z() << " }";
	return os;
}

/// @brief outputs a map
std::ostream & operator<<( std::ostream & os, std::map< char, core::Size > const & map )
{
	os << "{";
	std::map< char, core::Size >::const_iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) {
		if ( it != map.begin() ) {
			os << ", ";
		}
		os << " " << it->first << ":" << it->second;
	}
	os << "}";
	return os;
}

} // std
