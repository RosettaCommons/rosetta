// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/carbohydrates/pose_io.cc
/// @brief   Pose input/output function definitions for carbohydrate-specific data formats.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/io/carbohydrates/pose_io.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/carbohydrates/util.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <sstream>


// Construct tracer.
static thread_local basic::Tracer TR( "core.io.carbohydrates.pose_io" );


namespace core {
namespace io {
namespace carbohydrates {

using namespace std;
using namespace core;


// Output /////////////////////////////////////////////////////////////////////

// Return a GWS-formatted string for the given carbohydrate residue.
/// @details  The GWS format for (1->4)-beta-D-Galactopyranosyl, e.g., is "--4b1D-Gal,p".
/// (Question marks are used for cases where information is not known or not applicable.)
std::string
residue_gws_string( core::pose::Pose const & pose, core::uint const seqpos )
{
	using namespace conformation;
	using namespace chemical::carbohydrates;
	using namespace pose::carbohydrates;

	Residue const & res( pose.residue( seqpos ) );
	assert( res.is_carbohydrate() );

	stringstream gws_string( stringstream::out );

	CarbohydrateInfoCOP info( res.carbohydrate_info() );

	gws_string << "--";

	core::uint const parent_seqpos( find_seqpos_of_saccharides_parent_residue( res ) );
	if ( parent_seqpos ) {
		Residue const & parent( pose.residue( parent_seqpos ) );

		if ( ! parent.is_carbohydrate() ) {  // glycoconjugate case
			gws_string << '?';
		} else /* parent is cabohydrate */ {
			CarbohydrateInfoCOP parent_info( parent.carbohydrate_info() );
			if ( ! res.is_branch_lower_terminus() ) {
				gws_string << parent_info->mainchain_glycosidic_bond_acceptor();
			} else /* residue is branch lower terminus */ {
				Size const n_branches( parent_info->n_branches() );
				Size const n_non_branch_residue_connections( parent.n_residue_connections() - n_branches );
				for ( core::uint i( 1 ); i <= n_branches; ++i ) {
					// The below assumes that non-polymeric residue connections are in the same numerical order as
					// carbohydrate numbering schemes.  If this ever fails, something more complicated will need to be
					// figured out. ~Labonte
					core::uint const branch_res_conn_id( n_non_branch_residue_connections + i );
					if ( parent.residue_connection_partner( branch_res_conn_id ) == seqpos ) {
						gws_string << parent_info->branch_point( i );
					}
				}
			}
		}
	} else /* is lower terminus */ {
		gws_string << '?';
	}
	assert( gws_string.str().size() == 3 );  // If this fails, it probably indicates a design flaw; see above.

	gws_string << info->anomer()[ 0 ];
	gws_string << info->anomeric_carbon();
	gws_string << info->stereochem();
	gws_string << '-';
	gws_string << res.name3();

	// TODO: Add more sugar modifications as needed.
	if ( info->is_amino_sugar() && res.name3() != "Neu" ) {
		gws_string << 'N';
	}
	if ( info->is_acetylated() ) {
		gws_string << "Ac";
	}
	switch ( info->ring_size() ) {
		case 4:
			gws_string << ",o";
			break;
		case 5:
			gws_string << ",f";
			break;
		case 6:
			gws_string << ",p";
			break;
	}

	return gws_string.str();
}

// Return a GWS-formatted string for each carbohydrate residue in the given sequence range.
/// @details  This function is used recursively to handle branching.
/// Branching in GWS format uses parentheses following the residue to which the branch is attached.
std::string
residue_range_gws_string( core::pose::Pose const & pose, core::uint const begin, core::uint const end ) {
	stringstream gws_string( stringstream::out );

	for ( core::uint i( begin ); i <= end; ++i ) {
		// First, output the information for the ith residue.
		gws_string << residue_gws_string( pose, i );

		// Now, see if this residue has any branches off the main chain.
		conformation::Residue const & res_i( pose.residue( i ) );
		assert( res_i.is_carbohydrate() );
		Size const n_branches( res_i.carbohydrate_info()->n_branches() );
		Size const n_non_branch_residue_connections( res_i.n_residue_connections() - n_branches );

		// If so, begin nesting parentheses for the number of branches.
		for ( core::uint j( 1 ); j <= n_branches; ++j ) {
			gws_string << '(';
		}

		// Now, recursively call self for each branch.
		for ( core::uint j( 1 ); j <= n_branches; ++j ) {
			// Determine the start and stop of each branch.
			core::uint const branch_j_res_conn_id( n_non_branch_residue_connections + j );
			core::uint const branch_j_begin( res_i.residue_connection_partner( branch_j_res_conn_id ) );
			core::uint branch_j_end( branch_j_begin );
			while ( ! pose.residue( branch_j_end ).is_upper_terminus() ) {
				++branch_j_end;
			}

			// Recurse.
			gws_string << residue_range_gws_string( pose, branch_j_begin, branch_j_end );

			// ...And close parentheses for branch j.
			gws_string << ')';
		}
	}

	return gws_string.str();
}

// Return a GWS-formatted string for the given carbohydrate chain, including branches.
std::string
chain_gws_string( core::pose::Pose const & pose, core::uint const chain_id )
{
	using namespace conformation;
	using namespace pose::carbohydrates;

	Conformation const & conf( pose.conformation() );

	assert( chain_id <= conf.num_chains() );
	PyAssert( chain_id <= conf.num_chains(),
			"dump_gws_chain( core::pose::Pose const & pose, core::uint const chain_id, std::string const & filename ): "
			"variable chain_id is out of range!" );

	core::uint begin( conf.chain_begin( chain_id ) );
	core::uint const end( conf.chain_end( chain_id ) );

	if ( begin != end && pose.residue( begin ).type().is_lipid() ) {
		// If the glycan could be a glycolipid, start processing from the 2nd residue of the chain, as the 1st residue
		// will be the lipid head group.
		++begin;
	}

	Residue const & first_res( pose.residue( begin ) );
	if ( first_res.is_carbohydrate() ) {
		stringstream gws_string( stringstream::out );

		// First, deal with the reducing end.
		if ( first_res.is_lower_terminus() ) {
			gws_string << "redEnd";
		} else if ( false ) {
			; // TODO: Add code to check for glycosides and output the R group.
		} else {
			core::uint const parent_seqpos( find_seqpos_of_saccharides_parent_residue( first_res ) );
			Residue const & parent_res( pose.residue( parent_seqpos ) );
			if ( parent_res.is_carbohydrate() ) {
				gws_string << "redEnd";  // perhaps this should be "freeEnd"?
			} else {
				gws_string << parent_res.name3();
				gws_string << pose.pdb_info()->number( parent_seqpos );
				gws_string << "=0.0u";  // TODO: I could grab the mass, but that's probably overkill.
			}
		}

		// Then, loop through each residue in the chain, including branches.
		gws_string << residue_range_gws_string( pose, begin, end ) << ';';

		return gws_string.str();
	} else /* residue is not carbohydrate */ {
		TR.Warning << "Chain " << chain_id << " is not a glycan." << endl;
		return "";
	}
}


// Write the GlycoWorkbench structure file for the given pose chain to <filename>.
/// @details  If the given chain is a glycan, this function writes a GWS file describing the CFG "structure" of that
/// glycan.  If the chain contains branch-point residues, attached branches are also included in the output.
void
dump_gws_chain( core::pose::Pose const & pose, core::uint const chain_id, std::string const & filename )
{
	utility::io::ozstream file( filename.c_str() );
	if ( ! file ) {
		TR.Error << "Unable to write file: " << filename << endl;
		return;
	}
	file << chain_gws_string( pose, chain_id ) << endl;
	file.close();
}

// Write the GlycoWorkbench structure file for all carbohydrate chains of the given pose to <filename>.
/// @details  Chains that are branches are only output once, as a part of their parent chain.
void
dump_gws( core::pose::Pose const & pose, std::string const & filename )
{
	using namespace conformation;
	using namespace pose::carbohydrates;

	stringstream gws_string( stringstream::out );

	Conformation const & conf( pose.conformation() );
	Size const n_chains( conf.num_chains() );
	for ( core::uint i( 1 ); i <= n_chains; ++i ) {
		core::uint begin( conf.chain_begin( i ) );
		if ( ! ( pose.residue( begin ).is_carbohydrate() || pose.residue( begin ).type().is_lipid() ) ) {
			continue;  // Skip non-carbohydrate chains.  (Chains starting with a lipid MIGHT be a glycolipid.)
		}

		core::uint const end( conf.chain_end( i ) );
		if ( begin != end && pose.residue( begin ).type().is_lipid() ) {
			// If the chain could be a glycolipid, start processing from the 2nd residue of the chain, as the 1st
			// residue will be the lipid head group, (but not if the chain is just a single-residue lipid!)
			++begin;
		}

		Residue const & first_res( pose.residue( begin ) );
		if ( first_res.is_carbohydrate() ) {
			// If the first residue of this chain is a carbohydrate, check its parent to ascertain if it has already
			// been output.
			core::uint const parent_seqpos( find_seqpos_of_saccharides_parent_residue( first_res ) );
			if ( parent_seqpos ) {
				Residue const & parent( pose.residue( parent_seqpos ) );
				if ( parent.is_carbohydrate() ) {
					continue;  // Skip this chain, because we already dealt with it when we dealt with its parent.
				}
			}
		}
		gws_string << chain_gws_string( pose, i );
	}
	utility::io::ozstream file( filename.c_str() );
	if ( ! file ) {
		TR.Error << "Unable to write file: " << filename << endl;
		return;
	}
	file << gws_string.str() << endl;
	file.close();
}

}  // namespace carbohydrates
}  // namespace io
}  // namespace core
