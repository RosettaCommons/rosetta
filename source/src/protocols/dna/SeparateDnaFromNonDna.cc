// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author ashworth

// Unit headers
#include <protocols/dna/SeparateDnaFromNonDna.hh>
#include <protocols/dna/SeparateDnaFromNonDnaCreator.hh>

#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
// note: Pose.hh includes Conformation.hh, FoldTree.hh, Jump.hh, and ResidueType.hh
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace dna {

using utility::vector1;
using namespace core;
using namespace chemical;
using namespace kinematics;
using namespace conformation;

static thread_local basic::Tracer TR( "protocols.dna.SeparateDnaFromNonDna", basic::t_info );

std::string
SeparateDnaFromNonDnaCreator::keyname() const
{
	return SeparateDnaFromNonDnaCreator::mover_name();
}

protocols::moves::MoverOP
SeparateDnaFromNonDnaCreator::create_mover() const {
	return protocols::moves::MoverOP( new SeparateDnaFromNonDna );
}

std::string
SeparateDnaFromNonDnaCreator::mover_name()
{
	return "SeparateDnaFromNonDna";
}

SeparateDnaFromNonDna::SeparateDnaFromNonDna()
	: Mover( SeparateDnaFromNonDnaCreator::mover_name() ),
		translation_( 1000, 0, 0 )
{}

SeparateDnaFromNonDna::SeparateDnaFromNonDna( Real x, Real y, Real z )
	: Mover( SeparateDnaFromNonDnaCreator::mover_name() ),
		translation_( x, y, z )
{}

SeparateDnaFromNonDna::SeparateDnaFromNonDna( numeric::xyzVector< Real > const & xyz )
	: Mover( SeparateDnaFromNonDnaCreator::mover_name() ),
		translation_( xyz )
{}

SeparateDnaFromNonDna::SeparateDnaFromNonDna( SeparateDnaFromNonDna const & other ) :
	//utility::pointer::ReferenceCount(),
	Mover( other ),
	translation_( other.translation() )
{}

SeparateDnaFromNonDna::~SeparateDnaFromNonDna(){}

/// @brief set up an appropriate fold tree for a non-DNA/DNA interface, and simply pull apart non-DNA and DNA
/// @details addition fold-tree considerations will have to be made if there are any 'non-DNA' chains that should stick to the DNA instead of to the non-DNA group
void
SeparateDnaFromNonDna::apply( pose::Pose & pose )
{
	TR << "old fold tree:\n" << pose.fold_tree() << std::endl;

	// construct a new fold tree that links non-DNA chains together, then DNA chains together, then links non-DNA group to DNA group by a single jump
	FoldTree fold_tree( pose.total_residue() );

	// collect non-DNA and DNA chain indices
	vector1< Size > non_DNA_chain_indices, DNA_chain_indices;
	Conformation const & conf( pose.conformation() );
	for ( Size chain_index(1); chain_index <= conf.num_chains(); ++chain_index ) {
		if ( pose.residue_type( conf.chain_begin( chain_index ) ).is_DNA() ) {
			DNA_chain_indices.push_back( chain_index );
			TR << "chain " << chain_index << " is a DNA chain" << '\n';
		} else {
			non_DNA_chain_indices.push_back( chain_index );
			TR << "chain " << chain_index << " is a non-DNA chain" << '\n';
		}
		// manually add continuous segment edges
		// NO! FoldTree will attempt to do this on its own, and this will create duplications that will cause assertion/runtime failures
//		fold_tree.add_edge( conf.chain_begin( chain_index ), conf.chain_end( chain_index ), Edge::PEPTIDE );
	}

	// connect non-DNA chains
	Size const first_nonDNA_chain( non_DNA_chain_indices.front() );
	Size const nonDNA_anchor_index( conf.chain_begin( first_nonDNA_chain ) );
	for ( vector1< Size >::const_iterator chain_iter( non_DNA_chain_indices.begin()+1 );
	      chain_iter != non_DNA_chain_indices.end(); ++chain_iter ) {
		Size const chain_begin( conf.chain_begin( *chain_iter ) );

		TR << "linking non-DNA chains " << first_nonDNA_chain << " and " << *chain_iter
			<< " residues " << nonDNA_anchor_index << " (";
		if ( pose.pdb_info() ) {
			TR << pose.pdb_info()->number( nonDNA_anchor_index ) << ") and "
				<< chain_begin << " (" << pose.pdb_info()->number( chain_begin );
		} else {
			TR << nonDNA_anchor_index << ") and " << chain_begin << " (" << chain_begin;
		}
		TR << ") with a jump" << '\n';

		fold_tree.new_jump( nonDNA_anchor_index, chain_begin, chain_begin-1 );

	}

	// connect DNA chains
	Size const first_DNA_chain( DNA_chain_indices.front() );
	Size const DNA_anchor_index( conf.chain_begin( first_DNA_chain ) );
	for ( vector1< Size >::const_iterator chain_iter( DNA_chain_indices.begin()+1 );
	      chain_iter != DNA_chain_indices.end(); ++chain_iter ) {
		Size const chain_begin( conf.chain_begin( *chain_iter ) );

		TR << "linking DNA chains " << first_DNA_chain << " and " << *chain_iter << " residues " << DNA_anchor_index << " (";
		if ( pose.pdb_info() ) {
			TR << pose.pdb_info()->number( DNA_anchor_index ) << ") and "
				<< chain_begin << " (" << pose.pdb_info()->number( chain_begin );
		} else {
			TR << DNA_anchor_index << ") and " << chain_begin << " (" << chain_begin;
		}
		TR << ") with a jump" << '\n';

		fold_tree.new_jump( DNA_anchor_index, chain_begin, chain_begin-1 );
	}

	// connect non-DNA to DNA with a single jump
	TR << "linking non-DNA chain group to DNA chain group by jump between non-DNA chain "
		<< first_nonDNA_chain << " and DNA chain " << first_DNA_chain << " residues "
		<< nonDNA_anchor_index << " (";
	if ( pose.pdb_info() ) {
		TR << pose.pdb_info()->number( nonDNA_anchor_index ) << ") and "
			<< DNA_anchor_index << " (" << pose.pdb_info()->number( DNA_anchor_index );
	} else {
		TR << nonDNA_anchor_index << ") and " << DNA_anchor_index << " (" << DNA_anchor_index;
	}
	TR << ")" << '\n';

	// make sure that this jump goes in the 'forward' direction, in case the DNA chains came first
	Size const jump_start( std::min( nonDNA_anchor_index, DNA_anchor_index ) ),
	           jump_end( std::max( nonDNA_anchor_index, DNA_anchor_index ) );
	fold_tree.new_jump( jump_start, jump_end, jump_end-1 );

	// set the new fold tree
	pose.fold_tree( fold_tree );
	TR << "new fold tree:\n" << pose.fold_tree() << std::endl;

	// lengthen the final jump to isolate non-DNA from DNA
	Size const prot_DNA_jump_index( pose.num_jump() );
	Jump prot_DNA_jump( pose.jump( prot_DNA_jump_index ) );
	prot_DNA_jump.set_translation( translation_ );
	TR << std::fixed << std::setprecision(0);
	TR << "separating non-DNA chains from DNA chains by translation ("
	   << translation_.x() << "," << translation_.y() << "," << translation_.z() << ")" << '\n';
	pose.set_jump( prot_DNA_jump_index, prot_DNA_jump );
	TR.flush();
}

std::string
SeparateDnaFromNonDna::get_name() const {
	return SeparateDnaFromNonDnaCreator::mover_name();
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void SeparateDnaFromNonDna::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{}

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
SeparateDnaFromNonDna::fresh_instance() const
{
	return moves::MoverOP( new SeparateDnaFromNonDna );
}

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
SeparateDnaFromNonDna::clone() const
{
	return moves::MoverOP( new SeparateDnaFromNonDna( *this ) );
}

} // dna
} // protocols
