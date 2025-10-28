// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/constraints/FabConstraint.cc
/// @brief This class is specific to antibodies and penalizes presence of non-cdr residues
/// @brief at Antigen-Antibody interfaces. It has been ported from Fab constraint in rosetta++
/// @brief which uses a constant constraint score of 0.5 for flanking residues)
/// @author Krishna Kilambi (kkpraneeth@jhu.edu, April 2012)

#include <core/scoring/constraints/FabConstraint.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/func/ConstantFunc.hh>
#include <core/scoring/func/ConstantFunc.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>


#include <core/id/AtomID.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/chains_util.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <string>

static basic::Tracer TR( "core.scoring.constraints.FabConstraint" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
FabConstraint::FabConstraint():
	MultiConstraint()
{}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
FabConstraint::FabConstraint(ConstraintCOPs const & cst_in) :
	MultiConstraint(cst_in)
{}

///
ConstraintOP
FabConstraint::clone() const {
	return ConstraintOP( new FabConstraint( *this ) );
}

std::string
FabConstraint::type() const {
	return "FabConstraint";
}

void
FabConstraint::show(std::ostream& out) const
{
	out << "FabConstraint is an AmbiguousConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
	for ( auto const & cst_it : member_constraints() ) {
		cst_it->show( out );
	}

	out << " ...all member constraints of this FabConstraint shown." << std::endl;
}

void
FabConstraint::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	func::FuncFactory const & /*func_factory*/
) {
	utility::vector1<Size> res1;
	utility::vector1<Size> res2;
	utility::vector1<std::string> tempres1;
	utility::vector1<std::string> tempres2;
	std::string antchains;
	Size flag(1);


	std::string line;
	while ( getline(data, line) ) {
		std::string entry1, entry2, entry3, cstname;
		std::istringstream line_stream(line);
		if ( flag ) {
			line_stream >> entry1 >> entry2 >> entry3;
			flag = 0;
		} else {
			line_stream >> cstname >> entry1 >> entry2 >> entry3;
		}
		TR.Info << "Entry1:" << entry1 << " Entry2:" << entry2 << " Entry3:" << entry3 << std::endl;
		tempres1.push_back(entry1);
		tempres2.push_back(entry2);
		antchains = entry3;
	}

	for ( Size n=1; n<= tempres1.size(); ++n ) {
		res1.push_back(pose_res_no(pose, tempres1[n]));
		res2.push_back(pose_res_no(pose, tempres2[n]));
	}

	TR.Info << "Penalizing residues which are not in range "
		<< res1[1] << "-" << res2[1] << ", "
		<< res1[2] << "-" << res2[2] << ", "
		<< res1[3] << "-" << res2[3] << ", "
		<< res1[4] << "-" << res2[4] << ", "
		<< res1[5] << "-" << res2[5] << ", "
		<< res1[6] << "-" << res2[6]
		<< " and at interface with antigen chains " << antchains << std::endl;


	// The current input format only supports single letter chains.
	utility::vector1<std::string> antchains_vec;
	for ( char chain: antchains ) {
		antchains_vec.push_back( std::string{chain} );
	}
	setup_csts(pose, res1, res2, antchains_vec);

	if ( data.good() ) {
		//chu skip the rest of line since this is a single line defintion.
		while ( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) data.setstate( std::ios_base::eofbit );
	}
} // read_def

//return pose residue no
Size
FabConstraint::pose_res_no(
	core::pose::Pose const & pose,
	std::string const & res_designation
) {
	Size first_nondigit = res_designation.find_first_not_of("-1234567890");

	int resnum = std::stoi( res_designation.substr(0,first_nondigit) );

	char ins_code = res_designation[first_nondigit];
	std::string chain = res_designation.substr(first_nondigit+1);

	Size pose_resnum = pose.pdb_info()->pdb2pose(chain,resnum,ins_code);
	if ( pose_resnum != 0 ) { // We have the with-insertion code version
		return pose_resnum;
	} else {
		chain = res_designation.substr(first_nondigit);
		return pose.pdb_info()->pdb2pose(chain,resnum);
	}
}

//Build a vector of the associated penalty scores for each antibody residue in the sequence
//with antibody residue pose numbers as indices
utility::vector1<Real>
FabConstraint::calc_penalty_vec(
	Size start_res,
	Size stop_res,
	utility::vector1<Size> res1,
	utility::vector1<Size> res2
) {
	utility::vector1<Real> penalty;
	Size n = 1;
	for ( Size m = 1 ; m <= stop_res ; ++m ) {
		if ( m >= start_res && m <= stop_res ) {
			penalty.push_back(1.5);
			//if you hit the end of the flank region at c-term end of a cdr, go back and reassign the correct penalties
			//for cdrs and cdr flanking regions
			if ( m == res2[n]+2 ) {
				penalty[res1[n]-2] = 0.5;
				penalty[res1[n]-1] = 0.5;
				penalty[res2[n]+1] = 0.5;
				penalty[res2[n]+2] = 0.5;
				for ( Size p = res1[n] ; p <= res2[n] ; ++p ) {
					penalty[p] = 0.0;
				}
				if ( n < res2.size() ) n++;
			}
		} else {
			penalty.push_back(0);
		}
	}
	return penalty;
}


void
FabConstraint::setup_csts(
	core::pose::Pose const & pose,
	utility::vector1<Size> const & res1,
	utility::vector1<Size> const & res2,
	utility::vector1<std::string> const & antchains
) {
	utility::vector1<Real> abpenalty;
	func::ConstantFuncOP flankpenaltyfunc( new func::ConstantFunc(0.5) );
	func::ConstantFuncOP noncdrpenaltyfunc( new func::ConstantFunc(1.5) );

	//set up antigen and antibody chain limits
	Size ant_start_chain = pose::get_chain_id_from_chain(antchains[1], pose);
	Size ant_stop_chain = pose::get_chain_id_from_chain(antchains[antchains.size()], pose);
	Size ab_start_chain = pose.chain(res1[1]);
	Size ab_stop_chain = pose.chain(res1[res1.size()]);

	Size ant_start_res = pose.conformation().chain_begin(ant_start_chain);
	Size ant_stop_res = pose.conformation().chain_end(ant_stop_chain);
	Size ab_start_res = pose.conformation().chain_begin(ab_start_chain);
	Size ab_stop_res = pose.conformation().chain_end(ab_stop_chain);

	abpenalty = calc_penalty_vec(ab_start_res, ab_stop_res, res1, res2);

	/*    TR.Info << "Abpenalty Vector: ";
	for (Size p = 1 ; p <= abpenalty.size() ; ++p){
	TR.Info << abpenalty[p] << " ";
	}
	TR.Info << std::endl;*/

	//Find residues at interface and setup constraints
	//loop over all antibody residues
	for ( Size i = ab_start_res ; i <= ab_stop_res ; ++i ) {
		id::AtomID atom1(pose.residue_type(i).atom_index("CA"),i);
		//now loop over all antigen residues
		for ( Size j = ant_start_res ; j <= ant_stop_res ; ++j ) {
			//check for residues at interface
			if ( pose.residue(i).xyz("CA").distance(pose.residue(j).xyz("CA")) < 8.0 ) {
				TR.Info << "Residue " << pose.residue(i).name3() << " " << pose.pdb_info()->pose2pdb(i) << " is at interface" << std::endl;
				id::AtomID atom2(pose.residue_type(j).atom_index("CA"),j);
				runtime_assert(atom1.valid() && atom2.valid());
				if ( abpenalty[i] == 1.5 ) {
					add_individual_constraint( utility::pointer::make_shared< AtomPairConstraint >(atom1,atom2, noncdrpenaltyfunc));
				} else if ( abpenalty[i] == 0.5 ) {
					add_individual_constraint( utility::pointer::make_shared< AtomPairConstraint >(atom1,atom2, flankpenaltyfunc));
				}
				break;
			}
		}

	}

} // setup_csts

bool
FabConstraint::operator==( Constraint const & rhs ) const
{
	// The base class will ensure both this and rhs are of the same type
	return MultiConstraint::operator== ( rhs );
}

bool FabConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< FabConstraint const * > ( &other );
}

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::FabConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< MultiConstraint >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::FabConstraint::load( Archive & arc ) {
	arc( cereal::base_class< MultiConstraint >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::FabConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::FabConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_FabConstraint )
#endif // SERIALIZATION
