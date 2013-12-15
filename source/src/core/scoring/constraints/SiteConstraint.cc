// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/SiteConstraint.cc
/// @brief This class is an AmbiguousConstraint in which the set is comprised of AtomPairConstraints
/// @brief of an atom of interest in one chain versus the CA of all residues in another chain.
/// @author Brian Weitzner (brian.weitzner@jhu.edu, May 2011)

#include <core/scoring/constraints/SiteConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
// AUTO-REMOVED #include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <core/id/AtomID.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>

static basic::Tracer TR("core.scoring.constraints.SiteConstraint");

namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
SiteConstraint::SiteConstraint():
AmbiguousConstraint()
{}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
SiteConstraint::SiteConstraint( ConstraintCOPs & cst_in ):
AmbiguousConstraint( cst_in )
{}

void
SiteConstraint::show( std::ostream& out) const
{
	/// APL -- you cannot show the active constraint in the absence of a Pose.
	//out << "AmbiguousConstraint Active constraint:" << std::endl;
	//active_constraint()->show(out);
	out << "SiteConstraint is an AmbiguousConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
	for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
		(*cst_it)->show(out);
	}

	out << " ...all member constraints of this SiteConstraint shown." << std::endl;
}


void
SiteConstraint::read_def(
    std::istream & data,
    core::pose::Pose const & pose,
   func::FuncFactory const & func_factory
) {
    TR.Debug << "ConstraintIO::read_site_cst" << std::endl;
    Size res;
    std::string tempres;
    std::string name;
    std::string chain;
    std::string func_type;
    std::string type;
    data >> name >> tempres >> chain >> func_type;

    ConstraintIO::parse_residue( pose, tempres, res );
    TR.Info << "read: " << name << " " << res << " " << "constrain to chain: " << chain << " func: " << func_type << std::endl;



    func::FuncOP aFunc = func_factory.new_func( func_type );
    aFunc->read_data( data );

    if ( TR.Debug.visible() ) {
        aFunc->show_definition( TR.Debug ); TR.Debug<<std::endl;
    }
    setup_csts( res, name, chain, pose, aFunc );

    if ( data.good() ) {
        //chu skip the rest of line since this is a single line defintion.
		while( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) data.setstate( std::ios_base::eofbit );
	}

	if ( TR.Debug.visible() ) {
		aFunc->show_definition( TR.Debug );
		TR.Debug << std::endl;
	}

} // read_def
void
SiteConstraint::setup_csts(
    Size res,
    std::string name,
    std::string chain,
    core::pose::Pose const & pose,
    func::FuncOP const & func
) {
    id::AtomID target_atom( pose.residue_type( res ).atom_index( name ), res );
    //Size target_chain = pose.chain( res );
    Size constraint_chain = pose::get_chain_id_from_chain( chain, pose );

    Size start_res = pose.conformation().chain_begin( constraint_chain );
    Size end_res = pose.conformation().chain_end( constraint_chain );

    for (Size j = start_res ; j < end_res ; ++j ) {
                id::AtomID atom2( pose.residue_type( j ).atom_index( "CA" ), j );
                runtime_assert( target_atom.valid() && atom2.valid() );
                add_individual_constraint( new AtomPairConstraint( target_atom, atom2, func ) );
    }

} // setup_csts

} // constraints
} // scoring
} // core
