// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/SiteConstraintResidues.cc
/// @brief This class is an AmbiguousConstraint in which the set is comprised of AtomPairConstraints
/// @brief of an atom of interest in one chain versus the CA of of another sets of residues
/// @author Lei Shi (shilei@uw.edu)

#include <core/scoring/constraints/SiteConstraintResidues.hh>
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

static thread_local basic::Tracer TR( "core.scoring.constraints.SiteConstraintResidues" );

namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
SiteConstraintResidues::SiteConstraintResidues():
AmbiguousConstraint()
{}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
SiteConstraintResidues::SiteConstraintResidues( ConstraintCOPs & cst_in ):
AmbiguousConstraint( cst_in )
{}

void
SiteConstraintResidues::show( std::ostream& out) const
{
	/// APL -- you cannot show the active constraint in the absence of a Pose.
	//out << "AmbiguousConstraint Active constraint:" << std::endl;
	//active_constraint()->show(out);
	out << "SiteConstraintResidues is an AmbiguousConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
	for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
		(*cst_it)->show(out);
	}

	out << " ...all member constraints of this SiteConstraintResidues shown." << std::endl;
}


void
SiteConstraintResidues::read_def(
    std::istream & data,
    core::pose::Pose const & pose,
   func::FuncFactory const & func_factory
) {
    TR.Debug << "read_site_cst" << std::endl;
    std::string tempres;
    Size res1;
    std::string name;
    Size res2;
    Size res3;
    std::string func_type;
    std::string type;
    data >> res1 >> name >> res2 >> res3 >> func_type;
    TR.Info << "read: " << res1 << " "<< name << " constrain to residues " << res2 << ":" << res3  << " func: " << func_type << std::endl;
    func::FuncOP aFunc = func_factory.new_func( func_type );
    aFunc->read_data( data );

    if ( TR.Debug.visible() ) {
        aFunc->show_definition( TR.Debug ); TR.Debug<<std::endl;
    }
    setup_csts( res1, name, res2,res3, pose, aFunc );

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
SiteConstraintResidues::setup_csts(
    Size res1,
    std::string name,
    Size res2,
    Size res3,
    core::pose::Pose const & pose,
    func::FuncOP const & func
) {
    id::AtomID target_atom( pose.residue_type( res1 ).atom_index( name ), res1 );
    //Size target_chain = pose.chain( res );
    for (Size j = res2 ; j < res3 ; ++j ) {
                id::AtomID atom2( pose.residue_type( j ).atom_index( "CA" ), j );
                runtime_assert( target_atom.valid() && atom2.valid() );
                add_individual_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( target_atom, atom2, func ) ) ) );
    }

} // setup_csts

} // constraints
} // scoring
} // core
