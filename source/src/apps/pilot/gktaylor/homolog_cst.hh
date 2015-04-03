// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/RichSequenceMapping.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/MixtureFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/id/AtomID.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>

#include <utility/io/izstream.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/constraints.OptionKeys.gen.hh>


using namespace core;

void add_constraints (
	core::pose::Pose *pose,
	std::string fn
) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;

	static bool init( false );

	static ConstraintSetOP cst_set( new ConstraintSet() );

	if ( !init ) {
		utility::io::izstream data( fn.c_str() );
		std::string line;
		if ( !data ) {
			std::cerr << "ERROR:: Unable to open constraints file: "
								<< fn << '\n';
			std::exit( 1 );
		}

		core::Real max_dist = option[ constraints::max_cst_dist ]();

		getline(data,line); // header line
		while( getline( data, line ) ) {
			// line format:
			// distance
			// gaussian_param
			// exp_param
			// mixture_param
			// bg_anchor
			// bg_sd
			std::istringstream line_stream( line );
			// std::cout << "line = " << line << '\n';
			std::string query_id, template_id;

			int resi, resj;
			float dist, a, b, c, fit_quality;

			line_stream
				>> query_id  >> template_id
				>> resi >> resj
				>> dist
				>> a >> b >> c
				>> fit_quality;

			if ( dist > max_dist ) continue; // skip distances that are too big.

			if ( std::abs( resi - resj ) < 10 ) {
				continue;
			}
			if ( dist > 18 ) {
				continue;
			}

			if ( c < 0 ) {
				c = 0;
			}
			if ( c > 1 ) {
				c = 1;
			}

			cst_set->add_constraint(
				new core::scoring::constraints::AtomPairConstraint(
					core::id::AtomID(2,resi),
					core::id::AtomID(2,resj),
					new core::scoring::constraints::MixtureFunc( dist, a, b, c, 14, 5 )
				)
			);
		} // while ( getline(data,line) )

		init = true;
	} // if ( !init )


	ConstraintSetOP cst_setOP = cst_set->clone();

	pose->constraint_set( cst_setOP );
} // add_constraints

/// @brief reads a .homolog_stats file with information from templates, creates a starting structure
/// with appropriate constraints added. Assumes as input an extended Pose with the proper sequence.
/// Backbone torsions from the template structure are copied into into aligned residues within the provided
/// pose, and sidechain torsions are copied if the amino acids at this position are identical.
void create_starting_template (
	core::pose::Pose & pose,
	std::string fn
) {
	utility::io::izstream data( fn.c_str() );
	std::string line;
	if ( !data ) {
		std::cerr << "ERROR:: Unable to open file: "
							<< fn << '\n';
		std::exit( 1 );
	}

	getline(data,line); // header line
	core::pose::PoseOP template_pose( new core::pose::Pose );
	core::sequence::RichSequenceMapping mapping;
	core::scoring::constraints::ConstraintSetOP cst_set( new core::scoring::constraints::ConstraintSet() );

	std::string query_id, template_id;
	while( getline(data,line) ) { // data line
		std::istringstream line_stream( line );

		int resi, resj, t_resi, t_resj;
		float dist, a, b, c, fit_quality;

		line_stream
			>> query_id  >> template_id
			>> resi >> resj
			>> t_resi >> t_resj
			>> dist
			>> a >> b >> c
			>> fit_quality;

		if ( dist > 8 ) continue; // skip distances that are too big.

		if ( std::abs( resi - resj ) < 10 ) {
			continue;
		}

		if ( c < 0 ) {
			c = 0;
		}
		if ( c > 1 ) {
			c = 1;
		}

		cst_set->add_constraint(
			new core::scoring::constraints::AtomPairConstraint(
				core::id::AtomID(2,resi),
				core::id::AtomID(2,resj),
				new core::scoring::constraints::MixtureFunc( dist, a, b, c, 14, 5 )
			)
		);

		// std::cout << "resi = " << resi << " resj = " << resj << " dist = " << dist << '\n';

		mapping.insert_aligned_residue_safe( resi, t_resi );
		mapping.insert_aligned_residue_safe( resj, t_resj );
	}

	utility::vector1< core::scoring::constraints::ConstraintOP > cst_list = cst_set->get_all_constraints();
	std::cout << "added " << cst_list.size() << " constraints." << '\n';
	pose.constraint_set( cst_set );

	mapping.show();
	core::import_pose::pose_from_pdb( *template_pose, template_id );
	template_pose->dump_pdb("template.pdb");

	// thread the input sequence onto the template_pose
	utility::vector1< core::Real > torsions;
	utility::vector1< core::Real > chis;
	utility::vector1< bool > changed_torsions( pose.total_residue(), false );
	torsions.reserve( 3 );
	chis    .reserve( 4 ); // good guess for upper bound on number of chi angles
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( mapping.size1() < i ) break; // off the end of the mapping!

		Size template_pos = mapping[ i ];
		if ( template_pos != 0 ) {
			// std::cout << "looking at (" << i << "," << template_pos << ")" << '\n';
			// copy backbone torsions
			torsions = template_pose->residue(i).mainchain_torsions();
			pose.set_phi  ( i, torsions[1] );
			pose.set_psi  ( i, torsions[2] );
			pose.set_omega( i, torsions[3] );

			// if ResidueTypes are the same, copy side chain torsions

			// if ( pose.residue(i).name1() == template_pose->residue(template_pos).name1() ) {
			// 				// std::cout << "types are (" << pose.residue(i).name1() << ","
			// 				// 					<< template_pose->residue(template_pos).name1() << ")" << '\n';
			// 				chis = template_pose->residue(i).chi();
			//
			// 				for ( Size chino = 1; chino <= chis.size(); ++chino ) {
			// 					pose.set_chi( i, chino, chis[chino] );
			// 				} // for ( Size i = 1; i <= pose.total_residue(); ++i )
			// 			} // if ( pose.residue(i).name1() == template_pose->residue(template_pos).name1() )
			changed_torsions[ i ] = true;
		} // if ( template_pos != 0 )
	} // for ( Size i = 1; i <= pose.total_residue(); ++i )
}

// burial is represented as number of number of neighbor atoms within 8 angstroms
utility::vector1< int > calculate_burial(
	core::pose::Pose mypose
) {

	utility::vector1< int > burial;
	burial.resize( mypose.total_residue() );

	for ( unsigned int i = 1; i <= mypose.total_residue(); ++i ) {
		for ( unsigned int j = i + 1; j <= mypose.total_residue(); ++j ) {
			core::conformation::Residue resi = mypose.residue(i);
			core::conformation::Residue resj = mypose.residue(j);

			if ( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) < 8 ) {
				burial[ i ]++;
				burial[ j ]++;
			}
		}
	}

	return burial;
}
