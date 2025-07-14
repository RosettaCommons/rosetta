// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IO-functionality for Constraints
/// @brief
/// @author Oliver Lange olange@u.washington.edu

// Unit headers
#include <core/scoring/constraints/ConstraintIO.hh>

// Package headers
#include <core/scoring/constraints/Constraint.hh>
//#include <core/scoring/constraints/ConstraintForest.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
//#include <core/scoring/constraints/BindingSiteConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/util.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

#include <core/scoring/constraints/ConstraintFactory.hh> // AUTO IWYU For ConstraintFactory


// Boost headers

static basic::Tracer tr( "core.scoring.constraints.ConstraintsIO" );

namespace core {
namespace scoring {
namespace constraints {

func::FuncFactory & ConstraintIO::get_func_factory() {
	return get_instance()->func_factory_;
}

ConstraintFactory & ConstraintIO::get_cst_factory() {
	return * ConstraintFactory::get_instance();
}

void
ConstraintIO::read_cst_atom_pairs(
	std::istream & data,
	std::string & next_section_name,
	ConstraintSet & cst_set,
	pose::Pose const & pose,
	bool const force_pdb_info_mapping
) {
	tr.Debug << "ConstraintIO::read_cst_atom_pairs" << std::endl;
	std::string line;
	while ( getline( data, line ) ) {
		Size res1, res2;
		std::string tempres1, tempres2;
		std::string name1, name2;
		std::string func_type;
		std::istringstream line_stream( line );
		line_stream
			>> name1 >> tempres1
			>> name2 >> tempres2
			>> func_type;

		// backwards compatibility with old RNA/DNA atom names:
		name1 = utility::replace_in( name1, "*", "'" );
		name2 = utility::replace_in( name2, "*", "'" );

		tr.Debug << "\tabout to parse residues" << std::endl;
		parse_residue( pose, tempres1, res1, force_pdb_info_mapping );
		parse_residue( pose, tempres2, res2, force_pdb_info_mapping );

		if ( name1.find("[" )!=std::string::npos ) { //end of this section
			tr.Debug << "section end detected in line " << line << std::endl;
			next_section_name = line;
			return;
		}
		tr.Debug  << "read: " << name1 << " " << name2 << " "
			<< res1 << " " << res2 << " func: " << func_type
			<< std::endl;

		if ( res1==0 || res2==0 || res1>pose.size() || res2> pose.size() ||
				( !pose.residue_type( res1 ).has( name1 ) ) ||
				( !pose.residue_type( res2 ).has( name2 ) )
				) {

			tr.Error << "error in constraint (no such atom in pose!)"
				<< name1 << " " << name2 << " " << res1 << " " << res2 << " func: " << func_type << std::endl;
			utility_exit_with_message ("Constraint data referred to atom which is not present in pose");
		}

		id::AtomID atom1( pose.residue_type( res1 ).atom_index( name1 ), res1 );
		id::AtomID atom2( pose.residue_type( res2 ).atom_index( name2 ), res2 );

		func::FuncOP aFunc = get_instance()->func_factory_.new_func( func_type );
		aFunc->read_data( line_stream );

		if ( tr.Debug.visible() ) {
			aFunc->show_definition( tr.Debug ); tr.Debug << std::endl;
		}

		cst_set.add_constraint( utility::pointer::make_shared< AtomPairConstraint >( atom1, atom2, aFunc ) );
	} // while getline
	tr.Debug << "end of file reached" << std::endl;
	next_section_name = "";
}

ConstraintOP ConstraintIO::parse_atom_pair_constraint(
	std::istream & data,
	core::pose::Pose pose,
	bool const force_pdb_info_mapping
) {
	Size res1, res2;
	std::string tempres1, tempres2;
	std::string name1, name2;
	std::string func_type;

	data
		>> name1 >> tempres1
		>> name2 >> tempres2
		>> func_type;

	parse_residue( pose, tempres1, res1, force_pdb_info_mapping );
	parse_residue( pose, tempres2, res2, force_pdb_info_mapping  );

	tr.Info  << "read: " << name1 << " " << name2 << " "
		<< res1 << " " << res2 << " func: " << func_type
		<< std::endl;
	if ( res1>pose.size() || res2> pose.size() ) {
		tr.Warning  << "ignored constraint (no such atom in pose!)"
			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		return nullptr;
	}

	id::AtomID atom1( pose.residue_type( res1 ).atom_index( name1 ), res1 );
	id::AtomID atom2( pose.residue_type( res2 ).atom_index( name2 ), res2 );

	func::FuncOP aFunc = get_instance()->func_factory_.new_func( func_type );
	aFunc->read_data( data );

	if ( tr.Debug.visible() ) {
		aFunc->show_definition( tr.Debug ); tr.Debug<<std::endl;
	}

	ConstraintOP cst_op( utility::pointer::make_shared< AtomPairConstraint >( atom1, atom2, aFunc ) );
	return cst_op;
} // parse_atom_pair_constraint

void
ConstraintIO::read_cst_coordinates(
	std::istream & data,
	std::string & next_section_name,
	ConstraintSet & cst_set,
	pose::Pose const & pose,
	bool const force_pdb_info_mapping
) {

	tr.Debug << "ConstraintIO::read_cst_coordinate" << std::endl;
	std::string line;
	while ( getline( data, line ) ) {

		Real x, y, z;
		Size res1, res2;
		std::string name1, name2;
		std::string func_type;

		std::istringstream line_stream( line );

		line_stream
			>> name1;

		if ( name1.find("[" )!=std::string::npos ) { //end of this section
			tr.Debug << "section end detected in line " << line << std::endl;
			next_section_name = line;
			return;
		}

		line_stream >> res1
			>> name2 >> res2
			>> x >> y >> z
			>> func_type;

		// backwards compatibility with old RNA/DNA atom names:
		name1 = utility::replace_in( name1, "*", "'" );
		name2 = utility::replace_in( name2, "*", "'" );

		res1 = parse_residue( pose, res1, 0/*chain*/, force_pdb_info_mapping );
		res2 = parse_residue( pose, res2, 0/*chain*/, force_pdb_info_mapping );

		tr.Debug  << "read: " << name1 << " " << name2 << " "
			<< res1 << " " << res2 << " func: " << func_type
			<< std::endl;

		if ( res1>pose.size() || res2> pose.size() ||
				( !pose.residue_type( res1 ).has( name1 ) ) ||
				( !pose.residue_type( res2 ).has( name2 ) )
				) {
			tr.Error << "error in constraint (no such atom in pose!)"
				<< name1 << " " << name2 << " " << res1 << " " << res2 << " func: " << func_type << std::endl;
			utility_exit_with_message ("Constraint data referred to atom which is not present in pose");
		}


		id::AtomID atom1( pose.residue_type( res1 ).atom_index( name1 ), res1 );
		id::AtomID atom2( pose.residue_type( res2 ).atom_index( name2 ), res2 );

		func::FuncOP aFunc = get_instance()->func_factory_.new_func( func_type );
		aFunc->read_data( line_stream );

		//  if ( tr.Debug.visible() ) {
		aFunc->show_definition( std::cout ); std::cout<<std::endl;
		//  }

		Vector transform( x, y, z );
		cst_set.add_constraint( utility::pointer::make_shared< CoordinateConstraint >( atom1, atom2, transform, aFunc ) );

	} // while getline
	tr.Debug << "end of file reached" << std::endl;
	next_section_name = "";
}

ConstraintOP ConstraintIO::parse_coordinate_constraint(
	std::istream & data,
	core::pose::Pose pose,
	bool const force_pdb_info_mapping
) {
	Real x, y, z;
	Size fixed_res, other_res;
	std::string tempfixed_res, tempother_res;
	std::string fixed_res_name, other_res_name;
	std::string func_type;

	data
		>> fixed_res_name >> tempfixed_res
		>> other_res_name >> tempother_res
		>> x >> y >> z
		>> func_type;

	parse_residue( pose, tempfixed_res, fixed_res, force_pdb_info_mapping );
	parse_residue( pose, tempother_res, other_res, force_pdb_info_mapping );

	tr.Debug  << "read: " << fixed_res_name << " " << other_res_name << " "
		<< fixed_res << " " << other_res << " func: " << func_type << std::endl;
	if ( fixed_res > pose.size() || other_res > pose.size() ) {
		tr.Warning  << "ignored constraint (no such atom in pose!)"
			<< fixed_res_name << " " << other_res_name << " "
			<< fixed_res << " " << other_res << std::endl;
		return nullptr;
	}

	id::AtomID atom1( pose.residue_type( fixed_res ).atom_index( fixed_res_name ), fixed_res );
	id::AtomID atom2( pose.residue_type( other_res ).atom_index( other_res_name ), other_res );

	func::FuncOP aFunc = get_instance()->func_factory_.new_func( func_type );
	aFunc->read_data( data );

	if ( tr.Debug.visible() ) {
		aFunc->show_definition( tr.Debug ); tr.Debug << std::endl;
	}

	Vector transform( x, y, z );
	ConstraintOP cst_op( utility::pointer::make_shared< CoordinateConstraint >( atom1, atom2, transform, aFunc ) );
	return cst_op;
} // parse_coordinate_constraint


void
ConstraintIO::read_cst_angles(
	std::istream & data,
	std::string & next_section_name,
	ConstraintSet & cst_set,
	pose::Pose const & pose,
	bool const force_pdb_info_mapping
) {
	tr.Debug << "ConstraintIO::read_cst_angles" << std::endl;
	std::string line;
	while ( getline( data, line ) ) {
		Size res1, res2, res3;
		std::string tempres1, tempres2, tempres3;
		std::string name1, name2, name3;
		std::string func_type;
		std::istringstream line_stream( line );
		line_stream
			>> name1 >> tempres1
			>> name2 >> tempres2
			>> name3 >> tempres3
			>> func_type;

		parse_residue( pose, tempres1, res1, force_pdb_info_mapping );
		parse_residue( pose, tempres2, res2, force_pdb_info_mapping );
		parse_residue( pose, tempres3, res3, force_pdb_info_mapping );

		if ( name1.find("[" )!=std::string::npos ) { //end of this section
			tr.Debug << "section end detected in line " << line << std::endl;
			next_section_name = line;
			return;
		}
		id::AtomID atom1( pose.residue_type( res1 ).atom_index( name1 ), res1 );
		id::AtomID atom2( pose.residue_type( res2 ).atom_index( name2 ), res2 );
		id::AtomID atom3( pose.residue_type( res3 ).atom_index( name3 ), res3 );

		tr.Debug << "read: " << name1 << " " << res1 << " "  <<  name2 << " "
			<< res2 << " " << name3 << " " << res3 << std::endl;
		func::FuncOP aFunc = get_instance()->func_factory_.new_func( func_type );
		aFunc->read_data( line_stream );

		if ( tr.Debug.visible() ) {
			aFunc->show_definition( tr.Debug ); tr.Debug<<std::endl;
		}

		if ( res1 > pose.size() || res2 > pose.size() || res3 > pose.size() ) {
			tr.Warning << "ignored constraint" << name1 << " " << name2 << " "
				<< res1 << " " << res2 << std::endl;
			continue;
		}

		cst_set.add_constraint( utility::pointer::make_shared< AngleConstraint >( atom1, atom2, atom3, aFunc ) );
	}// while getline
	tr.Debug << "end of file reached" << std::endl;
	next_section_name = "";
}

std::string
get_section_name ( std::string line ) {
	if ( line.size() == 0 ) return line;
	std::istringstream line_stream( line );
	std::string tok;
	line_stream >> tok;

	if ( tok == "[" ) {
		line_stream >> tok;
	} else {
		std::string::size_type start = tok.find("[");
		if ( start != 0 ) return "NO_SECTION";
	}

	std::string::size_type loc = tok.find("]");
	if ( loc != std::string::npos ) {
		return tok;
	} else {
		int start = 0;
		return tok.substr(start,loc);
	}
}

void
ConstraintIO::read_cst_bindingsites(
	std::istream & data,
	std::string & next_section_name,
	ConstraintSet & cst_set,
	pose::Pose const & pose,
	bool const //force_pdb_info_mapping // = false
) {
	tr.Debug << "ConstraintIO::read_cst_angles" << std::endl;
	std::string line;
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
		if ( line.find("[" )!=std::string::npos ) { //end of this section
			tr.Debug << "section end detected in line " << line << std::endl;
			next_section_name = line;
			return;
		}

		ConstraintOP bsc = ConstraintFactory::get_instance()->newConstraint( "BindingSite" );
		bsc->read_def( line_stream, pose, get_func_factory() );
		cst_set.add_constraint( bsc );
	} // while getline
	tr.Debug << "end of file reached" << std::endl;
	next_section_name = "";
}

ConstraintSetOP ConstraintIO::read_constraints(
	std::istream & data,
	ConstraintSetOP cset,
	pose::Pose const& pose,
	bool const force_pdb_info_mapping// = false
) {
	std::string line;
	std::streampos original_pos = data.tellg();
	while ( line.size() == 0 ) getline(data,line); // header line
	std::string section = get_section_name( line );
	std::string pre_read;
	while ( section.size() ) {
		tr.Info << "read constraints section --" << section << "---" << std::endl;
		if ( section ==  "atompairs" ) {
			read_cst_atom_pairs( data, pre_read, *cset, pose, force_pdb_info_mapping );
		} else if ( section == "coordinates" ) {
			read_cst_coordinates( data, pre_read, *cset, pose, force_pdb_info_mapping );
		} else if ( section == "angles" ) {
			read_cst_angles( data, pre_read, *cset, pose, force_pdb_info_mapping );
		} else if ( section == "bindingsites" ) {
			read_cst_bindingsites( data, pre_read, *cset, pose, force_pdb_info_mapping );
		} else if ( section == "NO_SECTION" ) {
			tr.Info << " no section header [ xxx ] found, try reading line-based format... DON'T MIX"
				<< std::endl;
			/// izstreams cannot be rewound with seekg()
			/// never call this function if you have constructed an izstream
			/// and you haven't deteremined that indeed the file format its representing
			/// is the old style (as opposed to the new style) constraint format.
			debug_assert( dynamic_cast< zlib_stream::zip_istream * > ( &data ) == nullptr );

			data.seekg( original_pos );
			return read_constraints_new( data, cset, pose, force_pdb_info_mapping );
		} else { //section header, but unknown name
			utility_exit_with_message(
				"constraint-file: section " + section + " not recognized!"
			);
		}
		tr.Trace << "pre_read: " << pre_read << std::endl;
		section = get_section_name( pre_read );
	}
	// pose.constraint_set( cset );
	return cset;
}

/// @details All the heavy lifting is done by read_constraints( istream &, ConstraintSetOP, Pose const & ), or
/// by read_constraints_new( istream &, ConstraintSetOP, Pose const & ), but the logic for deciding which of
/// two execution paths to follow that lives inside read_constraints( isteam &, ... ) will not work if given
/// an izstream constructed from a zipped file.  SO instead, we check the file format of the input constraint
/// file here and then rewind to the beginning of the file using the izstream seek_beg() function.
ConstraintSetOP
ConstraintIO::read_constraints(
	std::string const & fname,
	ConstraintSetOP cset,
	pose::Pose const & pose,
	bool const force_pdb_info_mapping // = false
) {
	utility::io::izstream data( fname.c_str() );
	tr.Info << "read constraints from " << fname << std::endl;
	if ( !data ) {
		utility_exit_with_message( "[ERROR] Unable to open constraints file: "+ fname );
	}

	std::string line;
	getline(data,line); // header line
	std::string section = get_section_name( line );
	data.seek_beg();
	if ( section == "NO_SECTION" ) {
		return read_constraints_new( data, cset, pose, force_pdb_info_mapping );
	} else {
		return read_constraints( data, cset, pose, force_pdb_info_mapping );
	}

	return read_constraints( data, cset, pose, force_pdb_info_mapping );
} // read_constraints

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details read in constraints with new format from a file
ConstraintSetOP
ConstraintIO::read_constraints_new(
	std::string const & fname,
	ConstraintSetOP cset,
	pose::Pose const & pose,
	bool const force_pdb_info_mapping
) {
	utility::io::izstream data( fname.c_str() );
	tr.Info << "read constraints from " << fname << std::endl;
	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open constraints file: " + fname );
	}
	return read_constraints_new( data, cset, pose, force_pdb_info_mapping );
}

ConstraintSetOP
ConstraintIO::read_constraints_new(
	std::istream & data,
	ConstraintSetOP cset,
	pose::Pose const & pose,
	bool const force_pdb_info_mapping
) {
	Size count_constraints(0);
	while ( data.good() ) { // check if we reach the end of file or not
		// read in each constraint and add it constraint_set
		ConstraintOP cst_op( read_individual_constraint_new( data, pose, get_func_factory(), force_pdb_info_mapping ) );
		if ( cst_op != nullptr ) {
			++count_constraints;
			cset->add_constraint( cst_op );
		} else if ( !data.eof() ) { // not end of line
			tr.Error << "reading constraints from file" << std::endl;
			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			if ( option[ OptionKeys::constraints::exit_on_bad_read ]() ) {
				utility_exit_with_message( "ERROR: reading constraints from file"  );
			}
			break;
		}
	} // while
	tr.Info << "Read in " << count_constraints << " constraints" << std::endl;
	return cset;
} // read_constraints_new

ConstraintOP
ConstraintIO::read_individual_constraint_new(
	std::istream & data,
	pose::Pose const& pose,
	func::FuncFactory const & func_factory,
	std::string tag,
	bool const //force_pdb_info_mapping // = false
)
{
	ConstraintOP cst_op;
	cst_op = get_cst_factory().newConstraint( tag );

	std::string error_msg("");
	// bool error_seen( false );

	if ( cst_op != nullptr ) {
		//  try {
		cst_op->read_def( data, pose, func_factory );
		// } catch (utility::excn::Exception &excn  ) {
		//    tr.Error << "reading of " + tag + " failed.\n" << excn << std::endl;
		//    cst_op = NULL;
		//}
		if ( !data.good() && !data.eof() ) {
			error_msg += "reading of " + tag + " failed.\n";
			tr.Error << error_msg << std::endl;
			cst_op = nullptr;
		}
	} else {
		error_msg += "constraint type " + tag + " not known.\n";
		tr.Error << error_msg << std::endl;
		cst_op = nullptr;
	}
	if ( !cst_op ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if ( option[ OptionKeys::constraints::exit_on_bad_read ]() ) { //OL Aug 12 2010,
			//changed option default from 'true' to 'false'. This was way confusing!
			// hard exits can be rather annoying... after waiting for 24h in the queue of a cluster.
			//changed back to 'true' -- found several 'silent' cases where I wish
			// there had been a hard exit-- rhiju
			utility_exit_with_message( "[ ERROR ]" + error_msg );
		}
	}
	return cst_op;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details read in each individual constraint. The constraints could be single-line ones such as AtomPair,
/// Angle and Diehedral, or multi-line ones such as Multiconstraint and Ambiguous. Return owing pointer of the
/// constraint if read in successfully, otherwise return NULL. Identify the type of constraint and call
/// each constraint's read_def function to finish reading. Skip lines beginning with '#' or '\n'
/// May be called recursively in the case of Multiconstraint and AmbiguousConstraint.
/// @note the istream data should point to the beginning of a line when this function is called.
ConstraintOP
ConstraintIO::read_individual_constraint_new(
	std::istream & data,
	pose::Pose const& pose,
	func::FuncFactory const & func_factory,
	bool const force_pdb_info_mapping // = false
)
{
	std::string tag;
	// get the ConstraintType tag
	while ( true ) {
		char const c = data.peek(); // get first char of the line but not moving istream pointer
		if ( c == '#' || c == '\n' ) { //ignore # comment line and empty line
			while ( data.good() && (data.get() != '\n') ) {}
			continue;
		}
		if ( data.eof() ) return nullptr;
		data >> tag;
		if ( data.fail() ) {
			tr.Error << "can't read constraint type" << std::endl;
			return nullptr;
		}
		break;
	}
	if ( ( tag.substr(0,3) == "END" )||( tag.substr(0,3) == "End" ) ) return nullptr; // stopper for MultiConstraint or AmbiguousConstraint
	return read_individual_constraint_new( data, pose, func_factory, tag, force_pdb_info_mapping );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void
ConstraintIO::write_constraints( std::ostream& out, ConstraintSet const& cst_set, pose::Pose const& pose ) {
	cst_set.show_definition( out, pose );
}

void
ConstraintIO::write_constraints( std::string const& filename, ConstraintSet const& cst_set, pose::Pose const& pose ) {
	utility::io::ozstream dump_cst( filename );
	write_constraints( dump_cst, cst_set, pose );
}

void
ConstraintIO::parse_residue( pose::Pose const& pose, std::string const & residue_string, Size & residue_num, bool const force_pdb_info_mapping )
{
	std::stringstream data;
	std::string chain;
	int resnum;

	data.str( residue_string );

	data >> resnum;

	if ( (data >> chain).fail() ) chain = "";

	residue_num = parse_residue( pose, resnum, chain, force_pdb_info_mapping );
}


Size
ConstraintIO::parse_residue( pose::Pose const& pose, int const resnum, std::string const & chain /* = i"" */, bool const force_pdb_info_mapping_in )
{
	// this option is a vector1< bool > for pretty arcane reasons -- rosetta does not provide a set default option for bool, but does so for vector< bool>.
	using namespace basic::options;
	using namespace core::pose::full_model_info;
	// Presence of FullModelInfo overrides... in which case, this might indicate a residue that's not yet in the pose.
	if ( /*chain != 0 &&*/ full_model_info_defined( pose ) ) {
		// AMW TODO: presence of segid in constraints here.
		auto const & foo = const_full_model_info( pose ).full_to_sub();
		if ( foo.find( const_full_model_info( pose ).full_model_parameters()->conventional_to_full( resnum, chain, "    " ) ) != foo.end() ) {
			//tr << "I am mapping apparent resnum " << resnum << " chain \'" << chain << "\' to "
			// << const_full_model_info( pose ).full_to_sub()[ const_full_model_info( pose ).full_model_parameters()->conventional_to_full( resnum, chain, "    " ) ] << std::endl;
			//tr << const_full_model_info( pose ).full_to_sub() << std::endl;
			return foo.at( const_full_model_info( pose ).full_model_parameters()->conventional_to_full( resnum, chain, "    " ) );
		}
	}
	if ( core::conformation::is_chain_valid(chain) && pose.pdb_info() ) {
		// Enter this option if user provided a PDB chain ID and a residue number in the cst file
		Size resnum_out( pose.pdb_info()->pdb2pose( chain, resnum ) );
		// If resnum_out == 0, then user did not specify a Pose residue that (currently) exists
		// Inform user of the error via a Warning
		// Note: a constraint to a residue that does not exist may exist later during Stepwise RNA protocol
		// Each Constraint method should then Warn user that it cannot add the specified constraint
		if ( resnum_out == 0 ) {
			tr.Warning << "Residue specified by constraint file does not exist. "
				"Does residue number " + utility::to_string(resnum) +
				" chain " + chain + " exist in both your Pose and constraint file?" << std::endl;
		}
		return resnum_out;
	}
	bool force_pdb_info_mapping = pose.pdb_info() &&
		( force_pdb_info_mapping_in || ( option[ OptionKeys::constraints::force_pdb_info_mapping ]().size() ?
		option[ OptionKeys::constraints::force_pdb_info_mapping ]()[1] : false ) );
	if ( force_pdb_info_mapping && pose.pdb_info() ) {
		Size resnum_out = pose.pdb_info()->pdb2pose( "A", resnum );
		if ( resnum_out > 0 ) return resnum_out;
		// some legacy PDB's have ' ' instead of 'A' remaining as default for chains...
		return pose.pdb_info()->pdb2pose( " ", resnum );
	}
	// Return the input resnum if user provided a Pose residue number in the cst file
	// If resnum > pose.size(), then user did not specify a Pose residue that (currently) exists
	// Inform user of the error via a Warning
	// Note: a constraint to a residue that does not exist may exist later during Stepwise RNA protocol
	// Each Constraint method should then Warn user that it cannot add the specified constraint
	if ( Size(resnum) > pose.size() ) {
		tr.Warning << "Residue specified by constraint file does not exist. "
			"Does Pose residue number " + utility::to_string(resnum) +
			" exist in both your Pose and constraint file?" << std::endl;
	}
	// However, if user specified to constrain to Pose residue 0, exit as that does not make sense
	if ( Size(resnum) == 0 ) {
		utility_exit_with_message("Constraint file specified to constrain using Pose residue 0."
			" That cannot be possible! Rosetta numbering starts at 1");
	}

	return Size( resnum );
}

} //constraints
} //scoring
} //core
