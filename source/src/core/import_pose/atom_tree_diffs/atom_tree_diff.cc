// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/silent/atom_tree_diff.cc
///
/// @brief  Silent-file format based on "diffs" of AtomTree DOFs
/// @author Ian W. Davis


#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <algorithm>
#include <set>
#include <sstream>

namespace core {
namespace import_pose {
namespace atom_tree_diffs {


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Provides a StrictWeakOrdering comparator for sorting elements of a ScoresMap by one particular score type.
/// @details Obviously, all entries in the map must have that score type present, or this dies a horrible death.
class ScoreLessThanComparator {
public:
	typedef ScoresPairList::value_type value_type;
	ScoreLessThanComparator(std::string const & score_name, bool reverse=false): score_name_(score_name), reverse_(reverse) {}
	~ScoreLessThanComparator() {}
	bool operator() (value_type a, value_type b)
	{ return reverse_ ^ (a.second[score_name_] < b.second[score_name_]); }
private:
	std::string score_name_;
	bool reverse_;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
AtomTreeDiff::AtomTreeDiff():
	scores_(),
	tag_idx_(),
	offsets_(),
	ref_poses_(),
	in_(),
	file_read_(false)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
AtomTreeDiff::AtomTreeDiff(std::string filename):
	scores_(),
	tag_idx_(),
	offsets_(),
	ref_poses_(),
	in_(),
	file_read_(false)
{
	read_file(filename);
}

void AtomTreeDiff::read_file(std::string filename){
	debug_assert( ! file_read_ ); // only read file once, by constructor or elsewhere. can't read multiple files yet
	in_.open(filename.c_str());
	core::pose::PoseOP ref_pose;
	std::set<std::string> used_tags;
	while ( true ) {
		std::string tag;
		std::map< std::string, core::Real > scores;
		// in_ is now positioned just after the SCORES line of the *previous* structure
		// If we use this position, we must call header_from_atom_tree_diff() when re-reading.
		long curr_pos = in_.tellg();
		if ( ! header_from_atom_tree_diff(in_, tag, scores) ) break;
		// in_ is now positioned just after the SCORES line of *this* structure
		// If we use this position, we must NOT call header_from_atom_tree_diff() when re-reading.
		//long curr_pos = in_.tellg();

		if ( used_tags.count(tag) > 0 ) {
			basic::Error() << "Tag " << tag << " appears at least twice in the atom_tree_diff file!  Discarding structure..." << std::endl;
			continue;
		}
		core::Size const end= tag.find_last_of('_');
		if ( end == std::string::npos ) {
			utility_exit_with_message(tag+" doesn't end with a 4-digit code");
		}
		std::string ref_tag= tag.substr(0, end);
		if ( scores.find("is_reference_pose") != scores.end() ) {
			core::pose::Pose empty_pose;
			ref_pose = core::pose::PoseOP( new core::pose::Pose() );
			core::import_pose::atom_tree_diffs::pose_from_atom_tree_diff(in_, empty_pose, *ref_pose);
			unique_ref_poses_.push_back(ref_pose);
			core::Size const start= ref_tag.find_first_of('_');
			ref_tag= ref_tag.substr(start+1);
			if ( ref_tags_.find(ref_tag) != ref_tags_.end() ) {
				basic::Error() << "Tag " << ref_tag << " appears at least twice in the atom_tree_diff file!  Discarding reference structure..." << std::endl;
				continue;
			}
			ref_tags_[ref_tag]= 0;
			ref_poses_[ref_tag] = ref_pose;
		} else { // not a reference pose, actual data
			scores_.push_back( std::make_pair(tag, scores) );
			tag_idx_[tag] = scores_.size();
			offsets_[tag] = curr_pos;
			++ref_tags_[ref_tag];
			ref_poses_[tag] = ref_pose;
			// Missing ref_pose will fail hard when trying to read_pose(), so no worries.
			// Some silent files will deliberately omit a reference pose to reduce size.
			//if( ref_pose() == NULL ) basic::Error() << "Tag " << tag << " has no reference pose!" << std::endl;
		}
	}
	// Once we've hit EOF, we need to reopen the file before we can successfully seek again.
	in_.close();
	in_.clear(); // this is required for Linux, though not OS X
	in_.open(filename.c_str());
	file_read_= true;
}

AtomTreeDiff::~AtomTreeDiff()
{
	in_.close();
}

bool AtomTreeDiff::has_ref_pose(std::string const & tag) const{
	return  ref_poses_.find(tag) != ref_poses_.end() ;
}

bool AtomTreeDiff::has_tag(std::string const & tag) const{
	return tag_idx_.find(tag) != tag_idx_.end();
}

bool AtomTreeDiff::has_ref_tag(std::string const & tag) const{
	return ref_poses_.find(tag) != ref_poses_.end();
}

//static
void AtomTreeDiff::sort_by(
	std::string const & score_name,
	ScoresPairList & scores,
	bool descending /*=false*/
)
{
	std::stable_sort(scores.begin(), scores.end(), ScoreLessThanComparator(score_name, descending));
}


void AtomTreeDiff::read_pose(std::string const & tag, core::pose::Pose & pose_out)
{
	if ( ref_poses_.find(tag) == ref_poses_.end() || ref_poses_[tag].get() == NULL ) {
		utility_exit_with_message("No reference pose available for "+tag);
	}
	read_pose( tag, pose_out, *(ref_poses_[tag]) );
}

void AtomTreeDiff::read_pose(std::string const & tag, core::pose::Pose & pose_out, core::pose::Pose const & ref_pose)
{
	//std::cout << "Seeking to " << offsets_[tag] << std::endl;
	in_.seekg( offsets_[tag] );
	//std::cout << "Now at " << in_.tellg() << std::endl;
	if ( offsets_[tag] != (long)in_.tellg() ) {
		utility_exit_with_message("Unable to seek in input file!");
	}

	std::string tag_reread;
	Scores scores_reread;
	if ( !core::import_pose::atom_tree_diffs::header_from_atom_tree_diff(in_, tag_reread, scores_reread) || tag != tag_reread ) {
		utility_exit_with_message("Seek took us to the wrong entry!");
	}
	pose_from_atom_tree_diff( in_, ref_pose, pose_out );
}


// For reasons I don't understand, GCC will not compile this if I mark this "const":
core::pose::PoseCOP AtomTreeDiff::ref_pose_for(std::string const & tag)
{
	if ( ref_poses_.find(tag) == ref_poses_.end() || ref_poses_[tag].get() == NULL ) {
		utility_exit_with_message("No reference pose available for "+tag);
	}
	return ref_poses_[tag];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dump_score_line(
	std::ostream & out,
	std::string const & pose_tag,
	std::map< std::string, core::Real > const & scores
)
{
	// Repeating the pose tag in the score line is a great convenience
	// for grepping out score data and later correlating it with a model.
	out << "SCORES " << pose_tag;
	for ( auto const & pair : scores ) {
		// Scores that are very near zero have artificially high precision
		// when rendered in scientific notation, leading to numerical instability in benchmarks.
		// Stupid C++ doesn't seem to have a fixed-point output mode
		// that doesn't also append extra (unneeded) zeros to the end of every float.
		core::Real val = pair.second;
		if ( std::abs(val) < 1e-8 ) val = 0.0;
		out << ' ' << pair.first << ' ' << val;
	}
	out << '\n';
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Writes out the given pose in PDB format (plus fold tree).
/// When reading back in, this pose will be used as the reference point for all
/// subsequent poses, until another reference pose is reached.
void dump_reference_pose(
	std::ostream & out,
	std::string const & pose_tag,
	Scores const & scores,
	core::pose::Pose const & pose
)
{
	out << "POSE_TAG " << pose_tag << '\n';

	Scores mod_scores = scores; // a mutable copy
	mod_scores["is_reference_pose"] = 1;
	dump_score_line(out, pose_tag, mod_scores);

	out << "BEGIN_PDB_FORMAT " << pose_tag << '\n';
	io::pdb::dump_pdb( pose, out );
	out << "END_PDB_FORMAT " << pose_tag << '\n';

	// Some variants introduce extra DOFs but do not change the three-letter name.
	// Useless MUTATE entries will be ignored when reading in.
	for ( Size rsd = 1, rsd_end = pose.size(); rsd <= rsd_end; ++rsd ) {
		out << "MUTATE " << rsd << " " << pose.residue_type(rsd).name() << "\n";
	}

	core::kinematics::FoldTree const & foldtree = pose.fold_tree();

	//basic::OWrapperStream ws(out);
	out << foldtree << std::endl; // writes FOLD_TREE

	out << "END_POSE_TAG " << pose_tag << '\n';
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Writes out the given pose in a atom_tree_diff-type format by only
/// recording AtomTree DOFs that differ from ref_pose.
/// Thus pose and ref_pose should have the same sequence and AtomTree topology.
/// The entry is labeled with the given textual tag and 0+ named numeric scores.
/// The number of digits of precision after the decimal point can be specified
/// independently for backbone and sidechain.
/// More precision reduces propagated coordinate error but increases both the
/// number of characters used per value and the total number of values written.
/// The default values should be very good for most applications.
void dump_atom_tree_diff(
	std::ostream & out,
	std::string const & pose_tag,
	std::map< std::string, core::Real > const & scores,
	core::pose::Pose const & ref_pose_in,
	core::pose::Pose const & pose,
	int bb_precision /*= 6*/,
	int sc_precision /*= 4*/,
	int bondlen_precision /*= 2*/
)
{
	using core::Size;
	using core::Real;
	using basic::Warning;
	using namespace core::id;
	using namespace core::scoring;

	if ( bb_precision < sc_precision ) Warning() << "Crazy fool, bb_precision should be >= sc_precision!" << std::endl;
	Real bb_tol = 1.0, sc_tol = 1.0, bondlen_tol = 1.0;
	for ( int i = 0; i < bb_precision; ++i ) bb_tol /= 10.0;
	for ( int i = 0; i < sc_precision; ++i ) sc_tol /= 10.0;
	for ( int i = 0; i < bondlen_precision; ++i ) bondlen_tol /= 10.0;

	core::kinematics::AtomTree const & atom_tree = pose.atom_tree();
	core::kinematics::FoldTree const & foldtree = pose.fold_tree();

	out << "POSE_TAG " << pose_tag << '\n';
	dump_score_line(out, pose_tag, scores);

	// Also diff sequence, by transforming copy of ref_pose
	// to match seq. of pose, then computing atom_tree diffs.
	core::pose::Pose ref_pose(ref_pose_in);
	for ( Size rsd = 1, rsd_end = pose.size(); rsd <= rsd_end; ++rsd ) {
		if ( pose.residue_type(rsd).name() == ref_pose.residue_type(rsd).name() ) continue;
		using namespace core::conformation;
		ResidueOP newres = ResidueFactory::create_residue(pose.residue_type(rsd), ref_pose.residue(rsd), ref_pose.conformation());
		ref_pose.replace_residue(rsd, *newres, true /*orient backbone*/);
		out << "MUTATE " << rsd << " " << pose.residue_type(rsd).name() << "\n";
	}

	// Needed for custom foldtrees in docking.  Hasn't been tested yet...
	// Probably want to set the fold tree before setting any DOFs, so we put it near the top.
	out << foldtree << std::endl; // writes FOLD_TREE

	// Save stream state
	// We use fixed # of digits after the decimal because phi and theta
	// are in radians, so all DOF values fall within the range [-4,4].
	std::ios_base::fmtflags orig_flags = out.flags();
	std::streamsize orig_precision = out.precision();
	out.setf( std::ios_base::fixed );

	// DOFs for bonded atoms
	for ( Size rsd = 1, rsd_end = pose.size(); rsd <= rsd_end; ++rsd ) {
		bool const is_jump_residue = foldtree.is_jump_point(rsd);
		for ( Size atom = 1, atom_end = pose.residue(rsd).natoms(); atom <= atom_end; ++atom ) {
			AtomID aid(atom, rsd);
			DOF_ID dof_phi(aid, PHI), dof_theta(aid, THETA), dof_d(aid, D);
			if ( atom_tree.atom(aid).is_jump() ) {
				// Jump atoms have the rbN values all set to zero -- maybe vestigal part of atom tree?
			} else {
				// Backbone heavyatoms have to be very precise, or else a lot of error
				// can propagate to the end of the chain.
				// Sidechain atoms and hydrogens don't have to worry nearly so much.
				// (Except for jump residues, where the stubs on either end must be very precise!)
				bool is_backbone = is_jump_residue // i.e. things should be "backbone" if they're anchors for a jump!
					|| (pose.residue(rsd).atom_type(atom).is_heavyatom() && atom <= pose.residue(rsd).last_backbone_atom());
				int precision = (is_backbone ? bb_precision : sc_precision);
				Real tol = (is_backbone ? bb_tol : sc_tol);

				Real before_phi = ref_pose.dof(dof_phi), before_theta = ref_pose.dof(dof_theta), before_d = ref_pose.dof(dof_d);
				Real after_phi = pose.dof(dof_phi), after_theta = pose.dof(dof_theta), after_d = pose.dof(dof_d);
				// Format: resno atomno phi [theta [d]], only for atoms with changed DOFs
				// Phi comes first because dihedrals change most often.
				// Theta comes next because in most cases bond angles don't change, so we can omit it!
				// D comes last and requires least precision, because it doesn't control a lever-arm effect.
				bool const changed_phi = (std::abs(numeric::nearest_angle_radians(after_phi, before_phi) - before_phi) > tol); // otherwise get bogus diffs near pi
				bool const changed_theta = (std::abs(after_theta - before_theta) > tol); // rarely near pi, so I haven't put in the check
				bool const changed_d = (std::abs(after_d - before_d) > bondlen_tol);
				if ( changed_phi || changed_theta || changed_d ) {
					out.precision( precision );
					out << rsd << ' ' << atom;
					//out << ' ' << pose.residue(rsd).name() << ' ' << pose.residue(rsd).atom_name(atom); // debugging
					out << ' ' << after_phi;
					if ( changed_theta || changed_d ) {
						out << ' ' << after_theta;
						if ( changed_d ) {
							out.precision( bondlen_precision );
							out << ' ' << after_d;
						}
					}
					out << '\n';
				}
			}// end non-jump atom
		}// end loop over atoms
	}// end loop over residues

	// DOFs for jumps
	for ( int jump = 1, jump_end = pose.num_jump(); jump <= jump_end; ++jump ) {
		out << "JUMP " << jump;
		out.precision( 12 ); // matrix needs to be very exact
		for ( int i = 1; i <= 3; ++i ) {
			for ( int j = 1; j <= 3; ++j ) {
				out << ' ' << pose.jump(jump).get_rotation()(i,j);
			}
		}
		out.precision( 4 ); // distance needs one more decimal than PDB coords, max
		for ( int i = 1; i <= 3; ++i ) {
			out << ' ' << pose.jump(jump).get_translation()(i);
		}
		out << '\n';
	}

	// Explicit END_POSE_TAG is a good safeguard against file corruption, where
	// a process gets killed halfway through outputting a pose.
	// The missing END_POSE_TAG line(s) and/or POSE_TAG occuring in the middle of a line
	// are signs of this kind of file corruption (should be repairable by hand).
	out << "END_POSE_TAG " << pose_tag << '\n';

	// Restore stream state
	out.flags( orig_flags );
	out.precision( orig_precision );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Reads from the current stream position until it finds both a
/// pose_tag line and a scores line, then extracts the data and returns it.
/// Returns true on success and false on failure.
bool header_from_atom_tree_diff(
	std::istream & in,
	std::string & pose_tag_out,
	Scores & scores_out
)
{
	scores_out.clear();
	bool found_tag = false, found_scores = false;
	while ( in.good() && !( found_tag && found_scores ) ) {
		std::string line, key;
		core::Real value;
		getline(in, line);
		if ( in.fail() ) return false;
		std::istringstream is(line);
		is >> key;
		if ( key == "POSE_TAG" ) {
			found_tag = true;
			is >> pose_tag_out;
		} else if ( key == "SCORES" ) {
			is >> key; // discard redundant pose_tag field
			found_scores = true;
			while ( is.good() ) {
				is >> key >> value;
				if ( is.fail() ) break;
				scores_out[ key ] = value;
			}
		}
		// else just ignore and discard this line
	}
	return ( found_tag && found_scores );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Extracts a pose-diff and sets the DOFs of "pose" accordingly,
/// working from the current stream position.  Use header_from_atom_tree_diff()
/// to locate the pose you wish to extract prior to calling this function.
/// Returns false if there is some sort of IO or format error.
bool pose_from_atom_tree_diff(
	std::istream & in,
	core::pose::Pose const & ref_pose,
	core::pose::Pose & pose
)
{
	using core::kinematics::Jump;
	using namespace core::id;
	using namespace core::scoring;

	// start with pose as a copy of ref_pose:
	pose = ref_pose; // deep copy

	basic::Tracer TR( "core.import_pose.atom_tree_diffs.atom_tree_diff.pose_from_atom_tree_diff" );
	while ( in.good() ) {
		std::string line, key;
		getline(in, line);
		if ( in.fail() ) {
			TR << "getline() failed" << std::endl;
			return false;
		}
		//std::cout << "    " << line << std::endl; // for debugging
		std::istringstream is(line);
		// Start by assuming it will be an atom line, then backtrack if not
		core::Size rsd_no, atom_no;
		is >> rsd_no >> atom_no;
		if ( is.fail() ) {
			//TR << "oops, not an atom line!" << std::endl;
			//TR << ">> " << line << std::endl;
			is.clear();
			is.seekg(0, std::ios::beg);
			is >> key;
			if ( key == "END_POSE_TAG" ) return true;
			// I changed to uppercase for format consistency on 6 Dec 2007.
			// Check for the deprecated lowercase version could eventually be removed.
			else if ( key == "JUMP" || key == "jump" ) {
				core::Size jump_no;
				is >> jump_no;
				if ( is.fail() || jump_no < 1 || jump_no > pose.num_jump() ) {
					TR << "uh-oh, pose doesn't have a jump " << jump_no << std::endl;
				} else {
					Jump::Matrix mat;
					Jump::Vector vec;
					Jump jump;
					for ( int i = 1; i <= 3; ++i ) {
						for ( int j = 1; j <= 3; ++j ) {
							is >> mat(i,j);
						}
					}
					for ( int i = 1; i <= 3; ++i ) {
						is >> vec(i);
					}
					if ( is.fail() ) {
						TR << "error reading in jump data for jump " << jump_no << std::endl;
					} else {
						jump.set_rotation(mat);
						jump.set_translation(vec);
						pose.set_jump(jump_no, jump);
					}
				}
			} else if ( key == "MUTATE" ) {
				core::Size resnum;
				std::string resname;
				is >> resnum >> resname;
				resname = chemical::fixup_patches( resname );
				if ( is.fail() ) {
					TR << "error reading sequence mutation data" << std::endl;
				} else if ( resnum < 1 || resnum > pose.size() ) {
					TR << "d'oh, pose doesn't have a residue " << resnum << std::endl;
				} else if ( !pose.residue(resnum).residue_type_set()->has_name(resname) ) {
					TR << "unrecognized residue type name '" << resname << "'" << std::endl;
				} else if ( pose.residue_type(resnum).name() == resname ) {
					// These will happen routinely when reading a reference pose,
					// where the entire sequence is specified explicitly.
					TR.Debug << "ignoring no-op mutation: MUTATE " << resnum << " " << resname << std::endl;
				} else {
					using namespace core::conformation;
					ResidueOP newres = ResidueFactory::create_residue(
						// if can't find residue name, name_map() calls exit()
						pose.residue(resnum).residue_type_set()->name_map(resname),
						pose.residue(resnum), pose.conformation());
					pose.replace_residue(resnum, *newres, true /*orient backbone*/);
				}
			} else if ( key == "FOLD_TREE" ) {
				// This is typically present only when using embedded PDB format (see below).
				core::kinematics::FoldTree foldtree;
				is.clear();
				is.seekg(0, std::ios::beg);
				is >> foldtree;
				if ( !is.fail() && Size(foldtree.nres()) == pose.size() ) {
					pose.fold_tree(foldtree);
				} else {
					TR << "danger danger, error reading fold tree" << std::endl;
				}
			} else if ( key == "BEGIN_PDB_FORMAT" ) {
				// Just an embedded PDB-format file, separated by delimiters.
				// This is currently used for writing reference structures for the diffs,
				// but could also be used in cases where a diff just isn't efficient.
				using namespace core::chemical;
				using namespace core::import_pose;
				std::string const end_pdb_key = "END_PDB_FORMAT";
				utility::vector1< io::pdb::Record > pdb_data;
				while ( in.good() ) {
					getline(in, line);
					if ( in.fail() ) {
						TR << "getline() failed while reading embedded PDB" << std::endl;
						return false;
					}
					if ( line.size() >= end_pdb_key.size() && line.compare(0, end_pdb_key.size(), end_pdb_key) == 0 ) {
						break; // my kingdom for startswith()!

						// I can't figure out if getline() might leave behind a \n or \r\n,
						// but all real PDB records are longer than 2 characters anyway.
					} else if ( line.size() > 2 ) {
						pdb_data.push_back( io::pdb::create_record_from_pdb_line( line ) );
					}
				}
				core::io::StructFileRep sfr( core::io::pdb::create_sfr_from_pdb_records(pdb_data) );
				sfr.filename() = "atom_tree_diff.pdb"; // I'm afraid to leave this empty...
				pose.clear();
				core::import_pose::build_pose( sfr.clone(), pose, *(ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) ) );
			}
		} else { // it's an atom line...
			if ( rsd_no < 1 || rsd_no > pose.size() ) {
				TR << "uh-oh, pose doesn't have a residue " << rsd_no << std::endl;
			} else if ( atom_no < 1 || atom_no > pose.residue_type(rsd_no).natoms() ) {
				TR << "uh-oh, residue " << rsd_no << " doesn't have an atom " << atom_no << std::endl;
			} else {
				core::Real phi_value, theta_value, d_value;
				//TR << "setting dofs for " << rsd_no << " " << atom_no << " ...";
				AtomID aid(atom_no, rsd_no);
				DOF_ID dof_phi(aid, PHI), dof_theta(aid, THETA), dof_d(aid, D);
				is >> phi_value;
				if ( !is.fail() ) {
					//TR << " phi";
					pose.set_dof( dof_phi, phi_value );
					is >> theta_value;
					if ( !is.fail() ) {
						//TR << " theta";
						pose.set_dof( dof_theta, theta_value );
						is >> d_value;
						if ( !is.fail() ) {
							//TR << " d (!)";
							pose.set_dof( dof_d, d_value );
						}
					}
				} else TR << "somebody screwed up ... no phi to set for " << rsd_no << " " << atom_no << " !" << std::endl;
				//TR << std::endl;
			}
		}
	}
	TR << "beeg trouble for moose and squirrel -- reached EOF" << std::endl;
	return false;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Helper for dump_atom_tree_diff(), fills map with weighted score terms.
void map_of_weighted_scores(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	Scores & scores_out
)
{
	using namespace core::scoring;

	core::Real const tot_score = sfxn( pose );
	/// Now handled automatically.  sfxn.accumulate_residue_total_energies( pose );

	// Which score terms to use
	typedef utility::vector1<ScoreType> ScoreTypeVec;
	ScoreTypeVec score_types;
	for ( int i = 1; i <= n_score_types; ++i ) {
		ScoreType ii = ScoreType(i);
		if ( sfxn.has_nonzero_weight(ii) ) score_types.push_back(ii);
	}

	scores_out.clear();
	for ( auto const & score_type : score_types ) {
		scores_out[ name_from_score_type(score_type) ] = ( sfxn.get_weight(score_type) * pose.energies().total_energies()[ score_type ] );
	}
	scores_out[ name_from_score_type(core::scoring::total_score) ] = tot_score;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Copies the pose and adds random noise to the AtomTree DOFs
/// in the Nth digit (specified by bb_precision and sc_precision).
/// This simulates the roundoff error introduced by writing an atom_tree diff
/// to file and then reconstructing a pose from it.
/// Backbone needs to be more precise than sidechains because in typical foldtrees
/// it propagates errors over a much longer distance (i.e. the whole protein chain
/// rather than just to the end of the sidechain).
/// Using 200-500 residue single-chain globular proteins with one ligand each,
/// 6 digits of backbone precision and 4 digits of sidechain precision give
/// RMS coordinate errors of 0.0001 to 0.0003.  Since this is smaller than
/// the least significant digit in a PDB file, it should have an imperceptible
/// effect on the final structure.
/// Even in release mode, this function takes a long time:  30 sec to 1 minute.
void rms_error_with_noise(
	core::pose::Pose const & ref_pose,
	int bb_precision /*= 6*/,
	int sc_precision /*= 4*/
)
{
	using core::Size;
	using core::Real;
	using basic::T;
	using basic::Warning;

	using namespace core::id;
	using namespace core::scoring;

	if ( bb_precision < sc_precision ) Warning() << "Crazy fool, bb_precision should be >= sc_precision!" << std::endl;
	Real bb_tol = 1.0, sc_tol = 1.0;
	for ( int i = 0; i < bb_precision; ++i ) bb_tol /= 10.0;
	for ( int i = 0; i < sc_precision; ++i ) sc_tol /= 10.0;

	core::pose::Pose pose(ref_pose);
	core::kinematics::AtomTree const & atom_tree = pose.atom_tree();
	core::kinematics::FoldTree const & foldtree = pose.fold_tree();

	for ( Size rsd = 1, rsd_end = pose.size(); rsd <= rsd_end; ++rsd ) {
		bool const is_jump_residue = foldtree.is_jump_point(rsd);
		for ( Size atom = 1, atom_end = pose.residue(rsd).natoms(); atom <= atom_end; ++atom ) {
			AtomID aid(atom, rsd);
			DOF_ID dof_phi(aid, PHI), dof_theta(aid, THETA), dof_d(aid, D);
			if ( atom_tree.atom(aid).is_jump() ) {
				// Jump atoms have the rbN values all set to zero -- maybe vestigal part of atom tree?
			} else {
				// Backbone heavyatoms have to be very precise, or else a lot of error
				// can propagate to the end of the chain.
				// Sidechain atoms and hydrogens don't have to worry nearly so much.
				//bool is_backbone = pose.residue(rsd).atom_type(atom).is_heavyatom() && atom <= pose.residue(rsd).last_backbone_atom();
				bool is_backbone = is_jump_residue // i.e. things should be "backbone" if they're anchors for a jump!
					|| (pose.residue(rsd).atom_type(atom).is_heavyatom() && atom <= pose.residue(rsd).last_backbone_atom());
				//Real precision = (is_backbone ? bb_precision : sc_precision);
				Real tol = (is_backbone ? bb_tol : sc_tol);

				// Introduce random error into our copy pose, +/- tol/2
				pose.set_dof(dof_phi,   pose.dof(dof_phi)   + tol * (numeric::random::rg().uniform() - 0.5));
				pose.set_dof(dof_theta, pose.dof(dof_theta) + tol * (numeric::random::rg().uniform() - 0.5));
				pose.set_dof(dof_d,     pose.dof(dof_d)     + tol * (numeric::random::rg().uniform() - 0.5));
			}// end non-jump atom
		}// end loop over atoms
	}// end loop over residues

	T("rms_error_with_noise") << "bb=" << bb_precision << "," << bb_tol << "  sc=" << sc_precision << "," << sc_tol
		<< "  rms_no_super=" << rmsd_no_super(ref_pose, pose, is_heavyatom) << "  rms_with_super=" << rmsd_with_super(ref_pose, pose, is_heavyatom)
		<< "  polymer_rms_no_super=" << rmsd_no_super(ref_pose, pose, is_polymer_heavyatom) << "  polymer_rms_with_super=" << rmsd_with_super(ref_pose, pose, is_polymer_heavyatom)
		<< std::endl;
}

/// @brief Test if given file is an atom_tree_diff
bool file_is_atom_tree_diff( std::string const & filename ) {
	utility::io::izstream infile( filename );
	return file_is_atom_tree_diff( infile );
}

/// @brief Test if given stream is an atom_tree_diff
/// @details If everything goes right, after the call, the read position should be at the same place it was to start with
bool file_is_atom_tree_diff( std::istream & in ) {
	std::streampos startpos(in.tellg());
	std::string ignored_pose_tag;
	Scores ignored_scores;
	// There's probably a better heuristic
	bool retval = header_from_atom_tree_diff(in, ignored_pose_tag, ignored_scores);
	in.seekg(startpos);
	return retval;
}

} // atom_tree_diffs
} // import_pose
} // core
