// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @buf
/// @brief
/// @author will sheffler

#include <core/io/serialization/serialize_pose.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <numeric/xyzVector.hh>


#include <utility/vector1.hh>
#include <sstream>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

#define THROW_EXCEPTION(X) {std::cerr << "OVERFLOW ERROR: " <<  X;return;}

namespace core {
namespace io {
namespace serialization {

	//const std::string PSEUDORESIDUE_NAME = "VRT";

	bool
		is_pseudoresidue(const conformation::Residue & residue)
	{
		// check the residue name against the pseudoresidue name
		return (residue.aa() == core::chemical::aa_vrt);
	}


	/// Helper function to write all binary files as little endian.
	/// If the current architecture is big endian, swap the bytes.
#ifdef IS_BIG_ENDIAN
	static void
		swap_bytes(char * x, unsigned int n)
	{
		for (unsigned int i = 0; i < n / 2; ++ i) {
			char tmp = x[i];
			x[i] = x[n - i - 1];
			x[n - i - 1] = tmp;
		}
	}
#else
	static void swap_bytes(char *, unsigned int) {}
#endif


	/// Helper functions to read/write raw memory

	void
		write_bytes(char * x, unsigned int n, BUFFER & buf)
	{
		//std::cerr << "write bytes in : '" << x << "' " << n << std::endl;
		swap_bytes(x, n);
		if (buf.write(x, n) != 1) {
			THROW_EXCEPTION("Error writing\n");
		}
		//std::cerr << "write bytes out: '" << x << "'" << std::endl;
	}

	void
		read_bytes(char * x, unsigned int n, BUFFER & buf)
	{
		//std::cerr << "read bytes in   : '" << x << "' " << n << std::endl;
		if (buf.read(x, n) != 1) {
			THROW_EXCEPTION("Error reading\n");
		}
		swap_bytes(x, n);
		//std::cerr << "read bytes out  : '" << x << "'" << std::endl;
	}


	/// Helper functions to read/write any single primitive type

	template <typename T>
	static void
		write_bytes(T x, BUFFER & buf)
	{
		write_bytes((char *)&x, sizeof(T), buf);
	}

	template <typename T>
	static void
		read_bytes(T & x, BUFFER & buf)
	{
		read_bytes((char *)&x, sizeof(T), buf);
	}


	/// Read/write a single primitive type to a buf.

	void
		write_binary(char x, BUFFER & buf)
	{
		write_bytes<char>(x, buf);
	}

	void
		read_binary(char & x, BUFFER & buf)
	{
		read_bytes<char>(x, buf);
	}

	void
		write_binary(bool x, BUFFER & buf)
	{
		char x_char = x ? 0xFF : 0x00;
		write_bytes(x_char, buf);
	}

	void
		read_binary(bool & x, BUFFER & buf)
	{
		char x_char(0);
		read_bytes(x_char, buf);
		x = (x_char != 0x00);
	}

	void
		write_binary(unsigned int x, BUFFER & buf)
	{
		write_bytes(x, buf);
	}

	void
		read_binary(unsigned int & x, BUFFER & buf)
	{
		read_bytes(x, buf);
	}

	void
		write_binary(float x, BUFFER & buf)
	{
		write_bytes(x, buf);
	}

	void
		read_binary(float & x, BUFFER & buf)
	{
		read_bytes(x, buf);
	}


	void
		write_binary(double x, BUFFER & buf)
	{
		write_bytes(x, buf);
	}

	void
		read_binary(double & x, BUFFER & buf)
	{
		read_bytes(x, buf);
	}


	/// Read/write simple structure to a buf.

	void
		write_binary(const utility::vector1_bool & x, BUFFER & buf)
	{
		write_binary((unsigned int)(x.size()), buf);
		for (unsigned int ii = 1; ii <= x.size(); ++ ii) {
			write_binary(bool(x[ii]), buf);
		}
	}

	void
		read_binary(utility::vector1_bool & x, BUFFER & buf)
	{
		unsigned int size(0);
		read_binary(size, buf);

		x.resize(size);
		for (unsigned int ii = 1; ii <= x.size(); ++ ii) {
			bool x_ii;
			read_binary(x_ii, buf);
			x[ii] = x_ii;
		}
	}

	void
		write_binary(const std::vector<std::string> & x, BUFFER & buf)
	{
		write_binary((unsigned int)(x.size()), buf);
		for (unsigned int ii = 0; ii < x.size(); ++ ii) {
			write_binary(std::string(x[ii]), buf);
		}
	}

	void
		read_binary(std::vector<std::string> & x, BUFFER & buf)
	{
		unsigned int size(0);
		read_binary(size, buf);

		x.resize(size);
		for (unsigned int ii = 0; ii < x.size(); ++ ii) {
			std::string x_ii;
			read_binary(x_ii, buf);
			x[ii] = x_ii;
		}
	}

	void
		write_binary(const std::string & x, BUFFER & buf)
	{
		write_binary((unsigned int)(x.size()), buf);
		for (unsigned int ii = 0; ii < x.size(); ++ ii) {
			write_binary(x[ii], buf);
		}
	}

	void
		read_binary(std::string & x, BUFFER & buf)
	{
		unsigned int size(0);
		read_binary(size, buf);

		x.resize(size);
		for (unsigned int ii = 0; ii < x.size(); ++ ii) {
			char x_ii(0);
			read_binary(x_ii, buf);
			x[ii] = x_ii;
		}
	}

	void
		write_binary(const core::Vector & x, BUFFER & buf)
	{
		write_binary(x.x(), buf);
		write_binary(x.y(), buf);
		write_binary(x.z(), buf);
	}

	void
		read_binary(core::Vector & x, BUFFER & buf)
	{
		read_binary(x.x(), buf);
		read_binary(x.y(), buf);
		read_binary(x.z(), buf);
	}


	/// Utility read/write.

	void
		check_binary_unsigned_int(unsigned int x, BUFFER & buf)
	{
		unsigned int input = 0;
		read_binary(input, buf);
		if (input != x) {
			THROW_EXCEPTION("Error reading .\n");
		}
	}

	void
		write_binary_chars(const char *x, BUFFER & buf)
	{
		size_t len = strlen(x);

		for (unsigned int ii = 0; ii < len; ++ ii) {
			write_binary(x[ii], buf);
		}
	}

	void
		check_binary_chars(const char *x, BUFFER & buf)
	{
		size_t len = strlen(x);
		std::vector< char > input(len, 0);

		for (unsigned int ii = 0; ii < len; ++ ii) {
			read_binary(input[ii], buf);
		}

		if (memcmp(&input[0], x, len)) {
			THROW_EXCEPTION("Error reading .\n");
		}
	}


	void
		write_binary(const core::pose::Pose & pose, BUFFER & buf)
	{
		const unsigned int WRITE_VERSION = 2;
		write_binary(WRITE_VERSION, buf);

// residues
		write_binary((unsigned int)pose.total_residue(), buf);
		for (size_t j = 1; j <= pose.total_residue(); ++ j) {
			const core::conformation::Residue & residue = pose.residue(j);
		//	if (is_pseudoresidue(residue)) continue;
// residue name
			//write_binary(residue.name(), buf);
// residue type set
			write_binary(residue.type().residue_type_set().name(), buf);
// residue type
			write_binary(residue.type().name(), buf);
// residue chain
      write_binary((unsigned int)residue.chain(), buf);
// atoms xyz
			write_binary((unsigned int)residue.atoms().size(), buf);
			for (size_t k = 1; k <= residue.atoms().size(); ++ k) {
				write_binary(residue.atoms()[k].xyz(), buf);
			}
		}

// fold tree
    write_binary_chars("FOLDTREE", buf);
    std::ostringstream foldtree_outs;
    foldtree_outs << pose.fold_tree();
    write_binary( foldtree_outs.str(), buf);

// jumps
    write_binary_chars("JUMP", buf);
    write_binary((unsigned int)pose.num_jump(), buf);
    for (core::Size j=1; j <= pose.num_jump(); ++j) {
      std::ostringstream jump_outs;
      jump_outs << pose.jump(j);
      write_binary( jump_outs.str(), buf);
    }

// SS
		write_binary_chars("SSTR", buf);
		write_binary((unsigned int)pose.total_residue(), buf);
		for (size_t j = 1; j <= pose.total_residue(); ++ j) {
			write_binary(pose.secstruct(j), buf);
		}
	}


  void
    read_binary(core::pose::Pose & pose, BUFFER & buf)
  {

    using namespace core::chemical;
    using namespace core::conformation;
    unsigned int read_version = 0;
    read_binary(read_version, buf);

    switch (read_version) {
    case 2:
     {

				pose.clear();

				//typedef utility::vector1< std::pair< Size, conformation::ResidueOP > > ResidueVector;
				typedef utility::vector1< std::pair< id::AtomID, Vector > > AtomVector;
// total residue
				unsigned int total_residue = 0;
				read_binary( total_residue, buf);
				AtomVector atom_vec;
				atom_vec.reserve(total_residue*15);

				unsigned int prevchain = 1;
				for (size_t j = 1; j <= total_residue; ++ j) {
// type set
					std::string type_set;
					read_binary(type_set, buf);
// type name
					std::string type_name;
					read_binary(type_name, buf);
// chain
					unsigned int chainid(0);
					read_binary(chainid, buf);
					bool const is_lower_terminus( j==1 || chainid != prevchain );
// new residue
					ResidueTypeSetCOP target_residue_type_set( ChemicalManager::get_instance()->residue_type_set( type_set ) );
					ResidueOP new_rsd = ResidueFactory::create_residue( target_residue_type_set->name_map( type_name ) ); //, residue, pose.conformation() );
					new_rsd->chain(chainid);
					if (j > 1 && (is_lower_terminus || !new_rsd->is_polymer() || !pose.residue_type(j-1).is_polymer() )){
						pose.append_residue_by_jump( *new_rsd, 1 );
					} else {
						pose.append_residue_by_bond( *new_rsd ); // pose.append_residue_by_jump( *new_rsd, j-1,"", "", true );
					}
					prevchain = chainid;
					// update the pose-internal chain label if necessary
					if ( is_lower_terminus && pose.total_residue() > 1 ) {
						pose.conformation().insert_chain_ending( pose.total_residue() - 1 );
					}
// atoms xyz
					unsigned int natoms(0);
					read_binary(natoms, buf);
					for (size_t k = 1; k <= natoms; ++ k) {
						core::Vector xyz(0.0, 0.0, 0.0);
						read_binary(xyz, buf);
						atom_vec.push_back(std::make_pair(id::AtomID(k, j), xyz));
					}
				}
// fold tree
				check_binary_chars("FOLDTREE", buf);
				std::string foldtree_str;
				read_binary( foldtree_str, buf );
				std::istringstream foldtree_in;
				foldtree_in.str(foldtree_str);
				kinematics::FoldTree f;
				foldtree_in >> f;
				pose.fold_tree( f );
// jumps
				check_binary_chars("JUMP", buf);
				unsigned int jumps(0);
				read_binary(jumps, buf);
				utility::vector1< kinematics::Jump > jumpsv;
				for (size_t j = 1; j <= jumps; ++ j) {
					std::string jump_str;
					read_binary( jump_str, buf );
					std::istringstream jump_in;
					jump_in.str( jump_str );
					kinematics::Jump jump;
					jump_in >> jump;
					jumpsv.push_back( jump );
				}
// set jumps
				for ( Size nr = 1; nr <= f.num_jump(); nr++)  {
					pose.set_jump( nr, jumpsv[nr] );
				}

// SS
				check_binary_chars("SSTR", buf);
				check_binary_unsigned_int(pose.total_residue(), buf);
				for (size_t j = 1; j <= pose.total_residue(); ++ j) {
					char sstr(0);
					read_binary(sstr, buf);
					pose.set_secstruct(j, sstr);
				}

// set atom xyz
				for (AtomVector::iterator ii = atom_vec.begin(); ii != atom_vec.end(); ++ ii) {
					pose.set_xyz(ii->first, ii->second);
				}

			} break;

		default:
			THROW_EXCEPTION(std::string("Error reading ; bad pose version.\n"));
			break;
		}
	}

} // pose
} // io
} // core
