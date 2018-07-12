// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/util.cc
/// @brief: various supplemental declarations for network layer
///
/// @author Sergey Lyskov


#ifdef ZEROMQ

#include <protocols/network/util.hh>

#include <core/pose/Pose.hh>

#include <cereal/archives/binary.hpp>
#include <utility/io/zipstream.hpp>

#include <utility/pointer/memory.hh>

namespace protocols {
namespace network {


zmq::context_t &zmq_context()
{
	static ContextUP context( new zmq::context_t);
	return *context;
}


//  Receive ZeroMQ message as a string
std::string receive_message(zmq::socket_t & socket)
{
	zmq::message_t message;
	socket.recv(&message);

	return std::string( static_cast<char*>( message.data() ), message.size() );
}

//  Send ZeroMQ message as string
bool send_message(zmq::socket_t & socket, std::string const & string_message, int flags)
{
	zmq::message_t message( string_message.size() );
	memcpy(message.data(), string_message.data(), string_message.size());

	return socket.send(message, flags);
}


// bool send_message(zmq::socket_t & socket, void const *data, int size, int flags)
// {
// 	zmq::message_t message(size);
// 	memcpy(message.data(), data, ssize);
// 	return socket.send(message, flags);
// }



// serialize and compress Pose UI way
PoseBinary pose_to_bytes(core::pose::Pose const &pose)
{
	std::ostringstream buffer;
	auto arc = cereal::BinaryOutputArchive(buffer);
	pose.save_with_options(arc, /* save_observers = */ false);

	auto uncompressed_pose = buffer.str();

	return uncompressed_pose;

	//qDebug() << "Pose serialized into " << uncompressed_pose.size() << " bytes";
	/* // turn uncompressed_pose into compressed byte stream
	std::ostringstream compressed_buffer;
	zlib_stream::zip_ostream zos(compressed_buffer);
	int u_size = uncompressed_pose.size();
	zos.write(reinterpret_cast<char*>(&u_size), sizeof(u_size) );
	zos.write(uncompressed_pose.c_str(), u_size);
	//zos << uncompressed_pose.size() << uncompressed_pose;
	zos.zflush();  // zflush_finalize zflush

	return compressed_buffer.str();
	*/
	//qDebug() << "Pose serialized and compressed into " << compressed_pose.size() << " bytes";
}


// serialize and compress Pose UI way, if nullptr is given generate empty byte-string
PoseBinary pose_to_bytes(core::pose::PoseCOP const &pose)
{
	if(pose) return pose_to_bytes(*pose);
	else return PoseBinary();
}


// uncompress and deserialize aPose UI way
core::pose::PoseOP bytes_to_pose(PoseBinary const &compressed_pose)
{
	if( compressed_pose.empty() ) return core::pose::PoseOP();
	else {
		auto pose = utility::pointer::make_shared<core::pose::Pose>();

		/* // turn compressed_pose into decompressed pose
		   std::vector<char> uncompressed_pose2;
		   std::istringstream uncompressed_buffer2(compressed_pose);
		   zlib_stream::zip_istream zis(uncompressed_buffer2);
		   int c_size = 0;
		   zis.read(reinterpret_cast<char*>(&c_size), sizeof(c_size) );  //uncompressed_pose2;
		   uncompressed_pose2.resize(c_size);
		   zis.read(&uncompressed_pose2[0], uncompressed_pose2.size());

		   struct MemoryBuffer : std::streambuf
		   {
		   MemoryBuffer(std::vector<char> & source) { this->setg(source.data(), source.data(), source.data()+source.size() ); }
		   };

		   std::string uncompressed_pose2_str(uncompressed_pose2.begin(), uncompressed_pose2.end());

		   MemoryBuffer ibuff(uncompressed_pose2);
		   std::istream in(&ibuff);
		   auto darc = cereal::BinaryInputArchive(in);

		   pose->load(darc);
		*/

		// struct MemoryBuffer : std::streambuf
		// {
		// 	MemoryBuffer(std::string const & source) { this->setg(const_cast<char*>( source.c_str() ), const_cast<char*>( source.c_str() ) , const_cast<char*>( source.c_str()+source.size() ) ); }
		// };

		std::vector<char> uncompressed_pose2(compressed_pose.begin(), compressed_pose.end());
		struct MemoryBuffer : std::streambuf
		{
			MemoryBuffer(std::vector<char> & source) { this->setg(source.data(), source.data(), source.data()+source.size() ); }
		};
		MemoryBuffer ibuff(uncompressed_pose2);
		std::istream in(&ibuff);
		auto darc = cereal::BinaryInputArchive(in);

		pose->load_with_options(darc, /* load_observers = */ false);


		//qDebug() << "Pose decompressed into " << c_size << " bytes equal:" << (uncompressed_pose2_str == uncompressed_pose);
		//qDebug() << "Pose decompressed into " << zis.get_out_size() << " _ " << zis.get_in_size() <<" bytes";

		return pose;
	}
}


} // namespace network
} // namespace protocols

#endif // ZEROMQ
