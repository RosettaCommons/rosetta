// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/PyMolMover.hh
/// @brief  Send infromation to PyMol. Contain classes PyMolMover, PyMolObserver and helper classes.
/// @author Sergey Lyskov

#ifndef INCLUDED_protocols_moves_PyMolMover_hh
#define INCLUDED_protocols_moves_PyMolMover_hh

// unit headers
#include <protocols/moves/PyMolMover.fwd.hh>

// protocol headers
#include <protocols/moves/Mover.hh>

// core headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/datacache/CacheableObserver.hh>

#include <core/scoring/ScoreType.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/signals/Link.hh>

// c++ headers
#include <string>

// REQUIRED FOR WINDOWS
#ifndef __native_client__
#ifdef WIN32
#include <winsock.h>
#undef min
#undef max
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif
#endif

#ifdef __native_client__
  typedef int sockaddr_in;
#endif

//#ifdef WIN_PYROSETTA
//  typedef int sockaddr_in;
//#endif

namespace protocols {
namespace moves {

enum X11Colors {
	XC_first_color = 0,
	XC_black = XC_first_color, // 0 0 0
	XC_AntiqueWhite = 1,  // 250 235 215
	XC_BlanchedAlmond = 2,  // 255 235 205
	XC_BlueViolet = 3,  // 138 43 226
	XC_CadetBlue = 4,  // 95 158 160
	XC_CornflowerBlue = 5,  // 100 149 237
	XC_DarkBlue = 6,  // 0 0 139
	XC_DarkCyan = 7,  // 0 139 139
	XC_DarkGoldenrod = 8,  // 184 134 11
	XC_DarkGray = 9,  // 169 169 169
	XC_DarkGreen = 10,  // 0 100 0
	XC_DarkGrey = 11,  // 169 169 169
	XC_DarkKhaki = 12,  // 189 183 107
	XC_DarkMagenta = 13,  // 139 0 139
	XC_DarkOliveGreen = 14,  // 85 107 47
	XC_DarkOrange = 15,  // 255 140 0
	XC_DarkOrchid = 16,  // 153 50 204
	XC_DarkRed = 17,  // 139 0 0
	XC_DarkSalmon = 18,  // 233 150 122
	XC_DarkSeaGreen = 19,  // 143 188 143
	XC_DarkSlateBlue = 20,  // 72 61 139
	XC_DarkSlateGray = 21,  // 47 79 79
	XC_DarkSlateGrey = 22,  // 47 79 79
	XC_DarkTurquoise = 23,  // 0 206 209
	XC_DarkViolet = 24,  // 148 0 211
	XC_DebianRed = 25,  // 215 7 81
	XC_DeepPink = 26,  // 255 20 147
	XC_DeepSkyBlue = 27,  // 0 191 255
	XC_DimGray = 28,  // 105 105 105
	XC_DimGrey = 29,  // 105 105 105
	XC_DodgerBlue = 30,  // 30 144 255
	XC_FloralWhite = 31,  // 255 250 240
	XC_ForestGreen = 32,  // 34 139 34
	XC_GhostWhite = 33,  // 248 248 255
	XC_GreenYellow = 34,  // 173 255 47
	XC_HotPink = 35,  // 255 105 180
	XC_IndianRed = 36,  // 205 92 92
	XC_LavenderBlush = 37,  // 255 240 245
	XC_LawnGreen = 38,  // 124 252 0
	XC_LemonChiffon = 39,  // 255 250 205
	XC_LightBlue = 40,  // 173 216 230
	XC_LightCoral = 41,  // 240 128 128
	XC_LightCyan = 42,  // 224 255 255
	XC_LightGoldenrod = 43,  // 238 221 130
	XC_LightGoldenrodYellow = 44,  // 250 250 210
	XC_LightGray = 45,  // 211 211 211
	XC_LightGreen = 46,  // 144 238 144
	XC_LightGrey = 47,  // 211 211 211
	XC_LightPink = 48,  // 255 182 193
	XC_LightSalmon = 49,  // 255 160 122
	XC_LightSeaGreen = 50,  // 32 178 170
	XC_LightSkyBlue = 51,  // 135 206 250
	XC_LightSlateBlue = 52,  // 132 112 255
	XC_LightSlateGray = 53,  // 119 136 153
	XC_LightSlateGrey = 54,  // 119 136 153
	XC_LightSteelBlue = 55,  // 176 196 222
	XC_LightYellow = 56,  // 255 255 224
	XC_LimeGreen = 57,  // 50 205 50
	XC_MediumAquamarine = 58,  // 102 205 170
	XC_MediumBlue = 59,  // 0 0 205
	XC_MediumOrchid = 60,  // 186 85 211
	XC_MediumPurple = 61,  // 147 112 219
	XC_MediumSeaGreen = 62,  // 60 179 113
	XC_MediumSlateBlue = 63,  // 123 104 238
	XC_MediumSpringGreen = 64,  // 0 250 154
	XC_MediumTurquoise = 65,  // 72 209 204
	XC_MediumVioletRed = 66,  // 199 21 133
	XC_MidnightBlue = 67,  // 25 25 112
	XC_MintCream = 68,  // 245 255 250
	XC_MistyRose = 69,  // 255 228 225
	XC_NavajoWhite = 70,  // 255 222 173
	XC_NavyBlue = 71,  // 0 0 128
	XC_OldLace = 72,  // 253 245 230
	XC_OliveDrab = 73,  // 107 142 35
	XC_OrangeRed = 74,  // 255 69 0
	XC_PaleGoldenrod = 75,  // 238 232 170
	XC_PaleGreen = 76,  // 152 251 152
	XC_PaleTurquoise = 77,  // 175 238 238
	XC_PaleVioletRed = 78,  // 219 112 147
	XC_PapayaWhip = 79,  // 255 239 213
	XC_PeachPuff = 80,  // 255 218 185
	XC_PowderBlue = 81,  // 176 224 230
	XC_RosyBrown = 82,  // 188 143 143
	XC_RoyalBlue = 83,  // 65 105 225
	XC_SaddleBrown = 84,  // 139 69 19
	XC_SandyBrown = 85,  // 244 164 96
	XC_SeaGreen = 86,  // 46 139 87
	XC_SkyBlue = 87,  // 135 206 235
	XC_SlateBlue = 88,  // 106 90 205
	XC_SlateGray = 89,  // 112 128 144
	XC_SlateGrey = 90,  // 112 128 144
	XC_SpringGreen = 91,  // 0 255 127
	XC_SteelBlue = 92,  // 70 130 180
	XC_VioletRed = 93,  // 208 32 144
	XC_WhiteSmoke = 94,  // 245 245 245
	XC_YellowGreen = 95,  // 154 205 50
	XC_aquamarine = 96,  // 127 255 212
	XC_azure = 97,  // 240 255 255
	XC_beige = 98,  // 245 245 220
	XC_bisque = 99,  // 255 228 196
	XC_AliceBlue = 100,  // 240 248 255
	XC_blue = 101,  // 0 0 255
	XC_blue1 = 102,  // 0 0 255
	XC_blue2 = 103,  // 0 0 238
	XC_blue3 = 104,  // 0 0 205
	XC_blue4 = 105,  // 0 0 139
	XC_brown = 106,  // 165 42 42
	XC_burlywood = 107,  // 222 184 135
	XC_chartreuse = 108,  // 127 255 0
	XC_chocolate = 109,  // 210 105 30
	XC_coral = 110,  // 255 127 80
	XC_cornsilk = 111,  // 255 248 220
	XC_cyan = 112,  // 0 255 255
	XC_firebrick = 113,  // 178 34 34
	XC_gainsboro = 114,  // 220 220 220
	XC_gold = 115,  // 255 215 0
	XC_goldenrod = 116,  // 218 165 32
	XC_gray = 117,  // 190 190 190
	XC_gray0 = 118,  // 0 0 0
	XC_gray10 = 119,  // 26 26 26
	XC_gray100 = 120,  // 255 255 255
	XC_gray20 = 121,  // 51 51 51
	XC_gray30 = 122,  // 77 77 77
	XC_gray40 = 123,  // 102 102 102
	XC_gray50 = 124,  // 127 127 127
	XC_gray60 = 125,  // 153 153 153
	XC_gray70 = 126,  // 179 179 179
	XC_gray80 = 127,  // 204 204 204
	XC_gray90 = 128,  // 229 229 229
	XC_green = 129,  // 0 255 0
	XC_green1 = 130,  // 0 255 0
	XC_green2 = 131,  // 0 238 0
	XC_green3 = 132,  // 0 205 0
	XC_green4 = 133,  // 0 139 0
	XC_honeydew = 134,  // 240 255 240
	XC_ivory = 135,  // 255 255 240
	XC_khaki = 136,  // 240 230 140
	XC_lavender = 137,  // 230 230 250
	XC_linen = 138,  // 250 240 230
	XC_magenta = 139,  // 255 0 255
	XC_maroon = 140,  // 176 48 96
	XC_moccasin = 141,  // 255 228 181
	XC_navy = 142,  // 0 0 128
	XC_orange = 143,  // 255 165 0
	XC_orchid = 144,  // 218 112 214
	XC_peru = 145,  // 205 133 63
	XC_pink = 146,  // 255 192 203
	XC_plum = 147,  // 221 160 221
	XC_purple = 148,  // 160 32 240
	XC_red = 149,  // 255 0 0
	XC_red1 = 150,  // 255 0 0
	XC_red2 = 151,  // 238 0 0
	XC_red3 = 152,  // 205 0 0
	XC_red4 = 153,  // 139 0 0
	XC_salmon = 154,  // 250 128 114
	XC_seashell = 155,  // 255 245 238
	XC_sienna = 156,  // 160 82 45
	XC_snow = 157,  // 255 250 250
	XC_snow1 = 158,  // 255 250 250
	XC_snow2 = 159,  // 238 233 233
	XC_snow3 = 160,  // 205 201 201
	XC_snow4 = 161,  // 139 137 137
	XC_tan = 162,  // 210 180 140
	XC_thistle = 163,  // 216 191 216
	XC_tomato = 164,  // 255 99 71
	XC_turquoise = 165,  // 64 224 208
	XC_violet = 166,  // 238 130 238
	XC_wheat = 167,  // 245 222 179
	XC_white = 168,  // 255 255 255
	XC_yellow = 169, // 255 255 0
	XC_last_color
};


/// @brief PyMolMover helper class. Handle low level UDP transactions stuff.
///        This is a port of original Python version of UDP socket client written writen for PyRosetta
class UDPSocketClient
{
public:
	/// @brief ctor
	UDPSocketClient();

	/// @brief cctor
	UDPSocketClient( UDPSocketClient const & other );

	/// @brief dtor
	~UDPSocketClient();

	void sendMessage(std::string msg);

	void show(std::ostream & output) const;

private:

	void sendRAWMessage(int globalPacketID, int packetI, int packetCount, char * msg_begin, char *msg_end);

	/// @brief last know mximum size of suspenseful sended UDP packet.
	///        ~64k for local connection and ~10k for inet connection
	unsigned int max_packet_size_;  // last know mximum size of suspenseful sended UDP packet.

	/// @brief unique id of this socket client
	union UUID {
		char bytes_[16];
		short unsigned int shorts_[8];
		};

	/// @brief Almost real UUID, but for simplicity we just use random sequence
	UUID uuid_;

	/// @brief counter for number of packet already sent
	int sentCount_;

	/// @brief socket address and handle
	sockaddr_in socket_addr_;
	int socket_h_;
};

std::ostream &operator<< (std::ostream & output, UDPSocketClient const & client);

class PyMolMover : public protocols::moves::Mover
{
public:
	/// @brief ctor
	PyMolMover();

	/// @brief cctor
	PyMolMover( PyMolMover const & other );

	/// @brief dtor
	virtual ~PyMolMover();

	virtual std::string get_name() const;
	virtual void apply( Pose & );

	using protocols::moves::Mover::apply;

	/// @brief Actually our mover does not change the Pose object, so we have additional const version...
	virtual void apply( Pose const & );

	/// @brief Send message for PyMol to print
	void print(std::string const & message);

	/// @brief Send specified energy to PyMOL.
	void send_energy(Pose const &, core::scoring::ScoreType stype = core::scoring::total_score);

	/// @brief Send specified energy to PyMOL.
	void send_energy(Pose const &, std::string const & stype);

	/// @brief Send RAW energy array for coloring by PyMOL
	void send_RAW_Energies(Pose const &, std::string energyType, utility::vector1<int> const & energies);

	/// @brief Send Membrane Planes to PyMol
	/// @details If pose is a membrane pose
	/// pymol viewer will build CGO planes from points specified
	void send_membrane_planes( Pose const & );

	/// @brief Tell PyMOL to color protein with supplied custom colors
	virtual void send_colors(Pose const &, std::map<int, int> const & colors, X11Colors default_color=protocols::moves::XC_blue );

	/// @brief Specify PyMOL model name
	std::string pymol_name() { return pymol_name_; }
	void pymol_name(std::string new_pymol_name) { pymol_name_ = new_pymol_name; }


	/// @brief Flag that specify if PyMOL mover should send current Pose energy on every apply. If name set to empty string (default) then name derived from pdb_info will be used
	bool update_energy() { return update_energy_; }
	void update_energy(bool f) { update_energy_ = f; }


	core::scoring::ScoreType energy_type(void) { return energy_type_; }
	void energy_type(core::scoring::ScoreType t) { energy_type_ = t; }

	/// @brief Set the keep history flag. If set to True - PyMol will keep track of all frames that
	///        was sent.
	void keep_history(bool kh) { keep_history_ = kh; }


	/// @brief Return current keep_history flag.
	bool keep_history(void) { return keep_history_; }


	/// @brief Set current minimum update interval 't' in seconds. Mover will not send anyinformation to
	///         PyMOL before atleast 't' secons has passes since last packet was sent.
	///         default value is 0, - no packets are skipped.
	void update_interval(core::Real t) { update_interval_ = t; };


	/// @brief  Return current update interval.
	core::Real update_interval() { return update_interval_; };

	void show(std::ostream & output=std::cout) const;

	/// @brief Parses tag for rosetta scripts
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP clone() const;

	///@brief Set the pymol model name
	void
	set_PyMol_model_name( std::string name);

	std::string
	get_PyMol_model_name(Pose const & pose) const;

private:
	//void send_message(std::string const & message_type, bool keep_history, std::string const & name, std::string const &message);



	bool is_it_time();

	/// @brief Actual link object
	UDPSocketClient link_;

	bool update_energy_;

	core::scoring::ScoreType energy_type_;

	/// @brief If pose is a membrane pose, send planes
	bool update_membrane_;

	/// @brief Should PyMol keep history of all models that was sent? - Default is false.
	bool keep_history_;

	/// @brief  Update interval in seconds.
	core::Real update_interval_;

	/// @brief  Time stamp from where last packet was sent.
	core::Real last_packet_sent_time_;

	/// @brief  Name of model in pymol.
	std::string pymol_name_;
};  // class PyMolMover

// Insertion operator (overloaded so that PyMolMover can be "printed") in PyRosetta).
std::ostream &operator<< (std::ostream & output, PyMolMover const & mover);

class PyMolObserver : public core::pose::datacache::CacheableObserver
{
public:
	// This is set up to allow multiple settings with bit twiddling
	enum ObserverType {
		no_observer = 0,
		general_observer = 1,
		energy_observer = 2,
		conformation_observer = 4
	};

	PyMolObserver();
	PyMolObserver(PyMolObserver const & rval);
	virtual ~PyMolObserver();

	PyMolObserver &
	operator= (PyMolObserver const &rval);

	virtual
	core::pose::datacache::CacheableObserverOP
	clone();

	virtual
	core::pose::datacache::CacheableObserverOP
	create();

	void
	set_type( ObserverType setting);

	void
	add_type( ObserverType setting);

	ObserverType
	get_type() const { return type_; };

	virtual bool
	is_attached() const;

	virtual void generalEvent( core::pose::signals::GeneralEvent const & ev) {
		pymol_.apply( *ev.pose );
	};

	virtual void energyEvent( core::pose::signals::EnergyEvent const & ev) {
		pymol_.apply( *ev.pose );
	}

	virtual void conformationEvent( core::pose::signals::ConformationEvent const & ev) {
		pymol_.apply( *ev.pose);
	}

	PyMolMover & pymol() { return pymol_; };

	/// Attach observer to the pose object
	void attach(core::pose::Pose &p);

	/// Detach observer from the pose object
	void detach(core::pose::Pose &p);

protected:

	virtual void
	attach_impl(core::pose::Pose &pose);

	virtual void
	detach_impl();

	void
	update_links();

private:

	ObserverType type_;
	PyMolMover pymol_;

	utility::signals::Link general_event_link_;
	utility::signals::Link energy_event_link_;
	utility::signals::Link conformation_event_link_;
};

// Because C++ is silly about enum conversions
inline PyMolObserver::ObserverType operator| (PyMolObserver::ObserverType & l, PyMolObserver::ObserverType & r) {
	// We need to cast to int to avoid infinite loops
	return static_cast<PyMolObserver::ObserverType>( static_cast<int>(l) | static_cast<int>(r) );
}

/// @brief (Internal) helper function to create a PyMolObserver and add it to the given pose
/// NOTE: You NEED to adjust the observer type and call attach() on the return - by default a new PyMolObserver isn't attached/observing.
PyMolObserverOP
get_pymol_observer(core::pose::Pose & pose);

/// @brief Helper function that create PyMolObserver Object and add it to the give Pose.
///        This is the most likely the only function that you need to call...
PyMolObserverOP AddPyMolObserver(core::pose::Pose &p, bool keep_history=false, core::Real update_interval=0);

/// @brief Helper function that create PyMolObserver Object and add it to the give Pose energies object so pymol only updates on energy changes.
PyMolObserverOP AddPyMolObserver_to_energies(core::pose::Pose & p, bool keep_history=false, core::Real update_interval=0);

/// @brief Helper function that create PyMolObserver Object and add it to the give Pose conformation object so pymol only updates on conformation changes.
PyMolObserverOP AddPyMolObserver_to_conformation(core::pose::Pose & p, bool keep_history = false, core::Real update_interval = 0);

} // moves
} // protocols


#endif // INCLUDED_protocols_moves_PyMolMover_HH
