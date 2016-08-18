# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

"""
Brief:   The PyMOLMover class and associated methods.

Remarks: Use in conjunction with PyMOLPyRosettaServer.py.

Authors: Sergey Lyskov & Evan Baugh

Edits:   Labonte
"""

# Imports.
import uuid
import socket
import bz2
from array import array

import rosetta
import rosetta.core.pose
import rosetta.core.scoring.hbonds
import rosetta.core.scoring
import rosetta.core.scoring.dssp
import rosetta.protocols.moves


###############################################################################
# Constants.
X11_COLORS = {
    'black': (0, 0, 0, 0),
    'AntiqueWhite': (1, 250, 235, 215),
    'BlanchedAlmond': (2, 255, 235, 205),
    'BlueViolet': (3, 138, 43, 226),
    'CadetBlue': (4, 95, 158, 160),
    'CornflowerBlue': (5, 100, 149, 237),
    'DarkBlue': (6, 0, 0, 139),
    'DarkCyan': (7, 0, 139, 139),
    'DarkGoldenrod': (8, 184, 134, 11),
    'DarkGray': (9, 169, 169, 169),
    'DarkGreen': (10, 0, 100, 0),
    'DarkGrey': (11, 169, 169, 169),
    'DarkKhaki': (12, 189, 183, 107),
    'DarkMagenta': (13, 139, 0, 139),
    'DarkOliveGreen': (14, 85, 107, 47),
    'DarkOrange': (15, 255, 140, 0),
    'DarkOrchid': (16, 153, 50, 204),
    'DarkRed': (17, 139, 0, 0),
    'DarkSalmon': (18, 233, 150, 122),
    'DarkSeaGreen': (19, 143, 188, 143),
    'DarkSlateBlue': (20, 72, 61, 139),
    'DarkSlateGray': (21, 47, 79, 79),
    'DarkSlateGrey': (22, 47, 79, 79),
    'DarkTurquoise': (23, 0, 206, 209),
    'DarkViolet': (24, 148, 0, 211),
    'DebianRed': (25, 215, 7, 81),
    'DeepPink': (26, 255, 20, 147),
    'DeepSkyBlue': (27, 0, 191, 255),
    'DimGray': (28, 105, 105, 105),
    'DimGrey': (29, 105, 105, 105),
    'DodgerBlue': (30, 30, 144, 255),
    'FloralWhite': (31, 255, 250, 240),
    'ForestGreen': (32, 34, 139, 34),
    'GhostWhite': (33, 248, 248, 255),
    'GreenYellow': (34, 173, 255, 47),
    'HotPink': (35, 255, 105, 180),
    'IndianRed': (36, 205, 92, 92),
    'LavenderBlush': (37, 255, 240, 245),
    'LawnGreen': (38, 124, 252, 0),
    'LemonChiffon': (39, 255, 250, 205),
    'LightBlue': (40, 173, 216, 230),
    'LightCoral': (41, 240, 128, 128),
    'LightCyan': (42, 224, 255, 255),
    'LightGoldenrod': (43, 238, 221, 130),
    'LightGoldenrodYellow': (44, 250, 250, 210),
    'LightGray': (45, 211, 211, 211),
    'LightGreen': (46, 144, 238, 144),
    'LightGrey': (47, 211, 211, 211),
    'LightPink': (48, 255, 182, 193),
    'LightSalmon': (49, 255, 160, 122),
    'LightSeaGreen': (50, 32, 178, 170),
    'LightSkyBlue': (51, 135, 206, 250),
    'LightSlateBlue': (52, 132, 112, 255),
    'LightSlateGray': (53, 119, 136, 153),
    'LightSlateGrey': (54, 119, 136, 153),
    'LightSteelBlue': (55, 176, 196, 222),
    'LightYellow': (56, 255, 255, 224),
    'LimeGreen': (57, 50, 205, 50),
    'MediumAquamarine': (58, 102, 205, 170),
    'MediumBlue': (59, 0, 0, 205),
    'MediumOrchid': (60, 186, 85, 211),
    'MediumPurple': (61, 147, 112, 219),
    'MediumSeaGreen': (62, 60, 179, 113),
    'MediumSlateBlue': (63, 123, 104, 238),
    'MediumSpringGreen': (64, 0, 250, 154),
    'MediumTurquoise': (65, 72, 209, 204),
    'MediumVioletRed': (66, 199, 21, 133),
    'MidnightBlue': (67, 25, 25, 112),
    'MintCream': (68, 245, 255, 250),
    'MistyRose': (69, 255, 228, 225),
    'NavajoWhite': (70, 255, 222, 173),
    'NavyBlue': (71, 0, 0, 128),
    'OldLace': (72, 253, 245, 230),
    'OliveDrab': (73, 107, 142, 35),
    'OrangeRed': (74, 255, 69, 0),
    'PaleGoldenrod': (75, 238, 232, 170),
    'PaleGreen': (76, 152, 251, 152),
    'PaleTurquoise': (77, 175, 238, 238),
    'PaleVioletRed': (78, 219, 112, 147),
    'PapayaWhip': (79, 255, 239, 213),
    'PeachPuff': (80, 255, 218, 185),
    'PowderBlue': (81, 176, 224, 230),
    'RosyBrown': (82, 188, 143, 143),
    'RoyalBlue': (83, 65, 105, 225),
    'SaddleBrown': (84, 139, 69, 19),
    'SandyBrown': (85, 244, 164, 96),
    'SeaGreen': (86, 46, 139, 87),
    'SkyBlue': (87, 135, 206, 235),
    'SlateBlue': (88, 106, 90, 205),
    'SlateGray': (89, 112, 128, 144),
    'SlateGrey': (90, 112, 128, 144),
    'SpringGreen': (91, 0, 255, 127),
    'SteelBlue': (92, 70, 130, 180),
    'VioletRed': (93, 208, 32, 144),
    'WhiteSmoke': (94, 245, 245, 245),
    'YellowGreen': (95, 154, 205, 50),
    'aquamarine': (96, 127, 255, 212),
    'azure': (97, 240, 255, 255),
    'beige': (98, 245, 245, 220),
    'bisque': (99, 255, 228, 196),
    'AliceBlue': (100, 240, 248, 255),
    'blue': (101, 0, 0, 255),
    'blue1': (102, 0, 0, 255),
    'blue2': (103, 0, 0, 238),
    'blue3': (104, 0, 0, 205),
    'blue4': (105, 0, 0, 139),
    'brown': (106, 165, 42, 42),
    'burlywood': (107, 222, 184, 135),
    'chartreuse': (108, 127, 255, 0),
    'chocolate': (109, 210, 105, 30),
    'coral': (110, 255, 127, 80),
    'cornsilk': (111, 255, 248, 220),
    'cyan': (112, 0, 255, 255),
    'firebrick': (113, 178, 34, 34),
    'gainsboro': (114, 220, 220, 220),
    'gold': (115, 255, 215, 0),
    'goldenrod': (116, 218, 165, 32),
    'gray': (117, 190, 190, 190),
    'gray0': (118, 0, 0, 0),
    'gray10': (119, 26, 26, 26),
    'gray100': (120, 255, 255, 255),
    'gray20': (121, 51, 51, 51),
    'gray30': (122, 77, 77, 77),
    'gray40': (123, 102, 102, 102),
    'gray50': (124, 127, 127, 127),
    'gray60': (125, 153, 153, 153),
    'gray70': (126, 179, 179, 179),
    'gray80': (127, 204, 204, 204),
    'gray90': (128, 229, 229, 229),
    'green': (129, 0, 255, 0),
    'green1': (130, 0, 255, 0),
    'green2': (131, 0, 238, 0),
    'green3': (132, 0, 205, 0),
    'green4': (133, 0, 139, 0),
    'honeydew': (134, 240, 255, 240),
    'ivory': (135, 255, 255, 240),
    'khaki': (136, 240, 230, 140),
    'lavender': (137, 230, 230, 250),
    'linen': (138, 250, 240, 230),
    'magenta': (139, 255, 0, 255),
    'maroon': (140, 176, 48, 96),
    'moccasin': (141, 255, 228, 181),
    'navy': (142, 0, 0, 128),
    'orange': (143, 255, 165, 0),
    'orchid': (144, 218, 112, 214),
    'peru': (145, 205, 133, 63),
    'pink': (146, 255, 192, 203),
    'plum': (147, 221, 160, 221),
    'purple': (148, 160, 32, 240),
    'red': (149, 255, 0, 0),
    'red1': (150, 255, 0, 0),
    'red2': (151, 238, 0, 0),
    'red3': (152, 205, 0, 0),
    'red4': (153, 139, 0, 0),
    'salmon': (154, 250, 128, 114),
    'seashell': (155, 255, 245, 238),
    'sienna': (156, 160, 82, 45),
    'snow': (157, 255, 250, 250),
    'snow1': (158, 255, 250, 250),
    'snow2': (159, 238, 233, 233),
    'snow3': (160, 205, 201, 201),
    'snow4': (161, 139, 137, 137),
    'tan': (162, 210, 180, 140),
    'thistle': (163, 216, 191, 216),
    'tomato': (164, 255, 99, 71),
    'turquoise': (165, 64, 224, 208),
    'violet': (166, 238, 130, 238),
    'wheat': (167, 245, 222, 179),
    'white': (168, 255, 255, 255),
    'yellow': (169, 255, 255, 0),
}


###############################################################################
# Classes
class PySocketClient:
    def __init__(self, udp_port=65000, udp_ip='127.0.0.1'):
        """Create socket client targeting the given endpoint.

            udp_port - Target port.
            udp_ip - IP address or hostname.
        """
        # udp_ip provided as keyword argument is previous version,
        # maintaining name although a hostname may be provided.

        (hostname, aliaslist, ipaddrlist) = socket.gethostbyaddr(udp_ip)

        if udp_ip in ipaddrlist:
            # Explict address specified, use this address.
            self.udp_ip, self.udp_port = udp_ip, udp_port
        else:
            # Hostname specified, use resolved ip.
            self.udp_ip, self.udp_port = ipaddrlist[0], udp_port

        # Of course this will not work..., but it will be readjusted
        # automatically.
        self.last_accepted_packet_size = 1024*8  # ...well, maybe next time...
        self.uuid = uuid.uuid4()
        self.sentCount = array('H', [0, 0, 0])  # packet_id, N, count
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    def _send_raw_message(self, msg):
        buf = array('c', self.uuid.bytes)
        buf.extend(self.sentCount.tostring())
        #buf.extend(array('H', [1, 2]).tostring())
        buf.extend(msg)

        self.socket.sendto(buf, (self.udp_ip, self.udp_port))

    def send_message(self, msg):
        count = 1
        if len(msg) > self.last_accepted_packet_size:
            count = len(msg) / self.last_accepted_packet_size + 1

        #print 'Sending messgae in %s packets...' % count

        for i in range(count):
            self.sentCount[1] = i
            self.sentCount[2] = count
            m = msg[i * self.last_accepted_packet_size:
                    (i + 1) * self.last_accepted_packet_size]
            self._send_raw_message(m)

        self.sentCount[0] += 1


###############################################################################
class PyMOLMover(rosetta.protocols.moves.PyMolMover):
    def __init__(
            self,
            keep_history=False,
            update_energy=False,
            energy_type=rosetta.core.scoring.total_score,
            target_host="127.0.0.1",
            target_port=65000):
        rosetta.protocols.moves.PyMolMover.__init__(self)
        self.keep_history(keep_history)
        self.update_energy(update_energy)
        self.energy_type(energy_type)
        self.link = PySocketClient(target_port, target_host)

    # Private methods.
    def _get_pose_name(self, pose):
        #if 'name' in pose.__dict__:
        #    return pose.name[:255]
        if not pose.pdb_info():
            return 'pose'
        else:
            p1, p2, p3 = pose.pdb_info().name().rpartition('.pdb')
            name = p1 or p3

            # Workaround for annoying issue with paths as the name.
            if '/' in name:
                print 'The name "' + name + '" may cause problems for PyMOL.'
                print 'The displayed name will be different in PyMOL.'
                print 'You can change the Pose (PDBInfo) name using:'
                print 'pose.pdb_info().name(<new_name>)'

                name = name.split( '/' )[-1]

            return name[:255]

    def _find_variable_name(target, main_vars):
        """
        Returns all variables in workspace matching target.

        main_vars must be globals() or vars() when called.
        i.e., pose = Pose(); target = pose will return matches = ['pose']
        """
        # Search 'main_vars' for the object matching 'target'.
        # Return all matches.
        matches = [0]
        matches = [index for index, value in main_vars.items()
                                                            if value == target]
        if len(matches) > 1:
            print 'Consider reassignment;',
            print 'multiple variables have this object data.'
        return matches

    def _scale_list(self, input_list):
        """
        Scales an array from 0 to 255.

        Used for scaling a list of values, particularly for coloring.
        """
        mi = min(input_list)
        ma = max(input_list)
        if ma - mi < 1e-100:
            ma += 1e-100

        r = [int((i - mi) * 255. / (ma - mi)) for i in input_list]
        return r

    def _get_hex_per_residue_message(self, pose, per_residue_values,
                                     autoscale=True):
        """
        Converts a list to a hex per-residue message with appropriate residue
        names.
        """
        message = ''
        info = pose.pdb_info()

        # Handy scaling, since we are using hex.
        if autoscale:
            per_residue_values = self._scale_list(per_residue_values)

        # Check if the PDBInfo exists.
        if info is not None and info.nres() != 0:
            # The PDBInfo is defined.
            for i in xrange(len(per_residue_values)):
                chain = info.chain(i + 1)[0]
                res = info.number(i + 1)
                icode = info.icode(i + 1)
                message += '%s%4d%c%02x' % (chain, res, icode,
                                            per_residue_values[i])
        else:
            # PDBInfo is undefined.
            # Rather than seg fault, let's try pose numbering....
            for i in xrange(len(per_residue_values)):
                chain = ' '
                res = i + 1
                message += '%s%4d %02x' % (chain, res, per_residue_values[i])

        return message

    def _get_energies(self, pose, energy_type):
        """Returns a list of energies of type <energy_type> from the pose."""
        # Check to make sure the energies are available.
        if not pose.energies().energies_updated():
            raise IOError('PyMOL_Mover::send_specific_energy: ' +
                         'Energy is not updated; please score the pose first!')

        # Get the proper score type.
        if energy_type == rosetta.end_of_score_type_enumeration:
            # WORKAROUND: If one sets self.energy_type(total_score),
            # self.energy_type() returns instead end_of_score_type_enumeration.
            score_type = rosetta.core.scoring.total_score
        elif isinstance(energy_type, rosetta.core.scoring.ScoreType):
            score_type = energy_type
        elif isinstance(energy_type, str):
            score_type = rosetta.score_type_from_name(energy_type)
        else:
            raise TypeError('energy_type must be a string or ScoreType.')

        return [pose.energies().residue_total_energies(i + 1)[score_type]
                                         for i in xrange(pose.total_residue())]

    def _send_raw_energies(self, pose, energy_type, energies, autoscale=True):
        energy_type = str(energy_type)[:255]
        name = self._get_pose_name(pose)

        if autoscale:
            energies = self._scale_list(energies)

        e = self._get_hex_per_residue_message(pose, energies, autoscale)

        message = 'Ener.bz2' + chr(self.keep_history()) + chr(len(name)) + \
                   name + chr(len(energy_type)) + energy_type + bz2.compress(e)
        self.link.send_message(message)

    def _send_any(self, ptype, pose, data, size=6):
        """
        A generic message sender.

        ptype tags the data with how it is to be interpreted by
        PyMOLPyRosettaServer.py.  ptype MUST be the same on both sides -- the
        same here as in PyMOLPyRosettaServer.py.
        size sets the length of data units.  size MUST be a single digit!
        All data[i] MUST be at least size long!
        """
        to_send = str(size)
        info = pose.pdb_info()
        for i in xrange(1, pose.total_residue() + 1):
            if info is not None and info.nres() != 0:
                pdb = (str(info.number(i)) + info.icode(i) +
                                                        info.chain(i)).rjust(6)
            else:
                pdb = (str(i) + ' A').rjust(6)
            dat = str(data[i - 1])[0:size]
            if not len(dat) == size:
                raise IOError('Error: Not all data is the same size!')

            to_send += '%s%s' % (pdb, dat)

        name = self._get_pose_name(pose)

        message = ptype + chr(self.keep_history()) + chr(len(name)) + name + \
                                                          bz2.compress(to_send)

        self.link.send_message(message)


    ###########################################################################
    # Public methods.
    # Standard methods also available in some form in Rosetta 3's PyMolMover.
    # Send the structure.
    def apply(self, input_pose):
        pose = input_pose.get()  # Necessary for use in Mover containers.

        self.send_structure(pose)

        if self.update_energy():
            self.send_energy(pose, self.energy_type())

    def send_structure(self, pose):
        name = self._get_pose_name(pose)
        #print name

        # Creating message...
        os = rosetta.utility.OStringStream()
        pose.dump_pdb(os)
        #if rosetta._PLATFORM == "cygwin":
        #    message = 'PDB.gzip' + chr(len(name)) + name + \
        #                                               gzip.compress(os.str())
        #else:
        message = 'PDB.bz2 ' + chr(self.keep_history()) + chr(len(name)) + \
                                                  name + bz2.compress(os.str())
        self.link.send_message(message)

    # Energy output.
    def send_energy(self, input_pose, energy_type=rosetta.total_score,
                    label=False, sigs=6):
        """
        Sends cummulative energy score to PyMOL.

        This method will color a pose in PyMOL based on relative residue
        energies.  <energy_type> is a string representation of a specific
        scoring component.
        The <label> option displays <sigs> number of characters for each energy
        with labels on the CA of each residue.

        Examples:
            pymol.send_energy(pose, 'fa_atr')
            pymol.send_energy(pose, label=True)
        """
        energy_list = self._get_energies(input_pose, energy_type)
        self._send_raw_energies(input_pose, energy_type, energy_list)
        if label:
            self._send_any('lbE1.bz2', input_pose, energy_list, sigs)


    # Personalized coloring.
    def send_colors(self, pose, colors={}, default_color='blue'):
        """
        Colors protein, using a color dictionary as map.

        The color dictionary keys are pose residue numbers.  Colors are
        strings.  The default color for any residue not in the dictionary is
        "blue" but can be set to any desired color.

        Examples:
            pymol.send_colors(pose)  # Colors the protein blue.

            color_map = {1:"blue", 2:"red", 5:"purple"}
            pymol.send_colors(pose, color_map)
        """
        try:
            energies = [X11_COLORS[default_color][0]] * pose.total_residue()
        except KeyError:
            print default_color, "is not a valid color. ",
            print 'Type "list_valid_PyMOL_colors()" for a list of valid',
            print 'colors.'
            raise
        for r in colors:
            try:
                energies[r - 1] = X11_COLORS[colors[r]][0]
            except KeyError:
                print colors[r], "is not a valid color. ",
                print 'Type "list_valid_PyMOL_colors()" for a list of valid',
                print 'colors.'
                raise
        self._send_raw_energies(pose, 'X11Colors', energies, autoscale=False)


    ###########################################################################
    # "Luxury" methods.  (Many of which call _send_any() above.)
    # Labeling.
    def label_energy(self, input_pose, energy_type='total_score', sigs=6):
        """
        Displays <sigs> number of characters for each energy with labels on CA.
        """
        energy_list = self._get_energies(input_pose, energy_type)
        self._send_any('lbE1.bz2', input_pose, energy_list, sigs)

    # H-bonds.
    def send_hbonds(self, pose):
        """
        Sends list of hydrogen bonds and displays them in PyMOL.

        Makes use of PyMOL's "distance" function.
        """
        # Check that the energies are updated.
        if not pose.energies().energies_updated():
            raise IOError('PyMOL_Mover::send_hbonds: Energy is not updated;' +
                                                'please score the pose first!')

        # Get the H-bonds.
        hbset = rosetta.core.scoring.hbonds.HBondSet()
        rosetta.core.scoring.hbonds.fill_hbond_set(pose, False, hbset)

        to_send = str(hbset.nhbonds()).rjust(5)
        info = pose.pdb_info()

        # Get the energies.
        energy_list = [hbset.hbond(i + 1).energy()
                                              for i in xrange(hbset.nhbonds())]

        if energy_list:

            max_e = max(energy_list)
            min_e = min(energy_list)

            for i in xrange(1, hbset.nhbonds() + 1):
                hb = hbset.hbond(i)
                accatm = pose.residue(hb.acc_res()).atom_name(hb.acc_atm())
                donatm = pose.residue(hb.don_res()).atom_name(hb.don_hatm())

                # Each H-bond is 6 + 4 + 6 + 4 + 2 = 22 chars.
                if info is not None and info.nres() != 0:
                    acc_res = str(info.number(hb.acc_res()))
                    acc_icode = info.icode(hb.acc_res())
                    acc_chain = info.chain(hb.acc_res())
                    don_res = str(info.number(hb.don_res()))
                    don_icode = info.icode(hb.don_res())
                    don_chain = info.chain(hb.don_res())

                    to_send += (acc_res + acc_icode + acc_chain).rjust(6) + \
                               accatm + \
                               (don_res + don_icode + don_chain).rjust(6) + \
                               donatm + \
                               ('%02x' % int((energy_list[i - 1] - min_e) * 255.
                                                                / (max_e - min_e)))
                else:
                    to_send += str(hb.acc_res()).rjust(5) + \
                               accatm + \
                               str(hb.don_res()).rjust(5) + \
                               donatm + \
                               str(hb.energy).rjust(5)

            name = self._get_pose_name(pose)
            message = 'hbd.bz2 ' + chr(self.keep_history()) + chr(len(name)) + \
                                                       name + bz2.compress(to_send)
            self.link.send_message(message)
        else:
            print "No H-bonds could be determined for your pose!"

    # Secondary-structure assignment using DSSP.
    def send_ss(self, pose, ss=''):
        """
        Sends the DSSP assignment for pose to PyMOL and shows as a cartoon.

        Useful for when you are making moves to a pose that change secondary
        structure, and you wish for PyMOL to display the changes.
        """
        # Get ss.
        if not ss:
            dssp = rosetta.core.scoring.dssp.Dssp(pose)
            dssp.insert_ss_into_pose(pose)
            ss = pose.secstruct()

        # Send it!
        self._send_any(' ss.bz2 ', pose, ss, 1)

    # Polar identity per residue.
    def send_polars(self, pose):
        """
        Colors polar residues red and nonpolar residues blue.
        """
        # Send 1 or 0, if polar or not.
        data = [0] * pose.total_residue()  # Nonpolar by default.
        for i in xrange(1, pose.total_residue() + 1):
            data[i - 1] = int(pose.residue(i).is_polar())
        self._send_any('pol.bz2 ', pose, data, 1)

    # MoveMap DOF info per residue in pose.
    def send_movemap(self, pose, movemap):
        """
        Colors movable regions green and non-movable regions red.
        """
        data = [0] * pose.total_residue()
        # First digit is for bb, second for sc -- 1 = off, 2 = on.
        for i in xrange(1 , pose.total_residue() + 1):
            data[i - 1] = 11 + int(movemap.get_bb(i)) * 10 + \
                               int(movemap.get_chi(i))
        self._send_any('mm1.bz2 ', pose, data, 2)

    # Colors pdb by foldtree.
    def send_foldtree(self, pose, foldtree=''):
        """
        Colors the pose by fold tree information.

        Cutpoints (e.g., the C-termini of protein chains) are colored red.
        Jump points are colored orange. (Unfortunately, no indication of which
        jump point connects to which jump point is given at this time.)
        Loops are colored an assortment of colors other than red or orange.
        All other residues are colored gray.

        See also:
            PyMOL_Mover.view_foldtree_diagram()
        """
        # If not sent, use pose's FoldTree.
        if not foldtree:
            foldtree = pose.fold_tree()

        data = [0] * pose.total_residue()

        # Remove jump data for identification later.
        njump = foldtree.num_jump()
        starts = [0] * njump
        stops = [0] * njump

        in_loops = []

        # List of start and stop of jump edges, i.e., loops.
        for x in range(0, njump):
            loop = foldtree.jump_edge(x + 1)
            s1 = loop.start()
            s2 = loop.stop()
            if s1 < s2:
                starts[x] = s1
                stops[x] = s2
            else:
                starts[x] = s2
                stops[x] = s1

        # Each residue is either a jump point, a cutpoint, in a loop, or else
        # in a regular edge of the fold tree.
        for i in xrange(1, pose.total_residue() + 1):
            # Keeps the identity relative to the jump, not entry.
            if foldtree.is_jump_point(i):  # Jump point residue: color 1.
                data[i - 1] = 1

                # After this point, we will either be entering a loop or
                # leaving a loop.
                for j in range(0 , len(starts)):
                    if i == starts[j]:
                        in_loops.append(j + 1)
                    elif i == stops[j]:
                        in_loops.remove(j + 1)
            elif foldtree.is_cutpoint(i):  # Cutpoint residue: color 0
                data[i - 1] = 0
            elif in_loops:  # Residue in a loop: color varies.
                # Up to 7 loops can easily be viewed.
                # (PyMOLPyRosettaServer.py will need to be edited for more
                # colors for more loops.)
                data[i - 1] = 2 + max(in_loops)
            else:  # All other residues: color 2
                data[i - 1] = 2
        self._send_any('ft1.bz2 ', pose, data, 1)


    # #########################################################################
    # Graphing methods.
    # Send a foldtree diagram.
    def view_foldtree_diagram(self, pose, foldtree=None, scale=10, r=3):
        """
        Draws a 3-D fold tree diagram in the PyMOL viewer window.

        Chains are indicated with a straight line and colored by chain.
        Jumps are represented by "bridges" in unique colors that connect to the
        start and stop jump points on the line, and the bridge is labeled with
        start and stop residues and the jump number.
        Cutpoints are represented by a white and red line bisecting the chain
        and are labeled with the residue number.

        See also:
            PyMOL_Mover.send_foldtree()
        """
        # Default to the pose fold tree.
        if not foldtree:
            foldtree = pose.fold_tree()
        # Check the pose name.
        name = self._get_pose_name(pose)

        # Geometry
        to_send = str(scale).rjust(2)[:2] + str(r)[0]

        # Start with total length -- needed for scaling.
        to_send += str(pose.total_residue()).rjust(4)

        # First, send the chain info.
        # Get the unique chains and indices.
        chain_string = ''.join([pose.pdb_info().chain(i + 1)
                                      for i in xrange(pose.pdb_info().nres())])
        chains = {}
        for i in xrange(len(chain_string)):
            if not chain_string[i] in chains.keys():
                chains[chain_string[i]] = i + 1  # 1-indexed

        # Add to the message.
        # 2 + c*4, where the 2 tells c.
        to_send += str(len(chains)).rjust(2)  # How many chains
        for i in chains.keys():
            to_send += str(chains[i]).rjust(4)

        # Then, add the jump info.
        njump = foldtree.num_jump()
        to_send += str(njump).rjust(2)

        # Add the message.
        # 2 + j*9, where the 2 tells j.
        for j in xrange(njump):
            jump = foldtree.jump_edge(j + 1)
            jump_start = str(jump.start())
            jump_stop = str(jump.stop())
            jump_cut = str(foldtree.cutpoint_by_jump(j + 1))

            to_send += jump_start.rjust(3) + jump_stop.rjust(3) + \
                                                              jump_cut.rjust(3)

        # (2 + 1) + 4 + (2 + c*4) + (2 + j*9)
        message = 'ftd.bz2 ' + chr(self.keep_history()) + chr(len(name)) + \
                                                   name + bz2.compress(to_send)
        self.link.send_message(message)

    def plot_graph(self, name, connect, x_array, y_array, z_array=False,
                   scale=True, axis_color='blue', num=0):
        """
        Plots a list of points in PyMOL and draws axes.

        name is the name of the object in PyMOL.
        connect is as a color string or empty ('') for no color/connection.
        scale turns aces on or off.
        axis_color colors the axes.
        num is the number of internal numberings on the scale.

        The length of the given lists/arrays must be equal.

        Example:
            pymol.plot_graph("Line", "white", [0, 1, 2, 3, 4], [0, 2, 4, 6, 8])
        See also:
            PyMOL_Mover.send_point()
        """
        if not z_array:
            z_array = [0] * len(x_array)
        if not len(x_array) == len(y_array) == len(z_array):
            raise IOError('FAIL: arrays are not the same length!')
        to_send = connect.rjust(7) + str(int(scale)) + axis_color.rjust(7) + \
                                                              str(num).rjust(6)

        for i in range(0, len(x_array)):
            to_send += str(x_array[i])[:9].rjust(9)
            to_send += str(y_array[i])[:9].rjust(9)
            to_send += str(z_array[i])[:9].rjust(9)
        message = 'grp1.bz2' + chr(self.keep_history()) + chr(len(name)) + \
                                                   name + bz2.compress(to_send)
        self.link.send_message(message)


    def send_point(self, name, connect, x, y, z=False, rescale=True,
                   scale=True, axis_color='blue', num=0, banner=''):
        """
        Sends a point to add to a plot already drawn in PyMOL

        Be sure to call PyMOL_Mover.send_graph() first!
        Options are similar to those for PyMOL_Mover.send_graph(); see the
        documentation for that method.
        banner is the point label, if desired.
        rescale determines if the axes are redrawn due to the additional point.

        Example:
            pymol.plot_graph("Line", "white", [0, 1, 2, 3, 4], [0, 2, 4, 6, 8])
            pymol.send_point("Line", "white", 5, 10)
        See also:
            PyMOL_Mover.plot_graph()
        """
        if not z:
            z = 0
        to_send = str(x)[:9].rjust(9) + str(y)[:9].rjust(9) + \
                  str(z)[:9].rjust(9) + connect.rjust(7) + \
                  str(int(rescale)) + str(int(scale)) + axis_color.rjust(7) + \
                  str(num).rjust(6) + banner
        message = 'pnt.bz2 ' + chr(self.keep_history()) + chr(len(name)) + \
                                                   name + bz2.compress(to_send)
        self.link.send_message(message)


###############################################################################
class PyMOLObserver(rosetta.core.pose.PosePyObserver):
    """Responds to general events (changes of geometry and energies) to pose
    and sends updates to PyMOL.

    WARNING: This will slow up resources and cause PyMOL to crawl if enabled
    during protocols with large numbers of moves to pose.

    Usage:
    observer = PyMOL_Observer(pose)  # Construct observer and begin observing.
    observer.add_observer(pose2)  # Begin watching pose 2 also.
    observer.remove_observer(pose)  # Stop watching pose.

    """
    def __init__(self, pose_to_observe=None, keep_history=False,
                 update_energy=False,
                 energy_type=rosetta.core.scoring.total_score):
        rosetta.core.pose.PosePyObserver.__init__(self)
        self.pymol = PyMOL_Mover()
        self.pymol.keep_history(keep_history)
        self.pymol.update_energy(update_energy)
        self.pymol.energy_type(energy_type)
        if pose_to_observe is not None:
            self.add_observer(pose_to_observe)

    def generalEvent(self, event):
        #print 'PyMOL_Observer...'
        #print 'PyMOL_Observer:generalEvent', event.pose
        self.pymol.apply(event.getPose())


###############################################################################
# Helper methods.
def list_valid_PyMOL_colors():
    """
    A simple exposed helper method that prints a sorted list of valid color
    strings for use with PyMOL_Mover methods.
    """
    colors = X11_COLORS.keys()
    colors.sort()
    for color in colors:
        print color


# Providing old names for compatability
PyMOL_Mover = PyMOLMover
PyMOL_Observer = PyMOLObserver
