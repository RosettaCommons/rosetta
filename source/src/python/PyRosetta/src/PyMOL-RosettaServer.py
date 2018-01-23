#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license. The Rosetta software is developed by the contributing
# (c) members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington
# (c) CoMotion, email: license@uw.edu.

"""
@file    PyMOL-RosettaServer.py
@brief   Establishes a link between PyMOL and (Py)Rosetta
@author  Sergey Lyskov, Johns Hopkins University
@edits   Evan Baugh, Brian D. Weitzner, Rebecca Alford, Jason Labonte &
         Julia Koehler Leman
@details This is the script to run to let you view (Py)Rosetta runs inside PyMOL.
         To run it:
             1) Open PyMOL
             2) run /path/to/this/script/PyMOL-RosettaServer.py
             3) Instantiate a PyMOL_Mover in your Python code.
             4) Call PyMOL_Mover.apply(pose)
         See http://www.pyrosetta.org/pymol_mover-tutorial for more info.
         See https://doi.org/10.1371/journal.pone.0021931 for citation.
"""

from __future__ import print_function


# Imports.
import math
import time
import socket
import gzip
import bz2
import threading

# python 2/3 compatibility
try:
    from cStringIO import StringIO
    BytesIO = StringIO
except ImportError:
    from io import StringIO, BytesIO

from array import array

import pymol

from pymol.cgo import *


###############################################################################
# Constants.
COLOR_LIB = {
    'white': [1, 1, 1],
    'yellow': [1, 1, 0],
    'magenta': [1, 0, 1],
    'cyan': [0, 1, 1],
    'red': [1, 0, 0],
    'blue': [0, 0, 1],
    'green': [0, 1, 0],
    'black': [0, 0, 0]}

X11_COLORS = {
    'AliceBlue': (100, 240, 248, 255),
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
    'black': (0, 0, 0, 0),
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
# Exceptions.
class StartUpError(Exception):
    """Exception class for server start-up failures."""
    def __str__(self):
        return "FAILED TO START PyRosetta-PyMOL server." + \
               "\nDo you already have another instance of it running?"


# Other classes.
class PR_UDPServer:
    def __init__(self, udp_ip='127.0.0.1', udp_port=65000):
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        try:
            self.socket.bind((udp_ip, udp_port))
        except:
            raise StartUpError()
        self.buf = {}
        self.last_cleanup_time = time.time()

    def listen(self):
        data, addr = self.socket.recvfrom(1024 * 64)  # 64k buffer
        # print 'Got Message from %s' % str(addr), len(data)

        packet_id = data[:(16 + 2)]
        counts = array('H', data[18:22])  # should be 2 short integers

        # print 'Packet count info:', counts
        if counts[1] == 1:  # only one messgage in pack...
            return array('b', data[22:])  # bz2.decompress(data[22:])
        else:
            if packet_id not in self.buf:
                self.buf[packet_id] = [0., {}]

            c, d = self.buf[packet_id]
            d[counts[0]] = data[22:]
            self.buf[packet_id][0] = time.time()

            # Now, let's check if we can find all the pieces of the message....
            if len(d) == counts[1]:  # Yes, they are all here....
                # print 'Assembling message from %s pieces....' % counts[1]
                m = array('b')
                for i in range(counts[1]):
                    m.extend(array('b', d[i]))
                del self.buf[packet_id]

                # print 'Message is:', len(m), m, d
                # print 'Leftover buffer len:', len(self.buf)
                return m  # bz2.decompress(m)

            else:
                # There are no ready-to-return packets...; However, let's check
                # if the buffer can be cleaned up....
                # Anything older than 10 seconds should be discarded....
                current_time = time.time()
                if current_time - self.last_cleanup_time > 2.:
                    # Cleaning up every 2 s
                    for k in self.buf.keys():
                        if current_time - self.buf[k][0] > 10.0:
                            print('Buffer clean up: %s' % repr(k))
                            del self.buf[k]

                return None


class PR_PyMOLServer:
    def _color_model(self, name, etype, s):
        if etype == 'X11Colors':
            palette = 'X'
        else:
            palette = 'R'

        for i in xrange(0, len(s), 8):
            score = ('%s' % s[(i + 6):(i + 8)])
            color = palette + score
            target = '%s and chain %s and resi %s' % (name, s[i],
                                                      s[(i + 1):(i + 6)])
            # print 'Color: %s, target:%s' % (color, target)
            pymol.cmd.color(color, target)

    # Code for processing commands from PyRosetta goes here within the if
    # statement blocks below.
    def process_packet(self, msg):
        """ Format description:
            bytes 0:8     - packet type
            byte  8       - flags
            byte  9       - len(name)
            bytes 10:...  - name
            bytes ...:... - message itself
        """
        # print("the initial message is {0}".format(msg))
        ptype = msg[:8].tostring()
        flags = msg[8]
        name_len = msg[9]
        name = msg[10:(10 + name_len)].tostring()
        data = msg[10 + name_len:].tostring()
        # print 'Decoy type: %s, name: %s' % (ptype, name)

        # String is just text that we need to print.
        if ptype == 'Text    ':
            print(data.tostring())

        #######################################################################
        # PDB data.
        # String is just a pdb file, no compression.
        elif ptype == 'PDB     ':
            # print 'Getting PDB packet "%s"...' % name
            # print 'Processing pdb...'
            # pymol.cmd.delete(name)
            pymol.cmd.read_pdbstr(data, name, 1)
            # pymol.cmd.show("cartoon", name)
            # pymol.cmd.forward()
            # pymol.cmd.refresh()

        # String is a pdb file with compression.
        elif ptype.startswith(b'PDB.'):
            # print 'Getting PDB packet "%s"...' % name

            # Decompress.
            if ptype.endswith(b'.gzip'):
                s = gzip.GzipFile('', 'r', 0, BytesIO(data)).read()
            elif ptype.endswith(b'.bz2 '):
                s = bz2.decompress(data)

            # print 'Got compressed PDB:', name

            pymol.cmd.read_pdbstr(s, name, flags ^ 1)
            if flags:  # Go to the new frame.
                pymol.cmd.frame(pymol.cmd.count_frames())

        #######################################################################
        # Energy data.
        elif ptype.startswith(b'Ene'):
            # print 'Getting energy packet "%s"...' % name
            e_type_len = ord(data[0])
            e_type = data[1:(1 + e_type_len)].tostring()
            # print 'etype=%s  msg=%s' % (e_type, data)

            # Decompress.
            if ptype.endswith(b'.gzip'):
                s = gzip.GzipFile('', 'r', 0,
                                  StringIO(data[(1 + e_type_len):])).read()
            elif ptype.endswith(b'.bz2'):
                s = bz2.decompress(data[(1 + e_type_len):])
            # print 'Compression stats: %s-->%s' % \
            #       (len(data[(1 + e_type_len):]), len(s))
            # print("the decompressed data is {0} \n".format(s))
            try:
                # print 'Coloring model:', name
                self._color_model(name, e_type, s)
            except pymol.parsing.QuietException:
                print("Coloring failed...")
                print("Did you forget to send the pose geometry first?")

        # Label energy per residue.
        elif ptype == 'lbE1.bz2':
            data = bz2.decompress(data)
            # This int must be a single digit, the size of data pieces.
            size = int(data[0])
            data = data[1:]
            try:
                for i in range(0, len(data), size + 6):
                    dat = data[(i + 6):(i + 6 + size)]
                    sel = '%s and name ca and chain %s and resi %s' % \
                        (name, data[i + 5], data[i:(i + 5)].strip())
                    pymol.cmd.label(sel, dat)
            except pymol.parsing.QuietException:
                print("Commands failed...")
                print("Did you forget to send the pose geometry first?")

        #######################################################################
        # Display membrane planes
        elif ptype.startswith('Mem'):
            e_type_len = ord(data[0])
            e_type = data[1:(1 + e_type_len)].tostring()
            pymol.cmd.delete(name + 'membrane_planes')

            pymol.cmd.delete('membrane_planes')

            # Decompress data (currently only supporting .gzip from C++ end)
            if ptype.endswith('.gzip'):
                s = gzip.GzipFile('', 'r', 0,
                                  StringIO(data[(1 + e_type_len):])).read()

                # Read top points, bottom points, and normal coordinate
                mem_data = s.split(',')

                # Read in center position
                center = XYZCoord()
                center.x = float(mem_data[0])
                center.y = float(mem_data[1])
                center.z = float(mem_data[2])

                # Read in normal vector
                normal = XYZCoord()
                normal.x = float(mem_data[3])
                normal.y = float(mem_data[4])
                normal.z = float(mem_data[5])

                # Read in thickness
                thickness = float(mem_data[6])

                # Read in radius of gyration
                rg = float(mem_data[7])

                # Compute position of membrane planes
                planes = compute_plane_positions(center, normal, thickness, rg)

                # Draw membrane planes from points
                draw_membrane_planes(planes, normal)

        # #####################################################################
        # Display hydrogen bonds.
                # Energy data.
        elif ptype.startswith('hbd'):
            # print 'etype=%s  msg=%s' % (e_type, data)

            # Decompress.
            if ptype.endswith('.gzip'):
                data = gzip.GzipFile('', 'r', 0, StringIO(data)).read()
                # print("data is {0}\n".format(data))
            elif ptype == 'hbd.bz2 ':
                data = bz2.decompress(data)
            # First 5 characters are the # of H-bonds.
            nhbonds = data[:5]
            data = data[5:]  # 22 char per H-bond: 6+4 + 6+4 + 2
            try:
                pymol.cmd.delete(name + '_hbonds')
                for i in xrange(int(nhbonds)):
                    c = 22 * i
                    # Acceptor atom
                    acc_res = data[c:(c + 5)].strip()
                    acc_chain = data[c + 5]
                    acc_name = data[(c + 6):(c + 10)].strip()
                    # Donor atom
                    don_res = data[(c + 10):(c + 15)].strip()
                    don_chain = data[c + 15]
                    don_name = data[(c + 16):(c + 20)].strip()

                    # Make selection.
                    hbname = 'hb_' + acc_res + acc_chain + acc_name + '_' + \
                        don_res + don_chain + don_name + '_' + name
                    hbname = "hb_name"
                    pymol.cmd.distance(hbname,
                                       name + ' and chain ' + acc_chain +
                                       ' and res ' + acc_res + ' and name ' +
                                       acc_name,
                                       name + ' and chain ' + don_chain +
                                       ' and res ' + don_res + ' and name ' +
                                       don_name)
                    pymol.cmd.color('Hb%s' % data[(c + 20):(c + 22)], hbname)
                pymol.cmd.sync()  # Ensures above is complete before continuing
                pymol.cmd.hide('labels', 'hb_*_' + name)
                pymol.cmd.group(name + '_hbonds', 'hb_*_' + name)
            except pymol.parsing.QuietException:
                print("Commands failed...")
                print("Did you forget to send the pose geometry first?")

        #######################################################################
        # Update display of secondary structure.
        elif ptype.startswith(' ss'):
            # print 'etype=%s  msg=%s' % (e_type, data)

            # Decompress.
            if ptype.endswith('.gzip'):
                data = gzip.GzipFile('', 'r', 0, StringIO(data)).read()
            elif ptype == ' ss.bz2 ':
                data = bz2.decompress(data)

            size = int(data[0])
            data = data[1:]
            ss_map = {'H': 'H', 'E': 'S', 'L': 'L'}
            try:
                for i in xrange(0, len(data), size + 6):
                    dat = ss_map[data[(i + 6):(i + 6 + size)]]
                    sel = '%s and chain %s and resi %s' % \
                          (name, data[i + 5], data[i:(i + 5)].strip())
                    pymol.cmd.alter(sel, 'ss=\'' + dat + '\'')
                    pymol.cmd.show('cartoon', name)
            except pymol.parsing.QuietException:
                print("Commands failed...")
                print("Did you forget to send the pose geometry first?")

        #######################################################################
        # Color by a boolean value for polar residues.
        elif ptype.startswith('pol'):
            # print 'etype=%s  msg=%s' % (e_type, data)

            # Decompress.
            if ptype.endswith('.gzip'):
                data = gzip.GzipFile('', 'r', 0, StringIO(data)).read()
            elif ptype == 'pol.bz2 ':
                data = bz2.decompress(data)

            # print("in the pol the data is {0}\n".format(data))
            size = int(data[0])
            data = data[1:]
            try:
                for i in range(0, len(data), size + 6):
                    dat = data[(i + 6):(i + 6 + size)]
                    sel = '%s and chain %s and resi %s' % \
                          (name, data[i + 5], data[i:(i + 5)].strip())
                    # print("dat is {0}\nsel is {1}\n".format(dat, sel))
                    color = 'blue'
                    if int(dat):
                        color = 'red'
                    pymol.cmd.color(color, sel)
            except pymol.parsing.QuietException:
                print("Commands failed...")
                print("Did you forget to send the pose geometry first?")

        #######################################################################
        # Color by MoveMap DOF.
        elif ptype.startswith('mm'):
            # print 'etype=%s  msg=%s' % (e_type, data)

            # Decompress.
            if ptype.endswith('.gzip'):
                data = gzip.GzipFile('', 'r', 0, StringIO(data)).read()
            elif ptype == 'mm1.bz2 ':
                data = bz2.decompress(data)

            size = int(data[0])
            data = data[1:]
            try:
                pymol.cmd.remove('hydro')
                bb = '(name N or name CA or name C or name O)'
                for i in range(0, len(data), size + 6):
                    dat = data[(i + 6):(i + 6 + size)]
                    sel = '%s and chain %s and resi %s' % \
                          (name, data[i + 5], data[i:(i + 5)].strip())
                    # First digit is for bb, second for sc.
                    # 1 = off, 2 = on
                    # bb
                    color = 'red'
                    if int(dat[0]) - 1:
                        color = 'green'
                    pymol.cmd.color(color, sel + ' and ' + bb)
                    # sc
                    color = 'red'
                    if int(dat[1]) - 1:
                        color = 'green'
                    pymol.cmd.color(color, sel + ' and not ' + bb)
            except pymol.parsing.QuietException:
                print("Commands failed...")
                print("Did you forget to send the pose geometry first?")

        #######################################################################
        # Color by foldtree edges.
        elif ptype.startswith('ft1'):
            # print 'etype=%s  msg=%s' % (e_type, data)

            # Decompress.
            if ptype.endswith('.gzip'):
                data = gzip.GzipFile('', 'r', 0, StringIO(data)).read()
            elif ptype == 'ft1.bz2 ':
                data = bz2.decompress(data)

            size = int(data[0])
            data = data[1:]
            try:
                for i in range(0, len(data), size + 6):
                    dat = data[(i + 6):(i + 6 + size)]
                    sel = '%s and chain %s and resi %s' % \
                          (name, data[i + 5], data[i:(i + 5)].strip())
                    # Color is based on int sent.
                    #     Cutpoint residue = 0; color = 2000 (red)
                    #     Jump point residue = 1; color = 2100 (orange)
                    #     Normal edge residue = 2; color = 2200 (gray)
                    #     Residue in a loop = >2; color >= 2300 (varies)
                    pymol.cmd.color(str(int(dat) * 100 + 2000), sel)
            except pymol.parsing.QuietException:
                print("Commands failed...")
                print("Did you forget to send the pose geometry first?")

        #######################################################################
        # Foldtree diagram.
        elif ptype == 'ftd.bz2 ':
            data = bz2.decompress(data)

            # Delete previous fold trees and jumps.
            # (Later, make simultaneous viewing.)
            pymol.cmd.delete('jumps_' + name + ' or foldtree_' + name)

            # Get the scale.
            scale = int(data[:2])
            r = int(data[2])
            data = data[3:]

            # Process the chains.
            total = float(data[:4])
            data = data[4:]

            # Get the start points.
            nchains = int(data[:2])
            data = data[2:]
            chains = []
            for i in xrange(nchains):
                chains.append(int(data[:4]))
                data = data[4:]

            # Visualize the chains.
            chains.append(total)
            for i in xrange(len(chains) - 1):
                connect = str(25 + i)
                add_point('foldtree_' + name,
                          [chains[i] / total * scale, 0, 0], connect)
                add_point('foldtree_' + name,
                          [chains[i + 1] / total * scale, 0, 0], connect)

            # The number of jumps
            njump = int(data[:2])
            data = data[2:]

            # Use for size
            jumps = []
            for j in xrange(njump):
                # Store jump start, stop, and cut.
                jumps.append([int(data[:3]), int(data[3:6]), int(data[6:9])])
                data = data[9:]

            # Visualize the jumps and cuts.
            for j in range(njump):
                # Change colors later.
                connect = str(8 + j)  # Magic number for coloring

                # Used frequently below.
                th = 2 * math.pi/njump * j
                cylx = r * math.cos(th)
                cyly = r * math.sin(th)
                start = jumps[j][0] - 1
                stop = jumps[j][1] - 1
                cut = jumps[j][2] - 1

                jump = 'jump_' + str(j + 1) + '_' + name

                # Draw the jump as a bridge.
                # Start bridge at jump point.
                add_point(jump, [start/total * scale, 0, 0], connect)
                # Up one r.
                add_point(jump, [start/total * scale, cylx, cyly], connect,
                          False, False, '', 0, str(jumps[j][0]))
                # Over to center of bridge and name the jump.
                add_point(jump, [(start + stop)/total / 2 * scale, cylx, cyly],
                          connect, False, False, '', 0, 'j' + str(j + 1))
                # Over to connecting jump point.
                add_point(jump, [stop/total * scale, cylx, cyly], connect,
                          False, False, '', 0, str(jumps[j][1]))
                # Down one r to stop bridge.
                add_point(jump, [stop/total * scale, 0, 0], connect)

                # Draw the cutpoint.
                # Up 1/2 r from cutpoint.
                add_point(jump, [cut/total * scale, cylx/2, cyly/2], '', False,
                          False, '', 0, str(jumps[j][2]))
                # Down 1 r.
                add_point(jump, [cut/total * scale, -cylx/2, -cyly/2], 'red')

            # Group all jumps, view centered on them.
            if njump:
                pymol.cmd.group('jumps_' + name, 'jump_*_' + name)
                pymol.cmd.label('jumps_' + name, 'resn')
                pymol.cmd.center('jumps_' + name)

        #######################################################################
        # Generate a graph from data sent.
        elif ptype == 'grp1.bz2':
            data = bz2.decompress(data)
            # First 21 data char are options.
            options = data[:21]
            data = data[21:]
            connect = options[:7].strip()
            scale = bool(int(options[7:8]))
            axis_color = options[8:15].strip()
            num = int(options[15:])

            x_array = [0] * (len(data)/27)
            y_array = [0] * (len(data)/27)
            z_array = [0] * (len(data)/27)
            for i in range(0, len(data), 27):
                x_array[i/27] = float(data[i:(i + 9)].strip())
                y_array[i/27] = float(data[(i + 9):(i + 18)].strip())
                z_array[i/27] = float(data[(i + 18):(i + 27)].strip())
            plot3d(name, x_array, y_array, z_array, connect, scale,
                   axis_color, num)

        # Add a point to graph data.
        elif ptype == 'pnt.bz2 ':
            data = bz2.decompress(data)
            # First 27 are data; last 22+ are options and banner.
            add_point(name,
                      [float(data[0:9].strip()), float(data[9:18].strip()),
                       float(data[18:27].strip())],
                      data[27:34].strip(), bool(int(data[34:35])),
                      bool(int(data[35:36])), data[36:43].strip(),
                      int(data[43:49]), data[49:])

        #######################################################################
        # Template for translating new packers
        # ptype tags the data with how it is to be interpreted here.
        # ptype MUST be the same on both sides -- the same here as in
        # the PyMOL_Mover code.
        # size sets the length of data units.  size MUST be a single digit!
        # All data[i] MUST be at least size long!
        elif ptype == 'temp.bz2':
            data = bz2.decompress(data)
            size = int(data[0])
            data = data[1:]
            try:
                for i in range(0, len(data), size + 6):
                    dat = data[(i + 6):(i + 6 + size)]
                    sel = '%s and chain %s and resi %s' % \
                          (name, data[i + 5], data[i:(i + 5)].strip())
                    # dat-dependent commands to be performed on sel go here.
            except pymol.parsing.QuietException:
                print("Commands failed...")
                print("Did you forget to send the pose geometry first?")

        #######################################################################
        else:
            print('Unknown packet type: %s, - ignoring...' % ptype)


###############################################################################
# Membrane
class XYZCoord:
    """
    Class for storing xyz coord or just a triple of real values
    """

    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "(%8.3f, %8.3f, %8.3f)" % (self.x, self.y, self.z)


class PlanePoints:
    """
    Class for storing position of membrane planes
    """

    def __init__(self):
        self.upper_points = []
        self.lower_points = []


class Matrix:
    """
    Class for storing a 3x3 matrix (used for rotation matrices
    in membranes)
    """

    def __init__(self, row1=[], row2=[], row3=[]):
        self.row1 = row1
        self.row2 = row2
        self.row3 = row3


def add_vectors(v1, v2):
    """
    Add a pair of vectors represented as XYZcoord classes (see above)
    """
    x = v1.x + v2.x
    y = v1.y + v2.y
    z = v1.z + v2.z

    vf = XYZCoord(x, y, z)
    return vf


def scalar_multiply(v, n):
    """
    Multiply a vector by some scalar factor
    """
    return XYZCoord(v.x * n, v.y * n, v.z * n)


def normalize(normal):
    """
    Ensure that the provided normal is a unit normal_vector
    """

    mag = math.sqrt(normal.x * normal.x + normal.y * normal.y +
                    normal.z * normal.z)
    x = normal.x / mag
    y = normal.y / mag
    z = normal.z / mag

    return XYZCoord(x, y, z)


def rotation_matrix_radians(axis, theta):
    """
    Compute and return the rotation matrix associated with rotation about
    the given axis by theta radians
    """
    x = math.sqrt(math.pow(axis.x, 2) + math.pow(axis.y, 2) +
                  math.pow(axis.z, 2))

    a = math.cos(theta/2.0)
    b = -axis.x/x * math.sin(theta/2.0)
    c = -axis.y/x * math.sin(theta/2.0)
    d = -axis.z/x * math.sin(theta/2.0)

    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

    row1 = [aa + bb - cc - dd, 2 * (bc - ad), 2 * (bd + ac)]
    row2 = [2 * (bc + ad), aa+cc-bb-dd, 2*(cd-ab)]
    row3 = [2 * (bd - ac), 2*(cd+ab), aa+dd-bb-cc]

    return Matrix(row1, row2, row3)


def multiply_matrix(matrix, vector):
    """
    multiply matrix by 3D vector
    """

    x = matrix.row1[0] * vector.x + matrix.row1[1] * vector.y + \
        matrix.row1[2] * vector.z

    y = matrix.row2[0] * vector.x + matrix.row2[1] * vector.y + \
        matrix.row2[2] * vector.z

    z = matrix.row3[0] * vector.x + matrix.row3[1] * vector.y + \
        matrix.row3[2] * vector.z

    vf = XYZCoord(x, y, z)
    return vf


def rotate_vector(vector, theta, axis):
    """
    Rotate some vector about a specified axis by some angle theta
    """

    rot_matrix = rotation_matrix_radians(axis, theta)
    rotated_vec = multiply_matrix(rot_matrix, vector)
    return rotated_vec


def compute_plane_positions(center, normal, thickness, rg, npoints=4):
    """
    Compute the position of the membrane planes based on provided
    normal, center, thickness, and radius of gyration
    """
    # Create a new class to store plane points
    plane_position = PlanePoints()

    # Pick a set of angles defining the planes based on no. of desired points
    angles = []
    for i in range(npoints):
        angle = ((i+1)*2*math.pi)/npoints
        angles.append(angle)

    # Pick an arbitrary orthogonal unit vector
    tolerance = math.pow(10, -7)
    p = XYZCoord(0, 0, 0)
    if (math.fabs(normal.x + normal.y) < tolerance):
        p.x = -normal.y - normal.z
        p.y = normal.x
        p.z = normal.x
    else:
        p.x = normal.z
        p.y = normal.z
        p.z = -normal.x - normal.y

    # Scale point by radius of gyration
#    p = XYZCoord(2*rg*p.x, 2*rg*p.y, 2*rg*p.z)
    p = XYZCoord(20 * 15 * p.x, 20 * 15 * p.y, 20 * 15 * p.z)

    # Project center onto the normal by thickness to compute upper and lower
    t = thickness
    upper = add_vectors(center, XYZCoord(normal.x * t, normal.y * t,
                                         normal.z * t))
    lower = add_vectors(center, XYZCoord(-normal.x * t, -normal.y * t,
                                         -normal.z * t))

    # For remaining angles, rotate p about the normal vector
    for i in range(len(angles)):
        new_point = add_vectors(rotate_vector(p, angles[i], normal), upper)
        plane_position.upper_points.append(new_point)

    for i in range(len(angles)):
        new_point = add_vectors(rotate_vector(p, angles[i], normal), lower)
        plane_position.lower_points.append(new_point)

    return plane_position


def draw_membrane_planes(plane_points, normal_vector):
    """
    Draw CGO Planes Representing the upper & lower membrane planes
    Adapted from Evan Baugh's Draw Object code by Rebecca Alford to work
    with the framework
    """
    # Settings
    name = 'membrane_planes'
    color = XYZCoord(0.50, 0.50, 0.50)

    normal = normalize(normal_vector)

    # Create the top plane
    top_plane = [
        BEGIN, TRIANGLE_FAN,
        ALPHA, 0.5,
        COLOR, color.x, color.y, color.z,
        NORMAL, normal.x, normal.y, normal.z,
        ]

    for i in plane_points.upper_points:

        top_plane.append(VERTEX)
        top_plane.append(i.x)
        top_plane.append(i.y)
        top_plane.append(i.z)

    top_plane.append(END)

    # Create the bottom plane
    bottom_plane = [
        BEGIN, TRIANGLE_FAN,
        ALPHA, 0.5,
        COLOR, color.x, color.y, color.z,
        NORMAL, normal.x, normal.y, normal.z,
        ]

    for i in plane_points.lower_points:
        bottom_plane.append(VERTEX)
        bottom_plane.append(i.x)
        bottom_plane.append(i.y)
        bottom_plane.append(i.z)

    bottom_plane.append(END)

    # Display the upper & lower planes
    cmd.load_cgo(top_plane + bottom_plane, 'membrane_planes')


###############################################################################
# Graphing methods.
def make_axis(name, ends, dr, scale=False, axis_color='', num=0):
    """
    Make an axis 'name' from 'ends[0]' to 'ends[1]' on the 'dr' (direction).
    'scale' bool tells to label axis or not; 'axis_color' determines the axis
    color; 'num' allows intermediate scale points.
    """
    # Quit if bad data was fed.
    if ends[0] == ends[1]:
        return

    # Ends must be [max, min].
    d_max = ends[0]
    d_min = ends[1]

    # Create or alter axis object.
    if name in pymol.cmd.get_names():
        pymol.cmd.remove(name)
    else:
        pymol.cmd.create(name, None)

    # Add max and min points; connect them.
    pymol.cmd.pseudoatom(name,
                         name='max',
                         resn=str(d_max)[:5],
                         pos=[dr[0]*d_max, dr[1]*d_max, dr[2]*d_max],
                         color=axis_color)
    pymol.cmd.pseudoatom(name,
                         name='min',
                         resn=str(d_min)[:5],
                         pos=[dr[0]*d_min, dr[1]*d_min, dr[2]*d_min],
                         color=axis_color)
    pymol.cmd.bond(name + ' and name max', name + ' and name min')

    # Add number of intermediate atoms.
    for i in range(1, num + 1):
        val = (d_max - d_min)*float(i)/(num+1) + d_min
        pymol.cmd.pseudoatom(name,
                             name=str(i),
                             resn=str(val)[:5],
                             pos=[dr[0]*val, dr[1]*val, dr[2]*val],
                             color=axis_color)

    if scale:
        scale_axes(name)


def scale_axes(name='x_axis or y_axis or z_axis'):
    """
    Labels the selection name.
    Used here for convenience: axis point value str in resn.
    """
    pymol.cmd.hide('label', name)
    pymol.cmd.label(name, 'resn')


def get_ends(data):
    """
    Determines the ends for axes from list.
    Defaults to min, max, but if not pos/neg.
    Determines based on closest point.
    """
    d_max = max(data)
    if d_max > 0:
        d_max = d_max * 1.1
    else:
        d_max = abs(d_max * .1)
    d_min = min(data)
    if d_min < 0:
        d_min = d_min * 1.1
    else:
        d_min = -abs(d_min * .1)
    return [d_max, d_min]


def extract_coords(name):
    """Returns the coords of a data container as three arrays."""
    data = pymol.cmd.get_model(name)
    size = pymol.cmd.count_atoms(name)
    x_array = [0] * size
    y_array = [0] * size
    z_array = [0] * size
    for atom in data.atom:
        i = atom.index - 1
        x_array[i] = atom.coord[0]
        y_array[i] = atom.coord[1]
        z_array[i] = atom.coord[2]
    return x_array, y_array, z_array


def make_data(name, x_array, y_array, z_array, data_color='red'):
    """Generate a data object from three arrays."""
    # Erase old name, create new.
    pymol.cmd.delete(name)
    pymol.cmd.create(name, None)

    # Produce points.
    for i in range(0, len(x_array)):
        pymol.cmd.pseudoatom(name,
                             name=str(i + 1),
                             pos=[x_array[i], y_array[i], z_array[i]],
                             color=data_color)


def add_point(name, point, connect='', rescale=False, scale=False,
              axis_color='', num=0, banner=''):
    """
    Adds a point to existing data, reconnecting points optionally.
    Adds 'point' to 'name' and connects it with color 'connect' (empty for no
    connection).
    'rescale' will rescale the axes; 'scale' will label the scales;
    'axis_color' determines the axis color; 'num' the intermediate scales,
    'banner' is a resn name at 'point'.
    """
    # If 'name' doesn't exist, create it!
    created = False
    if name not in pymol.cmd.get_names():
        pymol.cmd.create(name, None)
        created = True
    ind = pymol.cmd.count_atoms(name) + 1
    pymol.cmd.pseudoatom(name,
                         name=str(ind),
                         resn=banner,
                         pos=point,
                         color=connect)

    if connect:
        pymol.cmd.bond(name + ' and name ' + str(ind),
                       name + ' and name ' + str(ind - 1))

    if rescale:
        rescale_cartesian(name, scale, axis_color, num)

    if created and rescale:
        pymol.cmd.group(name + '_data', name + ' or axes_' + name)


def rescale_cartesian(data, scale, axis_color, num):
    """
    Rescales the 'data' axes (x, y, z) based on 'data' data:
        'num', 'axis_color', 'scale'
    """
    x_array, y_array, z_array = extract_coords(data)
    make_axis('x_axis_' + str(data), get_ends(x_array), [1, 0, 0], scale,
              axis_color, num)
    make_axis('y_axis_' + str(data), get_ends(y_array), [0, 1, 0], scale,
              axis_color, num)
    make_axis('z_axis_' + str(data), get_ends(z_array), [0, 0, 1], scale,
              axis_color, num)
    pymol.cmd.group('axes_' + str(data), 'x_axis_' + str(data) +
                    ' or y_axis_' + str(data) + ' or z_axis_' + str(data))
    pymol.cmd.reset()


def connect_points(name, color='red'):
    """Draws colored lines between consecutive points."""
    for i in range(2, pymol.cmd.count_atoms(name) + 1):
        pymol.cmd.bond(name + ' and name ' + str(i),
                       name + ' and name ' + str(i - 1))
    pymol.cmd.color(color, name)


def plot3d(name, x_array, y_array, z_array=[], connect='red', scale=True,
           axis_color='blue', num=0):
    """
    Default plotting tool with optional z points, connection, and axis color.
    """
    if not z_array:
        z_array = [0] * len(x_array)
    if not len(x_array) == len(y_array) == len(z_array):
        raise IOError('FAIL: Arrays not the same length.')

    # Make the data.
    make_data(name, x_array, y_array, z_array, connect)

    # Make axes.
    rescale_cartesian(name, scale, axis_color, num)

    pymol.cmd.group(name + '_data', name + ' or axes_' + name)

    if connect:
        connect_points(name, connect)


###############################################################################
# PyMOL spectrum coloring.
def set_spectrum(low='blue', high='red'):
    """
    Sets the energy coloring spectrum using the above COLOR_LIB dictionary.
    """
    print('New spectrum:', low, '==>', high)
    low = COLOR_LIB[low]
    high = COLOR_LIB[high]
    # cust = [0] * 3
    # for j in range(0, 3):
    #    cust[j] = high[j] - low[j]
    cust = [high[j] - low[j] for j in range(3)]

    def x(i, el):
        cust[el] * (i - 255*(cust[el] < 0))/255.

    for i in range(256):
        pymol.cmd.set_color('R%02x' % i, [x(i, 0), x(i, 1), x(i, 2)])


###############################################################################
# Main PyMOLPyRosettaServer.py routines.
def main(ip, port):
    print('PyMOL <---> PyRosetta link started!')
    print('at', ip, 'port', port)

    udp_serv = PR_UDPServer(ip, port)
    PS = PR_PyMOLServer()
    while True:
        r = udp_serv.listen()
        if r:
            # print len(r)
            PS.process_packet(r)

    s.close()


def start_rosetta_server(ip='', port=65000):
    if not ip:
        ip = socket.gethostbyname(socket.gethostname())
        if ip == '127.0.0.1':
            print("Unable to automatically determine your IP address.",
                  end=' ')
            print("Please specify it manually.", end=' ')
            print("e.g., start_rosetta_server 192.168.0.1")
            return

    thread = threading.Thread(target=main, args=[ip, port])
    thread.setDaemon(1)
    thread.start()


###############################################################################
# Create our own color spectrum for PyMOL.
for i in range(256):
    pymol.cmd.set_color('R%02x' % i, [i/255., 0, 1 - i/255.])
for j in range(256):
    pymol.cmd.set_color('Hb%02x' % j, [j/255., 1, 0])
for i in X11_COLORS:
    pymol.cmd.set_color('X%02x' % X11_COLORS[i][0], X11_COLORS[i][1:])

# Add commands to PyMOL.
pymol.cmd.extend('scale_axes', scale_axes)
pymol.cmd.extend('extract_coords', extract_coords)
pymol.cmd.extend('make_data', make_data)
pymol.cmd.extend('add_point', add_point)
pymol.cmd.extend('plot3d', plot3d)

pymol.cmd.extend('set_spectrum', set_spectrum)
pymol.cmd.extend('start_rosetta_server', start_rosetta_server)

# Begin server connection.
start_rosetta_server('127.0.0.1', 65000)

# To use PyMOLPyRosettaServer.py over a network, uncomment the line below and
# set the first argument to your IP address:
# start_rosetta_server('192.168.0.1', 65000)
