#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   PyMOLPyRosettaServer.py
## @brief
## @author Sergey Lyskov, Johns Hopkins University


# This is the script to run to let you want Rosetta runs inside PyMOL.  To run it:
# 1) Open pymol
# 2) run /path/to/this/script/PyMOLPyRosettaServer.py
# 3) Start a Rosetta run with the flag -show_simulation_in_pymol X, where X is how many seconds pass between PyMOL updates.  5 is default; low values may overload your computer.

import time, socket, gzip, bz2, threading
from cStringIO import StringIO
from array import array


import pymol

#from NetLink import PR_UDPServer
# ^^^ this does not work on CygWin PyMOL so we just add our code here...

class PR_UDPServer:
    def __init__(self, udp_ip = '127.0.0.1', udp_port=65000):
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        try:
            self.socket.bind( (udp_ip, udp_port) )
        except:
            print 'FAILED TO START PyRosetta-PyMOL server, do you already have another instance of it running?'
            raise 'FAILED TO START PyRosetta-PyMOL server, do you already have another instance of it running?'
        self.buf = {}
        self.last_cleanup_time = time.time()



    def listen(self):
        data, addr = self.socket.recvfrom(1024*64)  # 64k buffer
        #print 'Got Message from %s' % str(addr), len(data)

        packet_id = data[:16+2]
        counts = array('H', data[18:22])  # should be 2 short integer

        #print 'Packet count info:', counts
        if counts[1] == 1:  # only one messgage in pack...
            return array('c', data[22:]) #bz2.decompress(data[22:])
        else:
            if packet_id not in self.buf: self.buf[packet_id] = [0., {}]

            c, d = self.buf[packet_id]
            d[counts[0]] = data[22:]
            self.buf[packet_id][0] = time.time()

            # now, lets check if we can find all the pieces of this message...
            if len(d) == counts[1]:  # yes, they all here...
                #print 'Asseblinbg messge from %s pieces...' % counts[1]
                m = array('c')
                for i in range(counts[1]):
                    m.extend( array('c', d[i]) )
                del self.buf[packet_id]

                #print 'Messge is:', len(m), m, d
                #print 'leftover buffer len:', len(self.buf)
                return m #bz2.decompress(m)

            else:
                # there is no ready-to-return packets... however lets check if buffer can be cleaned up...
                # anything older then 10second should be discarded...
                current_time = time.time()
                if current_time - self.last_cleanup_time > 2.: # cleaning up every 2s

                    for k in self.buf.keys():
                        if current_time - self.buf[k][0] > 10.0 :
                            print 'Buffer clean up: %s' % repr(k)
                            del self.buf[k]

                return None



class PR_PyMOLServer:

    def color_model(self, name, etype, s):
        if etype == 'X11Colors': palette = 'X'
        else: palette = 'R'

        for i in range(0, len(s), 7):
            pymol.cmd.color( palette +  ('%s' % s[i+5:i+7]), '%s and chain %s and resi %s' % (name, s[i], s[i+1:i+5]))


    def processPacket(self, msg):
        ''' Format description:
            bytes 0:8   - packet type
            byte  8     - flags
            byte  9     - len(name)
            bytes 10:... - name
            bytes ...:... - mesage its self
        '''
        ptype = msg[:8].tostring()
        flags = ord(msg[8])
        name_len = ord(msg[9])
        #print 'name_len:', name_len
        name = msg[10: 10+name_len].tostring()
        data = msg[10+name_len:]  #.tostring()
        #print 'Decoy type: %s, name: %s' % (ptype, name)

        if ptype == 'PDB     ':   # string is just a pdb file, no compression
            #print 'Getting PDB packet "%s"...' % name
            #print 'Processing pdb...'
            #pymol.cmd.delete(name)
            pymol.cmd.read_pdbstr(data, name, 1)
            #pymol.cmd.show("cartoon", name)
            #pymol.cmd.forward()
            #pymol.cmd.refresh()

        if ptype == 'Text    ':   # string is just text string that we need to print
            print data.tostring()

        elif ptype == 'PDB.gzip':   # string is just a pdb file, gzip compression
            #print 'Getting PDB.gzip packet "%s"...' % name
            pymol.cmd.read_pdbstr(gzip.GzipFile('', 'r', 0, StringIO(data)).read(), name, flags ^ 1)
            if flags: pymol.cmd.frame( pymol.cmd.count_frames() )


        elif ptype == 'PDB.bz2 ':   # string is just a pdb file, bz2 compression
            pymol.cmd.read_pdbstr(bz2.decompress(data), name, flags ^ 1)
            if flags: pymol.cmd.frame( pymol.cmd.count_frames() )


        elif ptype == 'Ener.bz2':   # energies info, bz2 compression
            #print 'Getting Ene2.bz2 packet...'
            e_type_len = ord(data[0])
            e_type = data[1: 1 + e_type_len].tostring()
            s = bz2.decompress( data[1+e_type_len:] )
            #print 'Compression stats: %s-->%s' % (len(data[1+e_type_len:]), len(s) )
            try:
                #for i in range(0, len(s), 7):
                #    pymol.cmd.color('R%s' % s[i+5:i+7], '%s and chain %s and resi %s' % (name, s[i], s[i+1:i+5]))
                self.color_model(name, e_type, s)

            except pymol.parsing.QuietException:
                print 'Coloring failed... did you forget to send geometry first?'


        elif ptype == 'Ene.gzip':   # energies info, gzip compression
            #print 'Getting Ene.gzip packet "%s"...' % name
            e_type_len = ord(data[0])
            e_type = data[1: 1 + e_type_len].tostring()
            s = gzip.GzipFile('', 'r', 0, StringIO(data[1+e_type_len:])).read()
            #print 'etype=%s  msg=%s' % (e_type, s)
            #print 'Compression stats: %s-->%s' % (len(data[1+e_type_len:]), len(s) )
            try:
                #for i in range(0, len(s), 7):
                #    pymol.cmd.color('R%s' % s[i+5:i+7], '%s and chain %s and resi %s' % (name, s[i], s[i+1:i+5]))
                self.color_model(name, e_type, s)

            except pymol.parsing.QuietException:
                print 'Coloring failed... did you forget to send geometry first?'

#############################################################################################
#################################

        # color by polar resi bool
        elif ptype == 'pol.bz2 ':
            data = bz2.decompress(data)
            size = int(data[0])
            data = data[1:]
            try:
                for i in range(0,len(data),size+6):
                    dat = data[i+6:i+6+size]
                    sel = '%s and chain %s and resi %s' % (name,data[i+5],data[i:i+5].strip())
                    color = 'blue'
                    if int(dat):
                        color = 'red'
                    pymol.cmd.color(color , sel)
            except pymol.parsing.QuietException:
                print 'Coloring failed... did you forget to send the geometry?'


        # color by movemap dof
        elif ptype == 'mm1.bz2 ':
            data = bz2.decompress(data)
            size = int(data[0])
            data = data[1:]
            try:
                cmd.remove('hydro')
                bb = '(name N or name CA or name C or name O)'
                for i in range(0,len(data),size+6):
                    dat = data[i+6:i+6+size]
                    sel = '%s and chain %s and resi %s' % (name,data[i+5],data[i:i+5].strip())
                    # bb
                    color = 'blue'
                    if int(dat[0])-1:
                        color = 'red'
                    pymol.cmd.color(color , sel+' and '+bb)
                    # sc
                    color = 'blue'
                    if int(dat[1])-1:
                        color = 'red'
                    pymol.cmd.color(color , sel+' and not '+bb)
            except pymol.parsing.QuietException:
                print 'Coloring failed... did you forget to send the geometry?'

        # color by foldtree edges
        elif ptype == 'ft1.bz2 ':
            data = bz2.decompress(data)
            size = int(data[0])
            data = data[1:]
            try:
                for i in range(0,len(data),size+6):
                    dat = data[i+6:i+6+size]
                    sel = '%s and chain %s and resi %s' % (name,data[i+5],data[i:i+5].strip())
                    # color based on int sent...
                    pymol.cmd.color(str(int(dat)*100+2000),sel)
            except pymol.parsing.QuietException:
                print 'Commands failed... did you forget to send the geometry?'


        # foldtree diagram
        elif ptype == 'ftd.bz2 ':
            data = bz2.decompress(data)
            # the number of jumps
            njump = int(data[:2])
            data = data[2:]
            # delete previous foldtrees and jumps, later make simultaneous viewing
            pymol.cmd.delete('jumps_'+name+' or foldtree_'+name)
            # use for size
            scale = .1    # arbitrary, get into viewing window easily
            jumps = []    # jumps assigned to index jump-1 with start(0) and stop(1)
            for j in range( 0 , njump ):
                # store jump start,stop,cut
                jumps.append([int(data[:3].strip()),int(data[3:6].strip()),int(data[6:9].strip())])
                data = data[9:]
            # add edges, color cutpoints red
            for i in range( 0 , len(data) ):
                if not int(data[i]):
                    connect = 'red'
                else:
                    # color by chain
                    connect = str(25+int(data[i]))    # magic number for coloring
                add_point( 'foldtree_'+name , [i*scale,0,0] , connect )
            # hacky string for selecting all jumps
            jumpsall = 'None'    # jumpsall will be used to select all the jumps
            r = 3    # arbitrary, distance jumps away, may want inside loop, different r per loop
            for j in range( 1 , njump + 1 ):
                # change color later
                connect = str(7+j)    # magic number for coloring
                # used frequently below
                th = 2*math.pi/njump*(j-1)
                cylx = r*math.cos(th)*scale
                cyly = r*math.sin(th)*scale
                start = (jumps[j-1][0]-1)*scale
                stop = (jumps[j-1][1]-1)*scale
                cut = (jumps[j-1][2]-1)*scale

                jump = 'jump_'+str(j)+'_'+name
                jumpsall = jumpsall + ' or ' +jump
                add_point( jump , [start,0,0] , connect )    # start bridge
                add_point( jump , [start,cylx,cyly] , connect ,False,False,'',0,str(jumps[j-1][0]))    # up one r
                add_point( jump , [(start+stop)/2,cylx,cyly] , connect ,False,False,'',0,'j'+str(j))    # name the jump
                add_point( jump , [stop,cylx,cyly] , connect ,False,False,'',0,str(jumps[j-1][1]))    # up one r
                add_point( jump , [stop,0,0] , connect )    # stop bridge
                add_point( jump , [cut,cylx,cyly] , '',False,False,'',0,str(jumps[j-1][2]))
                add_point( jump , [cut,0,0] , 'red' )
            # group all jumps, view centered on them
            if njump:
                pymol.cmd.group('jumps_'+name,jumpsall)
                pymol.cmd.label('jumps_'+name,'resn')
                pymol.cmd.center('jumps_'+name)


        # display hydrogen bonds
        elif ptype == 'hbd.bz2 ':
            data = bz2.decompress(data)
            # first 5 char is # of hbonds
            nhbonds = data[:5]
            data = data[5:] # 22 char per hbond 6+4 + 6+4 + 2
            try:
                pymol.cmd.delete(name+'_hbonds')
                for i in range(0,int(nhbonds)):
                    c = 22*i
                    # acceptor atom
                    acc_res = data[c:c+5].strip()
                    acc_chain = data[c+5]
                    acc_name = data[c+6:c+10].strip()
                    # donor atom
                    don_res = data[c+10:c+15].strip()
                    don_chain = data[c+15]
                    don_name = data[c+16:c+20].strip()
                    # make selection
                    hbname = 'hb_'+ acc_res+acc_chain+acc_name + '_' + don_res+don_chain+don_name + '_' + name
                    pymol.cmd.distance( hbname , name+' and chain '+acc_chain+' and res '+acc_res+' and name '+acc_name , name+' and chain '+don_chain+' and res '+don_res+' and name '+don_name )
                    pymol.cmd.color('Hb%s' % data[c+20:c+22], hbname)
                pymol.cmd.hide('labels','hb_*_'+name)
                pymol.cmd.group(name+'_hbonds','hb_*_'+name)
            except:
                print 'Commands failed... did you forget to send the geometry?'


        # generate a graph from data sent
        elif ptype == 'grp1.bz2':
            data = bz2.decompress(data)
            # first 21 data char are options

            options =data[:21]
            data = data[21:]
            connect = options[:7].strip()
            scale = bool(int(options[7:8]))
            axis_color = options[8:15].strip()
            num = int(options[15:])

            x_array = [0]*(len(data)/27)
            y_array = [0]*(len(data)/27)
            z_array = [0]*(len(data)/27)
            for i in range(0,len(data),27):
                x_array[i/27] = float(data[i:i+9].strip())
                y_array[i/27] = float(data[i+9:i+18].strip())
                z_array[i/27] = float(data[i+18:i+27].strip())
            plot3d(name,x_array,y_array,z_array,connect,scale,axis_color,num)



        # add a point to graph data
        elif ptype == 'pnt.bz2 ':
            data = bz2.decompress(data)
            # first 27 are data...last 22+ are options and banner
            add_point(name,[float(data[0:9].strip()),float(data[9:18].strip()),float(data[18:27].strip())],data[27:34].strip(),bool(int(data[34:35])),bool(int(data[35:36])),data[36:43].strip(),int(data[43:49]),data[49:])


        # label energy per residue
        elif ptype == 'lbE1.bz2':
            data = bz2.decompress(data)
            size = int(data[0])    # this int must be single digit, the size of data pieces
            data = data[1:]
            try:
                for i in range(0,len(data),size+6):
                    dat = data[i+6:i+6+size]
                    sel = '%s and name ca and chain %s and resi %s' % (name,data[i+5],data[i:i+5].strip())    # keep it for protein now
                    # dat dependent commands on sel
                    pymol.cmd.label(sel,dat)
            except pymol.parsing.QuietException:
                print 'Commands failed... did you forget to send the geometry?'




        # template
        elif ptype == 'temp.bz2':    # ptype must be the same on both sides
            data = bz2.decompress(data)
            size = int(data[0])    # this int must be single digit, the size of data pieces
            data = data[1:]
            try:
                for i in range(0,len(data),size+6):
                    dat = data[i+6:i+6+size]
                    sel = '%s and chain %s and resi %s' % (name,data[i+5],data[i:i+5].strip())
                    # dat dependent commands on sel
            except pymol.parsing.QuietException:
                print 'Commands failed... did you forget to send the geometry?'



#################################
#############################################################################################

        else:
            print 'Unknow packet type: %s, - ignoring...' % ptype



#############################################################################################
#################################
# graphs


def make_axis( name , ends , dr , scale = False , axis_color = '' , num = 0 ):
# make an axis 'name' from 'ends[0]' to 'ends[1]' on the 'dr' (direction)
#    'scale' bool tells to label axis or not, 'axis_color' determines the axis color
#    'num' allows intermediate scale points

    # quit if bad data was fed
    if ends[0]==ends[1]:
        return

    # ends must be [max,min]
    d_max = ends[0]
    d_min = ends[1]

    # create or alter axis object
    if name in pymol.cmd.get_names():
        pymol.cmd.remove(name)
    else:
        pymol.cmd.create(name,None)

    # add max and min points, connect them
    pymol.cmd.pseudoatom(name,
        name='max',
        resn=str(d_max)[:5],
        pos=[dr[0]*d_max,dr[1]*d_max,dr[2]*d_max],
        color=axis_color)
    pymol.cmd.pseudoatom(name,
        name='min',
        resn=str(d_min)[:5],
        pos=[dr[0]*d_min,dr[1]*d_min,dr[2]*d_min],
        color=axis_color)
    pymol.cmd.bond(name+' and name max',name+' and name min')


    # add num intermediate atoms
    for i in range( 1 , num + 1 ):
        val = (d_max - d_min)*float(i)/(num+1) + d_min
        pymol.cmd.pseudoatom(name,
            name=str(i),
            resn=str(val)[:5],
            pos=[dr[0]*val,dr[1]*val,dr[2]*val],
            color=axis_color)

    if scale:
        scale_axes(name)


# label the selection name
#    used here for convenience, axis point value str in resn
def scale_axes( name = 'x_axis or y_axis or z_axis' ):
    pymol.cmd.hide('label',name)
    pymol.cmd.label(name,'resn')
pymol.cmd.extend('scale_axes',scale_axes)


# determine the ends for axes from list
#    defaults to min, max but if not pos/neg
#    determines based on closest point
#    this is silly...
def get_ends( data ):
    d_max = max(data)
    if d_max > 0:
        d_max = d_max*1.1
    else:
        d_max = abs(d_max*.1)
    d_min = min(data)
    if d_min < 0:
        d_min = d_min*1.1
    else:
        d_min = -abs(d_min*.1)
    return [d_max,d_min]


# returns the coords of a data container as three arrays
def extract_coords( name ):
    data = pymol.cmd.get_model(name)
    size = pymol.cmd.count_atoms(name)
    x_array = [0]*size
    y_array = [0]*size
    z_array = [0]*size
    for atom in data.atom:
        i = atom.index-1
        x_array[i] = atom.coord[0]
        y_array[i] = atom.coord[1]
        z_array[i] = atom.coord[2]
    return x_array,y_array,z_array
pymol.cmd.extend('extract_coords',extract_coords)


# generate a data object from three arrays
def make_data( name , x_array , y_array , z_array , data_color = 'red' ):
    # erase old name, create new
    pymol.cmd.delete(name)
    pymol.cmd.create(name,None)

    # produce points
    for i in range( 0 , len(x_array) ):
        pymol.cmd.pseudoatom(name,
            name=str(i+1),
            pos=[x_array[i],y_array[i],z_array[i]],
            color=data_color)
pymol.cmd.extend('make_data',make_data)


# adds a point to existing data, reconnects points optionally
def add_point( name , point , connect = '' , rescale = False , scale = False , axis_color = '' , num = 0 , banner = '' ):
    # adds 'point' to 'name' and connects it with color 'connect' (empty for no connect)
    #    'rescale' will rescale the axes, 'scale' will label the scales,
    #    'axis_color' determines the axis color, 'num' intermediate scales,
    #    'banner' is a resn name at 'point'
    # if 'name' doesn't exist, create it!
    created = False
    if not name in pymol.cmd.get_names():
        pymol.cmd.create(name,None)
        created = True
    ind = pymol.cmd.count_atoms(name) + 1
    pymol.cmd.pseudoatom(name,
        name=str(ind),
        resn=banner,
        pos=point,
        color=connect)

    if connect:
        pymol.cmd.bond(name+' and name '+str(ind),name+' and name '+str(ind-1))

    if rescale:
        rescale_cartesian(name,scale,axis_color,num)

    if created and rescale:
        pymol.cmd.group(name+'_data',name+' or axes_'+name)
pymol.cmd.extend('add_point',add_point)


# creates axes (x,y,z) based on data
def rescale_cartesian( data , scale , axis_color , num ):
    # rescales the 'data' axes based on 'data' data,
    #    'num','axis_color', 'scale'
    x_array,y_array,z_array = extract_coords(data)
    make_axis('x_axis_'+str(data),get_ends(x_array),[1,0,0],scale,axis_color,num)
    make_axis('y_axis_'+str(data),get_ends(y_array),[0,1,0],scale,axis_color,num)
    make_axis('z_axis_'+str(data),get_ends(z_array),[0,0,1],scale,axis_color,num)
    pymol.cmd.group('axes_'+str(data),'x_axis_'+str(data)+' or y_axis_'+str(data)+' or z_axis_'+str(data))
    pymol.cmd.reset()


# draws colored lines between consecutive points
def connect_points( name , color = 'red' ):
    for i in range( 2 , pymol.cmd.count_atoms(name) + 1 ):
        pymol.cmd.bond(name+' and name '+str(i),name+' and name '+str(i-1))
    pymol.cmd.color(color,name)


# default plotting tool, optional z points, connection, axis color
def plot3d( name , x_array , y_array , z_array = [] , connect = 'red' , scale = True , axis_color = 'blue' , num = 0 ):
    if not z_array:
        z_array = [0]*len(x_array)
    if not len(x_array) == len(y_array) == len(z_array):
        raise IOError('FAIL: arrays not the same length')

    # make the data
    make_data(name,x_array,y_array,z_array,connect)

    # make axes
    rescale_cartesian(name,scale,axis_color,num)

    pymol.cmd.group(name+'_data',name+' or axes_'+name)

    if connect:
        connect_points(name,connect)
pymol.cmd.extend('plot3d',plot3d)


# end of graphs
#################################
#############################################################################################


color_lib={
'white':[1,1,1],
'yellow':[1,1,0],
'magenta':[1,0,1],
'cyan':[0,1,1],
'red':[1,0,0],
'blue':[0,0,1],
'green':[0,1,0],
'black':[0,0,0]}


X11Colors = {
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

# Creating our own color spectrum
for i in range(256): pymol.cmd.set_color('R%02x' % i, [i/255., 0, 1-i/255.])
for j in range(256): pymol.cmd.set_color('Hb%02x' % j, [j/255., 1, 0])

for i in X11Colors: pymol.cmd.set_color('X%02x' % X11Colors[i][0], X11Colors[i][1:])

# sets the energy coloring spectrumusing the above color_lib
def set_spectrum( low = 'blue' , high = 'red' ):
    print 'new spectrum:',low,'==>',high
    low = color_lib[low]
    high = color_lib[high]
    cust = [0]*3
    for j in range(0,3):
        cust[j]=high[j]-low[j]
    x=lambda i,el:cust[el]*( i - 255*(cust[el]<0) )/255.
    for i in range(256):
        pymol.cmd.set_color('R%02x' % i, [x(i,0), x(i,1), x(i,2)])

pymol.cmd.extend('set_spectrum',set_spectrum)



def main(ip, port):
    print 'PyMOL <---> PyRosetta link started!'
    print 'at',ip,'port',port

    udp_serv = PR_UDPServer(ip, port)
    PS = PR_PyMOLServer()
    while True:
        r = udp_serv.listen()
        if r:
            #print len(r)
            PS.processPacket(r)

    s.close()


def start_rosetta_server(ip='', port=65000):
   if not ip:
       ip = socket.gethostbyname(socket.gethostname())
       if ip == '127.0.0.1':
           print "Unable to automatically determine your IP address.  Please specify it manually. e.g. start_rosetta_server 192.168.0.1"
           return

   thread = threading.Thread(target=main, args=[ip,port])
   thread.setDaemon(1)
   thread.start()

pymol.cmd.extend('start_rosetta_server', start_rosetta_server)

start_rosetta_server('127.0.0.1', 65000)

#### To use PyMOLPyRosettaServer over a network, uncomment the line below and set the first argument to your IP address
#start_rosetta_server('192.168.0.1', 65000)
