#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import Tkinter, tkFileDialog, tkSimpleDialog


from pymol import cmd


def __init__( self ):
    self.menuBar.addmenuitem( 'Plugin', 'command', "RMSD_score_plotter",
                              label='RMSD Score Plotter',
                              command= lambda s=self: RMSD_score_plotter(s))


class RMSD_score_plotter:

    def __init__( self,  app ):


        self.structure_filenames = []
        self.loaded_structures = []
        self.score = []
        self.rmsd = []


        self.parent = app.root


        self.main = Tkinter.Toplevel( self.parent )


        self.main.title("RMSD Score Plotter")

        self.width = 600
        self.height = 600
        self.pad = 35
        self.major_mark_length = 20
        self.mark_color = 'black'

        self.canvas = Tkinter.Canvas( self.main, width=self.width, height=self.height )


        self.import_data()
        self.display_data()

        self.canvas.bind("<Double-Button-1>", self.event_double_click )
        self.canvas.bind("<Button-1>", self.event_one_click )


        self.canvas.pack(side="top")

        print "single click point to get name of structure"
        print "double click point to open structure in pymol"


    def import_data( self ):
        datafile = tkFileDialog.askopenfile()

        for line in datafile.read().split("\n")[1:]:

            split_line = line.split(" ")

            if len( split_line ) < 3: continue

            self.structure_filenames.append( split_line[0] )
            print "structure_filename:", split_line[0]

            self.score.append( float( split_line[1]) )
            print "score:", split_line[1]


            self.rmsd.append( float( split_line[2] ) )
            print "rmsd:", split_line[2]

            if len( split_line ) < 4: continue
            self.sasa_pack.append( float( split_line[3] ) )
            print "sasa pack:", split_line[3]


            print ""

        datafile.close()

    def draw_major_x_mark( self, point, text ):

        x, y = point
        self.canvas.create_line( x, y, x, y + self.major_mark_length, fill=self.mark_color )
        self.canvas.create_text( x + 3, y + 3,   text=text, fill=self.mark_color, anchor="nw" )

    def draw_major_y_mark( self, point, text ):
        x, y = point
        self.canvas.create_line( x - self.major_mark_length, y, x, y, fill=self.mark_color )
        self.canvas.create_text( x - 3, y + 3,   text=text, fill=self.mark_color, anchor="ne" )


    def draw_axes( self ):

        ul_corner = self.pad, self.pad
        ll_corner = self.pad, self.height-self.pad
        lr_corner = self.width-self.pad, self.height-self.pad
        lmid = float( self.width ) / 2, self.height-self.pad
        midl = self.pad, float( self.height ) / 2

        self.canvas.create_line( ll_corner[ 0 ], ll_corner[ 1 ], ul_corner[ 0 ], ul_corner[ 1 ], fill=self.mark_color )
        self.canvas.create_line( ll_corner[ 0 ], ll_corner[ 1 ], lr_corner[ 0 ], lr_corner[ 1 ], fill=self.mark_color )

        self.draw_major_x_mark( ll_corner, str( round( self.x_min, 1 ) ) )
        self.draw_major_x_mark( lmid,      str( round( float( self.x_min + self.x_max ) / 2, 1 ) ) )
        self.draw_major_x_mark( lr_corner, str( round( self.x_max, 1 ) ) )

        self.draw_major_y_mark( ll_corner, str( round( self.y_min, 1 ) ) )
        self.draw_major_y_mark( midl,      str( round( float( self.y_min + self.y_max ) / 2, 1 ) ) )
        self.draw_major_y_mark( ul_corner, str( round( self.y_max, 1 ) ) )


    def display_data( self ):
        self.x_min = min( self.rmsd )
        self.x_max = max( self.rmsd )

        x_scale = ( self.width - 2 * self.pad ) / ( self.x_max - self.x_min  )


        self.y_min = min( self.score )
        self.y_max = max( self.score )

        y_scale = ( self.height - 2 * self.pad ) / ( self.y_max - self.y_min )


        x_offset = self.pad - self.x_min*x_scale
        y_offset = self.pad - self.y_min*y_scale

        point_size = 8
        point_color = 'blue'


        print "x_min", self.x_min
        print "x_max", self.x_max
        print "x_scale", x_scale
        print "y_min", self.y_min
        print "y_max", self.y_max
        print "y_scale", y_scale
        print "x_offset", x_offset
        print "y_offset", y_offset


        self.draw_axes()


        for i in range( len( self.score ) ):


            x = self.rmsd[ i ]
            y = self.score[  i ]

            x = x_scale * x + x_offset
            y = self.height - (y_scale * y + y_offset)


            print "displaying point (rmsd, score)", self.rmsd[ i ], self.score[ i], "\tat \t",x, y


            x0 = x - point_size/2.
            x1 = x + point_size/2.

            y0 = y - point_size/2.
            y1 = y + point_size/2.

            self.canvas.create_oval(x0, y0, x1, y1, fill=point_color, tag=str( i ) )


    def event_double_click( self, event ):

        if event.widget == self.canvas:

            x = self.canvas.canvasx( event.x )
            y = self.canvas.canvasy( event.y )


            try:
                i =  int( self.canvas.gettags( self.canvas.find_closest( x, y ) )[0] )

                self.display_structure( self.structure_filenames[ i ] )
            except:
                pass

    def event_one_click( self, event ):

        if event.widget == self.canvas:

            x = self.canvas.canvasx( event.x )
            y = self.canvas.canvasy( event.y )

            try:
                i =  int( self.canvas.gettags( self.canvas.find_closest( x, y ) )[0] )

                structure_filename = self.structure_filenames[ i ]

                print "This is", structure_filename, self.rmsd[ i ], self.score[ i ]
            except:
                pass


    def display_structure( self, structure_filename ):
        if structure_filename in self.loaded_structures:
            print "Structure", structure_filename, "has already been loaded."
            return
        else:
            self.loaded_structures.append(structure_filename)

        print "opening structure file", structure_filename


        structure_file = open( structure_filename )


        pdb = structure_file.read()

        #if the structure_filename is
        #/path/to/pdb/af1prb0001.pdb
        #make the pdb_name: af1prb0001

        pdb_name = structure_filename.split("/")[-1]
        pdb_name = ".".join(pdb_name.split(".")[0:-1])

        print "pdb has ", len( pdb.split("\n") ), " lines."

        cmd.load_raw(content=pdb, format='pdb', object=pdb_name )

        cmd.hide("everything", pdb_name)
        cmd.show("cartoon", pdb_name)


