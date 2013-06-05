import blargs

import Tkinter, tkFileDialog, tkSimpleDialog
import numpy
import inspect, warnings
import math

from cvxopt import lapack, matrix

def fit_polynomial(X, Y, d):
    m = len(X)
    assert(len(Y) == m)

    A = matrix( [[X**k] for k in xrange(d)])

   # make a deep copy of Y
    xls = +Y

    # general least-squares: minimize ||A*x - y||_2
    lapack.gels(+A,xls)
    xls = xls[:d]

    return xls

def fit_polynomial_one_constraint(X, Y, d, x0, y0):
    m = len(X)
    assert(len(Y) == m)

    A = matrix( [[X**k] for k in xrange(d)])

    G = matrix(0.0, (2,d))
    G[0, range(d)] = matrix([[x0**k] for k in xrange(d)])
    G[1, range(1,d)] = matrix([[k*x0**(k-1)] for k in xrange(1,d)])

    # LS fit
    #
    #     minimize    (1/2) * || A*x - Y ||_2^2
    #     subject to  G*x = h
    #
    # Solve as a linear equation
    #
    #     [ A'*A  G' ] [ x ]   [ A'*Y ]
    #     [ G     0  ] [ Y ] = [ 0    ].

    K = matrix(0.0, (d+2,d+2))
    K[:d,:d] = A.T * A
    K[d:,:d] = G
    xls = matrix(0.0, (d+2,1))
    xls[:d] = A.T * Y
    xls[d] = y0
    lapack.sysv(K, xls)
    xls = xls[:d]
    return xls

# x0 is a list of xvals for zero-derivative constraints.
# The fit polynomial must have a derivative of zero
# at these points.  y0 is a list of yvals for x0
# x1 is a list of xvals for constraints where the
# derivative is unconstrained.  y1 is a list of yvals
# for x1.
# x0 must have the same length as y0.
# x1 must have the same length as y1
def fit_polynomial_multiple_constraints(X, Y, d, x0, y0, x1 = None, y1 = None):
    m = len(X)
    assert(len(Y) == m)

    A = matrix( [[X**k] for k in xrange(d)])

    ncst0 = len(x0)
    assert( ncst0 == len(y0))

    ncst1 = 0
    if x1 : ncst1 = len(x1)

    G = matrix(0.0, (2*ncst0+ncst1,d))
    for ii in xrange(ncst0) :
        G[2*ii+0, range(d)]   = matrix([[x0[ii]**k] for k in xrange(d)])
        G[2*ii+1, range(1,d)] = matrix([[k*x0[ii]**(k-1)] for k in xrange(1,d)])
    if x1 :
        for ii in xrange(ncst1) :
            G[ii+2*ncst0,range(d)] = matrix( [[x1[ii]**k] for k in xrange(d)] )

    # LS fit
    #
    #     minimize    (1/2) * || A*x - Y ||_2^2
    #     subject to  G*x = h
    #
    # Solve as a linear equation
    #
    #     [ A'*A  G' ] [ x ]   [ A'*Y ]
    #     [ G     0  ] [ Y ] = [ 0    ].

    K = matrix(0.0, (d+2*ncst0+ncst1,d+2*ncst0+ncst1))
    K[:d,:d] = A.T * A
    K[d:,:d] = G
    xls = matrix(0.0, (d+2*ncst0+ncst1,1))
    xls[:d] = A.T * Y
    for ii in xrange( ncst0 ) :
        xls[d+ii*2] = y0[ii]
    if x1 :
        for ii in xrange( ncst1 ) :
            xls[d+2*ncst0+ii] = y1[ii]
    lapack.sysv(K, xls)
    xls = xls[:d]
    return xls

    #coefs = numpy.array(xls.T)[0]
    #coefs = coefs[::-1]
    #return coefs


def old_fit_polynomial_multiple_constraints(X, Y, d, x0, y0):
    m = len(X)
    assert(len(Y) == m)

    A = matrix( [[X**k] for k in xrange(d)])

    ncst = len(x0)
    assert( ncst == len(y0))

    G = matrix(0.0, (2*ncst,d))
    for ii in xrange(ncst) :
        G[2*ii+0, range(d)]   = matrix([[x0[ii]**k] for k in xrange(d)])
        G[2*ii+1, range(1,d)] = matrix([[k*x0[ii]**(k-1)] for k in xrange(1,d)])

    # LS fit
    #
    #     minimize    (1/2) * || A*x - Y ||_2^2
    #     subject to  G*x = h
    #
    # Solve as a linear equation
    #
    #     [ A'*A  G' ] [ x ]   [ A'*Y ]
    #     [ G     0  ] [ Y ] = [ 0    ].

    K = matrix(0.0, (d+2*ncst,d+2*ncst))
    K[:d,:d] = A.T * A
    K[d:,:d] = G
    xls = matrix(0.0, (d+2*ncst,1))
    xls[:d] = A.T * Y
    for ii in xrange( ncst ) :
        xls[d+ii*2] = y0[ii]
    lapack.sysv(K, xls)
    xls = xls[:d]
    return xls

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno

class PolyFit:

    def __init__(self, initial_control_points):

        self.main = Tkinter.Tk()
        self.main.title("Polynomial Fitting")

        self.verbose = False
        self.canvas_width =  750
        self.canvas_height = 750
        self.plot_width = 680
        self.plot_height = 680
        self.left_gutter = 10
        self.top_gutter = 10
        self.fit_normal = True
        self.fit_deg_as_cos = False
        self.fit_deg_as_rad = False

        self.xaxis_tick_height = 5
        self.yaxis_tick_width = 5

        self.axis_color = "white"


        self.polynomial_dimensions = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
        self.nbins = 200

        self.coefficient_precision = 8

        self.dot_color = 'blue'
        self.dot_width = 6

        self.polynomial_colors = {2:"#D53E4F", 3:"#FC8D59", 4:"#FEE08B", 5:"#FFFFBF", 6:"#E6F598", 7:"#99D594", 8:"#3288BD", 9:"#FFFFFF", 10:"#000000",11:"#D53E4F", 12:"#FC8D59", 13:"#FEE08B", 14:"#FFFFBF", 15:"#E6F598", 16:"#99D594"}

        self.panel_frame = Tkinter.Frame(self.main)

        self.build_poly_frame(self.panel_frame)
        # inserting them in reverse order because it works...
        self.build_zoom_frame(self.panel_frame)
        self.build_constraints_frame(self.panel_frame)
        self.build_minima_frame(self.panel_frame)
#        self.build_reflect_frame(self.panel_frame)
        self.build_reference_polys_frame(self.panel_frame)


        self.panel_frame.pack(side="right")

        self.canvas = Tkinter.Canvas(\
            self.main,height=self.canvas_height, width=self.canvas_width)
        self.canvas.pack(side="top")

        self.btn_makeart = Tkinter.Button(\
            self.main, text="clear", command=self.clear_callback )
        self.btn_makeart.pack(side="top")

        self.hide_polys_checkbutton = Tkinter.Checkbutton(
            self.main, text="hide polys", command=self.toggle_hide_polys_callback )
        self.hide_polys_checkbutton.pack( side="top" )

        self.main.bind("<Button-1>", self.event_lclick )
        self.main.bind("<Button-3>", self.event_rclick )


        self.minima = [None, None] # a list of length two of x and y minima
        self.constraints = [None, None] # a list of length two of x and y values that the polynomials must pass through

        self.control_points = initial_control_points
        self.control_point_widgets = []

        self.polynomial_coefficients = {}
        self.polynomial_points = {}
        self.hide_polynomials = False

        self.ref_polynomial_coefficients= []
        self.ref_polynomial_points = {}

        self.xmin, self.xmax, self.ymin, self.ymax = -1, 1, -1, 1
        self.zoom_xmin_var.set(-1)
        self.zoom_xmax_var.set(1)
        self.zoom_ymin_var.set(-1)
        self.zoom_ymax_var.set(1)

        self.redraw_everything()


#        self.test_diagnostics()

    def go( self ) :
        self.main.mainloop()

    def build_poly_frame(self, parent):
        self.poly_text_frame = Tkinter.Frame(parent)
        self.poly_label = Tkinter.Label(\
            self.poly_text_frame, text="Polynomial Coefficients: largest coefficient first")
        self.poly_label.pack(side="top")
        self.poly_text = Tkinter.Text(\
            self.poly_text_frame, width = 50, background="lightgray")
        self.poly_text.pack(side="bottom")
        self.poly_text_frame.pack(side="top")

    def build_reflect_frame(self, parent):
        self.reflect_frame = Tkinter.Frame(parent)
        self.reflect_check_label= Tkinter.Label(self.reflect_frame, text="reflect y-axis:")
        self.reflect_check_label.pack(side="left")
        self.reflect_check_var = Tkinter.IntVar()
        self.reflect_check = Tkinter.Checkbutton(self.reflect_frame, variable=self.reflect_check_var, command=self.redraw_everything)
        self.reflect_check.pack(side="left")
        self.reflect_var = Tkinter.StringVar()
        self.reflect_var.trace("w", self.reflect_var_callback)
        self.reflect_value_label = Tkinter.Label(self.reflect_frame, text="y=")
        self.reflect_value_label.pack(side="left")
        self.reflect_value = Tkinter.Entry(self.reflect_frame, textvariable=self.reflect_var)
        self.reflect_value.pack(side="right", expand=1, fill=Tkinter.X)
        self.reflect_frame.pack(side="bottom", expand=1, fill=Tkinter.X)

    def build_minima_frame(self, parent):
        self.minima_frame = Tkinter.Frame(parent)
        self.minima_x_label= Tkinter.Label(self.minima_frame, text="Minima x:")
        self.minima_x_label.pack(side="left")
        self.minima_x_var = Tkinter.StringVar()
        self.minima_x_var.trace("w", self.minima_x_var_callback)
        self.minima_x_value = Tkinter.Entry(self.minima_frame, textvariable=self.minima_x_var)
        self.minima_x_value.pack(side="left", expand=1, fill=Tkinter.X)

        self.minima_y_label = Tkinter.Label(self.minima_frame, text="y:")
        self.minima_y_label.pack(side="left")
        self.minima_y_var = Tkinter.StringVar()
        self.minima_y_var.trace("w", self.minima_y_var_callback)
        self.minima_y_value = Tkinter.Entry(self.minima_frame, textvariable=self.minima_y_var)
        self.minima_y_value.pack(side="left", expand=1, fill=Tkinter.X)

        self.minima_frame.pack(side="bottom", expand=1, fill=Tkinter.X)

    def build_constraints_frame( self, parent ) :
        self.constraints_frame = Tkinter.Frame(parent)
        self.constraints_x_label= Tkinter.Label(self.constraints_frame, text="Constraints x:")
        self.constraints_x_label.pack(side="left")
        self.constraints_x_var = Tkinter.StringVar()
        self.constraints_x_var.trace("w", self.constraints_x_var_callback)
        self.constraints_x_value = Tkinter.Entry(self.constraints_frame, textvariable=self.constraints_x_var)
        self.constraints_x_value.pack(side="left", expand=1, fill=Tkinter.X)

        self.constraints_y_label = Tkinter.Label(self.constraints_frame, text="y:")
        self.constraints_y_label.pack(side="left")
        self.constraints_y_var = Tkinter.StringVar()
        self.constraints_y_var.trace("w", self.constraints_y_var_callback)
        self.constraints_y_value = Tkinter.Entry(self.constraints_frame, textvariable=self.constraints_y_var)
        self.constraints_y_value.pack(side="left", expand=1, fill=Tkinter.X)

        self.constraints_frame.pack(side="bottom", expand=1, fill=Tkinter.X)

    def build_reference_polys_frame(self, parent):
        self.reference_polys_frame = Tkinter.Frame(parent)
        self.reference_polys_label = Tkinter.Label(self.reference_polys_frame, text="Add reference polynomials here:")
        self.reference_polys_label.pack(side="top")
        self.reference_polys_text = Tkinter.Text(
            self.reference_polys_frame, width=50, background="lightgray")
#        self.reference_polys_text.tag_add("update", 1.0, Tkinter.END)
        self.reference_polys_text.bind("<Key>", self.reference_polys_text_callback)
        self.reference_polys_text.pack(side="bottom")
        self.reference_polys_frame.pack(side="bottom")

    def build_zoom_frame(self, parent, width=6):

        self.zoom_frame = Tkinter.Frame(parent)
        self.zoom_xmin_label = Tkinter.Label(self.zoom_frame, text="xmin:")
        self.zoom_xmin_label.pack(side="left")
        self.zoom_xmin_var = Tkinter.StringVar()
        self.zoom_xmin_var.trace("w", self.zoom_xmin_var_callback)
        self.zoom_xmin_value = Tkinter.Entry(self.zoom_frame, width=width, textvariable=self.zoom_xmin_var)
        self.zoom_xmin_value.pack(side="left", expand=1, fill=Tkinter.X)

        self.zoom_xmax_label = Tkinter.Label(self.zoom_frame, text="xmax:")
        self.zoom_xmax_label.pack(side="left")
        self.zoom_xmax_var = Tkinter.StringVar()
        self.zoom_xmax_var.trace("w", self.zoom_xmax_var_callback)
        self.zoom_xmax_value = Tkinter.Entry(self.zoom_frame, width=width, textvariable=self.zoom_xmax_var)
        self.zoom_xmax_value.pack(side="left", expand=1, fill=Tkinter.X)

        self.zoom_ymin_label = Tkinter.Label(self.zoom_frame, text="ymin:")
        self.zoom_ymin_label.pack(side="left")
        self.zoom_ymin_var = Tkinter.StringVar()
        self.zoom_ymin_var.trace("w", self.zoom_ymin_var_callback)
        self.zoom_ymin_value = Tkinter.Entry(self.zoom_frame, width=width, textvariable=self.zoom_ymin_var)
        self.zoom_ymin_value.pack(side="left", expand=1, fill=Tkinter.X)

        self.zoom_ymax_label = Tkinter.Label(self.zoom_frame, text="ymax:")
        self.zoom_ymax_label.pack(side="left")
        self.zoom_ymax_var = Tkinter.StringVar()
        self.zoom_ymax_var.trace("w", self.zoom_ymax_var_callback)
        self.zoom_ymax_value = Tkinter.Entry(self.zoom_frame, width=width, textvariable=self.zoom_ymax_var)
        self.zoom_ymax_value.pack(side="left", expand=1, fill=Tkinter.X)
        self.zoom_frame.pack(side="bottom", expand=1, fill=Tkinter.X)



    def reference_polys_text_callback(self, event):
        self.redraw_everything()

    def reflect_var_callback(self, name, index, mode):
        if self.verbose: print lineno(), "setting reflection to y=%s" % self.reflect_var.get()
        self.redraw_everything()

    def minima_x_var_callback(self, name, index, mode):
        minima_string = self.minima_x_var.get();
        if self.verbose : print "minima_x_var_callback",minima_string
        minima_strings = minima_string.split(",")
        if self.verbose : print minima_strings
        if len(minima_strings) == 0 :
            self.minima[ 0 ] = None
        else :
            try:
                if self.verbose: print "trying to read minima x values"
                self.minima[0] = []
                for ii in xrange( len( minima_strings ) ) :
                    if self.verbose: print "trying to extract minima", ii, minima_strings[ ii ]
                    iival = float( minima_strings[ ii ] )
                    if self.verbose: print iival
                    self.minima[0].append( iival )
                    if self.verbose: print lineno(), "setting minima x value #%d to x=%f" % (ii+1, self.minima[0][ii])
            except:
                self.minima[0] = None

        self.redraw_everything()

    def minima_y_var_callback(self, name, index, mode):
        minima_string = self.minima_y_var.get();
        if self.verbose: print "minima_y_var_callback",minima_string
        minima_strings = minima_string.split(",")
        if len(minima_strings) == 0 :
            self.minima[ 1 ] = None
        else :
            try:
                self.minima[1] = []
                for ii in xrange( len( minima_strings ) ) :
                    self.minima[1].append( float( minima_strings[ ii ] ))
                    if self.verbose: print lineno(), "setting minima y value #%d to x=%f" % (ii+1, self.minima[1][ii])
            except:
                self.minima[1] = None

        self.redraw_everything()
    def constraints_x_var_callback(self, name, index, mode):
        constraints_string = self.constraints_x_var.get();
        if self.verbose: print "constraints_x_var_callback",constraints_string
        constraints_strings = constraints_string.split(",")
        if self.verbose: print constraints_strings
        if len(constraints_strings) == 0 :
            self.constraints[ 0 ] = None
        else :
            try:
                if self.verbose: print "trying to read constraints x values"
                self.constraints[0] = []
                for ii in xrange( len( constraints_strings ) ) :
                    if self.verbose: print "trying to extract constraints", ii, constraints_strings[ ii ]
                    iival = float( constraints_strings[ ii ] )
                    if self.verbose: print iival
                    self.constraints[0].append( iival )
                    if self.verbose: print lineno(), "setting constraints x value #%d to x=%f" % (ii+1, self.constraints[0][ii])
            except:
                self.constraints[0] = None

        self.redraw_everything()

    def constraints_y_var_callback(self, name, index, mode):
        constraints_string = self.constraints_y_var.get();
        if self.verbose: print "constraints_y_var_callback",constraints_string
        constraints_strings = constraints_string.split(",")
        if len(constraints_strings) == 0 :
            self.constraints[ 1 ] = None
        else :
            try:
                self.constraints[1] = []
                for ii in xrange( len( constraints_strings ) ) :
                    self.constraints[1].append( float( constraints_strings[ ii ] ))
                    if self.verbose: print lineno(), "setting constraints y value #%d to x=%f" % (ii+1, self.constraints[1][ii])
            except:
                self.constraints[1] = None

        self.redraw_everything()

    def zoom_xmin_var_callback(self, name, index, mode):
        try:
            new_xmin = float(self.zoom_xmin_var.get())
        except:
            return

        if new_xmin > self.xmax:
            if self.verbose: print lineno(), "Attempting to set xmin=%s to a value greater than xmax=%s." % (new_xmin, self.xmax)
        else:
            self.xmin = new_xmin
            self.redraw_everything()

    def zoom_xmax_var_callback(self, name, index, mode):
        try:
            new_xmax = float(self.zoom_xmax_var.get())
        except:
            return

        if self.xmin > new_xmax:
            if self.verbose: print lineno(), "Attempting to set xmax=%s to a value less than xmin=%s." % (new_xmax, self.xmin)
        else:
            self.xmax = new_xmax
            self.redraw_everything()

    def zoom_ymin_var_callback(self, name, index, mode):
        try:
            new_ymin = float(self.zoom_ymin_var.get())
        except:
            return

        if new_ymin > self.ymax:
            if self.verbose: print lineno(), "Attempting to set ymin=%s to a value greater than ymax=%s." % (new_ymin, self.ymax)
        else:
            self.ymin = new_ymin
            self.redraw_everything()

    def zoom_ymax_var_callback(self, name, index, mode):
        try:
            new_ymax = float(self.zoom_ymax_var.get())
        except:
            return

        if self.ymin > new_ymax:
            if self.verbose: print lineno(), "Attempting to set ymax=%s to a value less than ymin=%s." % (new_ymax, self.ymin)
        else:
            self.ymax = new_ymax
            self.redraw_everything()

    def clear_callback(self):
        self.control_points = []
        self.redraw_everything()

    def toggle_hide_polys_callback(self) :
        self.hide_polynomials = False if self.hide_polynomials else True
        self.redraw_everything()

    def test_diagnostics(self):
        self.control_points = [(0,0), (1,1), (2,0)]
        self.compute_polynomial_fits()
        self.redraw_everything()


    def p2canv( self, x, y):
        px = (float(x)-self.xmin)*self.plot_width/(self.xmax-self.xmin) + \
            self.left_gutter
        py = -1*(float(y)-self.ymin)*self.plot_height/(self.ymax-self.ymin) + \
            self.top_gutter + self.plot_width

        return (px, py)

    def canv2p(self, px, py):
        x = (float(px) - self.left_gutter) * \
            (self.xmax - self.xmin)/self.plot_width + self.xmin
        y = -1*(float(py) - self.top_gutter - self.plot_width) * \
            (self.ymax - self.ymin)/self.plot_height + self.ymin

        return (x, y)

    def compute_zoom_and_axes(self):
        self.xaxis_ticks = [i * (self.xmax - self.xmin)/10 + self.xmin for i in range(11)]
        self.yaxis_ticks = [i * (self.ymax - self.ymin)/10 + self.ymin for i in range(11)]
        self.xaxis = (self.ymax-self.ymin)/2 + self.ymin
        self.yaxis = (self.xmax-self.xmin)/2 + self.xmin

    def clear( self ):
        for t in self.canvas.find_all():
            self.canvas.delete(t)

        #xaxis
        px0, py0 = self.p2canv(self.xmin, self.xaxis)
        px1, py1 = self.p2canv(self.xmax, self.xaxis)
        self.canvas.create_line( px0, py0, px1, py1, fill=self.axis_color)

        #yaxis
        px0, py0 = self.p2canv(self.yaxis, self.ymin)
        px1, py1 = self.p2canv(self.yaxis, self.ymax)
        self.canvas.create_line( px0, py0, px1, py1, fill=self.axis_color)

        for x in self.xaxis_ticks:
            px, py = self.p2canv(x,(self.ymax+self.ymin)/2)
            px0,px1 = px, px
            py0,py1 = py, py + self.xaxis_tick_height
            self.canvas.create_line( px0, py0, px1, py1, fill=self.axis_color)
            self.canvas.create_text( px1, py1+5, text=str(x), font=("helvetica", 7))

        for y in self.yaxis_ticks:
            px, py = self.p2canv((self.xmax+self.xmin)/2,y)
            px0,px1 = px, px + self.yaxis_tick_width
            py0,py1 = py, py
            self.canvas.create_line( px0, py0, px1, py1, fill=self.axis_color)
            self.canvas.create_text( px1+5, py1, text=str(y), font=("helvetica", 7))


    def redraw_everything( self ):

        self.compute_zoom_and_axes()
        self.clear()
        self.parse_reference_polynomials()
        self.compute_polynomial_fits()
        self.compute_reference_polynomial_points()
        self.draw_demo()
        self.print_coefficients()
        if self.verbose: print "Control point coordinates:"
        if self.verbose: print self.control_points
        if self.verbose: print ""

    def add_reference_polynomial( self, pstring ) :
        self.reference_polys_text.insert(Tkinter.END, pstring + "\n\n")

    def set_minima_x0s( self, x0 ):
        self.minima_x_var.set( x0 )
    def set_minima_y0s( self, y0 ):
        self.minima_y_var.set( y0 )

    def set_constraints_x1s( self, x1 ) :
        self.constraints_x_var.set( x1 )
    def set_constraints_y1s( self, y1 ) :
        self.constraints_y_var.set( y1 )

    def set_range_x( self, range_string ) :
        cols = range_string.split(",")
        self.zoom_xmin_var.set( cols[0] )
        self.zoom_xmax_var.set( cols[1] )

    def set_range_y( self, range_string ) :
        cols = range_string.split(",")
        self.zoom_ymin_var.set( cols[0] )
        self.zoom_ymax_var.set( cols[1] )

    def parse_reference_polynomials(self):
        contents = self.reference_polys_text.get(1.0, Tkinter.END)
        self.ref_polynomial_coefficients = []
        linenum = 1
        for line in contents.split("\n"):
            if line == "" or line[0] == "#": continue
            coefs = line.split(",")
            if len(coefs) > 11:
                if self.verbose: print "Reference polynomial has to many coefficients on line", linenum
                if self.verbose: print "line:'%s'" % line
                continue
            if len(coefs) < 3:
                if self.verbose: print "Reference polynomial does not have ateast 3 coefficients, use ',' to separate the coefficients on line", linenum
                if self.verbose: print "line:'%s'" % line
                continue
            try:
                coefs = map(float, coefs)
            except:
                if self.verbose: print "Reference polynomial has non-numeric coefficients, on line", linenum
                if self.verbose: print "line:'%s'" % line
                continue

            self.ref_polynomial_coefficients.append(coefs)


    def print_coefficients(self):
        self.poly_text.delete(1.0, Tkinter.END)
        for d, coefs in self.polynomial_coefficients.iteritems():
            if len(self.control_points) < d: break
            self.poly_text.tag_delete(str(d))

        lines = []
        for d, coefs in self.polynomial_coefficients.iteritems():
            if len(self.control_points) < d: break
            lines.append("#Dim=%s"%d)
            s = []
            for v in coefs:
                s.append( ("%." + str(self.coefficient_precision) + "f") %v)
            lines.append(",".join(s))
            lines.append("")
        self.poly_text.insert(Tkinter.END, "\n".join(lines))

        for d, coefs in self.polynomial_coefficients.iteritems():
            if len(self.control_points) < d: break
            tag_begin = str(((d-2)*3 + 2)) + ".0"
            tag_end = str(((d-2)*3 + 3)) + ".end"
            self.poly_text.tag_add(str(d), tag_begin, tag_end)
            self.poly_text.tag_config(str(d), foreground=self.polynomial_colors[d])

    def polynomial_x_points(self):
        x_points = [float(x)*(self.xmax - self.xmin)/self.nbins + self.xmin for x in range(self.nbins)]
        return x_points


    def fit_polynomials_using_polyfit(self, xs, ys, d):
        # silence RankWarning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            coefs = numpy.polyfit(xs, ys, deg=d-1)
        return coefs


    def fit_polynomials_using_cvxopt(self, xs, ys, d):
        coefs = fit_polynomial(matrix(xs), matrix(ys), d)
        coefs = numpy.array(coefs.T)[0]
        coefs = coefs[::-1]
        return coefs

    def fit_polynomials_constrained_using_cvxopt(self, xs, ys, d, minima, constraints):
        if len( minima[0] ) != len( minima[1]) or len(constraints[0]) != len(constraints[1]) :
            if self.verbose: print "Minima lists are not of the same length"
            coefs = fit_polynomial( matrix(xs), matrix(ys), d )
        #elif len( minima[0] ) == 1 :
        #    print "Fitting with a single constraint"
        #    coefs = fit_polynomial_one_constraint(matrix(xs), matrix(ys), d, minima[0][0], minima[1][0])
        else :
            try :
                coefs = fit_polynomial_multiple_constraints(matrix(xs), matrix(ys), d, minima[0], minima[1],constraints[0],constraints[1])
                if self.verbose: print "successfully fit multiply constrained polynomial of degree", d
            except :
                if self.verbose: print "Could not satisfy constraints"
                coefs = fit_polynomial(  matrix(xs), matrix(ys), d )
        coefs = numpy.array(coefs.T)[0]
        coefs = coefs[::-1]
        return coefs

    def get_x_polyval_points(self):
        x_points = self.polynomial_x_points()
        if self.fit_normal :
            x_polyval_points = x_points
        elif self.fit_deg_as_rad :
            x_polyval_points = [ x * math.pi / 180 for x in x_points ]
        elif self.fit_deg_as_cos :
            x_polyval_points = [ numpy.cos( math.pi/180 * ( 180 - x )) for x in x_points ]
        return x_polyval_points

    def compute_polynomial_fits(self):
        x_points = self.polynomial_x_points()
        x_polyval_points = self.get_x_polyval_points()

        self.polynomial_points = {}

        if len(self.control_points) == 0: return
        for d in self.polynomial_dimensions:
            if len(self.control_points) < d: break
            if self.fit_normal :
                xs = [ p[0] for p in self.control_points ]
            elif self.fit_deg_as_rad :
                xs = [ p[0] * math.pi / 180 for p in self.control_points ]
            elif self.fit_deg_as_cos :
                xs = [numpy.cos(math.pi/180*(180-p[0])) for p in self.control_points]
            ys = [p[1] for p in self.control_points]


# This is used to reflect the points across a line. The point of this
# was to control the derivative at a point Now just set the point of
# the minima.
#            try:
#                ry = float(self.reflect_var.get())
#            except:
#                ry = None
#            if self.reflect_check_var.get() and ry is not None:
#                xs.extend([2*ry - p[0] for p in self.control_points])
#                ys.extend(ys)


            #coefs = self.fit_polynomials_using_polyfit(xs, ys, d)
            if ( self.minima[0] is None or self.minima[1] is None ) and ( self.constraints[0] is None or self.constraints[1] is None ):
                if self.verbose: print "Fitting without constraints"
                coefs = self.fit_polynomials_using_cvxopt(xs, ys, d)
            else:
                minima = [ [], self.minima[1] if self.minima[1] is not None else [] ]
                if self.minima[0] is not None :
                    if self.fit_normal :
                        minima[0] = self.minima[0]
                    elif self.fit_deg_as_cos :
                        minima[0] = [ math.cos( math.pi - x*math.pi/180 ) for x in self.minima[0] ]
                    elif self.fit_deg_as_rad :
                        minima[0] = [ x*math.pi/180 for x in self.minima[0] ]
                                    
                constraints = [ [], self.constraints[1] if self.constraints[1] is not None else [] ]
                if self.constraints[0] is not None :
                    if self.fit_normal :
                        constraints[0] = self.constraints[0]
                    elif self.fit_deg_as_cos :
                        constraints[0] = [ math.cos( math.pi - x*math.pi/180 ) for x in self.constraints[0] ]
                    elif self.fit_deg_as_rad :
                        constraints[0] = [ x*math.pi/180 for x in self.constraints[0] ]
                                    
                coefs = self.fit_polynomials_constrained_using_cvxopt(xs, ys, d, minima, constraints)

            self.polynomial_coefficients[d] = coefs

            self.polynomial_points[d] = \
                (x_points, numpy.polyval(numpy.poly1d(coefs), x_polyval_points))

    def compute_reference_polynomial_points(self):
        self.ref_polynomial_points = {}

        x_points = self.polynomial_x_points()
        x_polyval_points = self.get_x_polyval_points()
        for i in range(len(self.ref_polynomial_coefficients)):
            ref_poly_coefs = self.ref_polynomial_coefficients[i]

            self.ref_polynomial_points[i] = \
                (x_points, numpy.polyval(numpy.poly1d(ref_poly_coefs), x_polyval_points))

    def draw_demo(self):
        self.draw_polynomial_fits(self.ref_polynomial_points, fixed_color="darkgray")
        if not self.hide_polynomials: self.draw_polynomial_fits(self.polynomial_points)
        self.draw_control_points()

    def draw_polynomial_fits(self, poly_points, fixed_color=None):
        canv_points = []
        for d, (xs, ys) in poly_points.iteritems():
            for i in range(len(xs)):
                px, py = self.p2canv(xs[i], ys[i])
                canv_points.append((d,px,py))

#        for i in range(len(canv_points)):
#            d,px,py = canv_points[i]
##            line_color = self.polynomial_colors[d]
#            px0,px1 = px - self.dot_width / 2, px + self.dot_width / 2
#            py0,py1 = py - self.dot_width / 2, py + self.dot_width / 2
#
#            if i == len(canv_points)-1:
#                self.canvas.create_text( px+10, py+10, text="dim=%s" %d)


        for i in range(1,len(canv_points)):
            last_i =  i-1

            ad, ax, ay = canv_points[last_i]
            d, px,py = canv_points[i]
            if ad != d: continue

            if fixed_color is None:
                line_color = self.polynomial_colors[d]
            else:
                line_color = fixed_color

            self.canvas.create_line( ax, ay, px, py, fill=line_color)


    def draw_control_points(self):

        self.control_point_widgets = []
        for x, y in self.control_points:
            px, py = self.p2canv(x, y)
            px0,px1 = px - self.dot_width / 2, px + self.dot_width / 2
            py0,py1 = py - self.dot_width / 2, py + self.dot_width / 2
            t = self.canvas.create_oval( px0, py0, px1, py1, fill=self.dot_color )
            self.control_point_widgets.append(t)



    def event_lclick(self, event):
        if event.widget == self.canvas:
            px = self.canvas.canvasx( event.x )
            py = self.canvas.canvasy( event.y )

            self.control_points.append(self.canv2p(px,py))

            self.redraw_everything()

    def event_rclick(self, event):
        if event.widget != self.canvas:
            if self.verbose: print "not canvas event"
            return


        if len(self.control_points) == 0:
            if self.verbose: print "no points to delete"
            return

        px = self.canvas.canvasx( event.x )
        py = self.canvas.canvasy( event.y )

        if self.verbose: print px, py, lineno()
        t = self.canvas.find_closest(px, py)
        if self.verbose: print lineno(), "found closest widget", t
        if self.verbose: print lineno(), self.control_point_widgets
        if self.verbose: print lineno(), self.control_points
        for i in range(len(self.control_point_widgets)):
            if t[0] == self.control_point_widgets[i]:
                if self.verbose: print lineno(), "deleting point", t, self.control_points[i]
                self.canvas.delete(t)

                self.control_points = \
                    self.control_points[:i] + self.control_points[i+1:]

        self.redraw_everything()

    def append_to_title( self, title ) :
        self.main.title( "Polynomial Fitting -- " + title )

    def set_fit_deg_as_rad( self ) :
        self.fit_normal = False
        self.fit_deg_as_rad = True
        self.fit_deg_as_cos = False
        self.main.title("Polynomial Fitting -- plotting degrees, fitting in radians")

    def set_fit_deg_as_cos( self ) :
        self.fit_normal = False
        self.fit_deg_as_rad = False
        self.fit_deg_as_cos = True
        self.main.title("Polynomial Fitting -- plotting degrees, fitting in cosine of the exterior angle")

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p:
        p.flag("fit_degrees_as_radians").conflicts( p.flag("fit_degrees_as_cosexterior"))
        p.str("reference_polynomial").shorthand("p").multiple()
        p.all_if_any( 
            p.str("x0").described_as("comma-separated list of floats for the x values where the fit polynomial should pass through a specific point and should have a derivative of zero"),
            p.str("y0").described_as("comma-separated list of floats for the y values that correspond to the x0s")
            )
        p.all_if_any(
            p.str("x1").described_as("comma-separated list of floats for the x values where the fit polynomial should pass through a specific point"),
            p.str("y1").described_as("comma-separated list of floats for the y values that correspond to the x1s")
            )
        p.str("range_x").described_as("comma-separated list of two floats for the lower- and upper-boundaries for the x axis")
        p.str("range_y").described_as("comma-separated list of two floats for the lower- and upper-boundaries for the y axis")
        p.str("control_points").described_as("comma-separated list of 2n floats representing the x and y coordinates of n control points: e.g. x1,y1,x2,y2,x3,y3")
        p.flag("verbose")
        p.str("title").shorthand( "t" )
        
    initial_control_points = []

    # read in any of the control points given on the command line
    if control_points :
        pts = control_points.split(",")
        if len(pts) % 2 != 0 : 
            print "Error the number of control points must be even"
            sys.exit(1)
        fpts = []
        for pt in pts :
            try :
                fpt = float( pt )
                fpts.append( fpt )
            except ValueError :
                print ( "Error, could not convert one of the control points to a floating point number: %s" % pt )
                sys.exit(1)
        for i in xrange(len(pts)/2) :
            initial_control_points.append( ( fpts[2*i], fpts[2*i+1] ))
        
    d = PolyFit(initial_control_points)
    if fit_degrees_as_radians :
        d.set_fit_deg_as_rad()
    elif fit_degrees_as_cosexterior :
        d.set_fit_deg_as_cos()

    if len( reference_polynomial ) != 0 :
        for pstring in reference_polynomial :
            d.add_reference_polynomial( pstring )

    if x0: d.set_minima_x0s( x0 )
    if y0: d.set_minima_y0s( y0 )

    if x1: d.set_constraints_x1s( x1 )
    if y1: d.set_constraints_y1s( y1 )

    if range_x: d.set_range_x( range_x )
    if range_y: d.set_range_y( range_y )

    if title : d.append_to_title( title )

    d.go()

