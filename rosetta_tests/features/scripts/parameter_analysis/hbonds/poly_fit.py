

import Tkinter, tkFileDialog, tkSimpleDialog
import numpy
import inspect, warnings

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno

class PolyFit:

    def __init__(self, initiali_control_points):

        self.main = Tkinter.Tk()
        self.main.title("Polynomial Fitting")
        
        self.canvas_width =  750
        self.canvas_height = 750
        self.plot_width = 680
        self.plot_height = 680
        self.left_gutter = 10
        self.top_gutter = 10

        self.xmin = -1
        self.xmax = 1

        self.ymin = -1
        self.ymax = 1

        self.xaxis = 0
        self.yaxis = 0

        self.xaxis_ticks = [-1, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1]
        self.xaxis_tick_height = 5
        self.yaxis_ticks = [-1, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1]
        self.yaxis_tick_width = 5

        self.axis_color = "white"


        self.polynomial_dimensions = [2,3,4,5,6,7,8]
        self.nbins = 200

        self.coefficient_precision = 3

        self.dot_color = 'blue'
        self.dot_width = 6

        self.polynomial_colors = {2:"#D53E4F", 3:"#FC8D59", 4:"#FEE08B", 5:"#FFFFBF", 6:"#E6F598", 7:"#99D594", 8:"#3288BD"}

        self.poly_text_frame = Tkinter.Frame(self.main)
        self.poly_label = Tkinter.Label(self.poly_text_frame, text="Polynomial Coefficients: smallest coefficient first")
        self.poly_label.pack(side="top")
        self.poly_text = Tkinter.Text(\
            self.poly_text_frame, width = 50, background="lightgray")
        self.poly_text.pack(side="bottom")
        self.poly_text_frame.pack(side="right")

        self.canvas = Tkinter.Canvas( self.main,height=self.canvas_height, width=self.canvas_width  )
        self.canvas.pack(side="top")
        
        self.btn_makeart = Tkinter.Button( self.main, text="clear", command=self.clear )
        self.btn_makeart.pack(side="top")

        self.main.bind("<Button-1>", self.event_lclick )
        self.main.bind("<Button-3>", self.event_rclick )


        self.control_points = initial_control_points

        self.control_point_widgets = []
        self.polynomial_points = {}
        self.polynomial_coefficients = {}
        self.redraw_everything()
    

#        self.test_diagnostics()

 
        self.main.mainloop()


    def test_diagnostics(self):

        self.control_points = [(0,0), (1,1), (2,0)]
        self.compute_polynomial_fits()

        self.redraw_everything()



    def p2canv( self, x, y):

        px = float(x)*self.plot_width/(self.xmax-self.xmin) + \
            self.left_gutter + self.plot_width/2
        py = -1*float(y)*self.plot_height/(self.ymax-self.ymin) + \
            self.top_gutter + self.plot_width/2

        return (px, py)

    def canv2p(self, px, py):
        x = float(px - self.plot_width/2 - self.left_gutter)
        x = x*(self.xmax-self.xmin)/self.plot_width
        y = float(py - self.plot_height/2 - self.top_gutter)
        y = -y*(self.ymax-self.ymin)/self.plot_height

        return (x, y)
    
    def clear( self ):
        for t in self.canvas.find_all():
            self.canvas.delete(t)

        #primal axes
        #xaxis
        px0, py0 = self.p2canv(self.xmin, self.xaxis)
        px1, py1 = self.p2canv(self.xmax, self.xaxis)
        self.canvas.create_line( px0, py0, px1, py1, fill=self.axis_color)

        #yaxis
        px0, py0 = self.p2canv(self.yaxis, self.ymin)
        px1, py1 = self.p2canv(self.yaxis, self.ymax)
        self.canvas.create_line( px0, py0, px1, py1, fill=self.axis_color)


        for x in self.xaxis_ticks:
            px, py = self.p2canv(x,0)
            px0,px1 = px, px 
            py0,py1 = py, py + self.xaxis_tick_height
            self.canvas.create_line( px0, py0, px1, py1, fill=self.axis_color)
            self.canvas.create_text( px1, py1+5, text=str(x), font=("helvetica", 7)) 
            
        for y in self.yaxis_ticks:
            px, py = self.p2canv(0,y)
            px0,px1 = px, px + self.yaxis_tick_width
            py0,py1 = py, py
            self.canvas.create_line( px0, py0, px1, py1, fill=self.axis_color)
            self.canvas.create_text( px1+5, py1, text=str(y), font=("helvetica", 7)) 


    def redraw_everything( self ):

        self.clear()
        self.compute_polynomial_fits()
        self.draw_demo()
        self.print_coefficients()
        print "Control point coordinates:"
        print self.control_points
        print ""
        
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
                s.append("%.4f" %v)
            s.reverse() # to get coefficients from low to high
            lines.append(",".join(s))
            lines.append("")
        self.poly_text.insert(Tkinter.END, "\n".join(lines))

        for d, coefs in self.polynomial_coefficients.iteritems():
            if len(self.control_points) < d: break
            tag_begin = str(((d-2)*3 + 2)) + ".0"
            tag_end = str(((d-2)*3 + 3)) + ".end"
            self.poly_text.tag_add(str(d), tag_begin, tag_end)
            self.poly_text.tag_config(str(d), foreground=self.polynomial_colors[d])


    def compute_polynomial_fits(self):

        x_points = [float(x)*(self.xmax - self.xmin)/self.nbins + self.xmin for x in range(self.nbins)]
        

        self.polynomial_points = {}

        if len(self.control_points) == 0: return
        for d in self.polynomial_dimensions:
            if len(self.control_points) < d: break
            xs = [p[0] for p in self.control_points]
            ys = [p[1] for p in self.control_points]

            # silence RankWarning
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                coefs = numpy.polyfit(xs,ys,deg=d)

            self.polynomial_coefficients[d] = coefs

            self.polynomial_points[d] = \
                (x_points, numpy.polyval(numpy.poly1d(coefs), x_points))

    def draw_demo(self):
        canv_points = []
        for d, (xs, ys) in self.polynomial_points.iteritems():
            for i in range(len(xs)):
                px, py = self.p2canv(xs[i], ys[i])
                canv_points.append((d,px,py))
        
        for i in range(len(canv_points)):
            d,px,py = canv_points[i]
            dot_color = self.polynomial_colors[d]
            px0,px1 = px - self.dot_width / 2, px + self.dot_width / 2
            py0,py1 = py - self.dot_width / 2, py + self.dot_width / 2
            
#            self.canvas.create_oval( px0, py0, px1, py1, fill=dot_color )
            if i == len(canv_points):
                self.canvas.create_text( px+10, py+10, text="dim=%s" %d)
            

        for i in range(1,len(canv_points)):
            last_i =  i-1
            
            ad, ax, ay = canv_points[last_i]            
            d, px,py = canv_points[i]
            if ad != d: continue

            dot_color = self.polynomial_colors[d]


            self.canvas.create_line( ax, ay, px, py, fill=dot_color)


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
            print "not canvas event"
            return
        

        if len(self.control_points) == 0:
            print "no points to delete"
            return

        px = self.canvas.canvasx( event.x )
        py = self.canvas.canvasy( event.y )

        print px, py, lineno()
        t = self.canvas.find_closest(px, py)
        print lineno(), "found closest widget", t
        print lineno(), self.control_point_widgets
        print lineno(), self.control_points
        for i in range(len(self.control_point_widgets)):
            if t[0] == self.control_point_widgets[i]:
                print lineno(), "deleting point", t, self.control_points[i]
                self.canvas.delete(t)

                self.control_points = \
                    self.control_points[:i] + self.control_points[i+1:]

        self.redraw_everything()


initial_control_points = [(0.17941176470588235, 0.21764705882352942), (-0.29999999999999999, 0.34705882352941175), (-0.49411764705882355, -0.055882352941176473), (-0.53823529411764703, -0.28235294117647058), (-0.30588235294117649, -0.35588235294117648)]


d = PolyFit(initial_control_points)


        
