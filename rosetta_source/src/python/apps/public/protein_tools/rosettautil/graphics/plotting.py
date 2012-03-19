import Gnuplot,Gnuplot.funcutils

def barplot(data_dict,title,outfile):
    plotter = Gnuplot.Gnuplot()
    plotter("set term postscript enhanced color")
    plotter("set output \""+outfile+ "\"")
    plotter("set boxwidth 1")
    plotter("unset xtics")
    plotter("set style fill solid border -1")
    plotter("set title \""+title+"\"")
    plotter("set yrange[0:"+str(data_dict[max(data_dict)])+"]")
    xtics = []
    datapoints = []
    counter = 1
    for key in data_dict:
        xtics.append("\""+key+"\" "+str(counter))
        height = float(data_dict[key])
        datapoints.append( Gnuplot.Data(counter,height,with="boxes") )
        counter += 1 

    plotter("set xtics(" + ",".join(xtics)  + ")")

    plotter.plot(*datapoints)


    
