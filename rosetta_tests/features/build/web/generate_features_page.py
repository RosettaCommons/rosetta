import glob


def create_features_page(output_fname, template_fname, plots_dir):
    plot_link_template =\
        '<a id="%s" href="%s" class="highslide" onclick="return hs.expand(this)"> %s </a><br/>'

    f = open(template_fname,'r')
    features_str = f.read()
    f.close()
    plot_fnames = glob.glob(plots_dir + "/*.png")
    plots_link_html = []
    for plot_fname in plot_fnames:
        name = plot_fname[len(plots_dir)+1:][:-4]
        plot_link_html = plot_link_template % (plot_fname, plot_fname, name)
        plots_link_html.append(plot_link_html)

    plots_link_html_str = "\n\t" + "\n\t".join(plots_link_html)
    features_str = features_str % {"plots" : plots_link_html_str}
    f = open(output_fname,'w')
    f.write(features_str)
    f.close()

create_features_page("features.html", "features.html.TEMPLATE", ".")


