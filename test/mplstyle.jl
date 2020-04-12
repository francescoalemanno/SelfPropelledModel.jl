using PyCall
using PyPlot
function matplotlibstyle()
    py"""
    from matplotlib import rcParams
    import numpy as np
    def mplfigsize(hscale, publicationcols=1,columnsinches=3.1):
        bonuscol=-0.12
        W=columnsinches*(publicationcols+bonuscol)
        return [W, 0.6180469715698392*hscale*W]
    def setfonts(fs):
        smaller=fs*9/11
        rcParams.update({
        "font.size": fs,
        'xtick.labelsize': smaller,
        'axes.labelsize': fs,
        'ytick.labelsize': smaller,
        'legend.fontsize': smaller,
        })
    def setstyle():
        rcParams.update({
        #    "pgf.texsystem": "pdflatex",
        #    "font.family": "serif",
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
        #    "text.usetex": True,
        #    "pgf.rcfonts": False,
            "figure.figsize": mplfigsize(1),
            "figure.autolayout": True,
            "figure.dpi": 300,
            'axes.spines.left': True,
            'axes.spines.bottom': True,
            'axes.spines.right': False,
            'axes.spines.top': False,
            'ytick.left': True,
            'xtick.bottom': True,
            'ytick.right': False,
            'xtick.top': False,
            'ytick.direction': 'in',
            'xtick.direction': 'in',
            'xtick.color': 'black',
            'ytick.color': 'black',
            'axes.edgecolor': 'black',
            'legend.frameon': False,
            'legend.columnspacing': 0.0,
            'legend.labelspacing': 0.2,
            'legend.borderaxespad': 0.5,
            'legend.borderpad': 0.0,
            'legend.handletextpad': 0.2,
            'legend.handlelength': 1.5,
            'lines.linewidth': 1.,
        })
        setfonts(9)
    setstyle()
    """
end
matplotlibstyle()
setstyle()=py"setstyle"()
mplfigsize(hscale, pubcols=1,colinches=3.1)=py"mplfigsize"(hscale, publicationcols=pubcols,columnsinches=colinches)
setfonts(x=9)=py"setfonts"(x)
