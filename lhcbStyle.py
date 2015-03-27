from ROOT import TStyle, gROOT, TLatex, TPaveText, TText

lhcbNames = []
lhcbStyle = None
lhcbLabel = None
lhcbLatex = None 

def setLHCbStyle():
    global lhcbStyle
    global lhcbText
    global lhcbLatex
    
    lhcbStyle = TStyle("lhcbStyle","Standard LHCb plots style")

    # use times new roman
    lhcbFont = 132
    # line thickness
    lhcbWidth = 2 
    lhcbTSize = 0.06

    #// use plain black on white colors
    lhcbStyle.SetFrameBorderMode(0)
    lhcbStyle.SetCanvasBorderMode(0)
    lhcbStyle.SetPadBorderMode(0)
    lhcbStyle.SetPadColor(0)
    lhcbStyle.SetCanvasColor(0)
    lhcbStyle.SetStatColor(0)
    lhcbStyle.SetPalette(55)
    
    lhcbStyle.SetLegendBorderSize(0)
    lhcbStyle.SetLegendFont(132)
    lhcbStyle.SetFillColor(1);
    lhcbStyle.SetFillStyle(1001)
    
    # set the paper & margin sizes
    lhcbStyle.SetPaperSize(20,26)
    lhcbStyle.SetPadTopMargin(0.15)
    lhcbStyle.SetPadRightMargin(0.15) 
    lhcbStyle.SetPadBottomMargin(0.15)
    lhcbStyle.SetPadLeftMargin(0.15)
    
    # use large fonts
    lhcbStyle.SetTextFont(lhcbFont)
    lhcbStyle.SetTextSize(lhcbTSize)
    #  lhcbStyle.SetTextSize(0.08)
    lhcbStyle.SetLabelFont(lhcbFont,"x")
    lhcbStyle.SetLabelFont(lhcbFont,"y")
    lhcbStyle.SetLabelFont(lhcbFont,"z")
    lhcbStyle.SetLabelSize(lhcbTSize,"x")
    lhcbStyle.SetLabelSize(lhcbTSize,"y")
    lhcbStyle.SetLabelSize(lhcbTSize,"z")
    lhcbStyle.SetTitleFont(lhcbFont)
    lhcbStyle.SetTitleFont(lhcbFont,"x");
    lhcbStyle.SetTitleFont(lhcbFont,"y");
    lhcbStyle.SetTitleFont(lhcbFont,"z");
    lhcbStyle.SetTitleSize(1.2*lhcbTSize,"x")
    lhcbStyle.SetTitleSize(1.2*lhcbTSize,"y")
    lhcbStyle.SetTitleSize(1.2*lhcbTSize,"z")
    
    # use bold lines and markers
    lhcbStyle.SetLineWidth(lhcbWidth)
    lhcbStyle.SetFrameLineWidth(lhcbWidth)
    lhcbStyle.SetHistLineWidth(lhcbWidth)
    lhcbStyle.SetFuncWidth(lhcbWidth)
    lhcbStyle.SetGridWidth(lhcbWidth)
    lhcbStyle.SetLineStyleString(2,"[12 12]") 
    lhcbStyle.SetMarkerStyle(20)
    lhcbStyle.SetMarkerSize(1.0)
    
    # label offsets
    lhcbStyle.SetLabelOffset(0.010)


    #titles
    lhcbStyle.SetTitleOffset(0.95,"X")
    lhcbStyle.SetTitleOffset(0.95,"Y")
    lhcbStyle.SetTitleOffset(1.2,"Z")
    lhcbStyle.SetTitleFillColor(0)
    lhcbStyle.SetTitleStyle(0)
    lhcbStyle.SetTitleBorderSize(0)
    lhcbStyle.SetTitleFont(lhcbFont,"title")
    lhcbStyle.SetTitleX(0.0)
    lhcbStyle.SetTitleY(1.0) 
    lhcbStyle.SetTitleW(1.0)
    lhcbStyle.SetTitleH(0.05)
    
    # by default, do not display histogram decorations:
    lhcbStyle.SetOptStat(0)  
    #lhcbStyle.SetOptStat("emr")     # show only nent -e , mean - m , rms -r
    #lhcbStyle.SetStatFormat("6.3g") # specified as c printf options
    lhcbStyle.SetOptTitle(0)
    lhcbStyle.SetOptFit(0)
    #lhcbStyle.SetOptFit(1011) # order is probability, Chi2, errors, parameters

    # look of the statistics box:
    lhcbStyle.SetStatBorderSize(0)
    lhcbStyle.SetStatFont(lhcbFont)
    lhcbStyle.SetStatFontSize(0.05)
    lhcbStyle.SetStatX(0.9)
    lhcbStyle.SetStatY(0.9)
    lhcbStyle.SetStatW(0.25)
    lhcbStyle.SetStatH(0.15)

    # put tick marks on top and RHS of plots
    lhcbStyle.SetPadTickX(1)
    lhcbStyle.SetPadTickY(1)

    # histogram divisions: only 5 in x to avoid label overlaps
    lhcbStyle.SetNdivisions(505,"x")
    lhcbStyle.SetNdivisions(505,"y")
    lhcbStyle.SetNdivisions(505,"z")


    # define style for text
    lhcbLabel =  TText()
    lhcbLabel.SetTextFont(lhcbFont)
    lhcbLabel.SetTextColor(1)
    lhcbLabel.SetTextSize(0.04)
    lhcbLabel.SetTextAlign(12)

    # define style of latex text
    lhcbLatex = TLatex()
    lhcbLatex.SetTextFont(lhcbFont)
    lhcbLatex.SetTextColor(1)
    lhcbLatex.SetTextSize(0.04)
    lhcbLatex.SetTextAlign(12)

    # set this style
    gROOT.SetStyle("lhcbStyle")
    gROOT.ForceStyle()
    return 


def printLHCb( optLR = 'L', isPrelim = False , optText = ''):
    global lhcbStyle
    global lhcbNames
    
    lhcbName = None
    
    if ( optLR is 'R' ):
        lhcbName = TPaveText( 0.70 - lhcbStyle.GetPadRightMargin(),
                              0.85 - lhcbStyle.GetPadTopMargin(),
                              0.95 - lhcbStyle.GetPadRightMargin(),
                              0.95 - lhcbStyle.GetPadTopMargin(),
                              "BRNDC" )
    elif ( optLR is 'L' ):
        lhcbName = TPaveText(lhcbStyle.GetPadLeftMargin() + 0.05,
                             0.85 - lhcbStyle.GetPadTopMargin(),
                             lhcbStyle.GetPadLeftMargin() + 0.30,
                             0.95 - lhcbStyle.GetPadTopMargin(),
                             "BRNDC")
    elif ( optLR is 'BR' ):
        lhcbName = TPaveText( 0.70 - lhcbStyle.GetPadRightMargin(),
                              0.05 + lhcbStyle.GetPadBottomMargin(),
                              0.95 - lhcbStyle.GetPadRightMargin(),
                              0.15 + lhcbStyle.GetPadBottomMargin(),
                              "BRNDC" )

    if ( isPrelim  ):
        lhcbName.AddText('#splitline{LHCb}{#scale[1.0]{Preliminary}}')
    else:
        lhcbName.AddText('LHCb')
        
    
    lhcbName.SetFillColor(0)
    lhcbName.SetTextAlign(12)
    lhcbName.SetBorderSize(0)
    lhcbName.Draw()


    lhcbNames += [ lhcbName ]
    return 

