#!/usr/bin/env python3

"""
Author:
    Name:           Cody Alexander Ramirez
    Email:          cody.ramirez@stjude.org
    Affiliation:    St. Jude Children's Research Hospital, Memphis, TN
    Date:           July 16th, 2021

Creates a scatter plot between the test sample vs. matched PDX tumor log2(CPM)

Usage:
    python3 ScatterPlotCorrelation.py <INPUT>

    Example:
        python3 ScatterPlotCorrelation.py TextFile.txt
            OR
        python3 ScatterPlotCorrelation.py DirectoryPath

    File header must contain the following colums:
    1) geneID
    2) geneSymbol
    3) Two columns starting with log2cpm_ for both samples

Myogenic gene markers are highlighted (23 genes that are differentially expressed in 1 of the 3 cell populations identified in this tumor type).
    Each subset of genes are highlighted for different cell populations (3 colors total)
        1) Mesoderm is Orange
        2) Myoblast is Cyan
        3) Myocytes is Lime
    Two genes are plotted to highlight the two major subtypes within this tumor group
        1) HMGA2 is for ERMS
        2) NOS1 is for ARMS

The coefficient of determination (R^2) is the proportion of the variance in the dependent variable that is predictable from the independent variable(s).
    1) Black includes all data points
    2) Orange includes Mesoderm genes
    3) Cyan includes Myoblast genes
    4) Lime includes Myocyte genes

References:
    1) Law CW, Alhamdoosh M, Su S, et al. RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR. F1000research. 2016;5. DOI: 10.12688/f1000research.9005.3.
    2) Zhao, Y., Li, MC., Konaté, M.M. et al. TPM, FPKM, or Normalized Counts? A Comparative Study of Quantification Measures for the Analysis of RNA-seq Data from the NCI Patient-Derived Models Repository. J Transl Med 19, 269 (2021). https://doi.org/10.1186/s12967-021-02936-w

"""





#Importing all necessary libraries
import os, sys, glob, fnmatch
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
from scipy import stats





#This function creates the scatter plot
def PlotData(CurrentFile):

    #To contain all GeneSymbols within the file
    AllGene = []
    #To contain all of the log2(CPM) values for the first sample
    AllX = []
    #To contain all of the log2(CPM) values for the second sample
    AllY = []

    #To contain all of the Mesoderm GeneSymbols in the order listed within the file
    MesoGene = []
    #To contain all of the log2(CPM) values for the first sample
    MesoX = []
    #To contain all of the log2(CPM) values for the second sample
    MesoY = []

    #To contain all of the Myoblast GeneSymbols in the order listed within the file
    MyoBGene = []
    #To contain all of the log2(CPM) values for the first sample
    MyoBX = []
    #To contain all of the log2(CPM) values for the second sample
    MyoBY = []

    #To contain all of the Myocyte GeneSymbols in the order listed within the file
    MyoCGene = []
    #To contain all of the log2(CPM) values for the first sample
    MyoCX = []
    #To contain all of the log2(CPM) values for the second sample
    MyoCY = []

    #Gene set markers
    MesodermGenes = ["MEOX2", "PAX3", "EGFR", "CD44", "DCN", "POSTN"]
    MyoblastGenes = ["PAX7", "MYF5", "MSC", "GPC3", "VIM"]
    MyocyteGenes = ["MYOD1", "MYOG", "MEF2A", "MEF2A", "MEF2C", "TTN", "NCAM1", "MYH3", "NEB", "CDH15", "NCOA1"]
    ERMS_Subtype = ["HMGA2"]
    ARMS_Subtype = ["NOS1"]

    #Reading through every line within the current file
    for line in open(CurrentFile, "r"):
        
        #Identifying the header line within the current file
        if line.startswith("geneID"):
            
            #Identifying the column number for geneSymbol
            GeneSymbolIndex = line.split().index("geneSymbol")
            
            #Identifying the first samples ID
            XAxisID = fnmatch.filter(line.split(), 'log2cpm*')[0].split('_')[-1]
            #Identifying the log2(CPM) data for the first sample
            XAxisIDIndex = line.split().index(fnmatch.filter(line.split(), 'log2cpm*')[0])
            
            #Identifying the second samples ID
            YAxisID = fnmatch.filter(line.split(), 'log2cpm*')[1].split('_')[-1]
            #Identifying the log2(CPM) data for the second sample
            YAxisIDIndex = line.split().index(fnmatch.filter(line.split(), 'log2cpm*')[1])
        
        #Gathering data from the remaining portion of the file
        else:
            
            #Making a variable to contain the GeneSymbol
            CurrentGene = line.split()[GeneSymbolIndex]
            #Making a variable to contain the first sample's corresponding log2(CPM) value
            CurrentXAxisValue = float(line.split()[XAxisIDIndex])
            #Making a variable to contain the second sample's corresponding log2(CPM) value
            CurrentYAxisValue = float(line.split()[YAxisIDIndex])
            
            #Collecting all GeneSymbol and log2(CPM) values within the entire file
            AllGene.append(CurrentGene)
            AllX.append(CurrentXAxisValue)
            AllY.append(CurrentYAxisValue)
            
            #Collecting all GeneSymbol and log2(CPM) values within the Mesoderm gene list
            if CurrentGene in MesodermGenes:
                MesoGene.append(CurrentGene)
                MesoX.append(CurrentXAxisValue)
                MesoY.append(CurrentYAxisValue)
            
            #Collecting all GeneSymbol and log2(CPM) values within the Myoblast gene list
            if CurrentGene in MyoblastGenes:
                MyoBGene.append(CurrentGene)
                MyoBX.append(CurrentXAxisValue)
                MyoBY.append(CurrentYAxisValue)
            
            #Collecting all GeneSymbol and log2(CPM) values within the Myocyte gene list
            if CurrentGene in MyocyteGenes:
                MyoCGene.append(CurrentGene)
                MyoCX.append(CurrentXAxisValue)
                MyoCY.append(CurrentYAxisValue)
            
            #Collecting the gene and log2(CPM) values for the ERMS subtype
            if CurrentGene in ERMS_Subtype:
                ERMS_Subtype = ERMS_Subtype + [CurrentXAxisValue, CurrentYAxisValue]
            
            #Collecting the gene and log2(CPM) values for the ARMS subtype
            if CurrentGene in ARMS_Subtype:
                ARMS_Subtype = ARMS_Subtype + [CurrentXAxisValue, CurrentYAxisValue]

    #Setting the default figure size to 10x10 inches
    plt.rcParams['figure.figsize'] = [10, 10]
    #Setting the default figure dots per inches (dpi) to 1,200
    plt.rcParams['figure.dpi'] = 1200 # 1200 e.g. is really fine, but slower
    #Creating a new figure
    fig = plt.figure()

    #Generating a scatter plot of all the data within the file colored black
    plt.scatter(AllX, AllY, s=1, c="Black")
    #Generating a linear regression line based off all the data points
    slope, intercept, TotalR, p, std_err = stats.linregress(AllX, AllY)
    def myfunc(AllX):
      return slope * AllX + intercept
    #Returning a list of Y values after applying the linear regression algorithm to all X values
    mymodel = list(map(myfunc, AllX))
    #Plotting the linear regression line
    plt.plot(AllX, mymodel, color="Black", linewidth=0.5,zorder=0)

    set_difference = set(MesodermGenes).symmetric_difference(set(MesoGene))
    list_difference = list(set_difference)
    if list_difference != []:
        print("The following Mesoderm gene(s) are missing: " + str(list_difference))
    #Generating a scatter plot of all the data within the Mesoderm gene list colored orange
    plt.scatter(MesoX, MesoY, s=20, color="Orange", label="Mesoderm", edgecolors='Black', marker="D", zorder=10)
    #Generating a linear regression line based off all the data points from the Mesoderm gene list
    slope, intercept, MesoR, p, std_err = stats.linregress(MesoX, MesoY)

    set_difference = set(MyoblastGenes).symmetric_difference(set(MyoBGene))
    list_difference = list(set_difference)
    if list_difference != []:
        print("The following Myoblast gene(s) are missing: " + str(list_difference))
    #Generating a scatter plot of all the data within the Myoblast gene list colored cyan
    plt.scatter(MyoBX, MyoBY, s=20, color="Cyan", label="Myoblast", edgecolors='Black', marker="^", zorder=10)
    #Generating a linear regression line based off all the data points from the Myoblast gene list
    slope, intercept, MyoBR, p, std_err = stats.linregress(MyoBX, MyoBY)

    set_difference = set(MyocyteGenes).symmetric_difference(set(MyoCGene))
    list_difference = list(set_difference)
    if list_difference != []:
        print("The following Myocytes gene(s) are missing: " + str(list_difference))
    #Generating a scatter plot of all the data within the Myocyte gene list colored lime
    plt.scatter(MyoCX, MyoCY, s=20, color="Lime", label="Myocytes", edgecolors='Black', zorder=10)
    #Generating a linear regression line based off all the data points from the Myocyte gene list
    slope, intercept, MyoCR, p, std_err = stats.linregress(MyoCX, MyoCY)

    #Creating a custom legend box to contain the R^2 values calculated above
    ax = plt.gca()
    texts = ["$\mathregular{R^{2}}$" + " = {:.2f}".format(TotalR**2),
            "$\mathregular{R^{2}}$" + " = {:.2f}".format(MesoR**2),
            "$\mathregular{R^{2}}$" + " = {:.2f}".format(MyoBR**2),
            "$\mathregular{R^{2}}$" + " = {:.2f}".format(MyoCR**2)]
    colors = ["Black", "Orange", "Cyan", "Lime"]
    Texts = []
    for t,c in zip(texts,colors):
        Texts.append(TextArea(t, textprops=dict(color=c)))
        texts_vbox = VPacker(children=Texts, pad=0, sep=0)
        ann = AnnotationBbox(texts_vbox,(0.02, 0.95),
                                xycoords=ax.transAxes,
                                box_alignment=(0,.5),
                                bboxprops=dict(color='black', facecolor='Dimgrey', boxstyle='round'))

    if len(ERMS_Subtype) == 3:
        #Plotting the ERMS subtype gene on the plot
        plt.annotate(ERMS_Subtype[0], (ERMS_Subtype[1], ERMS_Subtype[2]), ha="center", color="Red", path_effects=[pe.withStroke(linewidth=3, foreground="White")])
    else:
        print("HMGA2 (ERMS subtype gene marker) is not present within the file.")

    if len(ARMS_Subtype) == 3:
        #Plotting the ARMS subtype gene on the plot
        plt.annotate(ARMS_Subtype[0], (ARMS_Subtype[1], ARMS_Subtype[2]), ha="center", color="Red", path_effects=[pe.withStroke(linewidth=3, foreground="White")])
    else:
        print("NOS1 (ARMS subtype gene marker) is not present within the file.")

    #Plotting the graph legend markers
    plt.legend(title="Gene set markers", loc = "lower right")

    #Plotting the X label with the first sample
    plt.xlabel(str(XAxisID) + " $\mathregular{log_{2}(CPM)}$", fontsize="large")
    #Plotting the Y label with the second sample
    plt.ylabel(str(YAxisID) + " $\mathregular{log_{2}(CPM)}$", fontsize="large")

    #Ensuring that the X and Y tick marks are squared up by using the min and max from all data
    plt.xticks(np.arange(round(min(AllX+AllY), 0), round(max(AllX+AllY), 0)+1, 1.0))
    plt.yticks(np.arange(round(min(AllX+AllY), 0), round(max(AllX+AllY), 0)+1, 1.0))

    #Adding the annotation box to the figure
    ann.set_figure(fig)
    fig.artists.append(ann)
    plt.tight_layout()

    plt.savefig(YAxisID +"_vs_"+ XAxisID +".png", format="png")
    print(YAxisID +"_vs_"+ XAxisID +".png produced")

    #Saving the figure using the first sample and second sample's specific IDs
    plt.savefig(YAxisID +"_vs_"+ XAxisID +".svg", format="svg")
    print(YAxisID +"_vs_"+ XAxisID +".svg produced")

    plt.close(fig)





#This is the main portion of the program
#Check that only 2 command line arguments were given.
if (len(sys.argv) != 2):
    #Otherwise; prints an error statement, the documentation and then exits the program.
    sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

#Determines if the iputted data is a text file
elif "txt" in sys.argv[1]:
    print("\nFile inputted")
    #Generate the scatter plot
    PlotData(sys.argv[1])

#If the inputted data isn't a text file, I assume it is a directory path
else:
    print("\nDirectory path inputted")
    #Ittirates through all txt files within the inputted directory
    for filename in glob.glob(os.path.join(sys.argv[1], '*.txt')):
        print()
        print(filename)
        #Generates a scatter plot for each file
        PlotData(filename)
        print()

#Notifies the user when the program has completed.
print("Figure generation completed.\n")




