
import dicom
import os
from os import listdir
from os.path import isfile, isdir, join
import numpy as np
from matplotlib import pyplot, cm
from matplotlib.widgets import Slider
from matplotlib.backends.backend_pdf import PdfPages


subjDirs = ["01","02","03","04","05","06","07"]
dirM = "PRE-T1"
ext1 = ".npy"
ext2 = ".dcm"
outputMain = "C:\\XXX\\Results\\"
outputMainDir1 = join(outputMain,"MLeigen\\")
outputMainDir2 = join(outputMain,"MRmaps\\")
namesOutput1 = ["TI"]
namesOutput2 = []#["LLcorrection"]
#["A","B","TI","dataT1star","dataXsquare","minIndexImg","statusImg","numIterImg"]

save_flag = False

defMinInt = 0
defMaxInt = 2000

def myPlotImg(ArrayDicom,title="", minSc=[], maxSc=[], flagShow=True):
    if not len(minSc):
        minSc = ArrayDicom.min()
    else:
        minSc = minSc[0]
    if not maxSc:
        maxSc = ArrayDicom.max()
    else:
        maxSc = maxSc[0]
    
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    #fig.subplots_adjust(left=0.25, bottom=0.25)

    im = ax.imshow(np.flipud(ArrayDicom),interpolation='none',origin='lower')
    fig.colorbar(im)
    ax.set_aspect('equal', 'datalim')
    im.set_cmap(pyplot.gray())
    im.set_clim([minSc,maxSc])
    #fig.tight_layout()

    axcolor = 'lightgoldenrodyellow'
    #axIntMin = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    #axIntMax = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    #sIntMin = Slider(axIntMin, 'Min', minSc, maxSc, valinit=minSc)
    #sIntMax = Slider(axIntMax, 'Max', minSc, maxSc, valinit=maxSc)
    def update(val):
        im.set_clim([sIntMin.val,sIntMax.val])
        fig.canvas.draw()
    #sIntMin.on_changed(update)
    #sIntMax.on_changed(update)
    pyplot.title(title)

    if(flagShow):
        pyplot.show(block=False)
    return [fig,ax]


for subj in subjDirs:
    pathC1 = join(outputMainDir1,subj,dirM)
    print(pathC1)
    if(namesOutput2):
        pathC2 = join(outputMainDir2,subj,dirM,namesOutput2[0])
    else:
        pathC2 = join(outputMainDir2,subj,dirM)
    print(pathC2)

    print("sort out files")
    extf1 = namesOutput1[0] + ext1
    extf2 = namesOutput1[0] + ext2
    onlyfiles = [f for f in listdir(pathC1) if (isfile(join(pathC1, f)) and f.endswith(namesOutput1[0]+ext1))]

    if(save_flag):
        pdf = PdfPages(join(outputMain,(dirM+"_"+namesOutput1[0]+"_"+subj+".pdf"))) 
    for cc, f in enumerate(onlyfiles):
        digs = list()
        for i,x in enumerate(f):
            if (f[i].isdigit()==True):
                digs.append(f[i])
            else:
                break
            
        seriesNo = "".join(digs)
        fileN1 = join(pathC1,f)
        print(fileN1)
        dataPlot1 = np.load(fileN1)
        
        fig1,ax1 = myPlotImg(dataPlot1,title=(namesOutput1[0]+"-"+seriesNo+" "+"myfitting"),minSc=[defMinInt],maxSc=[defMaxInt])
        if(save_flag):
            pdf.savefig(fig1)
                      
        #find corresponding filename on different folder
        if(isdir(pathC2)):
            onlyfile2 = [f2 for f2 in listdir(pathC2) if((isfile(join(pathC2, f2)) and f2.endswith(seriesNo+ext2)))]
        else:
            onlyfile2 = []            
        if(onlyfile2):
            headerInfo2 = dicom.read_file(join(pathC2,onlyfile2[0]))
            dataPlot2 = headerInfo2.pixel_array
            
            fig2,ax2 = myPlotImg(dataPlot2,title=(namesOutput1[0]+"-"+seriesNo+" "+"MRmaps"),minSc=[defMinInt],maxSc=[defMaxInt])
            if(save_flag):
                pdf.savefig(fig2)

        fig3 = pyplot.figure()
        flagNZ1 = np.logical_and(dataPlot1>defMinInt, dataPlot1<defMaxInt)
        flagNZ2 = np.logical_and(dataPlot2>defMinInt, dataPlot2<defMaxInt)
        indexNZ = np.where( np.logical_and(flagNZ1,flagNZ2) )
        pyplot.scatter(dataPlot1[indexNZ[0],indexNZ[1]],dataPlot2[indexNZ[0],indexNZ[1]],marker='o')
        x = dataPlot1[indexNZ[0],indexNZ[1]]
        y = dataPlot2[indexNZ[0],indexNZ[1]]
        r = np.corrcoef(x,y)
        pyplot.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),color="red")
        pyplot.xlabel('myfitting')
        pyplot.ylabel('MRmaps')
        pyplot.text(1400,2000,'r='+"{:10.2f}".format(r[0,1]), bbox=dict(facecolor='red', alpha=0.5) )
        #pyplot.show(block=False)
        if(save_flag):
            pdf.savefig(fig3)
        
        #check the points that are of diagonal    
        
        flagOutL1 = np.less_equal(dataPlot1,dataPlot2/2.0)
        indexOutL1 = np.where(np.logical_and(np.logical_and(flagNZ1,flagNZ2),flagOutL1) )
        pyplot.scatter(dataPlot1[indexOutL1[0],indexOutL1[1]],dataPlot2[indexOutL1[0],indexOutL1[1]],marker='+',color="green")

        flagOutL2 = np.less_equal(dataPlot2,dataPlot1/2.0)
        indexOutL2 = np.where(np.logical_and(np.logical_and(flagNZ1,flagNZ2),flagOutL2) )
        pyplot.scatter(dataPlot1[indexOutL2[0],indexOutL2[1]],dataPlot2[indexOutL2[0],indexOutL2[1]],marker='+',color="magenta")
        
        pyplot.show(block=False)

        dataSh = dataPlot1.shape
        pyplot.figure(fig1.number)
        ax1.plot(indexOutL2[1],dataSh[0]-1-indexOutL2[0],'m+')#marker='+',color="magenta")
        ax1.plot(indexOutL1[1],dataSh[0]-1-indexOutL1[0],'g+')#marker='+',color="green")
        pyplot.show(block=False)

        pyplot.figure(fig2.number)
        ax2.scatter(indexOutL2[1],dataSh[0]-1-indexOutL2[0],marker='+',color="magenta")
        ax2.scatter(indexOutL1[1],dataSh[0]-1-indexOutL1[0],marker='+',color="green")
        pyplot.show(block=False)
        
        
    if(save_flag):
        pdf.close()
   
