import dicom
import os
from os import listdir
from os.path import isfile, isdir, join
import numpy as np
from matplotlib import pyplot, cm


mainDir = "C:\\XXX\\"
subjDirs = ["01","02","03","04","05","06","07"]
seriesNo = "28"
dirM = "PRE-T1"
ext = "IMA"
ext1 = ".npy"
outputMain = "C:\\XXX\\Results\\"
outputMainDir1 = join(outputMain,"MLeigen\\")
namesOutput1 = ["A","B","TI","dataT1star","dataXsquare","minIndexImg","statusImg","numIterImg"]
#["A","B","TI","dataT1star","dataXsquare","minIndexImg","statusImg","numIterImg"]

saveFlag = False

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

def onclick(event, figF, A, B, headerInfo, invTime, minIndexImg ):
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    
    N = 100
    x = np.linspace(0, max(invTime), N)
    shInfo = A.shape
    X = shInfo[0]-1-int(event.ydata)
    Y = int(event.xdata)
    y = A[X,Y]-B[X,Y]*np.exp(-x/T1[X,Y])
    nS = int(minIndexImg[X,Y])
    print(nS)
    pyplot.figure(figF.number)
    pyplot.clf()
    pyplot.plot(x, y)
    valPixs = [float(headerInfo[i].pixel_array[X,Y]) for i in indS]
    pyplot.scatter(invTime,valPixs,marker="x",color="red")
    if nS>0:
        newPixs = [-(valPixs[i]) for i in range(0,nS)]
        print(newPixs)
        newInvT = [invTime[i] for i in range(0,nS)]
        print(newInvT)
        pyplot.scatter(newInvT,newPixs,marker="x",color="magenta")
    pyplot.show(block=False)



subj = subjDirs[0]
pathC = join(mainDir,subj,dirM)
print(pathC)

#input files
onlyfiles = [f for f in listdir(pathC) if (isfile(join(pathC, f)) and f.endswith(ext))]

headerInfo = list()
invTime = list()
DicomDataS = list()
for cc, f in enumerate(onlyfiles):
    fileN = join(pathC,f)
    tmpDcm = dicom.read_file(fileN)
    if( str(tmpDcm.SeriesNumber) == seriesNo ):
        print(fileN)
        headerInfo.append(tmpDcm)
        invTime.append(float(tmpDcm.InversionTime))

indS = sorted(range(len(invTime)), key=lambda k: invTime[k])

#output files
pathC1 = join(outputMainDir1,subj,dirM)
print(pathC1)

extf1 = namesOutput1[0] + ext1
onlyfiles = [f for f in listdir(pathC1) if (isfile(join(pathC1, f)) and f.endswith(ext1))]

fig = list()
ax = list()
for cc, f in enumerate(onlyfiles):
    digs = list()
    for i,x in enumerate(f):
        if (f[i].isdigit()==True):
            digs.append(f[i])
        else:
            break
        
    sNo = "".join(digs)
    if sNo==seriesNo:
        fileN1 = join(pathC1,f)
        print(fileN1)
        if f.endswith("A"+ext1):
            A = np.load(fileN1)
            dataPlot = A
            tmpfig,tmpax = myPlotImg(dataPlot,title=(f+" "+"myfitting"),minSc=[defMinInt],maxSc=[500])
            fig.append(tmpfig)
            ax.append(tmpax)
        elif f.endswith("B"+ext1):
            B = np.load(fileN1)
            dataPlot = B
            tmpfig,tmpax = myPlotImg(dataPlot,title=(f+" "+"myfitting"),minSc=[defMinInt],maxSc=[500])
            fig.append(tmpfig)
            ax.append(tmpax)
        elif f.endswith("TI"+ext1):
            T1 = np.load(fileN1)
            dataPlot = T1
            tmpfig,tmpax = myPlotImg(dataPlot,title=(f+" "+"myfitting"),minSc=[defMinInt],maxSc=[defMaxInt])
            fig.append(tmpfig)
            ax.append(tmpax)
        elif f.endswith("dataT1star"+ext1):
            T1c = np.load(fileN1)
            dataPlot = T1c
            tmpfig,tmpax = myPlotImg(dataPlot,title=(f+" "+"myfitting"),minSc=[defMinInt],maxSc=[defMaxInt])
            fig.append(tmpfig)
            ax.append(tmpax)
        elif f.endswith("dataXsquare"+ext1):
            chiSq = np.load(fileN1)
            dataPlot = chiSq
            tmpfig,tmpax = myPlotImg(dataPlot,title=(f+" "+"myfitting"),minSc=[0],maxSc=[0.5])
            fig.append(tmpfig)
            ax.append(tmpax)
        elif f.endswith("minIndexImg"+ext1):
            minIter = np.load(fileN1)
            dataPlot = minIter
            tmpfig,tmpax = myPlotImg(dataPlot,title=(f+" "+"myfitting"))
            fig.append(tmpfig)
            ax.append(tmpax)
        elif f.endswith("numIterImg"+ext1):
            numLMiter = np.load(fileN1)
            dataPlot = numLMiter
            tmpfig,tmpax = myPlotImg(dataPlot,title=(f+" "+"myfitting"))
            fig.append(tmpfig)
            ax.append(tmpax)
        elif f.endswith("statusImg"+ext1):
            statusLM = np.load(fileN1)
            dataPlot = statusLM
            tmpfig,tmpax = myPlotImg(dataPlot,title=(f+" "+"myfitting"))
            fig.append(tmpfig)
            ax.append(tmpax)

#plot exponential fitting
figF = pyplot.figure()



cid = list()
for ff in fig:
    cid.append(ff.canvas.mpl_connect('button_press_event', lambda event: onclick(event, figF, A, B, headerInfo, invTime, minIter)) )

    
