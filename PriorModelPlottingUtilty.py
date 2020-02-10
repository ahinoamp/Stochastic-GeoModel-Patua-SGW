import itertools
import k3d
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import meshio
import numpy as np
import pandas as pd
import pynoddy
from scipy.interpolate import griddata
import scipy.interpolate as interp
import sys, os
import subprocess
from shutil import copyfile
import vtk
from vtkplotter import *
import pynoddy.history
import pynoddy.events
import time
from scipy.interpolate import interpn
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.lines import Line2D
import warnings

def index_containing_substring(the_list, substring):
    for i, s in enumerate(the_list):
        if substring in s:
              return i
    return -1

def get_wellbore_voxels_from_paths2(LithBlock, xi,yi,zi, xlim, ylim, zlim, delV):

    x = np.linspace(xlim[0], xlim[1], np.shape(LithBlock)[0])
    y = np.linspace(ylim[0], ylim[1], np.shape(LithBlock)[1])
    z = np.linspace(zlim[0], zlim[1], np.shape(LithBlock)[2])

    Vi = interpn((x,y,z), LithBlock, np.array([xi,yi,zi]).T, method='nearest')

    return Vi
    
def importModulesAndSetup():
    import pynoddy
    
    ##############
    # Some file fixes
    ##############
    basepynoddyfile = pynoddy.__file__[:-11]+'experiment/__init__.py'
    # Read in the file
    with open(basepynoddyfile, 'r') as file :
      filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('from . import util.sampling as Sample', 'from .util import sampling as Sample')

    # Write the file out again
    with open(basepynoddyfile, 'w') as file:
      file.write(filedata)

    target = pynoddy.__file__[:-11]+'output.py'

    source = 'output.py'

    copyfile(source, target)

    target = pynoddy.__file__[:-11]+'events.py'

    source = 'events.py'

    copyfile(source, target)

    ##############
    # Changing exection permissions
    ##############
    folder = os.getcwd()
    noddyEXE = folder+'/noddy.exe'
    strV = 'chmod 777 '+noddyEXE
    os.system(strV)
    import pynoddy.experiment
    
def CalculateModel(modelfile, output_name, outputoption = 'ALL', cubesize = 250):
    folder = os.getcwd()
    noddyEXE = folder+'/noddy.exe'
    H1 = pynoddy.history.NoddyHistory(modelfile)
    H1.change_cube_size(cubesize)
    H1.write_history(modelfile)
    subprocess.Popen([noddyEXE, modelfile, output_name, outputoption], 
                           shell=False, stderr=subprocess.PIPE, 
                           stdout=subprocess.PIPE).stdout.read()
    print('Finished calculating model')
    
def getDXF_parsed_structure(output_name):
    filename = output_name + '.dxf'
#    doc = ezdxf.readfile(filename)
    cell_data = []
    xpoint = []
    ypoint = []
    zpoint = []
    with open(filename) as f:
        cntr=0
        faceCounter=0
        for line in f:
            if(cntr==(7+faceCounter*28)):
                cell_data.append(line)
                faceCounter=faceCounter+1
            elif(cntr==(9+(faceCounter-1)*28)):
                xpoint.append(float(line))
            elif(cntr==(11+(faceCounter-1)*28)):
                ypoint.append(float(line))
            elif(cntr==(13+(faceCounter-1)*28)):
                zpoint.append(float(line))

            elif(cntr==(15+(faceCounter-1)*28)):
                xpoint.append(float(line))
            elif(cntr==(17+(faceCounter-1)*28)):
                ypoint.append(float(line))
            elif(cntr==(19+(faceCounter-1)*28)):
                zpoint.append(float(line))

            elif(cntr==(21+(faceCounter-1)*28)):
                xpoint.append(float(line))
            elif(cntr==(23+(faceCounter-1)*28)):
                ypoint.append(float(line))
            elif(cntr==(25+(faceCounter-1)*28)):
                zpoint.append(float(line))

            cntr=cntr+1

    points = np.column_stack((np.asarray(xpoint, dtype=float),
                             np.asarray(ypoint, dtype=float),
                             np.asarray(zpoint, dtype=float)))
    cell_data.pop()
    cell_data = np.asarray(cell_data, dtype=object)
#    print('Finished reading model')
   
    return points, cell_data, faceCounter

def convertSurfaces2VTK(points, cell_data, faceCounter, outputOption = 1, fileprefix='Surface',  xy_origin=[0,0,0]):
    
    # Choose output option
    num3Dfaces=faceCounter
    print('The number of triangle elements (cells/faces) is: ' + str(num3Dfaces))


    #apply origin transformation
    points[:, 0] = points[:, 0]+xy_origin[0]
    points[:, 1] = points[:, 1]+xy_origin[1]
    points[:, 2] = points[:, 2]+xy_origin[2]
    
    cell_data = pd.Series(cell_data.reshape((-1, )))

    CatCodes = np.zeros((len(cell_data),))
    filterB = (cell_data.str.contains('B')) 
    filterS = (cell_data.str.contains('S')) 

    CatCodes[filterB]= cell_data.loc[filterB].str[:-20].astype('category').cat.codes
    CatCodes[filterS]= -1*(cell_data.loc[filterS].str[:-12].astype('category').cat.codes+1)

#    print(np.unique(CatCodes[filterB]))
#    print(np.unique(CatCodes[filterS]))

    for i in range(1, len(CatCodes)):
        if(CatCodes[i]==0):
            CatCodes[i]=CatCodes[i-1]
            if(CatCodes[i-1]==0):
                CatCodes[i]=CatCodes[np.nonzero(CatCodes)[0][0]]

#    CatCodes[filterB]= 0
#    CatCodes[filterS]= 1

    UniqueCodes = np.unique(CatCodes)
    nSurfaces = len(UniqueCodes)
#    print(np.unique(CatCodes[filterB]))
#    print(np.unique(CatCodes[filterS]))
#    print(UniqueCodes)

    if (outputOption==2): ## if you would like a single vtk file
        cells = np.zeros((num3Dfaces, 3),dtype ='int')
        i=0
        for f in range(num3Dfaces):
            cells[f,:]= [i, i+1, i+2]
            i=i+3
        meshio.write_points_cells(
            "Model.vtk",
            points,
            cells={'triangle':cells},
            cell_data= {'triangle': {'cat':CatCodes}}   
            )
    else: ## option 1: make a separate file for each surface
        for i in range(nSurfaces):
            filterPoints = CatCodes==UniqueCodes[i]
            nCells = np.sum(filterPoints)
            Cells_i = np.zeros((nCells, 3),dtype ='int')
            cntr = 0
            for j in range(nCells):
                Cells_i[j]=[cntr, cntr+1, cntr+2]
                cntr=cntr+3

            booleanFilter = np.repeat(filterPoints,3)
   
            meshio.write_points_cells(
                fileprefix+str(i)+".vtk",
                points[np.repeat(filterPoints,3), :],
                cells={'triangle':Cells_i}
                )
#            print('surface ' +str(i) + '  code:'+str(UniqueCodes[i]))

#    print('Finished converting dxf to vtk')
    
    
    return nSurfaces, points, CatCodes

def CalculatePlotStructure(H1, output_name, plot, includeGravityCalc=0, cubesize = 250,  xy_origin=[317883,4379646, 1200-4000], plotwells =1):
    
    newmodelfile ='temp_model'
    H1.write_history(newmodelfile)

    #Alter the mesh size if desiring to speed up the process. Recommended size is 100
#    output_name = 'temp_noddy_out'
    
    
    #output options
    #1. BLOCK       
    #GEOPHYSICS   
    #SURFACES
    #BLOCK_GEOPHYS
    #BLOCK_SURFACES
    #TOPOLOGY
    #ANOM_FROM_BLOCK
    #ALL 
    if(includeGravityCalc==0):
        outputoption = 'BLOCK_SURFACES'
    else:
        outputoption = 'ALL'

    start = time.time()
    CalculateModel(newmodelfile, output_name, outputoption, cubesize)
    end = time.time()
    print('Calculation time took '+str(end - start) + ' seconds')

    ## Now need to change the DXF file (mesh format) to VTK. This is slow unfortunately
    start = time.time()
    points, cell_data, faceCounter = getDXF_parsed_structure(output_name)
    end = time.time()
    print('Parsing time took '+str(end - start) + ' seconds')

    ## Make a vtk file for each surface (option 1) or make a single vtk file for all surfaces (option 2)
    outputOption = 1
    fileprefix = 'TempSurface'

    start = time.time()
    nSurfaces, points, CatCodes = convertSurfaces2VTK(points, cell_data, faceCounter, outputOption, fileprefix,  xy_origin=xy_origin)
    
    end = time.time()
    print('Convert 2 VTK time took '+str(end - start) + ' seconds')

    N1 = pynoddy.output.NoddyOutput(output_name)

    lithology = N1.block

    [maxX, maxY, maxZ] = np.max(points, axis=0)
    [minX, minY, minZ] = np.min(points, axis=0)
    minZ = xy_origin[2]
    x = np.linspace(minX, maxX, N1.nx, dtype=np.float32)
    y = np.linspace(minY, maxY, N1.ny, dtype=np.float32)
    z = np.linspace(xy_origin[2], maxZ, N1.nz, dtype=np.float32)
#    z = np.linspace(0, 4000, N1.nz, dtype=np.float32)

    delx = x[1]-x[0]
    dely = y[1]-y[0]
    delz = z[1]-z[0]
    
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

    CoordXYZ = np.concatenate((xx.reshape(-1,1),yy.reshape(-1,1),zz.reshape(-1,1)), axis=1)

    Lithology = griddata(CoordXYZ, np.transpose(lithology, axes =(2, 1, 0)).reshape(-1,), (xx, yy, zz), method='nearest')
    
    vol = Volume(Lithology, c='jet', spacing=[delx, dely,delz], origin =[xy_origin[0], xy_origin[1], xy_origin[2]])
    lego = vol.legosurface(-1, np.max(Lithology)*2).opacity(0.15).c('jet')
#    vol = vol.color(['red', 'violet', 'green'])
    plot += lego
#    plot += vol

    colors = pl.cm.jet(np.linspace(0,1,nSurfaces))

    #outputOption = 1

    if(outputOption==1):
        for i in range(nSurfaces):
            filename = 'TempSurface'+str(i)+'.vtk'
            e=load(filename).c(colors[i, 0:3])

            plot += e
    else:
        filename = 'TempModel.vtk'
        e=load(filename)
        plot += e
    
    if(plotwells==1):
        Data = pd.read_csv('Data/WellNamesPaths.csv')

        filterV = (Data['X_m']>minX+0.5*delx) & (Data['Y_m']>minY+0.5*dely) &  (Data['Z_m']>minZ+0.5*delz) & (Data['X_m']<maxX-0.5*delx) & (Data['Y_m']<maxY-0.5*dely) & (Data['Z_m']<maxZ-0.5*delz)

        Data = Data[filterV]
        WellboreColors = get_wellbore_voxels_from_paths2(lithology, Data['X_m'], Data['Y_m'], Data['Z_m'], [minX, maxX], [minY, maxY], [minZ, maxZ], [delx, dely, delz])


        cm = plt.get_cmap('jet')
        maxL_V = np.max(WellboreColors)
        minL_V = np.min(WellboreColors)
        normL = (WellboreColors-minL_V)/(maxL_V-minL_V)
        RGBA = cm(normL)

        pointsPlot = Points(Data[['X_m', 'Y_m', 'Z_m']].values, r=4, c=RGBA)

        # color vertices based on their scalar value with any matplotlib color map
    #    points.pointColors(WellboreColors, cmap='coolwarm')

        plot += pointsPlot

    return points

def sampleSinglePriorModelParameters(ModelParametersTable, baseInputFile, OutputfileName):
    H1 = pynoddy.history.NoddyHistory(baseInputFile)
    
#    print('length is ' +str(len(ModelParametersTable)))
    nParameters2Set = int(len(ModelParametersTable)/3)
    Events = ModelParametersTable['Event Number']
    Properties = ModelParametersTable['Property']

    #    Events = ModelParametersTable['Event Name']
#    Events.remove('Parameter')
#    Events.remove('Property')
#    nEvents = len(Events)
    
#    Properties= pd.unique(ModelParametersTable['Property'].values).tolist()
#    nProperties = len(Properties)
    
    rowN=0
    for i in range(nParameters2Set):
        event = Events[rowN]
        prop = Properties[rowN]
        priorDistType = ModelParametersTable.loc[rowN, 'Parameter Value']
        priorparam1 = float(ModelParametersTable.loc[rowN+1, 'Parameter Value'])
        priorparam2 = float(ModelParametersTable.loc[rowN+2, 'Parameter Value'])
        if(i==1):
            layername = ModelParametersTable.loc[rowN, 'Event Name']
            indexLayer = index_containing_substring(H1.events[1].layer_names, layername)
            mu, sigma = priorparam1, priorparam2 # mean and standard deviation
            paramVal = np.random.normal(mu, sigma, 1)
            H1.events[1].layers[indexLayer].properties[prop]= paramVal[0]
        else:
            # Take care of uniform
            if(priorDistType=='Uniform'):
                paramVal = np.random.uniform(priorparam1,priorparam2,1)
                if(prop=='XAxis'):
                    RadiusLength=paramVal
                elif(prop=='YAxis'):
                    paramVal = RadiusLength
                elif(prop=='ZAxis'):
                    paramVal = RadiusLength
                H1.events[event].properties[prop] = paramVal[0]
        rowN = rowN +3
    
    nFaults = len(H1.events)-3
    arr = np.arange(3, nFaults+3)
    np.random.shuffle(arr)
    dictV = {}
    for i in range(nFaults):
        dictV[i+3]=arr[i]
    H1.reorder_events(dictV)
    H1.write_history(OutputfileName)
    return H1

def sampleMultiplePriorModelParameters(ModelParametersTable, baseInputFile, nrealizations, folder):
    if not os.path.exists(folder):
        os.mkdir(folder)
    
    H1 = pynoddy.history.NoddyHistory(baseInputFile)
    
    for r in range(nrealizations):
        OutputfileName = folder+'/realization_'+str(r)
        sampleSinglePriorModelParameters(ModelParametersTable, baseInputFile, OutputfileName)

def calcMultipleHistoryFiles(nrealizations, folder, cubesize = 250, outputoption = 'BLOCK_SURFACES', xy_origin=[317883,4379646, 1200-4000]):

    for r in range(nrealizations):
        inputfile = folder+'/realization_'+str(r)
        output_name = inputfile+'_output'
        CalculateModel(inputfile, output_name, outputoption, cubesize)

        if((outputoption=='BLOCK')|(outputoption=='BLOCK_GEOPHYS')):
            pass
        else:
            ## Now need to change the DXF file (mesh format) to VTK. This is slow unfortunately
            points, cell_data, faceCounter = getDXF_parsed_structure(output_name)

            ## Make a vtk file for each surface (option 1) or make a single vtk file for all surfaces (option 2)
            outputOption = 1
            nSurfaces, points, CatCodes = convertSurfaces2VTK(points, cell_data, faceCounter, outputOption, fileprefix=inputfile+'_Surface', xy_origin=xy_origin)
            data = np.concatenate((points,np.repeat(CatCodes,3).reshape(-1,1)), axis=1)
            np.savetxt(inputfile+'pointdata.csv', data, delimiter=',')        
        os.remove(output_name+".dxf")

def PlotScenarios(ScenarioN, plot, plotwells =0):
    
    if(ScenarioN=='Choose a Scenario'):
        print('Please choose a scenario in the dropdown menu above')
        return
    
    xy_origin=[317883,4379646, 1200-4000]
    folder = 'Data/'+ScenarioN+'/'
    output_name = folder+'noddy_out'
    includeGravityCalc = 0

    SurfacesScenario = {'Scenario1_MedResolution': 9, 'Scenario1_HighResolution':9,'Scenario2_MedResolution': 29,'Scenario2_HighResolution': 29, 'Scenario3_MedResolution': 12, 'Scenario3_HighResolution':12}
    nSurfaces = SurfacesScenario[ScenarioN]

    CubeSizes = {'Scenario1_MedResolution': 100, 'Scenario1_HighResolution':50,'Scenario2_MedResolution': 100,'Scenario2_HighResolution': 50, 'Scenario3_MedResolution': 100, 'Scenario3_HighResolution':50}
    cubesize = CubeSizes[ScenarioN]
    cubesize = 250
    modelfile = folder+ScenarioN+'.his'
    H1 = pynoddy.history.NoddyHistory(modelfile)
    outputoption = 'BLOCK_SURFACES'

    start = time.time()
    CalculateModel(modelfile, output_name, outputoption, cubesize)
    end = time.time()
    print('Calculation time took '+str(end - start) + ' seconds')

    ## Now need to change the DXF file (mesh format) to VTK. This is slow unfortunately
    start = time.time()
    points, cell_data, faceCounter = getDXF_parsed_structure(output_name)
    end = time.time()
    print('Parsing time took '+str(end - start) + ' seconds')

    ## Make a vtk file for each surface (option 1) or make a single vtk file for all surfaces (option 2)
    outputOption = 1
    fileprefix = folder+'TempSurface'

    start = time.time()
    nSurfaces, points, CatCodes = convertSurfaces2VTK(points, cell_data, faceCounter, outputOption, fileprefix,  xy_origin=xy_origin)
    
    end = time.time()
    print('Convert 2 VTK time took '+str(end - start) + ' seconds')

    
    N1 = pynoddy.output.NoddyOutput(folder+'noddy_out')
#    points = np.loadtxt(folder+'Points.csv', delimiter=',')
    lithology = N1.block

    [maxX, maxY, maxZ] = np.max(points, axis=0)
    [minX, minY, minZ] = np.min(points, axis=0)
    minZ = xy_origin[2]
    x = np.linspace(minX, maxX, N1.nx, dtype=np.float32)
    y = np.linspace(minY, maxY, N1.ny, dtype=np.float32)
    z = np.linspace(xy_origin[2], maxZ, N1.nz, dtype=np.float32)
#    z = np.linspace(0, 4000, N1.nz, dtype=np.float32)

    delx = x[1]-x[0]
    dely = y[1]-y[0]
    delz = z[1]-z[0]
    
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

    CoordXYZ = np.concatenate((xx.reshape(-1,1),yy.reshape(-1,1),zz.reshape(-1,1)), axis=1)

    Lithology = griddata(CoordXYZ, np.transpose(lithology, axes =(2, 1, 0)).reshape(-1,), (xx, yy, zz), method='nearest')
    
    vol = Volume(Lithology, c='jet', spacing=[delx, dely,delz], origin =[xy_origin[0], xy_origin[1], xy_origin[2]])
    lego = vol.legosurface(-1, np.max(Lithology)*2).opacity(0.15).c('jet')
#    vol = vol.color(['red', 'violet', 'green'])
    plot += lego
#    plot += vol

    colors = pl.cm.terrain(np.linspace(0,1,nSurfaces))
    colors = np.flipud(colors)
    #outputOption = 1
    for i in range(nSurfaces):
        filename = folder+'TempSurface'+str(i)+'.vtk'
        e=load(filename).c(colors[i, 0:3])
        plot += e
    
    if(plotwells==1):
        Data = pd.read_csv('Data/WellNamesPaths.csv')

        filterV = (Data['X_m']>minX+0.5*delx) & (Data['Y_m']>minY+0.5*dely) &  (Data['Z_m']>minZ+0.5*delz) & (Data['X_m']<maxX-0.5*delx) & (Data['Y_m']<maxY-0.5*dely) & (Data['Z_m']<maxZ-0.5*delz)

        Data = Data[filterV]
        WellboreColors = get_wellbore_voxels_from_paths2(lithology, Data['X_m'], Data['Y_m'], Data['Z_m'], [minX, maxX], [minY, maxY], [minZ, maxZ], [delx, dely, delz])


        cm = plt.get_cmap('jet')
        maxL_V = np.max(WellboreColors)
        minL_V = np.min(WellboreColors)
        normL = (WellboreColors-minL_V)/(maxL_V-minL_V)
        RGBA = cm(normL)

        pointsPlot = Points(Data[['X_m', 'Y_m', 'Z_m']].values, r=4, c=RGBA)

        # color vertices based on their scalar value with any matplotlib color map
    #    points.pointColors(WellboreColors, cmap='coolwarm')

        plot += pointsPlot

    return points

def GenerateRealizations(RealizationName, nrealizations, cubesizetxt = '150(m) - LowRes', outputoption = 'ALL', foldersuffix=''):
    Dict = {'50(m) - HighRes':50, '100(m) - MedRes':100, '150(m) - LowRes':150}
    cubesize = Dict[cubesizetxt]
    TableFile = 'Data/PriorUncertaintyTables/PriorUncertaintyTable_'+ RealizationName+'.csv'
    PriorStructuralModelTable = pd.read_csv(TableFile)
    baseInputFile = 'Data/RealizationsData/'+RealizationName+'.his'
    folder = 'Data/RealizationsData/'+RealizationName+foldersuffix

    sampleMultiplePriorModelParameters(PriorStructuralModelTable, baseInputFile, nrealizations, folder)
    calcMultipleHistoryFiles(nrealizations, folder, outputoption = outputoption, cubesize = cubesize)

def plotRealizations(RealizationName, foldersuffix, nRealizations=4):
    warnings.filterwarnings('ignore')
    xy_origin=[317883,4379646, 1200-4000]
    #xy_origin=[0,0, 0]
    folder = 'Data/RealizationsData/'+RealizationName+foldersuffix

    foldersuffix = ''
    nRows=4
    cubesize =50
    nRealizations = 4
    RealizationNumbers= [0, 1,2,3]
    fig, axs = plt.subplots(nRows,nRealizations, figsize=(nRealizations*5,nRows*3), gridspec_kw={'height_ratios': [3.5,2.5,2.5,3.5]})
    #levels = np.arange(115, 120, 0.5)
    for r in range(nRealizations):
    #    ax = fig.add_subplot(1,nRealizations,r+1, projection='3d')
        print(r)
        inputfile = folder+'/realization_'+str(RealizationNumbers[r])
        data = np.loadtxt(inputfile+'pointdata.csv', delimiter=',')
        CatCodesCorrect = data[:,3]
        points = data[:,0:3]
        # Data for three-dimensional scattered points
        Surfaces = np.unique(CatCodesCorrect)
        nSurfaces = len(Surfaces)
        cmapNames = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                    'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                    'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

    #    for i in range(nSurfaces):
    #        filterS= CatCodesCorrect==Surfaces[i]
    #        ax.scatter3D(points[filterS,0], points[filterS,1], points[filterS,2], c=points[filterS,1], cmap=cmapNames[i]);
    #        ax.plot_trisurf(points[filterS,0], points[filterS,1], points[filterS,2],cmap=cmapNames[i]);
    #    ax.set_xlabel('x')
    #    ax.set_ylabel('y')
    #    ax.set_zlabel('z')
    #    ax.view_init(azim=90, elev=10)
        fileGeo = inputfile+'_output'
        geophys = pynoddy.output.NoddyGeophysics(fileGeo)

    #    if(r==0):
        [maxX, maxY, maxZ] = np.max(points, axis=0)
        [minX, minY, minZ] = np.min(points, axis=0)
        x = np.linspace(minX, maxX, np.shape(geophys.grv_data)[1], dtype=np.float32)
        y = np.linspace(maxY, minY, np.shape(geophys.grv_data)[0], dtype=np.float32)
        xx, yy = np.meshgrid(x, y, indexing='xy')

        ax = axs[0,r]
    #    ax.imshow(geophys.grv_data, cmap = 'jet')
        cf = ax.contourf(xx, yy, geophys.grv_data, cmap = 'jet', zorder=0)
    #    cbar = plt.colorbar(cf, orientation = 'horizontal')
        cbar = plt.colorbar(cf, orientation = 'vertical', ax=ax, panchor=[0.6, 0.5], anchor=[-0.1, 0.5], pad=0)
        filterTop = (points[:,2]>(maxZ-float(cubesize)*1.2)) & (CatCodesCorrect>0)
        CatCodesCorrectFiltered = CatCodesCorrect[filterTop]
        uniqueFaults = np.unique(CatCodesCorrectFiltered)
        nUniqueFaults = len(uniqueFaults)
        for f in range(nUniqueFaults):
            filterTopSurf = filterTop &(CatCodesCorrect==uniqueFaults[f])
            X = points[filterTopSurf,0]
            Y = points[filterTopSurf,1]
            if(len(X)>300):
                fittingP = 10
            elif(len(X)>200):
                fittingP = 6
            elif(len(X)>100):
                fittingP = 4
            else:            
                fittingP = 3
            if((np.max(X)-np.min(X))>(np.max(Y)-np.min(Y))):
                z = np.polyfit(X, Y, fittingP)
                p = np.poly1d(z)
                X=np.sort(X)
                newY = p(X)
                ax.plot(X,newY, zorder=100, c='k')
            else:
                z = np.polyfit(Y, X, fittingP)
                p = np.poly1d(z)
                Y=np.sort(Y)
                newX = p(Y)
                ax.plot(newX,Y, zorder=100, c='k')

    #    ax.scatter(points[filterTop,0], points[filterTop,1], zorder=100, s=2, c='k')  

        ax.set_title('Realization ' + str(r+1))
        ax.axis('square')
        ax.set_xlim([np.min(xx), np.max(xx)])
        ax.set_ylim([np.min(yy), np.max(yy)])
        if(r!=0):
            ax.get_yaxis().set_ticks([])
            ax.get_xaxis().set_ticks([])
        else:
            ax.get_xaxis().set_ticks([])
            ax.set_ylabel('North (m - UTM 11N)')
            ax.get_yaxis().set_ticks([])
        ax.set_xlabel('East (m)')
        ax.set_ylabel('North (9000 m extent)')
        ax.set_xlabel('East (9000 m extent)')
        #plt.gca().invert_yaxis()

        #rerun the geopart to make things pretty
        historyfile = inputfile
        cubesize=50
        tempoutputfile= inputfile+'_temp'
        H1 = pynoddy.history.NoddyHistory(historyfile)
        H1.change_cube_size(cubesize)
        H1.events[1].RefineLayers(H1.get_cube_size())

        H1.write_history(tempoutputfile)
        output_name = inputfile+'_outtemp'

        CalculateModel(tempoutputfile, output_name, outputoption = 'BLOCK', cubesize = cubesize)
        N1 = pynoddy.output.NoddyOutput(output_name)

        startReading = 0
        fh = open(folder+'/realization_'+str(r) +'_outtemp.g00')

        Layer = []
        Number = []
        for line in fh:
            if "NUM ROCK TYPES" in line:
                n_events = int(line.split("=")[1])
                startReading=1
            if(startReading==1):
                if "ROCK DEFINITION" in line:
                    info = line[15:]
                    List= info.split('=')
                    name=List[0]
                    Num = int(List[1].strip())
                    if('Intrusive' in name):
                        Layer.append('Intrusive')
                        Number.append(Num)
                    elif('Sedimentary' in name):
                        Layer.append('Sedimentary')
                        Number.append(Num)
                    elif('Volcanic Mafic' in name):
                        Layer.append('Volcanic Mafic')
                        Number.append(Num)
                    elif('Volcanic Felsic' in name):
                        Layer.append('Volcanic Felsic')
                        Number.append(Num)
        fh.close()

        Data = pd.DataFrame({'Layer': Layer, 'Number': Number})
        grouped = Data.groupby('Layer').agg(
                MaxNumer=pd.NamedAgg(column='Number', aggfunc='max'), 
                MinNumer=pd.NamedAgg(column='Number', aggfunc='min')).reset_index()

        grouped['ActualMin']=1000
        grouped['ActualMax']=-1000
        minBlock = np.min(N1.block)
        maxBlock = np.max(N1.block)
        for i in range(4):
            nameLithology = grouped.loc[i, 'Layer']
            minV = grouped.loc[i, 'MinNumer']
            maxV = grouped.loc[i, 'MaxNumer']
            ActualMin = np.max([minV, minBlock])
            ActualMax = np.min([maxV, maxBlock])
            grouped.loc[i, 'ActualMin'] = ActualMin
            grouped.loc[i, 'ActualMax'] = ActualMax

        grouped['Number_colors'] = grouped['ActualMax']-grouped['ActualMin']


        cmapNames = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                        'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

        def getcolormap(name):
            color='k'
            if(name=='Intrusive'):
                color='Reds'
            elif(name=='Volcanic Felsic'):
                color='Purples'
            elif(name=='Volcanic Mafic'):
                color='binary'
            elif(name=='Sedimentary'):
                color='YlOrBr'
            return color

        grouped = grouped.sort_values(['MinNumer'], ascending=True).reset_index()

        for i in range(4):
            nameLithology = grouped.loc[i, 'Layer']
            colormapname = getcolormap(nameLithology)
            colormap = matplotlib.cm.get_cmap(colormapname, 100)
            newcolors = colormap(np.linspace(0.4, 0.8, int(grouped.loc[i, 'Number_colors'])))
            if(i==0):
                colormaplist = newcolors
            else:
                colormaplist = np.vstack((colormaplist, newcolors))
    #    colormaplist = np.flipud(colormaplist)
        newcmp = ListedColormap(colormaplist)

        ax = axs[1,r]
        N1.plot_section_faults('y', ax = ax, cmap=newcmp, xy_origin=xy_origin)
        yPlot = (N1.ymax+N1.ymin)/2 +xy_origin[1]
        filterTop = (points[:,1]>(yPlot-float(cubesize)*1.2)) & (points[:,1]<(yPlot+float(cubesize)*1.2)) &(CatCodesCorrect>0)
    #    ax.scatter(points[filterTop,0], points[filterTop,2], zorder=100, s=2, c='k')
        CatCodesCorrectFiltered = CatCodesCorrect[filterTop]
        uniqueFaults = np.unique(CatCodesCorrectFiltered)
        nUniqueFaults = len(uniqueFaults)
        for f in range(nUniqueFaults):
            filterTopSurf = filterTop &(CatCodesCorrect==uniqueFaults[f])
            X = points[filterTopSurf,0]
            Y = points[filterTopSurf,2]
            if(len(X)>400):
                fittingP = 4
            elif(len(X)>200):
                fittingP = 3
            elif(len(X)>100):
                fittingP = 3
            else:            
                fittingP = 3
            if((np.max(X)-np.min(X))>(np.max(Y)-np.min(Y))):
                z = np.polyfit(X, Y, fittingP)
                p = np.poly1d(z)
                X=np.sort(X)
                newY = p(X)
                ax.plot(X,newY, zorder=100, c='k')
            else:
                z = np.polyfit(Y, X, fittingP)
                p = np.poly1d(z)
                Y=np.sort(Y)
                newX = p(Y)
                ax.plot(newX,Y, zorder=100, c='k')

        if(r==0):
    #        ax.set_title('Constant Y Section (Model Center)', loc='left')
            ax.set_title('')
        else:
            ax.set_title('')

        xextent = N1.xmax-N1.xmin
        ax.set_xlabel('East: {:.0f} m extent'.format(xextent))
        ax.get_xaxis().set_ticks([])
        if(r!=0):
            ax.get_yaxis().set_ticks([])
            ax.set_ylabel('')
        else:
            ax.set_ylabel('Depth (m)')

        ax = axs[2,r] 

        N1.plot_section_faults('x', ax = ax, cmap=newcmp, xy_origin=xy_origin)
        xPlot = (N1.xmax+N1.xmin)/2 +xy_origin[0]
        filterTop = (points[:,0]>(xPlot-float(cubesize)*1.2)) & (points[:,0]<(xPlot+float(cubesize)*1.2)) &(CatCodesCorrect>0)
        CatCodesCorrectFiltered = CatCodesCorrect[filterTop]
        uniqueFaults = np.unique(CatCodesCorrectFiltered)
        nUniqueFaults = len(uniqueFaults)
        for f in range(nUniqueFaults):
            filterTopSurf = filterTop &(CatCodesCorrect==uniqueFaults[f])
            X = points[filterTopSurf,1]
            Y = points[filterTopSurf,2]
            if(len(X)>300):
                fittingP = 4
            elif(len(X)>200):
                fittingP = 3
            elif(len(X)>100):
                fittingP = 3
            else:            
                fittingP = 3
            if((np.max(X)-np.min(X))>(np.max(Y)-np.min(Y))):
                z = np.polyfit(X, Y, fittingP)
                p = np.poly1d(z)
                X=np.sort(X)
                newY = p(X)
                ax.plot(X,newY, zorder=100, c='k')
            else:
                z = np.polyfit(Y, X, fittingP)
                p = np.poly1d(z)
                Y=np.sort(Y)
                newX = p(Y)
                ax.plot(newX,Y, zorder=100, c='k')

        def getcolor(name):
            color='k'
            if(name=='Granitic'):
                color='r'
            elif(name=='Volcanic - Felsic'):
                color=[129/255,126/255,185/255]
            elif(name=='Volcanic - Mafic'):
                color=[112/255,112/255,112/255]
            elif(name=='Sedimentary'):
                color=[244/255, 130/255, 29/255]
            return color

        if(r==(nRealizations-1)):
            listLithosNames = ['Granitic', 'Volcanic - Felsic', 'Volcanic - Mafic', 'Sedimentary']    
            custom_lines = [Line2D([0], [0], color=getcolor(i), lw=4) for i in listLithosNames]   
            ax.legend(custom_lines, listLithosNames, bbox_to_anchor=(1.04, 1), loc='upper left', borderaxespad=0., fontsize=14)

        if(r==0):
    #        ax.set_title('Constant X Section (Model Center)', loc='left')
            ax.set_title('')
        else:
            ax.set_title('')

        yextent = N1.ymax-N1.ymin
        ax.set_xlabel('North: {:.0f} m extent'.format(yextent))
        ax.get_xaxis().set_ticks([])
        if(r!=0):
            ax.get_yaxis().set_ticks([])
            ax.set_ylabel('')
        else:
            ax.set_ylabel('Depth (m)')

        ax = axs[3,r]
        N1.plot_section_faults('z', ax = ax, cmap=newcmp, xy_origin=xy_origin)
        zPlot = (N1.zmax+N1.zmin)/2 +xy_origin[2]
        filterTop = (points[:,2]>(zPlot-float(cubesize)*1.2)) & (points[:,2]<(zPlot+float(cubesize)*1.2)) &(CatCodesCorrect>0)
    #    ax.scatter(points[filterTop,0], points[filterTop,1], zorder=100, s=2, c='k')
        CatCodesCorrectFiltered = CatCodesCorrect[filterTop]
        uniqueFaults = np.unique(CatCodesCorrectFiltered)
        nUniqueFaults = len(uniqueFaults)
        for f in range(nUniqueFaults):
            filterTopSurf = filterTop &(CatCodesCorrect==uniqueFaults[f])
            X = points[filterTopSurf,0]
            Y = points[filterTopSurf,1]
            if(len(X)>300):
                fittingP = 10
            elif(len(X)>200):
                fittingP = 6
            elif(len(X)>100):
                fittingP = 4
            else:            
                fittingP = 3
            if((np.max(X)-np.min(X))>(np.max(Y)-np.min(Y))):
                z = np.polyfit(X, Y, fittingP)
                p = np.poly1d(z)
                X=np.sort(X)
                newY = p(X)
                ax.plot(X,newY, zorder=100, c='k')
            else:
                z = np.polyfit(Y, X, fittingP)
                p = np.poly1d(z)
                Y=np.sort(Y)
                newX = p(Y)
                ax.plot(newX,Y, zorder=100, c='k')

        if(r==0):
    #        ax.set_title('Constant Z Section (Model Center)', loc='left')
            ax.set_title('')
        else:
            ax.set_title('')
        ax.get_yaxis().set_ticks([])
        ax.get_xaxis().set_ticks([])
        ax.set_ylabel('North (9000 m extent)')
        ax.set_xlabel('East (9000 m extent)')

    plt.subplots_adjust(wspace = 0.1)
    plt.subplots_adjust(hspace = 0.13)
    fig.savefig(RealizationName+'.png', dpi = 300, bbox_inches="tight")