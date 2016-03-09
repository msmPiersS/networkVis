######################################################################################
## R script for generating 2d interactive network visualisation
## 
## Process:
## 1) Loads square matrix from delimited file
## 2) Cleans up and calculates node positions (using igraph) based on different algorithms- visualsies locally in ggplot
## 3) Exports json files for use in sigmajs based interactive vis
##
## March 2018 - ps
#######################################################################################


#######################################################################################
## Set up

  #locations
  homeDir = '/Users/piers.stobbs/Documents/piers/Box Sync/datascience/networkVis/'
  dataDir = '/Users/piers.stobbs/Documents/piers/Box Sync/datascience/networkVis/data/'
  visDir = '/Users/piers.stobbs/Documents/piers/Box Sync/datascience/networkVis/2dvis/'
  setwd(homeDir)
  getwd()
  
  #files
  #nodesInputFile = 'sampleNodes.csv'
  #edgesInputFile = 'sampleEdges.csv'
  nodesInputFile = 'Mapping_List.txt'
  edgesInputFile = 'Mapping_Square.txt'
  
  #libraries
  library(data.table) #for efficient manipulation of data
  library(fpc) #to generate density based clusters
  library(igraph) #to generate the graph layouts
  library(ggplot2)  #to plot hte layours
  library(plyr) # for data manipulation
  library(jsonlite) #to create json export
  #library(RColorBrewer) #for colour generation
  
## End Set up
#######################################################################################


  
  
#######################################################################################
## Load data and process
  
  #load raw files
  rawNodes = fread(paste(dataDir, nodesInputFile, sep=""))
  rawEdges = fread(paste(dataDir, edgesInputFile, sep=""))
  
  #check edges file for text in column 1
  edgeLabels = ""
  edgeTypes = summary(rawEdges)
  if (grep("character", edgeTypes[2,1])==1) {
    print("Looks like first column of edges file is a string- assume it is labels and remove")
    edgeLables = rawEdges[, 1, with=FALSE]
    rawEdges = rawEdges[, -1, with=FALSE]
  }
  
  #quick check of dimensions
  dimEdges = dim(rawEdges)
  dimNodes = dim(rawNodes)
  if (dimEdges[1] != dimEdges[2]) {print("ERROR - input Edges file must be square matrix")}
  if (dimEdges[1] != dimNodes[1]) {print("ERROR - Nodes and Edges input files not of same dimension")}
  
  # assume we have raw counts coming in- need to normalise these
  # an interaction that occurs 80 out of 100 times should show as strongly as one that occurs 800 out of 1000 times
  # make sure diagonal is zero, remove any NANs and divide by row totals
  
  rawEdges = as.matrix(rawEdges) #convert to matrix
  diag(rawEdges) = 0 #zero out diagonal
  rawEdges[is.na(rawEdges)] = 0 #replace any nas with 0
  
  rawEdges = rawEdges/rowSums(rawEdges) #divide by rowsums
  rawEdges[rawEdges<1e-6] = 0 #clean up any tiny values
  
  #quick look at distribution of scores
  exploreDistFlag = 0
  if (exploreDistFlag==1) {
    hist(rawEdges)
    summary(rawEdges)
    sum(rawEdges>0.01)
    sum(rawEdges>0.05)
    sum(rawEdges>0.1)
  }
  

  # assume that node names are in column 1 of nodes input file and relative size is in column 2
  nodes = copy(rawNodes)
  edges = copy(rawEdges)
  setnames(nodes, colnames(nodes)[1:2], c("name", "relativeSize"))
  nodes[, name:=tolower(name)]
  colnames(edges) = as.matrix(tolower(nodes[, name]))
  setkey(nodes, name)

  #nodes = nodes[1:20,]
  #edges = edges[1:20, 1:20]    

## End Load data and process
#######################################################################################  
  
  
  
  
#######################################################################################
## Generate network layout and plot  
  
  # network graphs work best when there are not too many connections
  # so we will try filtering on the edges, removing any below a given threshold
  # we can iterate here to find a setting that looks ok, but lets start by making sure roughtly 75% of edges are zero
  
  edgesSort = sort(edges)
  thresh50 = edgesSort[round(0.5*length(edgesSort))]
  thresh75 = edgesSort[round(0.75*length(edgesSort))]
  thresh80 = edgesSort[round(0.80*length(edgesSort))]
  thresh90 = edgesSort[round(0.90*length(edgesSort))]
  
  edgesIn = copy(edges)
  edgesIn[edgesIn<thresh90] = 0 #change thresh90 to alternative to change filter setting
  
  # setup igraph - create graph from edges matrix and define what we are going to do
  # in  this instance we are going to use a "directed" which means the connection between A and B
  # can be different to the connection between B and A
  g = graph.adjacency(edgesIn, mode="directed", weighted="strength", diag=FALSE)
  
  # initiate initial location of the nodes- to ensure repeatability
  set.seed(1234)
  layoutStart = matrix(runif(2*nrow(nodes)), nrow(nodes), 2)
  
  # generate layout using Fruhterman-Reingolg layout algorithm
  tmpLayout = layout_with_fr(g, coords=layoutStart) 
  
  # setup table of nodes with coordinates
  nodesDT = as.data.table(get.data.frame(g, what="vertices"))
  nodesDT[, Xorg:=tmpLayout[, 1]]
  nodesDT[, Yorg:=tmpLayout[, 2]]
  
  #join in node sizes
  setkey(nodesDT, name)
  nodesDT = nodes[nodesDT]
  nodesDT[, relativeSizeOrg:=relativeSize]
  
  # cap outliers at 2 sigma from means positions to make sure we can view the layout successfully
  nodesDT[, X:=Xorg]
  nodesDT[, Y:=Yorg]
  stds = apply(nodesDT[, list(X, Y)], 2, sd)
  medians = apply(nodesDT[, list(X, Y)], 2, median)
  means = apply(nodesDT[, list(X, Y)], 2, mean)
  mins = apply(nodesDT[, list(X, Y)], 2, min)
  maxs = apply(nodesDT[, list(X, Y)], 2, max)
  
  capMax = means + 2*stds
  capMin = means - 2*stds
  
  nodesDT[nodesDT[, X] > capMax[1], X:= capMax[1]] 
  nodesDT[nodesDT[, Y] > capMax[2], Y:= capMax[2]] 
  nodesDT[nodesDT[, X] < capMin[1], X:= capMin[1]] 
  nodesDT[nodesDT[, Y] < capMin[2], Y:= capMin[2]] 
  
  #nodesDT[nodesDT[,Xorg]!=nodesDT[,X],]
  #nodesDT[nodesDT[,Yorg]!=nodesDT[,Y],]
  
  #clustering using dbscan
  # may need to play around with the inputs to get a good result
  # the first variable after the data is the range in which the algorithm 
  # looks for nodes to cluster togeter, the smaller the figure, the more the clusters
  # the second variable is the minimum number of points allowed in a cluster
  ds <- dbscan(nodesDT[, list(X, Y)], mean(stds)/5, MinPts = max(2,floor(nrow(nodesDT)/500)))
  
  useClusters = 1
  # if you dont want to use clusters you can simply group them based on size
  if (useClusters != 1) {
    nodesDT[, typeId:=floor(log(1+max(relativeSize))) - floor(log(1+relativeSize))]
    
  } else {
    nodesDT[, typeId:=ds$cluster+1] #make sure cluster ids are greater than zero
  }

  #create list of colours
  colourPalette = rainbow(max(nodesDT[,typeId]), s = 0.6, v = 0.75)
  nodesDT[, colourId:=colourPalette[typeId]]
  
  # pull edges and set to and from coordinates based on node positions
  edgesDT = as.data.table(get.data.frame(g, what="edges"))
  setkey(edgesDT, from)
  setkey(nodesDT, name)
  
  edgesDT = nodesDT[edgesDT][, list(to, from=name, strength, X1=X, Y1=Y)]
  setkey(edgesDT, to)
  edgesDT = nodesDT[edgesDT][, list(to = name, from, strength, X1, Y1, X2=X, Y2=Y)]
  setkey(edgesDT, from)
  
  # different options for how to scale the nodes- use log scale or not
  defaultSize = 10
  useLog = 0
  if (useLog==1) {
    nodesDT[, relativeSize:=defaultSize*log(1+nodesDT[, relativeSizeOrg])/log(1+max(nodesDT[, relativeSizeOrg]))]
  } else {
    nodesDT[, relativeSize:=defaultSize*(1+nodesDT[, relativeSizeOrg])/(1+max(nodesDT[, relativeSizeOrg]))]
  }
  
  #plot locally using ggplot
  localPlot = 0
  
  if (localPlot == 1) {
    
    pnet <- ggplot() + geom_segment(data=edgesDT, aes(x=X1, y=Y1, xend = X2, yend = Y2, size = strength), 
                                    colour="grey", alpha=0.5)
    
    pnet <- pnet  + 
      geom_point(data=nodesDT, aes(X, Y), color=nodesDT[, colourId], size=nodesDT[, relativeSize], alpha=0.8) +
      geom_text(data=nodesDT, aes(X, Y, label = name), size = 3) +
      #scale_colour_manual(guide=FALSE, name="",  values = c("node"="purple", "FALSE"="dark green"))+
      scale_colour_brewer(type = "seq", palette = 1) +
      scale_x_continuous(breaks = NULL) + 
      scale_y_continuous(breaks = NULL) +
      # discard default grid + titles in ggplot2 
      theme(panel.background = element_blank()) + 
      theme(legend.position="none") +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank()) +
      theme(legend.background = element_rect(colour = NA)) + 
      theme(panel.background = element_rect(fill = "white", 
                                            colour = NA)) + 
      theme(panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank())
    plot(pnet) 
  }
  
## Generate network layout and plot  
#######################################################################################
  
  
  
  
#######################################################################################
## export json
  
  #set up json export for js vis
  edgesExport = copy(edgesDT[, list(from, to, strength)])
  setnames(edgesExport, c("from", "to", "strength"), c("source", "target", "size"))
  edgesExport[, id:=1:nrow(edgesExport)]
  edgesExport[, label:=""]
  edgesExport[, color:="rgb(234,246,249)"]
  # round everything
  edgesExport[, size := round(size*1000)/1000]
  #make everything characters
  fieldClasses = laply(edgesExport, class)
  edgesExport[, (colnames(edgesExport)[fieldClasses!="character"]) := lapply(.SD, as.character), .SDcols = colnames(edgesExport)[fieldClasses!="character"]] 
  
  edgesExportJ = jsonlite::toJSON(edgesExport, pretty=FALSE)
  
  nodesExport = copy(nodesDT[, list(id = name, label = name, size = relativeSize, typeId, x = X, y = Y, color = colourId)])
  #setnames(nodesExport, c("name", "X", "Y"), c("label", "x", "y"))
  
  # turn colours into rgb for javascript
  for (i in 1:nrow(nodesExport)) {
    #i=1
    nodesExport[i, color:=paste("rgb(",paste(col2rgb(color), collapse=","),")", sep="")]
  }
  
  setkey(nodesExport, id) #make sure ordered by id
  #round everything
  nodesExport[, size := round(size*1000)/1000]
  nodesExport[, x := round(x*1000)/1000]
  nodesExport[, y := round(y*1000)/1000]
  
  #make everything characters
  fieldClasses = laply(nodesExport, class)
  nodesExport[, (colnames(nodesExport)[fieldClasses!="character"]) := lapply(.SD, as.character), .SDcols = colnames(nodesExport)[fieldClasses!="character"]] 
  
  
  #nodesExportJ = paste(ldply(seq(1,nrow(nodesExport),1), function(x) toJSON(nodesExport[x,])), sep=",")
  nodesExportJ = jsonlite::toJSON(nodesExport, pretty=FALSE)
  
  #export data json
  sink(paste(visDir,"data.json", sep=""))
  cat('{"edges":',edgesExportJ,',"nodes":',nodesExportJ,'}', sep="")
  sink()
  
  
  ## config file for javascript vis
  # note can play with node size and permanent lables using
  # maxNodeSize and labelThreshold
  
  configText = '
  {
  "type": "network",
  "version": "1.0",
  "data": "data.json",
  "logo": {
  "text": "",
  "file": "",
  "link": ""
  },
  "text": {
  "title": "visTitle",
  "more": "visMore",
  "intro": "visIntro"
  },
  "legend": {
  "edgeLabel": "edgeLabelText",
  "colorLabel": "Page Volume",
  "nodeLabel": "nodeLabelText"
  },
  "features": {
  "search": true,
  "groupSelectorAttribute": true,
  "hoverBehavior": "default"
  },
  "informationPanel": {
  "imageAttribute": false,
  "groupByEdgeDirection": false
  },
  "sigma": {
  "graphProperties": {
  "minEdgeSize": 0.5,
  "maxNodeSize": 15,
  "maxEdgeSize": 3,
  "minNodeSize": 1
  },
  "drawingProperties": {
  "labelThreshold": 10,
  "hoverFontStyle": "bold",
  "defaultEdgeType": "curve",
  "defaultLabelColor": "#fff",
  "defaultLabelHoverColor": "#fff",
  "defaultLabelSize": 14,
  "activeFontStyle": "bold",
  "fontStyle": "bold",
  "defaultHoverLabelBGColor": "#002147",
  "defaultLabelBGColor": "#ddd"
  },
  "mouseProperties": {
  "minRatio": 0.75,
  "maxRatio": 100
  }
  }
  }
  '
  
  
  
  #export config json
  visTitle = "Brand affinitiy visualisation"
  visMore = "Network visualisation based on proportion of insurance searchers who click on mulitple brands"
  visIntro = "Layout based on Force Directed Graph algorithm- attraction between brands proportional to volume of brand interactions"
  edgeLabelText = "Strength of brand interaction"
  nodeLabelText = "Relative size of brands"       
  
  configExport = configText
  configExport = gsub("visTitle", visTitle, configExport)
  configExport = gsub("visMore", visMore, configExport)
  configExport = gsub("visIntro", visIntro, configExport)
  configExport = gsub("edgeLabelText", edgeLabelText, configExport)
  configExport = gsub("nodeLabelText", nodeLabelText, configExport)
  
  #export data json
  sink(paste(visDir,"config.json", sep=""))
  cat(configExport)
  sink()
  
  # set web server running - start local websever
  # python -m SimpleHTTPServer 8000
  # then access at http://0.0.0.0:8000/ in browser
  
  
  
  ## End export json 
  #######################################################################################