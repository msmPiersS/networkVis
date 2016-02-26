######################################################################################
## R script for generating 2d interactive network visualisation
## 
## Process:
## 1) Loads square matrix from delimited file
## 2) Cleans up and calculates node positions (using igraph) based on different algorithms- visualsies locally in ggplot
## 3) Exports json files for use in sigmajs based interactive vis
##
## Feb 2018 - ps
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
  nodesInputFile = 'sampleNodes.csv'
  edgesInputFile = 'sampleEdges.csv'
  
  #set preferences
  clearUpperTriangle = TRUE ## set to FALSE to clear the upper triangle of the square matrix, set to FALSE to add upper triangle to lower triangle
  
  #libraries
  library(data.table) #for efficient manipulation of data
  #library(bit64)
  library(fpc) #to generate density based clusters
  library(igraph) #to generate the graph layouts
  library(ggplot2)  #to plot hte layours
  library(plyr) # for data manipulation
  library(jsonlite) #to create json export
  library(RColorBrewer) #for colour generation
  

## End Set up
#######################################################################################


  
  
#######################################################################################
## Load data and process
  
  #load raw files
  rawNodes = fread(paste(dataDir, nodesInputFile, sep=""))
  rawEdges = fread(paste(dataDir, edgesInputFile, sep=""))
  
  #quick check of dimensions
  dimEdges = dim(rawEdges)
  dimNodes = dim(rawNodes)
  if (dimEdges[1] != dimEdges[2]) {print("ERROR - input Edges file must be square matrix")}
  if (dimEdges[1] != dimNodes[1]) {print("ERROR - Nodes and Edges input files not of same dimension")}
  
  # going to make an "undirected" graph- in other words doesnt matter order
  # check edges- if non-empty upper triangle, offer to merge (add) to lower triangle or disregard
  nonZeroUpperFlag = (sum(rawEdges*upper.tri(rawEdges),na.rm = TRUE)>0)
  if (nonZeroUpperFlag) {
    # if we have values in the upper triangle of the edges matrix we have to choose either to remove them or to add them to the lower triangle
    if (clearUpperTriangle) {
      # if preference is to clear the upper triangle, do so
      edges = as.matrix(rawEdges*lower.tri(rawEdges))
    } else {
      # otherwise add upper to lower and then clear upper
      edges = as.matrix(rawEdges*lower.tri(rawEdges) + t(rawEdges*upper.tri(rawEdges)))
    }
    
  } else {
    edges = as.matrix(rawEdges)
  }
  
  # assume that node names are in column 1 of nodes input file and relative size is in column 2
  nodes = rawNodes
  setnames(nodes, colnames(nodes)[1:2], c("name", "relativeSize"))
  setkey(nodes, name)
  colnames(edges) = as.matrix(nodes[, name])
  
## End Load data and process
#######################################################################################  
  
  
  
  
#######################################################################################
## Generate network layout and plot  
  
  # setup igraph - create graph from edges matrix and define what we are going to do 
  # (use undirected graph with edge weights determined by entries in edge matrix, excluding the diagonal)
  g = graph.adjacency(edges, mode="undirected", weighted="strength", diag=FALSE)
  
  #str(g)
  #get.data.frame(g, what="edges")
  
  # initiate initial location of the nodes- to ensure repeatability
  set.seed(1234)
  layoutStart = matrix(runif(2*nrow(rawNodes)), nrow(rawNodes), 2)
  
  # generate layout using Fruhterman-Reingolg layout algorithm
  tmpLayout = layout_with_fr(g, coords=layoutStart) 
  
  # setup table of nodes with coordinates
  nodesDT = as.data.table(get.data.frame(g, what="vertices"))
  nodesDT[, Xorg:=tmpLayout[, 1]]
  nodesDT[, Yorg:=tmpLayout[, 2]]
  
  #join in node sizes
  setkey(nodesDT, name)
  nodesDT = nodes[nodesDT]
  
  # cap outliers
  capX = 100
  sum(abs(nodesDT[, Xorg])>capX)/nrow(nodesDT)
  capY = 100
  sum(abs(nodesDT[, Yorg])>capY)/nrow(nodesDT)
  
  nodesDT[, X:=Xorg]
  nodesDT[abs(Xorg)>capX, X:= capX*sign(Xorg)]
  nodesDT[, Y:=Yorg]
  nodesDT[abs(Yorg)>capY, Y:= capY*sign(Yorg)]
  
  
  #clustering using dbscan
  stds = apply(nodesDT[, list(X, Y)], 2, sd)
  ds <- dbscan(nodesDT[, list(X, Y)], mean(stds), MinPts = floor(nrow(nodesDT)/500))
  #ds <- dbscan(nodesDT[, list(Xorg, Yorg)], mean(stds)/10, MinPts = floor(nrow(nodesDT)/20))
  #dbscan(nodesDT[, list(X, Y)], eps=5, MinPts = 5)
  
  useClusters = 1
  if (useClusters != 1) {
    nodesDT[, typeId:=floor(log(1+max(relativeSize))) - floor(log(1+relativeSize))]
    
  } else {
    nodesDT[, typeId:=ds$cluster]
    
  }
  
  
  
  edgesDT = as.data.table(get.data.frame(g, what="edges"))
  setkey(edgesDT, from)
  setkey(nodesDT, name)
  
  edgesDT = nodesDT[edgesDT][, list(to, from=name, strength, X1=X, Y1=Y)]
  setkey(edgesDT, to)
  edgesDT = nodesDT[edgesDT][, list(to = name, from, strength, X1, Y1, X2=X, Y2=Y)]
  setkey(edgesDT, from)
  
  #pull in node ids names
  #setkey(edgesDT, toName)    
  #edgesDT[, to:=nodesDT[edgesDT, id]]
  #setkey(edgesDT, fromName)
  #edgesDT[, from:=networkDT[edgesDT, id]]
  
  localPlot = 0
  
  if (localPlot == 1) {
    defaultSize = 10
    pnet <- ggplot() + geom_segment(data=edgesDT, aes(x=X1, y=Y1, xend = X2, yend = Y2,size = strength), 
                                    colour="grey", alpha=0.5)
    
    
    pnet <- pnet  + 
      geom_point(data=nodesDT, aes(X, Y), color=nodesDT[, typeId], size=nodesDT[, relativeSize], alpha=0.8) +
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
  
  nodesExport = copy(nodesDT[, list(id = name, label = name, size = relativeSize, typeId, x = X, y = Y)])
  #setnames(nodesExport, c("name", "X", "Y"), c("label", "x", "y"))
  
  tmpcolList = brewer.pal(9,"Set1")
  
  tmpcol = data.table(typeId = sort(unique(nodesExport[, typeId])), color =  tmpcolList[1:length(unique(nodesExport[, typeId]))])
  
  for (i in 1:nrow(tmpcol)) {
    #i=1
    tmpcol[i, color:=paste("rgb(",paste(col2rgb(tmpcol[i, color]), collapse=","),")", sep="")]
  }
  
  setkey(tmpcol, typeId)
  setkey(nodesExport, typeId)  
  nodesExport = nodesExport[tmpcol]
  setkey(nodesExport, id)
  # round everything
  #nodesExport[, size := round(log(1+size)*1000)/1000]
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
  "maxNodeSize": 20,
  "maxEdgeSize": 3,
  "minNodeSize": 1
  },
  "drawingProperties": {
  "labelThreshold": 20,
  "hoverFontStyle": "bold",
  "defaultEdgeType": "curve",
  "defaultLabelColor": "#000",
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
  #system(paste("cd ", visDir, sep=""))
  #system("python -m SimpleHTTPServer 8000")
  
  #note you then need to top from command line and find the python pid
  # then sudo kill pid
  
  #then access at http://0.0.0.0:8000/
  
  
  
  ## End export json 
  #######################################################################################