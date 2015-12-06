from igraph import *

#All the parameter graph in functions defined below is an object instance of Graph class in igraph
graph=Graph()

#get parameters for a given graph
#http://med.bioinf.mpi-inf.mpg.de/netanalyzer/help/2.7/#refDong
def get_parameters(graph,directed=False,unconn=True):
	param_dict=OrderedDict()
	#number of connected components
	#indicate the connectivity of network and a lower number of connected components suggests a stronger connectivity
	param_dict['connected_componet_number']=len(graph.clusters())
	#parameters related to shortest paths
	#network diameter, network radius, average shortest path length
	param_dict['diameter']=graph.diameter()
	param_dict['radius']=graph.radius()
	param_dict['average_shortest_path_length']=mean(map(lambda x:len(x),graph.shortest_paths()))
	#parameters related to neighborhood
	#average number of neighbors, network density, number of isolated nodes, network centralization, network heterogeneity, number of multi-edge node pairs
	param_dict['average_neighbor_number']=mean(map(lambda x:len(x),graph.neighborhood()))
	param_dict['density']=graph.density()
	param_dict['isolated_node_number']=len(graph.vs.select(_degree_eq=0))
	param_dict['centralization']=len(graph.vs)/(len(graph.vs)-2)*(max(map(lambda x:len(x),graph.neighborhood()))/(len(graph.vs)-1)-graph.density())
	param_dict['heterogeneity']=cmath.sqrt(variance(map(lambda x:len(x),graph.neighborhood()))).real/mean(map(lambda x:len(x),graph.neighborhood()))
	param_dict['multiple_edge_node_pair_number']=len(list(set([sorted((x['from'],x['to'])) for x in graph.es if graph.is_multiple(x)])))
	#clustering coefficient
	#a measure of the degree to which nodes in a graph tend to cluster together
	param_dict['clustering_coefficient']=clustering_coefficient(graph)
	return param_dict



#analyze the stress centrality for a given graph
#stress centrality of a node is the number of shortest paths passing through
def stress_centrality(graph):
	if isinstance(graph,Graph):
		pass
	else:
		raise Exception('Error: Invalid input graph!')
	paths={}
	vertex_pairs=list(itertools.combinations([x.index for x in graph.vs],2))
	for p in vertex_pairs:
		source,target=p
		paths[p]=graph.get_shortest_paths(source,to=graph.vs.select(lambda x:x.index==target))

	centrality_dict={}
	centrality_dict=centrality_dict.fromkeys([x.index for x in graph.vs],0)
	for v in paths.values():
		for x in v:
			centrality_dict[x]+=1

	score_list=sorted([(graph.vs[x[0]],x[1]) for x in centrality_dict.items()],key=lambda x:x[1],reverse=True)
	return score_list

#analyze the radiality centrality for a given graph
#the radiality is defined as sum(d_G + 1 - d(v, w))/(n - 1). where d(w, v) is the length of the shortest path from node w to node v, d_G is the diameter of the network, n is the size of the network.
def radiality_centrality(graph):
	if isinstance(graph,Graph):
		pass
	else:
		raise Exception('Error: Invalid input graph!')

	diameter=graph.diameter()
	size=len(graph.vs)
	score_dict={}
	for v in graph.vs:
		score_dict[v]=sum(map(lambda x:graph.shortest_paths(source=v.index,target=x.index)[0][0]==float('Inf') and diameter+1 or diameter+1-graph.shortest_paths(source=v.index,target=x.index)[0][0]/(size-1),graph.vs))
	score_list=sorted(centrality_dict.items(),key=lambda x:x[1],reverse=True)

#analyze the centroid centrality for a given graph
#http://www.cbmc.it/~scardonig/centiscape/CentiScaPefiles/CentralitiesTutorial.pdf
def centroid_centrality(graph):
	pathm=np.array([])
	adjm=graph.get_adjacency().data
	for v1 in xrange(0,len(graph.vs)):
		temp=[]
		for v2 in xrange(0,len(graph,vs)):
			if adjm[i][j]>0:
				temp.append(graph.shortest_paths(source=v1,target=v2)[0][0])
			else:
				temp.append(0)
		pathm=np.append(pathm,temp)

	fm=np.zeros((len(graph.vs),len(graph.vs)))
	for v1 in xrange(0,len(graph.vs)):
		for v2 in xrange(0,len(graph.vs)):
			fm[v1][v2]=len(filter(lambda x:x[0]<x[1],zip(pathm[v1],pathm[v2])))-len(filter(lambda x:x[0]>x[1],zip(pathm[v1],pathm[v2])))

	score_dict={}
	for v in graph.vs:
		score_dict[v]=min(fm[v.index])

	return sorted(score_dict.items(),key=lambda x:x[1],reverse=True)

#analyze the eccentricity centrality for a given graph
def eccentricity_centrality(graph):
	score_list=graph.eccentricity()
	return sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)

#analyze the centrality for a given graph
#centrality measurements included:
#	degree centrality
def degree_centrality(graph):
	score_list=graph.degree()
	return sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
#	betweenness centrality
def betweenness_centrality(graph):
	score_list=graph.betweenness()
	return sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
#	closeness centrality
def closeness_centrality(graph):
	score_list=graph.closeness()
	return sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
#	eigenvector centrality
def eigenvector_centrality(graph):
	score_list=graph.evcent()
	return sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
#	pagerank centrality
def pagerank_centrality(graph):
	score_list=graph.pagerank()
	return sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
#analyze the centrality for a given graph
def analyze_centrality(graph):
	centrality_dict=OrderedDict()

	print 'Analyzing degree centrality...'
	score_list=graph.degree()
	centrality_dict['degree']=sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
	print 'Done!'
	print

	print 'Analyzing betweenness centrality...'
	score_list=graph.betweenness()
	centrality_dict['betweenness']=sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
	print 'Done!'
	print

	print 'Analyzing  closeness centrality...'
	score_list=graph.closeness()
	centrality_dict['closeness']=sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
	print 'Done!'
	print

	print 'Analyzing eigenvector centrality...'
	score_list=graph.evcent()
	centrality_dict['eigenvector']=sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
	print 'Done!'
	print

	print 'Analyzing pagerank centrality...'
	score_list=graph.pagerank()
	centrality_dict['pagerank']=sorted([(graph.vs[i],score_list[i]) for i in xrange(0,len(score_list))],key=lambda x:x[1],reverse=True)
	print 'Done!'
	print

	return centrality_dict
