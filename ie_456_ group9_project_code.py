import random
import numpy as np
import sys
import time


#Since object oriented programming approach used classes are defined in a seperate py document which is also in the project folder.
#3 classes are generated for graphs,vertices and the edges. 
from graph_factory import *


random.seed(1304)
np.random.seed(1304)
sys.setrecursionlimit(500000)


#In the following three function different degree sequence in which degrees follows normal, binomial or power_law distribution is generated. 
def get_degree_sequence_with_normal(n):
    #generating degrees between 1 and n-1 since 1 vertex cannot have degree more than n-1. Degrees are chosen wrt. random distribution 
    degrees = [int(random.choice(range(1, n-1))) for i in range(n)]
    # the first condition in the if statement ensures that sum of the degrees in the sequence is even since otherwise it is not possible to generate graphs
    # second condition which corresponds to m<n-1(minimally connected condition) ensures having degree sequences that can generate graphs that can be converted to connected graphs.  
    if sum(degrees) % 2 != 0 or sum(degrees)<2*n-2: 
        return get_degree_sequence_with_normal(n)
    return sorted(degrees,reverse=True)

#Two function below has the same conditions, they are just generating sequences with another distribution. 
def get_degree_sequence_with_binomial(n,p,k):

    degrees = np.random.binomial(n-2, p, n) +1+k
    if sum(degrees) % 2 != 0 or sum(degrees)<2*n-2: 
        return get_degree_sequence_with_binomial(n,p,k)
    return sorted(degrees,reverse=True)


def get_degree_sequence_with_power_law(n,alpha,k_min):

    degrees = [int(np.random.pareto(alpha) + k_min) for _ in range(n)]
    if sum(degrees) % 2 != 0 or sum(degrees)<2*n-2:
        return get_degree_sequence_with_power_law(n,alpha,k_min)
    return sorted(degrees,reverse=True)

#This function checks if the degree sequence is graphical or not
def is_graphical(degrees):

    n=len(degrees)
    if n==0:
        return True
    #Condition for the corrected Durfee number is checked, if holds only first m inequalities of the erdos gallai theorem is used. 
    m=0
    while n>m and degrees[m]>=m:
        m+=1
    #implementation of the erdos gallai theorem 
    for k in range(1, m+1):
        left_sum = sum(degrees[:k])
        right_sum = (k)*(k-1) + sum([min(k, i) for i in degrees[k:]])
        if left_sum > right_sum:
            return False
    return True

#The below 3 functions generates graphical degree sequence by checking if the previously constructed degree sequences graphical or not with is_graphical() function 
def generate_graphical_random_degree_sequence_with_normal(n):

    degrees=get_degree_sequence_with_normal(n)
    while not is_graphical(degrees):
        degrees=get_degree_sequence_with_normal(n)
    return degrees


def generate_graphical_random_degree_sequence_with_binomial(n,p=0.5,k=0):

    degrees=get_degree_sequence_with_binomial(n,p,k)
    while not is_graphical(degrees):
        degrees=get_degree_sequence_with_binomial(n,p,k)
    return degrees


def generate_graphical_random_degree_sequence_with_power_law(n,alpha=2.5,k_min=1):

    degrees=get_degree_sequence_with_power_law(n,alpha,k_min)
    while not is_graphical(degrees):
        degrees=get_degree_sequence_with_power_law(n,alpha,k_min)
    return degrees

#in order to make observations about graphs density, below code is constructed 
def get_density(degrees):

    n=len(degrees)
    m=sum(degrees)
    density=int(m/(n*(n-1))*100)/100
    return density


#algorithms
#4 funtion in below are for hh algorithm. 
def graph_generation_HH_random_to_high(vertecies):
    
    #In this algorithm vertex is selected randomly from the given degree sequence(vertices).
    vertex=random.choice(vertecies)
    #selected vertex is deleted from the vertecies. 
    vertecies.remove(vertex)
    #res degree is a attribute in the vertex class. adjacent vertices are added to the vertex until res_degree(residual degree) becomes zero. 
    #Res_degree decreases every time add_adjacent function is called in the vertex class.The function generates edges by calling the edge class between the choosen vertex and selected adjacent. 
    for i in range(vertex.res_degree):
        #degree sequences are sorted, it is sufficient to add the i'th (first one) since the i'th vertex will have the highest degree .
        vertex.add_adjacent(vertecies[i])
    if vertecies:
        #algortihm is recursive until they are no vertex left in vertices, the function will be called.
        return graph_generation_HH_random_to_high(Graph.sort_graph(vertecies))


#Most procedures are similar for other hh algortihms hence only minor differences will be explained
def graph_generation_HH_high_to_high(vertecies):   
    #choosing the first element from vertices which corresponds to highest degree vertex 
    vertex=vertecies[0]
    vertecies.remove(vertex)
    for i in range(vertex.res_degree):
        vertex.add_adjacent(vertecies[i])  
    if vertecies:
        return graph_generation_HH_high_to_high(Graph.sort_graph(vertecies))


def graph_generation_HH_small_to_high(vertecies):
    #choosing the last element from vertices which corresponds to smallest degree vertex since degree sequences are sorted. 
    vertex=vertecies[-1]
    vertecies.pop()
    for i in range(vertex.res_degree):
        vertex.add_adjacent(vertecies[i])
    if vertecies:
        return graph_generation_HH_small_to_high(Graph.sort_graph(vertecies))


def graph_generation_HH_small_to_small(vertecies):
    #since small to high algortihm is problematic errors are  prevented by try and except statement 
    try:
        vertex=vertecies[-1]
        vertecies.remove(vertex)
        for i in range(vertex.res_degree):
            #the vertices are distributed to last element of the vertecies which corresponds to the smallest degree vertex 
            vertex.add_adjacent(vertecies[-1-i])
        if vertecies:
            return graph_generation_HH_small_to_small(Graph.sort_graph(vertecies))
    except:
        return 

#sequential algorithm 
#Since the function takes graph object as a parameter, the functions for graph class will be used. Functions will be explained in the graph factory py document. 
def graph_generation_sequentiol(graph):
    #algorithm runs until all residuals become zero. 
    while graph.sum_res_degrees()>0:
        #after getting the vertices from graph class, min degree among them is chosen 
        min_degree_v=vertex_with_min_degree(graph.get_vertecies())
        #candidate list is constructed. 
        candidates=get_candidate_list(min_degree_v,graph.get_vertecies(),graph)
        while min_degree_v.res_degree>0 and candidates:
            #this if statement is written to fasten the algorithm. If the # of residuals degree of the min vertex is equal to the length of the candidate list,
            #there is no need to call candidate list for each residual like in else statement. All candidates are added to the vertex as an adjacent vertex,when condition holds.
            if min_degree_v.res_degree==len(candidates):
                for c in candidates:
                    min_degree_v.add_adjacent(c)
            #choosing randomly from the candidate list with a probability with respect to their weights 
            else:
                c = random.choices(candidates, weights=[i.res_degree for i in candidates])[0]
                #adding the chosen vertex from candidate list to min degree vertex as an adjacent 
                min_degree_v.add_adjacent(c)
                #removing the chosen edge from the candidate list 
                candidates.remove(c)
                #constructing the new candidate list with the function we created called get_candidate_list()
                candidates=get_candidate_list(min_degree_v,candidates,graph)


#tools
#function to get the vertex with min degree number. 
def vertex_with_min_degree(vertecies):

        def x(vertex:Vertex):
            if vertex.res_degree>0:
                return vertex.res_degree
            return 99999999
        return sorted(vertecies,key=x)[0]

#function to get the vertex with max degree number.
def vertex_with_max_degree(vertecies):
        
        def x(vertex:Vertex):
            if vertex.res_degree>0:
                return vertex.res_degree
            return 0
        return sorted(vertecies,key=x)[-1]

#function to get the candidate list for the sequnetial algortihm. 
def get_candidate_list(vertex,vertecies,graph):

    candidates=[]
    for i in vertecies:
        if i.res_degree>0 and i!= vertex and i not in vertex.adjacents.keys():
            degrees=convert_to_degree_sequence(graph.vertecies)
            if i.res_degree>vertex.res_degree:
                degrees[degrees.index(vertex.res_degree)]=vertex.res_degree-1
                degrees[degrees.index(i.res_degree)]=i.res_degree-1
            else:
                degrees[degrees.index(i.res_degree)]=i.res_degree-1
                degrees[degrees.index(vertex.res_degree)]=vertex.res_degree-1
            while 0 in degrees:
                degrees.remove(0)
            if is_graphical(sorted(degrees,reverse=True)):
                candidates.append(i)
    return candidates
                
#this function converts the sequnece into degree sequence again by updating residuals and returning the sorted graph
def convert_to_degree_sequence(verticies):

    degrees=[]
    for i in verticies:
        degrees.append(i.res_degree)
    return sorted(degrees,reverse=True)


#gives the sum of degrees of a given vertices by fist converting it to a degree sequence.
def sum_of_degrees(verticies):

    degrees=convert_to_degree_sequence(verticies)
    return sum(degrees)

#algorithm generated for pairwise interchange 
def pairwaise_interchanges(graph):
    #if graph is disconnected subgraphs of the graph is found by using get_subgraphs() function from graph class 
    if not graph.is_connected():
        subgraphs=graph.get_subgraphs()
        #subgraphs are labeled as cyclic or noncyclic 
        cyc=[]
        noncyc=[]
        for s in subgraphs:
            e=has_cycle(s)
            if e :
                #if the graph has cycle in additon to labeling as cyc, the number of cycles are hold with e 
                cyc.append((s,e))
            else:
                noncyc.append(s)
        #if there is a noncyclic subgraph, the following code runs until no noncyclic graphs are present 
        if noncyc:
            while noncyc:
                graph.num_of_pairways+=2
                #first the subgraph that has the highest # of cycles is chosen 
                cyc=sort_max_cycle(cyc)
                #edges are selected and removed from both noncylic and cyclic graph. 
                edge1=choose_edge(noncyc[-1])
                v1,j1=edge1.delete()
                #for the cyclic graph an edge is chosen from a cycle 
                edge2=random.choice(cyc[0][1])
                cyc[0][1].remove(edge2)
                v2,j2=edge2.delete()
                #by changing the adjacent vertices edges are swaped 
                v1.add_adjacent(j2)
                v2.add_adjacent(j1)
                noncyc.pop()
            return pairwaise_interchanges(graph)
        v_list=[]
        j_list=[]
        #if no noncyclic graph left, choose and edge from a cycle of each subgraph 
        for s in cyc:
            graph.num_of_pairways+=1
            v,j=s[1][0].delete()
            v_list.append(v)
            j_list.append(j)
        #By respecting the index number(taking the one before(i-1)) a cycle is formed.
        for i,v in enumerate(v_list):
            v.add_adjacent(j_list[i-1])

#a function to choose an egde from the subgraphs of pairwise interchange algorithm 
def choose_edge(vertecies)-> Edge:

    def x(vertex):
        return vertex.degree
    vertex=max(vertecies,key=x)
    other_vertex=max(vertex.adjacents,key=x)
    return vertex.adjacents[other_vertex]

#a function that founds if there is a cycle in the subgraph. Dfs algorithm is used for this purpose and the function is in the graph class 
def has_cycle(subgraph):
        
        edges=[]
        def _dfs(vertex,visited,pre_vertex=None):
            visited[vertex]+=1
            for neighbor in vertex.adjacents.keys():
                if neighbor!=pre_vertex:
                    if visited[neighbor]>=1:
                        visited[neighbor]+=1
                        if vertex.adjacents[neighbor] not in edges:
                            edges.append(vertex.adjacents[neighbor])
                    else:
                        _dfs(neighbor,visited,vertex)
        start_vertex = subgraph[0]
        visited={}
        for i in subgraph:
            visited[i]=0
        _dfs(start_vertex, visited) 
        return edges

#for cyclic graphs, to find the subgraph with the max number of cycles easily, they are sorted with respect to number of cycles they have. 
def sort_max_cycle(cyc):

    def x(tuple):
            return len(tuple[1])
    return sorted(cyc,key=x,reverse=True)







































#distributionlar için input


i=0
nler=[10,20,40,80,100]
binomial=[]
normal=[]
power_law=[]
for n in nler:
    for _ in  range(5)  :
        i+=1
        degrees=generate_graphical_random_degree_sequence_with_binomial(n)
        while True:
            if  degrees in binomial:
                degrees=generate_graphical_random_degree_sequence_with_binomial(n)
            else:
                break
        binomial.append(degrees)
        m=int(sum(degrees)/2)
        with open(f"Binomial/Inputs/Group9-{n}-{m}-Input-{i}.txt","w") as f:
            for d in degrees:
                f.write(f"{d} ")

        degrees=generate_graphical_random_degree_sequence_with_normal(n)
        while True:
            if  degrees in normal:
                degrees=generate_graphical_random_degree_sequence_with_normal(n)
            else:
                break
        normal.append(degrees)
        m=int(sum(degrees)/2)
        with open(f"Normal/Inputs/Group9-{n}-{m}-Input-{i}.txt","w") as f:
            for d in degrees:
                f.write(f"{d} ")

        degrees=generate_graphical_random_degree_sequence_with_power_law(n)
        while True:
            if  degrees in power_law:
                degrees=generate_graphical_random_degree_sequence_with_power_law(n)
            else:
                break
        power_law.append(degrees)
        m=int(sum(degrees)/2)
        with open(f"Power_law/Inputs/Group9-{n}-{m}-Input-{i}.txt","w") as f:
            for d in degrees:
                f.write(f"{d} ")




#binomial için output


i=0
o=0
for d in binomial :
    i+=1
    n=len(d)
    m=int(sum(d)/2)
    g1=Graph(d)
    g2=Graph(d)
    g3=Graph(d)
    g4=Graph(d)
    g5=Graph(d)
    o+=1
    start_time=time.time()
    graph_generation_HH_random_to_high(g1.get_vertecies())
    connection=g1.is_connected()
    while not g1.is_connected():
        pairwaise_interchanges(g1)
    process_time=time.time()-start_time
    with open(f"Binomial/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_random_to_high.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g1.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g1.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_HH_high_to_high(g2.get_vertecies())
    connection=g2.is_connected()
    while not g2.is_connected():
        pairwaise_interchanges(g2)
    process_time=time.time()-start_time
    with open(f"Binomial/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_high_to_high.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g2.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g2.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_HH_small_to_high(g3.get_vertecies())
    connection=g3.is_connected()
    while not g3.is_connected():
        pairwaise_interchanges(g3)
    process_time=time.time()-start_time
    with open(f"Binomial/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_small_to_high.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g3.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g3.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_HH_small_to_small(g4.get_vertecies())
    connection=g4.is_connected()
    while not g4.is_connected():
        pairwaise_interchanges(g4)
    process_time=time.time()-start_time
    with open(f"Binomial/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_small_to_small.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g4.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g4.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_sequentiol(g5)
    connection=g5.is_connected()
    while not g5.is_connected():
        pairwaise_interchanges(g5)
    process_time=time.time()-start_time
    with open(f"Binomial/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-Sequential.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g5.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g5.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")





#power_law için output


i=0
o=0
for d in power_law :
    i+=1
    n=len(d)
    m=int(sum(d)/2)
    g1=Graph(d)
    g2=Graph(d)
    g3=Graph(d)
    g4=Graph(d)
    g5=Graph(d)
    o+=1
    start_time=time.time()
    graph_generation_HH_random_to_high(g1.get_vertecies())
    connection=g1.is_connected()
    
    while not g1.is_connected():
        pairwaise_interchanges(g1)
    
    process_time=time.time()-start_time
    with open(f"Power_law/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_random_to_high.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g1.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g1.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_HH_high_to_high(g2.get_vertecies())
    connection=g2.is_connected()
    while not g2.is_connected():
        pairwaise_interchanges(g2)
    process_time=time.time()-start_time
    with open(f"Power_law/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_high_to_high.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g2.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g2.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_HH_small_to_high(g3.get_vertecies())
    connection=g3.is_connected()
    while not g3.is_connected():
        pairwaise_interchanges(g3)
    process_time=time.time()-start_time
    with open(f"Power_law/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_small_to_high.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g3.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g3.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_HH_small_to_small(g4.get_vertecies())
    connection=g4.is_connected()
    while not g4.is_connected():
        pairwaise_interchanges(g4)
    process_time=time.time()-start_time
    with open(f"Power_law/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_small_to_small.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g4.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g4.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_sequentiol(g5)
    connection=g5.is_connected()
    while not g5.is_connected():
        pairwaise_interchanges(g5)
    process_time=time.time()-start_time
    with open(f"Power_law/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-Sequential.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g5.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g5.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")




#normal için output


i=0
o=0
for d in normal  :
    i+=1
    n=len(d)
    m=int(sum(d)/2)
    g1=Graph(d)
    g2=Graph(d)
    g3=Graph(d)
    g4=Graph(d)
    g5=Graph(d)
    o+=1
    start_time=time.time()
    graph_generation_HH_random_to_high(g1.get_vertecies())
    connection=g1.is_connected()
    while not g1.is_connected():
        pairwaise_interchanges(g1)
    process_time=time.time()-start_time
    with open(f"Normal/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_random_to_high.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g1.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g1.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_HH_high_to_high(g2.get_vertecies())
    connection=g2.is_connected()
    while not g2.is_connected():
        pairwaise_interchanges(g2)
    process_time=time.time()-start_time
    with open(f"Normal/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_high_to_high.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g2.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g2.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_HH_small_to_high(g3.get_vertecies())
    connection=g3.is_connected()
    while not g3.is_connected():
        pairwaise_interchanges(g3)
    process_time=time.time()-start_time
    with open(f"Normal/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_small_to_high.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g3.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g3.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_HH_small_to_small(g4.get_vertecies())
    connection=g4.is_connected()
    while not g4.is_connected():
        pairwaise_interchanges(g4)
    process_time=time.time()-start_time
    with open(f"Normal/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-HH_small_to_small.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g4.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g4.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
    o+=1
    start_time=time.time()
    graph_generation_sequentiol(g5)
    connection=g5.is_connected()
    while not g5.is_connected():
        pairwaise_interchanges(g5)
    process_time=time.time()-start_time
    with open(f"Normal/Outputs/Group9-{n}-{m}-Input-{i}-Output-{o}-Sequential.txt","w") as f:
            f.write(f"{connection} \n")
            f.write(f"{g5.num_of_pairways}\n")
            f.write(f"{process_time}\n")
            for v in g5.vertecies:
                for neighbor in v.get_adjacents():
                    f.write(f"{neighbor} ")
                f.write("\n")
























