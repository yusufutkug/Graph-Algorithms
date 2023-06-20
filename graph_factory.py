import random
#Each graph object is generated with their own attributes which are vertices, edges, residuals and number of pairwise interchange needed. 
class Graph:
    #initilization of the graph object 
    def __init__(self,degrees):

        self.vertecies=[]
        self.residual_vertecies=[]
        self.edges=[]
        self.num_of_pairways=0
        for i,j in enumerate(degrees):
            self.vertecies.append(Vertex(self,i,j))
            self.residual_vertecies.append(self.vertecies[-1])

    #functions that will be used with graph objects 
    @staticmethod
    #below function takes the vertices of the graph and sorts them according to their residual degree. 
    def sort_graph(vertecies):

        def x(vertex:Vertex):
            return vertex.res_degree
        return sorted(vertecies,key=x,reverse=True)

    #if one ones to know the graph object's veritices attribute, they shoul call the get_vertices function which will return the sorted version of the vertices (degree sequence)
    def get_vertecies(self):

        self.vertecies=Graph.sort_graph(self.vertecies)
        return Graph.sort_graph(self.vertecies)

    #returns the graphs total # of vertices
    def get_num_of_verticies(self):

        return len(self.vertecies)

    #returns if the graph is connected or not by using the _dfs function which also in the graph class
    def is_connected(self):

        if self.sum_res_degrees()==0:
            start_vertex = self.vertecies[0]
            visited = [False] * len(self.vertecies)
            self._dfs(start_vertex, visited)
            if all(visited):
                return 1
            return 0
        return "fail"

    # depth first search algortihm that checks whether a graph is connected or not 
    def _dfs(self, vertex, visited):

        visited[vertex.id] = True
        for neighbor in vertex.adjacents.keys():
            if not visited[neighbor.id]:
                self._dfs(neighbor, visited)
    
    #The functions in below are mainly to get the information about the graphs attributes and since the name of the functions are self explanatory, no further explanation is made.
    def get_adjecent_matrix(self):

        n=len(self.vertecies)
        matrix=[[0 for i in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(n):
                if self. vertecies[i] in self.vertecies[j].adjacents.keys():
                    matrix[i][j]=1
        return matrix


    def sum_res_degrees(self):

        sum=0
        for i in self.vertecies:
            sum+=i.res_degree
        return sum


    def get_subgraphs(self):

        def __subgraph(vertex,subgraph):
            for neighbor in vertex.adjacents:
                    if neighbor not in subgraph:
                        subgraph.append(neighbor)
                        __subgraph(neighbor, subgraph)
            return subgraph
        subgraphs=[]
        other_part=[]
        while True:
            residual_part=[]
            for i in self.vertecies:
                if i not in other_part:
                    residual_part.append(i)
            if len(residual_part)==0:
                break
            sub=[]
            sub=__subgraph(random.choice(residual_part),sub)
            subgraphs.append(sub)
            for i in sub:
                other_part.append(i)
        return subgraphs


    def get_degree_sequence(self):

        degrees=[]
        for i in self.vertecies:
            degrees.append(i.degree)
        return sorted(degrees,reverse=True)
    

    def get_residual_vertex(self):

        out=[]
        for i in self.residual_vertecies:
            if i.res_degree<=0:
                out.append(i)
        for i in out:
            self.residual_vertecies.remove(i)
        return self.sort_graph(self.residual_vertecies)




class Vertex:
    #initilization of a vertex object. Each vertex has a graph id, degree, residual degree and adjacent list. 
    def __init__(self,graph,id,degree):

        self.my_graph=graph
        self.id=id
        self.degree=degree
        self.res_degree=degree
        self.adjacents={}
        
    #when this function is called for a vertex, a vertix is added as adjacent which means that there must be an edge between the vertex and its adjacent vertex.
    #hence by calling the edge class an edge object is created. 
    def add_adjacent(self,vertex):

        if vertex in self.adjacents.keys():
            return
        edge=Edge(self,vertex)
        #necessary related updates on graph and the attributes
        self.my_graph.edges.append(edge)
        self.adjacents[vertex]=edge
        vertex.add_me(self,edge)
        self.res_degree-=1
    
    #after an edge is constructed between the vertex and its adjacent, the edge constructed is send to the adjacent matrix and residual degree of the adjacent vertix is also updated.  
    def add_me(self,vertex,edge):

        self.adjacents[vertex]=edge
        self.res_degree-=1

    #two function below are self explanatory
    def remove_adjacent(self,vertex):

        del self.adjacents[vertex]
        self.res_degree+=1


    def get_adjacents(self):

        result=[]
        for i in self.adjacents.keys():
            result.append(i.id)
        return result
    
    #function to use when the edge of a chosen vertex wanted to be deleted. my_graph is the graph that the vertex is belongs to, one can remove an edge from the edges list of a given graph object. 
    def delete_edge(self,edge):

        self.my_graph.edges.remove(edge)




class Edge:

    #each edge object has vertex i and j and they have weights according to their degree 
    def __init__(self,i,j) -> None:

        self.vertex_i=i
        self.vertex_j=j
        self.label=(i.id,j.id)
        self.weight=i.degree+j.degree
        
    #function to delete an edge object    
    def delete(self):

        self.vertex_j.remove_adjacent(self.vertex_i)
        self.vertex_i.remove_adjacent(self.vertex_j)
        self.vertex_i.delete_edge(self)
        return self.vertex_i, self.vertex_j
    
    