import sys

# --------------------------------------------------------------------------
#   calculate the distance between two nodes
#   enter : point - longitude and latitude coordinate with radian
#   sort : the distance of the two nodes, the unit is meter
# --------------------------------------------------------------------------
def distanceLongLat(pointA, pointB):
    longA = pointA[0]
    longB =pointB[0]
    latA = pointA[1]
    latB = pointB[1]
    deltaLong = longB - longA
    # distance angulaire en radian
    S_AB = math.acos(math.sin(latA)*math.sin(latB)+math.cos(latA)*math.cos(latB)*math.cos(deltaLong))
    # distance en metre= distance angulaire * rayon terre
    return (S_AB* 6378000)


# --------------------------------------------------------------------------
#   class Vertex
#   This class saves the information of a vertex: Id, coordinates, adjacent
#   vertex, distance, whether it is vid asitend the previous vertex
# --------------------------------------------------------------------------
class Vertex:
    
    def __init__(self, node, x, y):
        self.id = node
        self.x = x
        self.y = y
        
        #all the adjacent vertexes of the current
        self.adjacent = {}
        
        #the distance between two vertexes
        self.distance = sys.maxint-1
        
        #whether the vertex is vistied
        self.visited = False
        
        #the previous vertex or from which
        self.previous = None
        
        #whether it is in the feasible list
        self.is_in_list = 0
    
    #add neighbor to a vertex
    def add_neighbor(self, neighbor, weight, cost):
        self.adjacent[neighbor] = [weight,cost]

    #return the id of the vertex
    def get_id(self):
        return self.id
    
    def get_x(self):
        return self.x
    
    def get_y(self):
        return self.y

    #return the distance between current vertex to his neighbor
    def get_weight(self, neighbor):
        return self.adjacent[neighbor][0]
    
    #return the cost when passing from current vertex to his neighbor
    def get_cost(self, neighbor):
        return self.adjacent[neighbor][1]

    #set the distance between the start vertex and the current vertex
    def set_distance(self, dist):
        self.distance = dist

    #get the distance between the start vertex and the current vertex
    def get_distance(self):
        return self.distance

    #set the previous of the current vertex
    def set_previous(self, prev):
        self.previous = prev

    #the current is already visited and set this value true
    def set_visited(self):
        self.visited = True

    #set wehther the current vertex is in the feasible list, 0 is not in and 1 is in the list
    def set_is_in_list(self,flag):
        self.is_in_list = flag

    #return 1 if the current vertex is in the feasible list , if not return 0
    def get_is_in_list(self):
        return self.is_in_list


# --------------------------------------------------------------------------
#   class Vertex
#   This class saves the information of a graph: vertexs and pathes
# --------------------------------------------------------------------------
class Graph:
    def __init__(self):
        self.vert_dict = {}
        self.num_vertices = 0

    def __iter__(self):
        return iter(self.vert_dict.values())

    #add vertex to the graph
    def add_vertex(self, node, x, y):
        self.num_vertices = self.num_vertices + 1
        new_vertex = Vertex(node, x, y)
        self.vert_dict[node] = new_vertex
        return new_vertex

    #return the vertex of the certain index n
    def get_vertex(self, n):
        if n in self.vert_dict:
            return self.vert_dict[n]
        else:
            return None

    #add an edge to the graph
    def add_edge(self, frm, to, distance, cost):
        if frm not in self.vert_dict:
            self.add_vertex(frm)
        if to not in self.vert_dict:
            self.add_vertex(to)
        #add neignbors to the vertexes
        self.vert_dict[frm].add_neighbor(self.vert_dict[to], distance, cost)
        self.vert_dict[to].add_neighbor(self.vert_dict[frm], distance, cost)



#-------------------------------------------------------------------------
#  the function shortest(v, path)
#  enter: v - target vertex, path - the list only saved the target vertex
#  effect: save the shortest we fount into the path list
#-------------------------------------------------------------------------
def shortest(v, path):
    if v.previous:
        path.append(v.previous.get_id())
        shortest(v.previous, path)



#initialize the visited list, this list is used to save the vertexes which are
# already visited. At last we will plot these vertexes in yellow
visited_list = []


#-------------------------------------------------------------------------
#  the version V0 of  dijkstra's algo
#  enter: grahe, start vertex, target vertex
#  effrct:search the vertexs between the start vertex and the target vertex
#   and find the shortest path from start vertex to target vertex
#-------------------------------------------------------------------------

def dijkstraV0(aGraph, start, target):
    print '''Dijkstra's shortest path'''
    
    #initialization
    start.set_distance(0)
    unvisited_list = [v for v in aGraph]
    
    while len(unvisited_list):
        
        #search the vertex with the smallest distance
        tmp_dis = sys.maxint
        for i in range(len(unvisited_list)):
            if tmp_dis > unvisited_list[i].get_distance():
                tmp_dis,tmp_v,tmp_i = unvisited_list[i].get_distance(),unvisited_list[i],i
        current = tmp_v

        #update the distance of the adjacent vertexes of current vertex
        for next in current.adjacent:
            
            #if next is visited, do nothing
            if next.visited:
                continue
            
            #calculate the new distance
            new_dist = current.get_distance() + current.get_weight(next)
            
            #if the new distance is less than the old one, update it
            if new_dist < next.get_distance():
                next.set_distance(new_dist)
                next.set_previous(current)
                    
        #current vertex sets visited and add it into the visited list
        del unvisited_list[tmp_i]
        current.set_visited()
        visited_list.append(current)



# --------------------------------------------------------------------------
#   the version V1 of  dijkstra's algo
#   calculate the shortest path using dijkstra algo
#   enter:graph, start vertex, target vertex
#   effect:search the vertexs between the start vertex and the target vertex
#   and find the shortest path from start vertex to target vertex
# --------------------------------------------------------------------------

def dijkstraV1(aGraph, start, target):
    print '''Dijkstra's shortest path'''  
    start.set_distance(0)

    # we change the unvisited_list with all the unvisited vertexes to the candidate list
    # with the feasible vertexes of the current vertex
    feasible_list = [start]
    while len(feasible_list):
        
        #search for the vertex with the smallest distance in the feasible list
        tmp_dis = sys.maxint
        for i in range(len(feasible_list)):
            noeud = feasible_list[i]
            if tmp_dis > feasible_list[i].get_distance():
                tmp_dis,tmp_v,tmp_i = feasible_list[i].get_distance(),feasible_list[i],i
        current = tmp_v
        
        # we have already found the shortest path ----the stop condition
        if current == target:
            break
        
        #update the distance of all the adjacent vertexes of the current vertex
        for next in current.adjacent:
            if next.visited:
                continue
            new_dist = current.get_distance() + current.get_weight(next)
            if new_dist < next.get_distance():
                next.set_distance(new_dist)
                next.set_previous(current)
            if next not in feasible_list:
                feasible_list.append(next)

        #for the current, set it visited true, add it into visited list and delete it from feasible list
        current.set_visited()
        visited_list.append(current)
        del feasible_list[tmp_i]



# --------------------------------------------------------------------------
#   the function euclideanDistance(start, target)
#   calculate the euclidean distance between two vertexes
#   enter: start vertex, target vertex
#   quit: the euclidean distance between two vertexes
# --------------------------------------------------------------------------
def euclideanDistance(start, target):
    return ((start.get_x()-target.get_x()) ** 2 + (start.get_y()-target.get_y()) ** 2) ** 0.5




# --------------------------------------------------------------------------
#   the version V1 of  Astar's algo
#   calculate the shortest path using A star algorithm
#   enter: graph, start vertex, target vertex
#   effect: search the vertexs between the start vertex and the target vertex
#   and find the shortest path from start vertex to target vertex
# --------------------------------------------------------------------------
#@profile
def astarV1(aGraph, start, target):
    print '''AStar's shortest path'''
    
    #initialization
    start.set_distance(0)
    visited_count=0
    feasible_list = [start]
    
    while len(feasible_list):
        visited_count = visited_count + 1
        
        #search the vertex with the smallest total weight in the feasible
        tmp_dis = sys.maxint
        for i in range(len(feasible_list)):
            noeud = feasible_list[i]
            tmp_e_dis = feasible_list[i].get_distance() + euclideanDistance(noeud, target)
            if tmp_dis > tmp_e_dis:
                tmp_dis,tmp_v,tmp_i = tmp_e_dis,feasible_list[i],i
        current = tmp_v
        
        #if the current vertex wo found in feasible is the target, the algo stops
        if current == target:
            return visited_count

        #update the distance of all the adjacent vertexes of the current vertex
        for next in current.adjacent:
            if next.visited:
                continue
            new_dist = current.get_distance() + current.get_weight(next)
            if new_dist < next.get_distance():
                next.set_distance(new_dist)
                next.set_previous(current)
            if next not in feasible_list:
                feasible_list.append(next)

        #for the current, set it visited true, add it into visited list and delete it from feasible list
        current.set_visited()
        visited_list.append(current)
        del feasible_list[tmp_i]



# --------------------------------------------------------------------------
#   the version V2 of  Astar's algo
#   calculate the shortest path using A star algorithm
#   enter: graph, start vertex, target vertex
#   effect: search the vertexs between the start vertex and the target vertex
#   and find the shortest path from start vertex to target vertex
#   Compared to V1 of Astar's algo, there are two improvements:
#   1. we add a new list called ed_list to save the eculidien distance of each vertex
#   2. we add a new parameter(is_in_list) to the vertex class to show whether the
#      vertex is in the feasible list
# --------------------------------------------------------------------------
#@profile
def astarV2(aGraph, start, target):
    print '''AStar's shortest path'''
    
    #initialization
    start.set_distance(0)
    feasible_list=[]
    ed_list=[]
    visited_count = 0
    feasible_list.append(start)
    start.set_is_in_list(1)
    ed_list.append(euclideanDistance(start, target))
    
    while len(feasible_list):
        visited_count = visited_count + 1
        
        #search the vertex with the smallest total weight in feasible list
        tmp_dis = sys.maxint
        for i in range(len(feasible_list)):
            tmp_e_dis = ed_list[i] + feasible_list[i].get_distance()
            if tmp_dis > tmp_e_dis:
                tmp_dis,tmp_v,tmp_i = tmp_e_dis,feasible_list[i],i
        current = tmp_v
    
        #for the current, set it visited true, add it into visited list and delete it from feasible list
        current.set_visited()
        visited_list.append(current)
        feasible_list[tmp_i].set_is_in_list(0)
        del feasible_list[tmp_i]
        del ed_list[tmp_i]
        
        #if the current vertex wo found in feasible is the target, the algo stops
        if current == target:
            return visited_count
        
        #update the distance of all the adjacent vertexes of the current vertex
        for next in current.adjacent:
            if next.visited:
                continue
            new_dist = current.get_distance() + current.get_weight(next)
            if new_dist < next.get_distance():
                next.set_distance(new_dist)
                next.set_previous(current)
            if next.get_is_in_list() == 0:
                feasible_list.append(next)
                next.set_is_in_list(1)
                ed_list.append(euclideanDistance(next, target))





# --------------------------------------------------------------------------
#   the version V2 of  Dijkstra's algo
#   calculate the shortest path using Dijkstra algorithm
#   enter: graph, start vertex, target vertex, alpha in [0,1]
#   effect: search the vertexs between the start vertex and the target vertex
#   and find the shortest path from start vertex to target vertex
#   Compared to V1 of Dijkstra algo, this version has two arguments -- distance
#   and danger
#   We add a new parameter alpha to calculat the total weight of the vertex
# --------------------------------------------------------------------------
def dijkstra2(aGraph, start, target, alpha):
    print '''Two parameters' shortest path'''
    
    #initialization
    start.set_distance(0)
    feasible_list = [start]
    
    while len(feasible_list):
        
        #search the vertex with the smallest total weight in feasible list
        tmp_dis = sys.maxint
        for i in range(len(feasible_list)):
            if tmp_dis > feasible_list[i].get_distance():
                tmp_dis,tmp_v,tmp_i = feasible_list[i].get_distance(),feasible_list[i],i
        current = tmp_v
        
        #if the current vertex wo found in feasible is the target, the algo stops
        if current == target:
            break
        
        #update the distance of all the adjacent vertexes of the current vertex
        for next in current.adjacent:
            if next.visited:
                continue
            new_dist = current.get_distance() + alpha * current.get_weight(next) + (1 - alpha)*current.get_cost(next)
            if new_dist < next.get_distance():
                next.set_distance(new_dist)
                next.set_previous(current)
            if next not in feasible_list:
                feasible_list.append(next)

        #for the current, set it visited true, add it into visited list and delete it from feasible list
        current.set_visited()
        visited_list.append(current)
        del feasible_list[tmp_i]





# the main method
if __name__ == '__main__':

    import csv
    import matplotlib.pyplot as plt
    import math
    import time
    g = Graph()
    
    # read file paris_noueds.csv
    csvfile = file('paris_noeuds.csv','rU')
    reader = csv.reader(csvfile)
    list_node=[];list_lon=[];list_lat=[]
    count=0
    for line in reader:
        count=count+1
        tmp = "".join(line)
        node,longitude,latitude = tmp.split("\t")
        list_node.append(int(node))
        list_lon.append(float(longitude)/180*math.pi)
        list_lat.append(float(latitude)/180*math.pi)
    csvfile.close()


    # to transform the longitude and the latitude of the node to x, y.
    longitude_min=min(list_lon)
    latitude_min=min(list_lat)
    list_x=[];list_y=[]
    for i in range(count):
        x = distanceLongLat([list_lon[i],list_lat[i]],[longitude_min,list_lat[i]])
        y = distanceLongLat([list_lon[i],latitude_min],[list_lon[i],list_lat[i]])
        g.add_vertex(list_node[i],x,y)
        list_x.append(x)
        list_y.append(y)
        
    # read file paris_arcs.csv
    csvfile = file('paris_arcs.csv','rU')
    reader = csv.reader(csvfile)

    for line in reader:
        tmp = "".join(line)
        node1,node2,dis,cost = tmp.split("\t")
        g.add_edge(int(node1),int(node2),int(dis),int(cost))
    csvfile.close()

    time_start=time.time()
    #dijkstraV0(g, g.get_vertex(4), g.get_vertex(2109))
    #dijkstraV1(g, g.get_vertex(4), g.get_vertex(2109))
    #print astarV1(g, g.get_vertex(4), g.get_vertex(2109))
    #print astarV2(g, g.get_vertex(4), g.get_vertex(2109))
    alpha = 0.1
    dijkstra2(g, g.get_vertex(4), g.get_vertex(2109), alpha)
    time_finish=time.time()

    #show the time the algorithm used
    print 'spent:'
    print time_finish-time_start

    #find the shortest by using the function shortest()
    target = g.get_vertex(2109)
    path = [target.get_id()]
    shortest(target, path)
    print 'The shortest path : %s' %(path[::-1])

    #show thw distance of the shortest path
    print 'The distance is:'
    print target.get_distance()

    #print all the nodes in blue
    plt.plot(list_x,list_y,'b.')

    #print all the visited nodes in yellow
    list_x_visited=[];list_y_visited=[]
    for v in visited_list:
        list_x_visited.append(v.get_x())
        list_y_visited.append(v.get_y())
    plt.plot(list_x_visited,list_y_visited,'y.')

    #print the shortest path in red
    for v in g:
        if v.get_id() in path:
            plt.plot(v.get_x(), v.get_y(),'r.')
    plt.show()



