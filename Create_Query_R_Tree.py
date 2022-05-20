#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 10:36:04 2022

@author: felixrosenbergernew
"""

import sys
import math
import time

B = 20

def main():
    """This function reads in the dataset, initiates the R-Tree creation, executes range queries, and displays the results."""
    # read data in
    points = [] # keep track of number of points
    n = 0 # set counter
    with open("R_Tree_Data.txt", 'r') as dataset: 
        for data in dataset.readlines(): # each line in data set
            data = data.split() # split by space
            points.append({ 
                'id': int(data[0]),
                'x': int(data[1]),
                'y': int(data[2])
                   }) # create dictionary where id and coordinates are defined and each line is appended accordingly
            n += 1 # increment counter with each line

    # build R-Tree
    rtree = RTree() # initiate class

    print("Building the R-Tree: Please wait...\n")
 
    for point in points: # insert data points from the root one by one 
        rtree.insert(rtree.root, point) # invoke insert algorithm from rtree class with root node and each point as parameter
   
    print("R-Tree construction completed\n")

    # read queries in from text file, same logic as above
    queries = []
    with open("test_queries.txt", 'r') as querydata:
        for line in querydata.readlines():
            query = line.split()
            queries.append({
                'x1': int(query[0]),
                'x2': int(query[1]),
                'y1': int(query[2]),
                'y2': int(query[3])
            })
    print("The current queries are", *queries, sep="\n") # prints queries in table format
    
    # execute sequential query and store the number of points found
    seq_results =  []
    sequential_query_start = time.time() # store the number of seconds passed since epoch before query execution
    for query in queries: # for each query, iterate through all points and check whether they lie within the query
        count = 0 # set counter to zero
        for point in points:
            if query["x1"] <= point["x"] <= query["x2"] and query["y1"] <= point["y"] <= query["y2"]:
                count += 1 # increment counter if a point lies within query
        seq_results.append(count) # append the sum of all points found per query to results
        
    sequential_query_end = time.time() # store the number of seconds passed since epoch after query execution
    seq_query_processing_time = sequential_query_end - sequential_query_start # calculate time required to run all queries
    seq_query_avg_processing_time = seq_query_processing_time/len(seq_results) # calculate average processing time per query
    
    print("The results for the sequential query are:\n")
    seq_query_num = 1
    for result in seq_results:
        print("For query", seq_query_num, "there are", result, "data points included in the range query.\n")
        seq_query_num += 1

    print("The total sequential query time is", seq_query_processing_time, "seconds.\n")
    print("The average sequential query time is", seq_query_avg_processing_time, "seconds.\n")  
    
    # execute R-Tree query and store the number of points found
    results = []
    tree_query_start = time.time() # store the number of seconds passed since epoch before query execution
    for query in queries: # for each query, invoke query algorithm and append number of points to results list
        results.append(rtree.query(rtree.root, query)) 
    tree_query_end = time.time() # store the number of seconds passed since epoch after query execution
    tree_query_processing_time = tree_query_end - tree_query_start # calculate time required to run all queries
    tree_avg_processing_time = tree_query_processing_time/len(queries) # calculate average processing time per query
    
    print("The results for the R-Tree query are:\n")
    query_num = 1
    for result in results: # for each query print the number of points found
        print("For query", query_num, "there are", result, "data points included in the range query.\n")
        query_num += 1
        
    print("The total R-Tree query time is", tree_query_processing_time, "seconds.\n")
    print("The average R-Tree query time is", tree_avg_processing_time, "seconds.\n")
    efficiency = seq_query_processing_time/tree_query_processing_time # calculate how much faster i.e. more efficient the R-Tree query algorithm is compared to sequential query
    print("The R-Tree query algorithm is", efficiency, "times faster than the standard sequential query.\n")

class Node(object): # node class
    def __init__(self):
        """The init function is always executed when a class is being initiated. It is used to assign values to object properties or 
        other necessary operations."""
        self.id = 0 
        self.child_nodes = [] # each internal node can contain a list of child nodes
        self.data_points = [] # each leaf node can contain a list of data points
        self.parent = None
        self.MBR = {
            'x1': -1,
            'y1': -1,
            'x2': -1,
            'y2': -1,
        }
        
    def perimeter(self):
        """This function calculates the half perimeter of a node's MBR property. It is used in the split() function."""
        return (self.MBR['x2'] - self.MBR['x1']) + (self.MBR['y2'] - self.MBR['y1'])

    def is_overflowing(self):
        """This function checks for both leaf and child nodes if they are overflowing by comparing their length to the parameter B. It
        returns True if a node is overflowing and False if not. It is used in the insert() and handle_overflow() functions."""
        if self.is_leaf():
            if self.data_points.__len__() > B: # check violation of upper bound where B is the upper bound.
                return True
            else:
                return False
        else:
            if self.child_nodes.__len__() > B: # check violation of upper bound where B is the upper bound.
                return True
            else:
                return False
            
    def is_underflowing(self):
        """This function checks whether nodes are underflowing which could happen if data points are deleted from the R-Tree."""
        if self.is_leaf():
            if self.data_points.__len__() < 0.4*B: # check violation of lower bound 
                return True
            else:
                return False
        else:
            if self.child_nodes.__len__() < 0.4*B: # check violation of lower bound
                return True
            else:
                return False
        

    def is_root(self):
        """This function checks whether a node is a root node by checking the value for a node's parent property. It returns True 
        if there is no parent node related to the node being checked and False otherwise. It is used in the handle_overflow() function."""
        if self.parent is None:
            return True
        else:
            return False

    def is_leaf(self):
        """This function checks whether a node is a leaf node by checking the length of the list for a node's child property. 
        It returns True if the length is zero and False otherwise. It is used in the is_underflowing(), is_overflowing(), 
        insert(), choose_subtree(), split(), query(), update_mbr() functions."""
        if self.child_nodes.__len__() == 0:
            return True
        else:
            return False

class RTree(object): # R-tree class

    def __init__(self):
        """This function invokes a node object with all its properties with every creation of an R-Tree object."""
        self.root = Node() # Create a root
                  
    def insert(self, u, p): # insert p(data point) to u (MBR)
        """This function inserts a data point into the R-Tree including overflow handling. Its input parameters are a node and a data
        point. For leaf nodes, it adds the point to the node, checks for overflow and invokes the handle_overflow() function to resolve it.
        For internal nodes, it chooses the best subtree until a leaf node is found recursively and updates the node's MBR."""
        if u.is_leaf(): 
            self.add_data_point(u, p) # add data point and update corresponding MBR
            if u.is_overflowing():
                self.handle_overflow(u) # handle overflow for leaf nodes by invoking handle overflow algo
        else:
            v = self.choose_subtree(u, p) # choose best subtree according to perimeter sum increase (minimum)
            self.insert(v, p) # invoke insert algorithm recursively
            self.update_mbr(v) # update the MBR for inserting the data point, for leaf not needed as included in handle overflow algo

    def choose_subtree(self, u, p): 
        """This function chooses either a leaf where a point is inserted or the appropriate child node. For each child node it saves 
        the necessary perimeter increase (if the node is chosen) and compares the increase to the next nodes with each iteration and 
        changes only if smaller increases are found. In the end, the child node with the smallest perimeter increase is returned 
        which is in line with R-Tree conditions. Its input parameters are a node and a data point. The function is used in the insert()
        function to decide which subtree of an R-Tree is best to insert a new data point into."""
        if u.is_leaf(): # if u is a leaf, it is chosen directly as the node where a point is inserted
            return u
        else:
            min_increase = sys.maxsize # set initial value to highest possible so that first iteration will experience smaller perimeter value that can be saved
            best_child = None
            for child in u.child_nodes: # iterate through u's child nodes
                if min_increase > self.peri_increase(child, p): # if new perimeter < current minimum
                    min_increase = self.peri_increase(child, p) # add new perimeter increase as current minimum
                    best_child = child # save current child as new best
            return best_child

    def handle_overflow(self, u):
        """This function handles a node overflow by splitting the node into two sub nodes. It then adds the two sub nodes either to 
        a newly created root (if the overflowed node was a root) or to a newly created parent (if the overflowed node was not a sub 
        node). In the latter case it invokes the function recursively in case the newly created parent overflows too. Its input parameter
        is a single node and the function is used in the insert() function for the case when a leaf node overflows."""
        u1, u2 = self.split(u) # u1 u2 are the two splits returned by the split algo
        if u.is_root(): # if u is the root, create a new root with u1 and u2 being its child nodes
            new_root = Node() # initiate new root as the parent of the two splits
            self.add_child(new_root, u1) # add the previously calculated splits to the new root
            self.add_child(new_root, u2)
            self.root = new_root
            self.update_mbr(new_root) # update the new root's MBR after creation
        else: # if u is not the root, delete u, set u1 and u2 as the children of u's parents
            w = u.parent # assign new parent to u
            w.child_nodes.remove(u) # remove u as child node
            self.add_child(w, u1) # add the splits to the new parent and update MBR
            self.add_child(w, u2)
            if w.is_overflowing():
                self.handle_overflow(w) # invoke handle overflow algorithm recursively
                
    def peri_increase(self, node, p): # calculate the increase of the perimeter after inserting the new data point
        """This function calculates the perimeter increase of MBRs a subtree needs to be chosen. From either the MBRs or the new 
        data points x/y coordinate, it selects the respective max/min value so that the new data point is included if it is outside 
        the existing MBR. It then subtracts the max from the min values of each dimension and subtracts the original perimeter from 
        the result to get the increase. Its input parameters are a node and a data point. This function is used in the choose_subtree()
        function to decide whether a child node has the minimum parameter increase when inserting a data point."""
        origin_mbr = node.MBR
        x1, x2, y1, y2 = origin_mbr['x1'], origin_mbr['x2'], origin_mbr['y1'], origin_mbr['y2']
        increase = (max([x1, x2, p['x']]) - min([x1, x2, p['x']]) + # from either existing MBR or new point, extract max/min x coordinate and subtract to get length
                    max([y1, y2, p['y']]) - min([y1, y2, p['y']])) - node.perimeter() # same for y dimension; then subtract original perimeter to get increase
        return increase
            
    def split(self, u):
        """This function splits nodes by returning the best combination of splits with regards to the total perimeter sum of 
        the two splits. The final split will be determined by the split that minimises the total perimeter increase depending 
        on the x and y (leafs) or MBR (nodes) coordinates. Its input parameter is a node and it is used in the handle_overflow() function."""
        # split u into s1 and s2
        best_s1 = Node()
        best_s2 = Node()
        best_perimeter = sys.maxsize
        if u.is_leaf(): # check whether u is a leaf node
            m = u.data_points.__len__() # get number of data points
            # sort x and y dimension coordinates respectively
            divides = [sorted(u.data_points, key=lambda data_point: data_point['x']), # specify key with lambda function --> x dimension coordinates only
                       sorted(u.data_points, key=lambda data_point: data_point['y'])] # y dimension
            for divide in divides: # iterate through both x and y dimension separately
                for i in range(math.ceil(0.4 * B), m - math.ceil(0.4 * B) + 1): # iterate through the number of possible combinations depending on the number of data points in the leaf and withouth violating R-Tree conditions
                    s1 = Node() # create new node for first split
                    s1.data_points = divide[0: i] # for each iteration cut the first half of the dimension depending on i
                    self.update_mbr(s1) # update the node's mbr
                    s2 = Node() # create new node for second split
                    s2.data_points = divide[i: divide.__len__()] # cut the second half of the dimension from i until the end
                    self.update_mbr(s2) # update the node's mbr
                    if best_perimeter > s1.perimeter() + s2.perimeter(): # check total sum of perimeter
                        best_perimeter = s1.perimeter() + s2.perimeter() # update best perimeter
                        # update best splits
                        best_s1 = s1
                        best_s2 = s2
        else: # if u is not a leaf node
            m = u.child_nodes.__len__() # get number of child nodes
            # instead of x and y dimension, sort based on MBR dimensions (2*x, 2*y)
            divides = [sorted(u.child_nodes, key=lambda child_node: child_node.MBR['x1']),  # create sorted list of x1 coordinates of all child nodes of u
                       sorted(u.child_nodes, key=lambda child_node: child_node.MBR['x2']),
                       sorted(u.child_nodes, key=lambda child_node: child_node.MBR['y1']),
                       sorted(u.child_nodes, key=lambda child_node: child_node.MBR['y2'])]
            for divide in divides:
                for i in range(math.ceil(0.4 * B), m - math.ceil(0.4 * B) + 1): 
                    s1 = Node()
                    s1.child_nodes = divide[0: i]
                    self.update_mbr(s1)
                    s2 = Node()
                    s2.child_nodes = divide[i: divide.__len__()]
                    self.update_mbr(s2)
                    if best_perimeter > s1.perimeter() + s2.perimeter(): # check total sum of perimeter
                        best_perimeter = s1.perimeter() + s2.perimeter() # update best perimeter
                        # update best splits
                        best_s1 = s1
                        best_s2 = s2
        # assign best nodes as parent nodes of respective child nodes
        for child in best_s1.child_nodes:
            child.parent = best_s1
        for child in best_s2.child_nodes:
            child.parent = best_s2

        return best_s1, best_s2

    def add_child(self, node, child):
        """This function adds a child node to the current parent node and updates the MBR. For both the x and y dimension 
        the bounding area of the child node is compared to the bounding area of the parent node. If the bounding area of 
        the child node in any dimension is outside the bounding area of the parent node, the parent's MBR is extended accordingly 
        to contain the child. Its input parameters are a node and a child and it is used in the handle_overflow() function when a
        splitted node needs to be attached as a child node to a parent."""
        node.child_nodes.append(child) # adds the child node to child node property of the parent (node)
        child.parent = node # assign the node as being the parent
        # return the child whose MBR requires the minimum increase in perimeter to cover p
        if child.MBR['x1'] < node.MBR['x1']:
            node.MBR['x1'] = child.MBR['x1']
        if child.MBR['x2'] > node.MBR['x2']:
            node.MBR['x2'] = child.MBR['x2']
        if child.MBR['y1'] < node.MBR['y1']:
            node.MBR['y1'] = child.MBR['y1']
        if child.MBR['y2'] > node.MBR['y2']:
            node.MBR['y2'] = child.MBR['y2']

    def add_data_point(self, node, data_point): 
        """This function adds a data point to a node and updates the corresponding MBR. For both the x and y dimension it is 
        checked whether the data point is smaller or larger than the MBR's border and the MBR is then updated accordingly to 
        include the data point. Its input parameters are a node and a data point. It is used in the insert() function when a
        leaf node is found."""
        node.data_points.append(data_point) # add data point to data points list property of node
        if data_point['x'] < node.MBR['x1']: # if points x coordinate is smaller than MBR's x1 coordinate
            node.MBR['x1'] = data_point['x'] # extend MBR's x1 coordinate to include the data point
        if data_point['x'] > node.MBR['x2']:
            node.MBR['x2'] = data_point['x']
        if data_point['y'] < node.MBR['y1']:
            node.MBR['y1'] = data_point['y']
        if data_point['y'] > node.MBR['y2']:
            node.MBR['y2'] = data_point['y']


    def update_mbr(self, node): # necessary after insertion
        """This function updates the MBR for either a leaf node or a parent after insertion. It iterates through the data points /
        child nodes and saves the x and y coordinates to a respective list. It then extracts the min and max values of each 
        dimension to form a new MBR which is then assigned as the node's new MBR. Its input parameter is a node and the function 
        is used in the insert(), handle_overflow(), and split() functions as these require MBRs to change."""
        x_list = []
        y_list = []
        if node.is_leaf():
            x_list = [point['x'] for point in node.data_points] # points from node's x dimension
            y_list = [point['y'] for point in node.data_points] # points from node's y dimension
        else:
            x_list = [child.MBR['x1'] for child in node.child_nodes] + [child.MBR['x2'] for child in node.child_nodes] # save all x coordinates for each child's MBR
            y_list = [child.MBR['y1'] for child in node.child_nodes] + [child.MBR['y2'] for child in node.child_nodes] # save all y coordinates for each child's MBR
        new_mbr = { # extract minimum coordinates for each dimension to form new MBR
            'x1': min(x_list),
            'x2': max(x_list),
            'y1': min(y_list),
            'y2': max(y_list)
        }
        node.MBR = new_mbr # assign newly formed MBR to update the node's MBR
        
    def query(self, node, query): #run to answer the query
        """This function returns the number of data points that intersect with a range query's coordinates. For leaf nodes, it returns the
        points that lie within the range query. For internal nodes, it checks for each of the node's child nodes whether their MBR intersects
        with the range query and invokes the search algorithm recursively on each child node. Its input parameters are a node and a range
        query. It is used as the query algorithm for the R-Tree."""
        num = 0
        if node.is_leaf(): # check if the node is a leaf node
            for point in node.data_points: # for each point in the leaf node
                if self.is_covered(point, query): # check if point lies within the range query
                    num += 1 # if yes, increment the number of found data points by 1
            return num # return the number of data points that lie within the query
        else: # if node is not a leaf node
            for child in node.child_nodes: # for each child node 
                if self.is_intersect(child, query): # check if child's MBR intersects with range query 
                    num = num + self.query(child, query) # call search algorithm on child recursively and increment number of data points only for data points that lie within the range query in later found leaf nodes
            return num

    def is_covered(self, point, query):
        """This function returns True if the input data point lies within the query range. This is done by comparing the coordinates of the 
        query with the coordinates of the data point. Its input parameters are a data point and a range query. This function is used
        in the query algorithm to check whether points in a leaf node lie within the range query."""
        x1, x2, y1, y2 = query['x1'], query['x2'], query['y1'], query['y2'] # extract coordinates from query range
        if x1 <= point['x'] <= x2 and y1 <= point['y'] <= y2: # check if point's x/y coordinate lies between query'x x/y coordinates
            return True
        else:
            return False    

    def is_intersect(self, node, query): 
        """This function checks for each dimension whether the distance between the two center points of a rectangle is less than
        or equal to the sum of either the length/width of each rectangle divided by two. The sum of the two halfs can only be smaller than
        the distance between the center points if the rectangles don't intersect and the condition must hold for all dimensions.
        Its input parameters are a node and the range query. It returns either True or False.
        This function is used in the query algorithm to check whether a child node intersects with the range query."""
        center1_x = (node.MBR['x2'] + node.MBR['x1']) / 2 # calculate center point of node's x-dimension
        center1_y = (node.MBR['y2'] + node.MBR['y1']) / 2 # calculate center point of node's y-dimension
        length1 = node.MBR['x2'] - node.MBR['x1'] # calculate length of node's x-dimension
        width1 = node.MBR['y2'] - node.MBR['y1'] # calculate length of node's y-dimension
        center2_x = (query['x2'] + query['x1']) / 2 # calculate center point of query's x-dimension
        center2_y = (query['y2'] + query['y1']) / 2 # calculate center point of query's y-dimension
        length2 = query['x2'] - query['x1'] # calculate length of query's x-dimension
        width2 = query['y2'] - query['y1'] # calculate length of query's y-dimension
        # if node and query range intersect, then the difference between their center points in each dimension must be smaller than the sum of the distances from each rectangle's center point to the respective border in each dimension
        if abs(center1_x - center2_x) <= length1 / 2 + length2 / 2 and abs(center1_y - center2_y) <= width1 / 2 + width2 / 2:  
            return True
        else:
            return False  

if __name__ == '__main__':
    main()