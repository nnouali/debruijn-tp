#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    :raises ArgumentTypeError: If file doesn't exist
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                    "{0} -h"
                                    .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open(fastq_file, "r") as file:
        for line in file:
            yield(next(file).strip())
            next(file)
            next(file)

def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(0, len(read) - kmer_size + 1):
        kmer = read[i:i + kmer_size]
        yield kmer

def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    # store k-mers and their occurrences
    kmer_dict = {}  
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1  # +1 existing k-mer
            else:
                kmer_dict[kmer] = 1  # + new k-mer in dic

    return kmer_dict


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    digraph = nx.DiGraph()

    # Iterate through the k-mer dictionary
    for kmer, occurrence in kmer_dict.items():
        # Split the k-mer into prefix and suffix
        prefix = kmer[:-1]
        suffix = kmer[1:]

        # Check if the prefix and suffix nodes exist, and create them if not
        if not digraph.has_node(prefix):
            digraph.add_node(prefix)
        if not digraph.has_node(suffix):
            digraph.add_node(suffix)

        # Add an edge from the prefix to the suffix with the weight (occurrence)
        digraph.add_edge(prefix, suffix, weight=occurrence)

    return digraph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list :
        if not path :
            continue
        if delete_entry_node and delete_sink_node :
            graph.remove_nodes_from(path)
        elif delete_entry_node : 
            graph.remove_nodes_from(path[:-1])
        elif delete_sink_node :
            graph.remove_nodes_from(path[1:])          
        else :
            graph.remove_nodes_from(path[1:-1])
    return graph

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                    delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    if len(path_list) == 0:
        return graph 
    best_path = None

    weight_stddev = statistics.stdev(weight_avg_list)

    if weight_stddev > 0:
        best_path_index = weight_avg_list.index(max(weight_avg_list))
        best_path = path_list[best_path_index]

    else:
        length_stddev = statistics.stdev(path_length)

        if length_stddev > 0:
            best_path_index = path_length.index(max(path_length))
            best_path = path_list[best_path_index]

        else:
            rd_index = random.randint(0, len(path_list) - 1)
            best_path = path_list[rd_index]


    # Remove unwanted paths from the graph
    for path in path_list:
        if path != best_path:
            remove_path = path.copy()
            if not delete_entry_node:
                remove_path.pop(0)
            if not delete_sink_node:
                remove_path.pop(-1)
            graph.remove_nodes_from(remove_path)


    return graph



def path_average_weight(graph, path):
    """Compute the weight of a path


    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue


    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    # Find all simple paths between the ancestor and descendant node
    all_paths = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))
    
    path_lengths = [len(path) for path in all_paths]
    weight_avg_list = []
    
    for path in all_paths:
        weights = []
        for i in range(len(path) - 1):
            edge_data = graph[path[i]][path[i+1]]
            weight = edge_data.get('weight', 1)
            weights.append(weight)
        
        avg_weight = sum(weights) / len(weights) if weights else 1
        weight_avg_list.append(avg_weight)
    
    cleaned_graph = select_best_path(graph, all_paths, path_lengths, weight_avg_list)
    
    return cleaned_graph

def simplify_bubbles(graph):
    """Detect and explode bubbles


    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble_detected = False 


    # Parcourir chaque nœud du graphe
    for node in graph.nodes():
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            for i in range(len(predecessors)):
                for j in range(i+1, len(predecessors)):
                    common_ancestor = nx.lowest_common_ancestor(graph, predecessors[i], predecessors[j])
                    if common_ancestor:
                        bubble_detected = True
                        break
                if bubble_detected:
                    break


        if bubble_detected:
            break
    
    # Si une bulle est détectée, résolvez-la et continuez de manière récursive jusqu'à ce qu'il n'y ait plus de bulles
    if bubble_detected:
        graph = simplify_bubbles(solve_bubble(graph, common_ancestor, node))


    return graph


def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips


    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) List of starting nodes
    :return: (nx.DiGraph) A directed graph object without unwanted entry paths
    """




    for node in list(graph.nodes):
        all_paths_for_node = []
        predecessors = list(graph.predecessors(node))
        if len(predecessors)>1:
                for start_node in starting_nodes:
                    if node not in starting_nodes:
                        paths =  list(nx.all_simple_paths(graph, start_node, node))
                        all_paths_for_node.append(paths[0])#toujours un seul path
                if len(all_paths_for_node)>1 :
                    path_length = [len(path) for path in all_paths_for_node]
                    path_weigths = [path_average_weight(graph, all_paths_for_node[i])  if path_length[i] >1 else graph[paths[i][0]][paths[i][1]]["weight"] for i in range(len(all_paths_for_node)) ]
                    graph = select_best_path(graph, all_paths_for_node,path_length, path_weigths, delete_entry_node=True, delete_sink_node=False)
                    #graph = solve_entry_tips(graph, starting_nodes)
                    break
    return graph




def solve_out_tips(graph, ending_nodes):
    """Remove out tips


    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    for node in list(graph.nodes):
        all_paths_for_node = []
        successors = list(graph.successors(node))
        if len(successors)>1:
                for end_node in ending_nodes:
                    if node not in ending_nodes:
                        paths =  list(nx.all_simple_paths(graph, node, end_node))
                        all_paths_for_node.append(paths[0])#toujours un seul path
                if len(all_paths_for_node)>1 :
                    path_length = [len(path) for path in all_paths_for_node]
                    path_weigths = [path_average_weight(graph, all_paths_for_node[i])  if path_length[i] >1 else graph[paths[i][0]][paths[i][1]]["weight"] for i in range(len(all_paths_for_node)) ]
                    graph = select_best_path(graph, all_paths_for_node,path_length, path_weigths, delete_entry_node=False, delete_sink_node=True)
                    #graph = solve_out_tips(graph, ending_nodes)
                    break
    return graph



def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    starting_nodes = []

    for node in graph.nodes():
        if not any(graph.predecessors(node)):
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    sink_nodes = []
    for node in graph.nodes():
        if not any(graph.successors(node)):
            sink_nodes.append(node)
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []

    for start_node in starting_nodes:
        for end_node in ending_nodes:

            if nx.has_path(graph, start_node, end_node):
                paths = nx.all_simple_paths(graph, source=start_node, target=end_node)
                for path in paths:
                    contig = path[0]
                    for i in range(1, len(path)):
                        contig += path[i][-1]
                    contigs.append((contig, len(contig)))
    
    return contigs


def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, 'w') as output:
        for i, (contig, length) in enumerate(contigs_list):
            header = f">contig_{i} len={length}\n"
            formatted_contig = textwrap.fill(contig, width=80)
            output.write(header)
            output.write(formatted_contig + '\n')


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)"""
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                        edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    #Partie a
    #fastQ
    reads = read_fastq(args.fastq_file)
    #Kmers
    #kmer_size = args.kmer_size  
    #for read in reads:
    #    kmers = list(cut_kmer(read, args.kmer_size))
    #    #lecture et ses k-mères
    #    print("Lecture : ", read)
    #    print("K-mères : ", kmers)
    #Kmer_dict
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    #print("\nDictionnaire des k-mères et de leurs occurrences :")
    #print(kmer_dict)
    #better print
    #for kmer, occurrence in kmer_dict.items():
    #    print(f"K-mer : {kmer}, Occurrence : {occurrence}")
    
    #Partie b
    digraph = build_graph(kmer_dict)
    digraph = simplify_bubbles(digraph)

    # Identification des nœuds d'entrée et de sortie
    entrance_nodes = get_starting_nodes(digraph)

    #print("entrance_nodes : ", entrance_nodes)
    output_nodes = get_sink_nodes(digraph)

    digraph = solve_entry_tips(digraph, entrance_nodes)

    digraph = solve_out_tips(digraph, output_nodes)
    entrance_nodes = get_starting_nodes(digraph)

    #print("entrance_nodes : ", entrance_nodes)
    output_nodes = get_sink_nodes(digraph)
    #print("output_nodes : ", output_nodes)
    # Extraction des contigs
    contigs = get_contigs(digraph, entrance_nodes, output_nodes)
    
    #for contig, length in contigs:
    #    print(f"Contig : {contig}, Longueur : {length}")
    save_contigs(contigs, args.output_file)
    # Dessiner le graphe
    pos = nx.spring_layout(digraph)
    nx.draw(digraph, pos, with_labels=False, node_size=10, node_color='b', font_size=8)
    plt.savefig("graph.png")
    plt.close()
    
if __name__ == '__main__': # pragma: no cover
    main()
