#!/usr/bin/env python3
from ete4 import Tree
import sys
from ete4.smartview import Layout,  BASIC_LAYOUT
'''
(c) Code written by Nikita Kulikov(1), Iker Irisarri(2), Miguel Aranda(3); 
1) Leibniz Institute for the Analysis of Biodiversity Change (LIB), Hamburg
2) Museo Nacional de Ciencias Naturales-CSIC (MNCN-CSIC), Madrid

This script visualises extended newick trees distinguishing between speciations (blue nodes) and duplication (red nodes) events. 

usage: python3 visualise_newick.py  extended_newick_tree

2025
'''
tree_file = sys.argv[1]

def color_events(node):
    if not node.is_leaf:
        if 'evoltype' in node.props:
            if node.props["evoltype"] == 'S':
                yield {'dot': {'fill': 'blue', 'radius': 5}}
            elif node.props["evoltype"] == 'D':
                yield {'dot': {'fill': 'red', 'radius': 5}}


tree = Tree(tree_file, parser=1)

events_ly = Layout(name='evol_event', draw_node=color_events)
tree.explore(layouts=[BASIC_LAYOUT, events_ly])
input()