"""
Module for handling simple networks, abolishing reliance on networkx.
"""
from typing import Literal, Any, Optional, Set, Dict
from collections import defaultdict

# Classes --------------------------------------------------------------------------------------------------------------
class Edge:
    """
    General class to link objects together for graphs, supporting optional weights and attributes (such as strand).
    Objects must have an ``id`` attribute or be strings.
    Note, this class intentionally holds references to nodes rather than the instances themselves.
    Represents a directed connection from 'fr' to 'to'.
    """
    def __init__(self, fr: Any, to: Any, fr_attribute: Any = None, to_attribute: Any = None,
                 weight: Optional[float] = None): # Weight is optional
        # Ensure node IDs are strings
        self.fr = fr if isinstance(fr, str) else getattr(fr, 'id', str(fr)) # type: str
        self.to = to if isinstance(to, str) else getattr(to, 'id', str(to))   # type: str
        self.fr_attribute = fr_attribute
        self.to_attribute = to_attribute
        # Store the weight; it can be None if not provided
        self.weight = weight

    def __repr__(self):
        weight_str = f" (w={self.weight})" if self.weight is not None else ""
        # Represent the fundamental direction of the edge object
        return f"Edge({self.fr} -> {self.to}{weight_str})"

    def __str__(self):
        return self.__repr__()

    def __format__(self, __format_spec: Literal['gfa', 'tsv'] = ''):
        if __format_spec == '':
            return self.__repr__()
        elif __format_spec == 'gfa':
            fr_orient = self.fr_attribute if self.fr_attribute in ['+', '-'] else '+'
            to_orient = self.to_attribute if self.to_attribute in ['+', '-'] else '+'
            # Optional weight tag could be added: f" W:{self.weight}"? Standard is less clear.
            return f"L\t{self.fr}\t{fr_orient}\t{self.to}\t{to_orient}\t*\n"
        elif __format_spec == 'tsv':
            weight_str = f"\t{self.weight}" if self.weight is not None else "\t" # Add empty tab if no weight
            return f"{self.fr}\t{self.to}{weight_str}\n"
        else:
            raise NotImplementedError(f'Format "{__format_spec}" not supported')

    # Equality and hash methods include weight
    def __eq__(self, other):
        if not isinstance(other, Edge):
            return NotImplemented
        # Direction matters for edge equality itself
        return (self.fr == other.fr and self.to == other.to and
                self.fr_attribute == other.fr_attribute and
                self.to_attribute == other.to_attribute and
                self.weight == other.weight) # Compare weight too

    def __hash__(self):
        # Hash includes weight and direction
        return hash((self.fr, self.to, self.fr_attribute, self.to_attribute, self.weight))

    def reverse(self) -> 'Edge':
        """Returns a new Edge object representing the reverse direction."""
        return Edge(self.to, self.fr, self.to_attribute, self.fr_attribute, self.weight)


class GraphError(Exception):
    pass


class Graph:
    """
    Represents a graph that can be either directed or undirected.
    Nodes are identified by string IDs.
    Edges are stored, and connectivity is managed based on the 'directed' flag.
    """
    def __init__(self, *edges: Edge, directed: bool = True):
        """
        Initializes the graph.

        Args:
            *edges: Variable number of Edge objects to initialize the graph with.
            directed: If True (default), the graph is treated as directed.
                      If False, the graph is treated as undirected, meaning adding an
                      edge A->B also allows traversal B->A with the same weight.
        """
        # Adjacency list: maps node ID to a set of outgoing Edge objects *starting* from that node.
        # For undirected graphs, this will include edges representing reverse traversal.
        self.adj: Dict[str, Set[Edge]] = defaultdict(set)
        # Set of unique Edge objects fundamentally added to the graph.
        self.edges: Set[Edge] = set()
        self._node_ids: Set[str] = set()
        self.directed: bool = directed

        for edge in edges:
            self.add_edge(edge)
    
    def __repr__(self):
        # Note: len(self.edges) counts only the *unique* edge objects added,
        # not the total number of traversable connections in the undirected case.
        return f"{'Directed' if self.directed else 'Undirected'} Graph with {len(self._node_ids)} nodes and {len(self.edges)} defined edges"

    def __str__(self):
        return self.__repr__()

    def add_edge(self, edge: Edge):
        """
        Adds an edge to the graph.

        If the graph is undirected (self.directed=False), adding edge (fr -> to)
        will allow traversal in both directions (fr -> to and to -> fr) with the
        same weight and reversed attributes. The original edge object is added
        to self.edges, but the adjacency list (self.adj) reflects reachability
        in both directions.
        """
        # Add the nodes to the set of known node IDs
        self._node_ids.add(edge.fr)
        self._node_ids.add(edge.to)

        # This helps track the originally added edges vs implicit reverse ones
        if edge not in self.edges:  # is_new_edge
             self.edges.add(edge)  # Add the primary edge representation if it's new

        # Check if this specific edge object is already in the adjacency list for the 'fr' node
        if edge not in self.adj[edge.fr]:  # Add forward connectivity to the adjacency list
            self.adj[edge.fr].add(edge)

        # If the graph is undirected, add reverse connectivity as well
        if not self.directed:
            # Create a conceptual reverse edge for traversal and add to the adjacency list of the 'to' node
            # Check if this specific reverse edge object is already in the adj list for the 'to' node
            if (reverse_edge := edge.reverse()) not in self.adj[edge.to]:
                self.adj[edge.to].add(reverse_edge)
            # Note: We do not add the reverse_edge to self.edges unless it's explicitly added later by the user

    def get_neighbors(self, node_id: str) -> Set[Edge]:
        """
        Returns the set of outgoing edges for a given node ID, respecting
        graph directionality (for undirected graphs, this includes edges
        allowing traversal back along an added edge).
        """
        return self.adj.get(node_id, set())
