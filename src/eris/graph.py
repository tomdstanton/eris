"""
Copyright 2025 Tom Stanton (tomdstanton@gmail.com)
https://github.com/tomdstanton/eris

This file is part of eris. eris is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. eris is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with eris.
If not, see <https://www.gnu.org/licenses/>.
"""
from typing import Literal, Union, Any, List, Optional, Tuple, Set, Dict
from collections import defaultdict, deque
from heapq import heappush, heappop
from warnings import warn

from eris import ErisWarning

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

    def get_effective_weight(self) -> float:
        """Returns the edge weight, defaulting to 1.0 if None (for BFS-like behavior in Dijkstra)."""
        return 1.0 if self.weight is None else self.weight

    def reverse(self) -> 'Edge':
        """Returns a new Edge object representing the reverse direction."""
        return Edge(self.to, self.fr, self.to_attribute, self.fr_attribute, self.weight)

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


class GraphError(Exception):
    pass


class GraphWarning(ErisWarning):
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
        # TODO: Add dead end count
        # Set of unique Edge objects fundamentally added to the graph.
        self.edges: Set[Edge] = set()
        self._node_ids: Set[str] = set()
        self.directed: bool = directed

        for edge in edges:
            self.add_edge(edge)
    
    def __repr__(self):
        directionality = "Directed" if self.directed else "Undirected"
        # Note: len(self.edges) counts only the *unique* edge objects added,
        # not the total number of traversable connections in the undirected case.
        return f"{directionality} Graph with {len(self._node_ids)} nodes and {len(self.edges)} defined edges"

    def __str__(self):
        return self.__repr__()

    def add_node(self, node_id: str):
        """Adds a node to the graph."""
        self._node_ids.add(node_id)

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

        # Add the primary edge representation if it's new
        # This helps track the originally added edges vs implicit reverse ones
        is_new_edge = edge not in self.edges
        if is_new_edge:
             self.edges.add(edge)

        # Add forward connectivity to the adjacency list
        # Check if this specific edge object is already in the adjacency list for the 'fr' node
        if edge not in self.adj[edge.fr]:
            self.adj[edge.fr].add(edge)

        # If the graph is undirected, add reverse connectivity as well
        if not self.directed:
            # Create a conceptual reverse edge for traversal
            reverse_edge = edge.reverse()
            # Add the reverse edge to the adjacency list of the 'to' node
            # Check if this specific reverse edge object is already in the adj list for the 'to' node
            if reverse_edge not in self.adj[edge.to]:
                self.adj[edge.to].add(reverse_edge)
            # Note: We do not add the reverse_edge to self.edges unless it's explicitly added later by the user

    def get_neighbors(self, node_id: str) -> Set[Edge]:
        """
        Returns the set of outgoing edges for a given node ID, respecting
        graph directionality (for undirected graphs, this includes edges
        allowing traversal back along an added edge).
        """
        return self.adj.get(node_id, set())

    def get_nodes(self) -> Set[str]:
        """Returns a set of all unique node IDs in the graph."""
        return self._node_ids

    # --- Pathfinding Methods ---
    # These methods remain largely the same, as they rely on get_neighbors,
    # which now correctly provides traversable edges based on graph directionality.

    def _find_shortest_path_bfs(self, start_node: str, end_node: str) -> Optional[Tuple[List[str], int]]:
        """
        Finds the shortest path (fewest edges) using BFS.
        Returns the path list and the number of edges, or None.
        """
        if start_node not in self._node_ids or end_node not in self._node_ids: return None
        if start_node == end_node: return [start_node], 0

        queue = deque([(start_node, [start_node])]) # Store (node, path_so_far)
        visited = {start_node}

        while queue:
            current_node, path = queue.popleft()

            # Use get_neighbors which respects graph directionality
            for edge in self.get_neighbors(current_node):
                neighbor = edge.to # The destination node of this edge
                if neighbor == end_node:
                    final_path = path + [neighbor]
                    return final_path, len(final_path) - 1 # Return path and edge count

                if neighbor not in visited:
                    visited.add(neighbor)
                    new_path = list(path)
                    new_path.append(neighbor)
                    queue.append((neighbor, new_path))
        return None # Path not found

    def _find_shortest_path_dijkstra(self, start_node: str, end_node: str) -> Optional[Tuple[List[str], float]]:
        """
        Finds the shortest path (minimum total weight) using Dijkstra.
        Uses default edge weight of 1.0 if weight is None.
        Returns the path list and the total weight, or None.
        """
        if start_node not in self._node_ids or end_node not in self._node_ids: return None
        if start_node == end_node: return [start_node], 0.0

        distances: Dict[str, float] = {node: float('inf') for node in self._node_ids}
        previous_nodes: Dict[str, Optional[str]] = {node: None for node in self._node_ids}
        pq: List[Tuple[float, str]] = [] # Priority queue: (distance, node_id)

        distances[start_node] = 0.0
        heappush(pq, (0.0, start_node))

        path_found = False
        while pq:
            current_distance, current_node = heappop(pq)

            # Optimization: If we extract the end_node, we found the shortest path
            if current_node == end_node:
                path_found = True
                break # Found the target

            # Optimization: Skip if we found a shorter path already
            if current_distance > distances[current_node]: continue

            # Use get_neighbors which respects graph directionality
            for edge in self.get_neighbors(current_node):
                neighbor = edge.to # The destination node of this edge
                weight = edge.get_effective_weight() # Handles None weight -> 1.0
                distance_through_current = distances[current_node] + weight

                if distance_through_current < distances[neighbor]:
                    distances[neighbor] = distance_through_current
                    previous_nodes[neighbor] = current_node
                    heappush(pq, (distance_through_current, neighbor))

        # Reconstruct path if end_node was reached
        if not path_found or distances[end_node] == float('inf'):
            return None # Path not found

        path, node = [], end_node
        while node is not None:
            path.append(node)
            node = previous_nodes[node] # Move to the predecessor
        if not path or path[-1] != start_node:
             return None # Should not happen if path_found is True, but safety check

        return path[::-1], distances[end_node] # Return reversed path and total weight

    def shortest_path(self, node_ids: List[str], method: Literal['bfs', 'dijkstra'] = 'bfs',
                      max_cost: float = float('inf')) -> Optional[Tuple[List[str], float]]:
        """
        Calculates the shortest path sequentially connecting a list of node IDs,
        using either BFS (fewest edges) or Dijkstra (lowest weight).

        Args:
            node_ids: A list of node IDs (strings) in the desired order (at least 2).
            method: The algorithm to use: 'bfs' (default) or 'dijkstra'.
            max_cost: Returns None if max_cost is exceeded

        Returns:
            A tuple containing:
            - The full path as a list of node IDs.
            - The total cost (either edge count for BFS or weight sum for Dijkstra).
            Returns None if the input list has fewer than 2 nodes, if any node
            is not in the graph, if the method is invalid, or if a path segment
            does not exist between any consecutive pair of nodes.

        Raises:
            ValueError: If an invalid method is specified.
            GraphError: Can be raised internally if nodes aren't found (though typically returns None).
        """
        if not node_ids or len(node_ids) < 2:
             warn("shortest_path requires at least two nodes.", GraphWarning)
             return None

        if  method_lower := method.lower() not in ['bfs', 'dijkstra']:  # Validate method choice
            raise ValueError(f"Invalid method '{method}'. Choose 'bfs' or 'dijkstra'.")

        for node_id in node_ids:  # Check if all nodes exist in the graph first
            if node_id not in self._node_ids:
                raise GraphError(f"Node '{node_id}' not found in the graph.")

        full_path: List[str] = []
        total_cost: float = 0.0

        for i in range(len(node_ids) - 1):   # Iterate through consecutive pairs of nodes
            start_node, end_node = node_ids[i], node_ids[i+1]

            segment_result: Optional[Tuple[List[str], Union[int, float]]] = None

            # Call the chosen pathfinding algorithm
            if method_lower == 'bfs':
                segment_result = self._find_shortest_path_bfs(start_node, end_node)
            elif method_lower == 'dijkstra':
                segment_result = self._find_shortest_path_dijkstra(start_node, end_node)

            if segment_result is None:  # Check if path segment was found
                return None
                # raise GraphError(f"No path found between '{start_node}' and '{end_node}' using {method}.")

            segment_path, segment_cost = segment_result

            # Append the segment to the full path
            # Avoid duplicating the start node of the next segment if path isn't empty
            full_path.extend(segment_path if not full_path else segment_path[1:])  # Add segment path excluding start node

            # Add the cost (edge count or weight) to the total
            total_cost += float(segment_cost) # Cast int cost from BFS to float

            if total_cost > max_cost:
                return None

        return full_path, total_cost
