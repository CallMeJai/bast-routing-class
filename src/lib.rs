use osmpbfreader::{NodeId, OsmObj};
use itertools::Itertools;
use geographiclib_rs::{Geodesic, InverseGeodesic};
use std::{cmp::Reverse, collections::{BinaryHeap, HashMap, HashSet}, path::Path};
use rand::{rngs::ThreadRng, seq::IteratorRandom, seq::SliceRandom};

#[derive(Copy, Clone, Debug)]
struct Arc {
    to_node: usize,
    cost: u32,
    arc_flag: bool,
}

impl Arc {
    fn new(to_node: usize, cost: u32, arc_flag: bool) -> Arc {
        Arc{to_node, cost, arc_flag}
    }
}

#[derive(Copy, Clone, Debug)]
struct Node {
    node_id: NodeId,
    lat: f64,
    lon: f64,
}

impl Node {
    fn new(node_id: NodeId, lat: f64, lon: f64) -> Node {
        Node{node_id, lat, lon}
    }
}

pub struct RoadNetwork {
    index_map: HashMap<NodeId, usize>,
    graph: Vec<Vec<Arc>>,
    nodes: Vec<Node>
}

impl RoadNetwork {
    fn new() -> RoadNetwork {
        RoadNetwork {
            index_map: HashMap::new(), 
            graph: Vec::new(), 
            nodes: Vec::new()
        }
    }

    pub fn read_from_osm_file<P: AsRef<Path>>(path: P) -> Result<RoadNetwork, std::io::Error> {
        let f = std::fs::File::open(path).unwrap();
        let mut pbf = osmpbfreader::OsmPbfReader::new(f);
        let mut rn = RoadNetwork::new();
        // assumes well-formed OSM file where all nodes come before all ways
        pbf.iter().map(Result::unwrap).for_each(|obj| {
            match obj {
                OsmObj::Node(node) => {
                    rn.index_map.insert(node.id, rn.nodes.len());
                    rn.nodes.push(Node::new(node.id, node.lat(), node.lon()));
                    rn.graph.push(Vec::new());
                },
                OsmObj::Way(way) => {
                    if let Some(v) = way.tags.into_inner().get("highway") {
                        if let Some(road_type) = RoadNetwork::classify_road(v) {
                            let speed = RoadNetwork::speed(road_type);
                            for (node_1, node_2) in way.nodes.iter().tuple_windows() {
                                let n_1_i = *rn.index_map.get(node_1).unwrap();
                                let n_2_i = *rn.index_map.get(node_2).unwrap();
                                let cost = (rn.approx_distance(n_1_i, n_2_i) / speed) as u32;
                                rn.graph[n_1_i].push(Arc::new(n_2_i, cost, false));
                                rn.graph[n_2_i].push(Arc::new(n_1_i, cost, false));
                            }
                        }
                    } 
                },
                _ => ()
            }
        });
        Ok(rn)
    }

    fn classify_road(v: &str) -> Option<RoadType> {
        match v {
            "motorway" => Some(RoadType::Motorway),
            "trunk" => Some(RoadType::Trunk),
            "primary" => Some(RoadType::Primary),
            "secondary" => Some(RoadType::Secondary),
            "tertiary" => Some(RoadType::Tertiary),
            "unclassified" => Some(RoadType::Unclassified),
            "residential" => Some(RoadType::Residential),
            "motorway_link" => Some(RoadType::MotorwayLink),
            "trunk_link" => Some(RoadType::TrunkLink),
            "primary_link" => Some(RoadType::PrimaryLink),
            "secondary_link" => Some(RoadType::SecondaryLink),
            "living_street" => Some(RoadType::LivingStreet),
            "service" => Some(RoadType::Service),
            "unsurfaced" => Some(RoadType::Unsurfaced),
            "road" => Some(RoadType::Road),
            _ => None
        }
    }

    // distance in m
    pub fn distance(&self, id_1: usize, id_2: usize) -> f64 {
        let g = Geodesic::wgs84();
        let (n_1, n_2) = (self.nodes[id_1], self.nodes[id_2]);
        g.inverse(n_1.lat, n_1.lon, n_2.lat, n_2.lon)
    }

    // approximate distance relevant to germany
    fn approx_distance(&self, id_1: usize, id_2: usize) -> f64 {
        let (n_1, n_2) = (self.nodes[id_1], self.nodes[id_2]);
        let diff_lat = (n_1.lat - n_2.lat) * 111_229.0;
        let diff_lon = (n_1.lon - n_2.lon) * 71_695.0;
        f64::sqrt(diff_lat * diff_lat + diff_lon * diff_lon)
    }

    // speed in m/s
    fn speed(road_type: RoadType) -> f64 {
        match road_type {
            RoadType::Motorway | RoadType::Trunk => 110_000.0 / 3600.0,
            RoadType::Primary => 70_000.0 / 3600.0,
            RoadType::Secondary => 60_000.0 / 3600.0,
            RoadType::Tertiary | RoadType::MotorwayLink | RoadType::TrunkLink | RoadType::PrimaryLink | RoadType::SecondaryLink => 50_000.0 / 3600.0,
            RoadType::Road | RoadType::Unclassified => 40_000.0 / 3600.0,
            RoadType::Residential | RoadType::Unsurfaced => 30_000.0 / 3600.0,
            RoadType::LivingStreet => 10_000.0 / 3600.0,
            RoadType::Service => 5_000.0 / 3600.0,
        }
    }

    pub fn reduce_to_largest_connected_component(&mut self) {
        let mut d = DijkstrasAlgorithm::new(self);
        let mut num_nodes_left = self.nodes.len();
        let num_edges = self.graph.iter().map(|v| v.len()).sum::<usize>() / 2;
        for source_node in 0..self.nodes.len() {
            if d.visited_nodes[source_node] != usize::MAX {
                continue;
            }
            d.compute_shortest_path(source_node, None, Some(source_node));
            num_nodes_left -= d.num_settled_nodes;
            if d.num_settled_nodes > num_nodes_left || d.num_settled_nodes > num_edges / 2 {
                break;
            }
        }
        let (&lcc_node, _) = d.visited_nodes.iter().filter(|&&marker| marker != usize::MAX)
            .counts().into_iter().max_by_key(|&(_, count)| count).unwrap();
        let index_map = d.visited_nodes.iter().enumerate()
            .filter(|(_i, &n)| n == lcc_node).enumerate()
            .map(|(new_idx, (old_idx, _lcc_node))| (self.nodes[old_idx].node_id, new_idx)).collect();
        let old_to_new_index_map: HashMap<usize, usize> = d.visited_nodes.iter().enumerate()
            .filter(|(_i, &n)| n == lcc_node).enumerate()
            .map(|(new_idx, (old_idx, _lcc_node))| (old_idx, new_idx)).collect();
        let mut nodes = vec![Node::new(NodeId(0), 0.0, 0.0); old_to_new_index_map.len()];
        for (i, val) in self.nodes.iter().enumerate() {
            if let Some(&idx) = old_to_new_index_map.get(&i) {
                nodes[idx] = *val;
            }
        }
        let mut graph: Vec<Vec<Arc>> = vec![Vec::new(); nodes.len()];
        for (i, val) in self.graph.iter().enumerate() {
            if let Some(&idx) = old_to_new_index_map.get(&i) {
                graph[idx] = val.iter().map(
                    |Arc{to_node: i, cost, arc_flag: flag}| 
                    Arc::new(
                        *old_to_new_index_map.get(i).unwrap(), 
                        *cost, 
                        *flag
                    )
                ).collect();
            }
        }
        self.index_map = index_map;
        self.graph = graph;
        self.nodes = nodes;
    }
}

impl std::fmt::Display for RoadNetwork {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "# of Nodes: {}, # of Arcs: {}", self.nodes.len(), 
            self.graph.iter().map(|v| v.len()).sum::<usize>() / 2)
    }
}

enum RoadType {
    Motorway,
    Trunk,
    Primary,
    Secondary,
    Tertiary,
    MotorwayLink,
    TrunkLink,
    PrimaryLink,
    SecondaryLink,
    Road,
    Unclassified,
    Residential,
    Unsurfaced,
    LivingStreet,
    Service
}

pub struct DijkstrasAlgorithm<'a> {
    rn: &'a RoadNetwork,
    visited_nodes: Vec<usize>,
    num_settled_nodes: usize,
    heuristic: Vec<u64>,
    consider_arc_flags: bool,
    parents: Vec<Option<usize>>,
}

impl DijkstrasAlgorithm<'_> {
    pub fn new(rn: &RoadNetwork) -> DijkstrasAlgorithm {
        DijkstrasAlgorithm{ 
            rn, 
            visited_nodes : vec![usize::MAX; rn.nodes.len()], 
            num_settled_nodes: 0, 
            heuristic: vec![0; rn.nodes.len()],
            consider_arc_flags: false,
            parents: vec![None; rn.nodes.len()],
        }
    }

    // returns cost of shortest path to target if target exists.
    // marks visited nodes with marker if marker exists.
    pub fn compute_shortest_path(&mut self, source: usize, target: Option<usize>,
        marker: Option<usize>) -> Option<u64> {
        self.parents = vec![None; self.parents.len()];
        let mut settled_nodes = vec![false; self.rn.nodes.len()];
        let mut pq = BinaryHeap::new(); // defaults to max-heap
        let mut node_costs: Vec<u64> = vec![u64::MAX; self.rn.nodes.len()];
        let mut h = self.heuristic[source];
        pq.push((Reverse(h), source, None));
        node_costs[source] = 0;
        if let Some(marker) = marker {
            self.visited_nodes[source] = marker;
        }
        while let Some((_, closest_node, parent)) = pq.pop() {
            if settled_nodes[closest_node] {
                continue; // no point going back over a settled node
            }
            settled_nodes[closest_node] = true;
            self.parents[closest_node] = parent;
            if target.is_some_and(|target| closest_node == target) { // found target
                self.num_settled_nodes = settled_nodes.iter().filter(|&x| *x).count();
                return Some(node_costs[closest_node]);
            }
            let edges = &self.rn.graph[closest_node];
            for &Arc{to_node: dest, cost, arc_flag} in edges {
                if self.consider_arc_flags && !arc_flag {
                    continue; // not marked by arc flag
                }
                if settled_nodes[dest] {
                    continue; // no point touching a settled node
                }
                // marking visited nodes
                if let Some(marker) = marker {
                    self.visited_nodes[dest] = marker;
                }
                let cost_from_closest = node_costs[closest_node] + cost as u64;
                let current_best_cost = node_costs[dest];
                if current_best_cost <= cost_from_closest {
                    continue; // can't relax edge
                }
                h = self.heuristic[dest];
                pq.push((Reverse(cost_from_closest + h), dest, Some(closest_node)));
                node_costs[dest] = cost_from_closest;
            }
        }
        self.num_settled_nodes = settled_nodes.iter().filter(|&x| *x).count();
        None
    }

    pub fn set_heuristic(&mut self, h: Vec<u64>) {
        self.heuristic = h;
    }

    pub fn simple_heuristic(&self, target: usize) -> Vec<u64> {
        (0..self.rn.nodes.len())
        .map(
            |p|
            (self.rn.approx_distance(p, target) * 3600.0 / 110_000.0) as u64)
        .collect()
    }

    pub fn set_consider_arc_flags(&mut self, consider_arc_flags: bool) {
        self.consider_arc_flags = consider_arc_flags;
    } 
}

pub struct LandmarkAlgorithm<'a> {
    rn: &'a RoadNetwork,
    landmarks: Vec<usize>,
    costs: Vec<Vec<u64>>, // vec len # nodes of vec len # landmarks
    rng: ThreadRng,
}

impl LandmarkAlgorithm<'_> {
    pub fn new(rn: &RoadNetwork) -> LandmarkAlgorithm<'_> {
        LandmarkAlgorithm {
            rn,
            landmarks: Vec::new(),
            costs: vec![Vec::new(); rn.nodes.len()],
            rng: rand::thread_rng(),
        }
    }

    pub fn select_landmarks(&mut self, n: usize) {
        self.landmarks.clear();
        self.landmarks = (0..self.rn.nodes.len()).choose_multiple(&mut self.rng, n);
    }

    // simplified and modified to give all distances from a single source as a vec
    fn dijkstra(&self, source: usize) -> Vec<u64> {
        let mut settled_nodes = vec![false; self.rn.nodes.len()];
        let mut pq = BinaryHeap::new(); // defaults to max-heap
        let mut node_costs: Vec<u64> = vec![u64::MAX; self.rn.nodes.len()];
        pq.push((Reverse(0), source));
        node_costs[source] = 0;
        while let Some((_, closest_node)) = pq.pop() {
            if settled_nodes[closest_node] {
                continue; // no point going back over a settled node
            }
            settled_nodes[closest_node] = true;
            let edges = &self.rn.graph[closest_node];
            for &Arc{to_node: dest, cost, ..} in edges {
                if settled_nodes[dest] {
                    continue; // no point touching a settled node
                }
                let cost_from_closest = node_costs[closest_node] + cost as u64;
                let current_best_cost = node_costs[dest];
                if current_best_cost <= cost_from_closest {
                    continue; // can't relax edge
                }
                pq.push((Reverse(cost_from_closest), dest));
                node_costs[dest] = cost_from_closest;
            }
        }
        node_costs
    }

    pub fn precompute_landmark_distances(&mut self) {
        self.costs.fill(Vec::new());
        for l in self.landmarks.iter() {
            for (i, cost) in self.dijkstra(*l).iter().enumerate() {
                self.costs[i].push(*cost);
            }
        }
    }

    pub fn landmark_heuristic(&self, target: usize) -> Vec<u64> {
        let mut h = vec![0; self.rn.nodes.len()];
        let l_t_dists: &Vec<u64> = &self.costs[target];
        for (node, dists) in self.costs.iter().enumerate() {
            h[node] = dists.iter().zip(l_t_dists.iter()).map(|(x, y)| x.abs_diff(*y)).max().unwrap();
        }
        h
    }
}

pub struct ArcFlagsAlgorithm<'a> {
    rn: &'a mut RoadNetwork,
    lat_min: f64,
    lat_max: f64,
    lon_min: f64,
    lon_max: f64,
    num_settled_nodes: usize,
}

impl<'a> ArcFlagsAlgorithm<'a> {
    pub fn new(rn: &'a mut RoadNetwork, lat_min: f64, lat_max: f64, lon_min: f64, lon_max: f64) -> ArcFlagsAlgorithm<'a> {
        ArcFlagsAlgorithm{rn, lat_min, lat_max, lon_min, lon_max, num_settled_nodes: 0}
    }

    fn contains(&self, node: usize) -> bool {
        let Node{lat, lon, ..} = self.rn.nodes[node];
        lat >= self.lat_min && lat <= self.lat_max && lon >= self.lon_min && lon <= self.lon_max
    }

    fn compute_boundary_nodes(&self) -> Vec<usize> {
        let mut boundary_nodes = Vec::new();
        for (n, edges) in self.rn.graph.iter().enumerate() {
            if !self.contains(n) {
                continue; // n not in region
            }
            for Arc{to_node, ..} in edges {
                if !self.contains(*to_node) {
                    boundary_nodes.push(n); // n is a boundary node
                    break;
                }
            }
        }
        boundary_nodes
    }

    // simplified and modified to give vec of parent nodes
    fn dijkstra(&self, source: usize, target: Option<usize>) -> Vec<Option<usize>> {
        let mut parents = vec![None; self.rn.nodes.len()];
        let mut settled_nodes: Vec<bool> = vec![false; self.rn.nodes.len()];
        let mut pq = BinaryHeap::new(); // defaults to max-heap
        let mut node_costs: Vec<u64> = vec![u64::MAX; self.rn.nodes.len()];
        pq.push((Reverse(0), source, None));
        node_costs[source] = 0;
        while let Some((_, closest_node, parent)) = pq.pop() {
            if settled_nodes[closest_node] {
                continue; // no point going back over a settled node
            }
            settled_nodes[closest_node] = true;
            parents[closest_node] = parent;
            if target.is_some_and(|target| closest_node == target) { // found target
                return parents;
            }
            let edges = &self.rn.graph[closest_node];
            for &Arc{to_node: dest, cost, ..} in edges {
                if settled_nodes[dest] {
                    continue; // no point touching a settled node
                }
                let cost_from_closest = node_costs[closest_node] + cost as u64;
                let current_best_cost = node_costs[dest];
                if current_best_cost <= cost_from_closest {
                    continue; // can't relax edge
                }
                pq.push((Reverse(cost_from_closest), dest, Some(closest_node)));
                node_costs[dest] = cost_from_closest;
            }
        }
        parents
    }


    pub fn precompute_arc_flags(&mut self) {
        let mut indices: Vec<(usize, usize)> = Vec::new();
        for (n, edges) in self.rn.graph.iter().enumerate() {
            if !self.contains(n) {
                continue; // n not in region
            }
            for (i, Arc{to_node, ..}) in edges.iter().enumerate() {
                if self.contains(*to_node) {
                    indices.push((n, i));
                }
            }
        }
        for (i, j) in indices {
            self.rn.graph[i][j].arc_flag = true;
        }
        let boundary_nodes = self.compute_boundary_nodes();
        for &node in &boundary_nodes {
            let parents = self.dijkstra(node, None);
            for (mut node, mut parent_option) in parents.iter().enumerate() {
                while parent_option.is_some() {
                    let parent = parent_option.unwrap();
                    if self.rn.graph[node].iter().any(|arc| arc.to_node == parent && arc.arc_flag) {
                        break; // arc bw node and node's parent already marked so all ancestors should already be marked
                    }
                    for arc in self.rn.graph[node].iter_mut().filter(|arc| arc.to_node == parent) {
                        if !arc.arc_flag {
                            arc.arc_flag = true;
                        }
                    }
                    node = parent;
                    parent_option = &parents[parent];
                }
            }
        }
    }

    pub fn compute_shortest_path(&mut self, source: usize, target: usize) -> Result<u64, String> {
        if !self.contains(target) {
            Err(format!("Target node at {}, {} not in bounds [{}, {}; {}, {}]",
                self.rn.nodes[target].lat,
                self.rn.nodes[target].lon,
                self.lat_min,
                self.lon_min,
                self.lat_max,
                self.lon_max
            ))
        } else {
            let mut d = DijkstrasAlgorithm::new(self.rn);
            d.set_consider_arc_flags(true);
            let res = d.compute_shortest_path(source, Some(target), None).unwrap();
            self.num_settled_nodes = d.num_settled_nodes;
            Ok(res)
        }
    }
}

pub struct ContractionHierarchies<'a> {
    rn: &'a mut RoadNetwork,
    node_ordering: Vec<usize>,
    rng: ThreadRng,
    visited_nodes: Vec<usize>,
    settled_nodes: Vec<bool>,
    dists: Vec<u64>,
    cost_upper_bound: u64,
    max_num_settled_nodes: u64,
    order_of_node: HashMap<NodeId, usize>,
    num_settled_nodes: u64,
}

impl<'a> ContractionHierarchies<'a> {
    pub fn new(rn: &'a mut RoadNetwork) -> ContractionHierarchies<'a> {
        let num_nodes = rn.nodes.len();
        let node_ordering = (0..num_nodes).collect();
        for arc_vec in rn.graph.iter_mut() {
            for arc in arc_vec {
                arc.arc_flag = true;
            }
        }
        ContractionHierarchies{
            rn,
            node_ordering,
            rng: rand::thread_rng(),
            visited_nodes: Vec::new(),
            settled_nodes: vec![false; num_nodes],
            dists: vec![u64::MAX; num_nodes],
            cost_upper_bound: u64::MAX,
            max_num_settled_nodes: u64::MAX,
            order_of_node: HashMap::new(),
            num_settled_nodes: 0,
        }
    }

    pub fn compute_random_node_ordering(&mut self) {
        self.node_ordering.shuffle(&mut self.rng);
    }

    fn contract_node(&mut self, n: usize, dry_run: bool) -> (u32, i32) {
        let mut shortcuts = 0;
        let mut num_arcs_removed = 0;
        let mut removed_arcs_indices = Vec::new();
        let v = self.node_ordering[n];
        // finding in-/out-bound nodes
        let adjacent_nodes = self.rn.graph[v].iter().filter(|x| x.arc_flag).map(|x| x.to_node).collect::<HashSet<_>>().into_iter().collect::<Vec<_>>();
        // calculating D_ijs
        let mut d = vec![vec![0; adjacent_nodes.len()]; adjacent_nodes.len()];
        for (i, &u) in adjacent_nodes.iter().enumerate() {
            let u_cost = self.rn.graph[v].iter().filter(|arc| arc.arc_flag)
            .find(|arc| arc.to_node == u).unwrap().cost;
            for j in 0..d[i].len() {
                if j != i {
                    d[i][j] += u_cost;
                    d[j][i] += u_cost;
                }
            }
        }
        // removing arcs adjacent to v
        for (i, arc) in self.rn.graph[v].iter_mut().enumerate() {
            if arc.arc_flag {
                num_arcs_removed += 1;
                arc.arc_flag = false;
                if dry_run {
                    removed_arcs_indices.push((v, i));
                }
            }
        }
        for &n in adjacent_nodes.iter() {
            for (i, arc) in self.rn.graph[n].iter_mut().enumerate() {
                if arc.arc_flag && arc.to_node == v {
                    arc.arc_flag = false;
                    if dry_run {
                        removed_arcs_indices.push((n, i));
                    }
                }
            }
        }
        // check if and add if shortcut is necessary
        for (i, &u) in adjacent_nodes.iter().enumerate() {
            self.set_cost_upper_bound(*d[i].iter().max().unwrap() as u64);
            self.dijkstra(u);
            for (j, &w) in adjacent_nodes.iter().enumerate() {
                if self.dists[w] > d[i][j] as u64 {
                    shortcuts += 1;
                    if !dry_run {    
                        self.rn.graph[u].push(Arc{to_node: w, cost: d[i][j], arc_flag: true});
                        self.rn.graph[w].push(Arc{to_node: u, cost: d[i][j], arc_flag: true});
                    }
                }
            }
        }
        if dry_run {
            for (x, y) in removed_arcs_indices {
                self.rn.graph[x][y].arc_flag = true;
            }
        }
        (shortcuts as u32, shortcuts - num_arcs_removed)
    }

    fn dijkstra(&mut self, source: usize) {
        let mut num_settled_nodes = 0;
        let mut pq = BinaryHeap::new(); // defaults to max-heap
        for &n in self.visited_nodes.iter() {
            self.dists[n] = u64::MAX;
            self.settled_nodes[n] = false;
        }
        self.visited_nodes.clear();
        pq.push((Reverse(0), source));
        self.dists[source] = 0;
        self.visited_nodes.push(source);
        while let Some((Reverse(cost), closest_node)) = pq.pop() {
            if self.settled_nodes[closest_node] {
                continue; // no point going back over a settled node
            }
            self.settled_nodes[closest_node] = true;
            num_settled_nodes += 1;
            if cost > self.cost_upper_bound {
                break;
            }
            if num_settled_nodes >= self.max_num_settled_nodes {
                break;
            }
            let edges = &self.rn.graph[closest_node];
            for &Arc{to_node: dest, cost, arc_flag} in edges {
                if !arc_flag {
                    continue; // not marked by arc flag
                }
                if self.settled_nodes[dest] {
                    continue; // no point touching a settled node
                }
                // marking visited nodes
                self.visited_nodes.push(dest);
                let cost_from_closest = self.dists[closest_node] + cost as u64;
                let current_best_cost = self.dists[dest];
                if current_best_cost <= cost_from_closest {
                    continue; // can't relax edge
                }
                pq.push((Reverse(cost_from_closest), dest));
                self.dists[dest] = cost_from_closest;
            }
        }
        self.num_settled_nodes = num_settled_nodes;
    }


    fn set_cost_upper_bound(&mut self, cost: u64) {
        self.cost_upper_bound = cost;
    }

    pub fn set_max_num_settled_nodes(&mut self, num: u64) {
        self.max_num_settled_nodes = num;
    }

    pub fn precompute(&mut self) -> u32 {
        let mut total_shortcuts = 0;
        // compute initial edge differences
        let mut pq = BinaryHeap::new();
        for i in 0..self.node_ordering.len() {
            let (_, edge_difference) = self.contract_node(i, true);
            pq.push((Reverse(edge_difference), i));
        }
        // actually contract all nodes
        let mut contracted_nodes = HashSet::new();
        for i in 0..self.node_ordering.len() {
            if let Some((Reverse(mut prev_ed), mut node)) = pq.pop() {
                // check that hasn't already been contracted or edge difference has increased
                loop {
                    while contracted_nodes.contains(&node) {
                        (Reverse(prev_ed), node) = pq.pop().unwrap();
                    }
                    let (_, edge_difference) = self.contract_node(node, true);
                    if prev_ed < edge_difference {
                        pq.push((Reverse(edge_difference), node));
                        (Reverse(prev_ed), node) = pq.pop().unwrap();
                    } else {
                        break;
                    }
                }
                contracted_nodes.insert(node);
                self.order_of_node.insert(self.rn.nodes[self.node_ordering[node]].node_id, i);
                let (shortcuts, _) = self.contract_node(node, false);
                total_shortcuts += shortcuts;
            }
        }
        // re-flag edges for upward graph
        for i in 0..self.node_ordering.len() {
            let order = self.order_of_node.get(&self.rn.nodes[self.node_ordering[i]].node_id).unwrap();
            for edge in self.rn.graph[self.node_ordering[i]].iter_mut() {
                if self.order_of_node[&self.rn.nodes[edge.to_node].node_id] > *order && !edge.arc_flag {
                    edge.arc_flag = true;
                }
            }
        }
        total_shortcuts
    }
    
    pub fn compute_shortest_path(&mut self, source: usize, target: usize) -> Result<u64, String> {
        self.set_max_num_settled_nodes(u64::MAX);
        self.set_cost_upper_bound(u64::MAX);
        self.dijkstra(source);
        let num_settled_nodes = self.num_settled_nodes;
        let mut dists = vec![u64::MAX; self.rn.nodes.len()];
        for &n in self.visited_nodes.iter() {
            dists[n] = self.dists[n];
        }
        let visited_nodes = self.visited_nodes.clone();
        self.dijkstra(target);
        let common_nodes: Vec<usize> = visited_nodes.into_iter().filter(|n| self.visited_nodes.contains(n)).collect();
        let mut full_dists = vec![u64::MAX; dists.len()];
        for n in common_nodes.into_iter() {
            full_dists[n] = dists[n] + self.dists[n];
        }
        self.num_settled_nodes += num_settled_nodes;
        Ok(full_dists.into_iter().min().unwrap())
    }
}

#[cfg(test)]
mod tests {
    use crate::{ArcFlagsAlgorithm, ContractionHierarchies, DijkstrasAlgorithm, LandmarkAlgorithm, RoadNetwork, Node, Arc};
    use std::time::{Duration, Instant};
    use osmpbfreader::NodeId;
    use rand::seq::IteratorRandom;
    use itertools::Itertools;

    #[test]
    fn lecture_node_contraction() {
        let mut rn = RoadNetwork::new();
        for i in 0..14 {
            rn.nodes.push(Node::new(NodeId(i as i64), 0.0, 0.0));
            rn.index_map.insert(NodeId(i as i64), i);
            rn.graph.push(Vec::new())
        }

        rn.graph[1].push(Arc::new(2, 3, true));
        rn.graph[1].push(Arc::new(6, 7, true));
        rn.graph[1].push(Arc::new(13, 4, true));
        rn.graph[2].push(Arc::new(1, 3, true));
        rn.graph[2].push(Arc::new(8, 2, true));
        rn.graph[2].push(Arc::new(13, 5, true));
        rn.graph[3].push(Arc::new(4, 4, true));
        rn.graph[3].push(Arc::new(5, 5, true));
        rn.graph[3].push(Arc::new(12, 2, true));
        rn.graph[4].push(Arc::new(3, 4, true));
        rn.graph[4].push(Arc::new(7, 4, true));
        rn.graph[4].push(Arc::new(12, 3, true));
        rn.graph[5].push(Arc::new(3, 5, true));
        rn.graph[5].push(Arc::new(6, 6, true));
        rn.graph[5].push(Arc::new(11, 3, true));
        rn.graph[6].push(Arc::new(1, 7, true));
        rn.graph[6].push(Arc::new(5, 6, true));
        rn.graph[6].push(Arc::new(10, 4, true));
        rn.graph[7].push(Arc::new(4, 4, true));
        rn.graph[7].push(Arc::new(9, 7, true));
        rn.graph[7].push(Arc::new(11, 3, true));
        rn.graph[8].push(Arc::new(2, 2, true));
        rn.graph[8].push(Arc::new(9, 5, true));
        rn.graph[8].push(Arc::new(13, 2, true));
        rn.graph[9].push(Arc::new(7, 7, true));
        rn.graph[9].push(Arc::new(8, 5, true));
        rn.graph[9].push(Arc::new(10, 3, true));
        rn.graph[10].push(Arc::new(6, 4, true));
        rn.graph[10].push(Arc::new(9, 3, true));
        rn.graph[10].push(Arc::new(11, 1, true));
        rn.graph[10].push(Arc::new(13, 1, true));
        rn.graph[11].push(Arc::new(5, 3, true));
        rn.graph[11].push(Arc::new(7, 3, true));
        rn.graph[11].push(Arc::new(10, 1, true));
        rn.graph[11].push(Arc::new(12, 1, true));
        rn.graph[12].push(Arc::new(3, 2, true));
        rn.graph[12].push(Arc::new(4, 3, true));
        rn.graph[12].push(Arc::new(11, 1, true));
        rn.graph[13].push(Arc::new(1, 4, true));
        rn.graph[13].push(Arc::new(2, 5, true));
        rn.graph[13].push(Arc::new(8, 2, true));
        rn.graph[13].push(Arc::new(10, 1, true));

        let mut ch = ContractionHierarchies::new(&mut rn);
        for i in 1..14 {
            let (shortcuts, _) = ch.contract_node(i, false);
            if i != 10 && i != 11 {
                assert_eq!(shortcuts, 0);
            }
            else {
                assert_eq!(shortcuts, 1);
            }
        }
        assert_eq!(rn.graph.iter().map(|edges| edges.len()).sum::<usize>(), 46);
        assert_eq!(rn.graph[11].iter().find(|arc| arc.to_node == 13).unwrap().cost, 2);
        assert_eq!(rn.graph[13].iter().find(|arc| arc.to_node == 11).unwrap().cost, 2);
        assert_eq!(rn.graph[12].iter().find(|arc| arc.to_node == 13).unwrap().cost, 3);
        assert_eq!(rn.graph[13].iter().find(|arc| arc.to_node == 12).unwrap().cost, 3);
        println!("----- Node Contraction (L) -----");
    }

    #[test]
    fn saarland_osm() {
        let now = Instant::now();
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
        let elapsed_time = now.elapsed();
        println!("Reading Saarland OSM took {} s", elapsed_time.as_secs_f32());
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 1_119_289);
        assert_eq!(rn.graph.iter().map(|v| v.len()).sum::<usize>() / 2, 227_826);
        let now = Instant::now();
        rn.reduce_to_largest_connected_component();
        let elapsed_time = now.elapsed();
        println!("Saarland LCC reduction took {} s", elapsed_time.as_secs_f32());
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 213_567);
        assert_eq!(rn.graph.iter().map(|v| v.len()).sum::<usize>() / 2, 225_506);
    }

    #[test]
    fn saarland_dijkstra() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time: Duration = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut d = DijkstrasAlgorithm::new(&rn);
        for (src, dst) in (0..rn.nodes.len()).choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            total_cost += d.compute_shortest_path(*src, Some(*dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("---------- Dijkstra (S) -----------");
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }

    #[test]
    fn saarland_a_star() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut total_heuristic_calc_time = Duration::ZERO;
        let mut d = DijkstrasAlgorithm::new(&rn);
        for (src, dst) in (0..rn.nodes.len()).choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            d.set_heuristic(d.simple_heuristic(*dst));
            total_heuristic_calc_time += now.elapsed();
            let now = Instant::now();
            total_cost += d.compute_shortest_path(*src, Some(*dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("----------  A-star (S) -----------");
        println!("Average heuristic calculation time: {} s", total_heuristic_calc_time.as_secs_f32() / 100.0);
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }

    #[test]
    fn saarland_a_star_landmarks() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut total_heuristic_calc_time = Duration::ZERO;
        let mut landmarks_precompute_time = Duration::ZERO;
        let mut d = DijkstrasAlgorithm::new(&rn);
        let mut l = LandmarkAlgorithm::new(&rn);
        l.select_landmarks(42);
        let now = Instant::now();
        l.precompute_landmark_distances();
        landmarks_precompute_time += now.elapsed();
        for (src, dst) in (0..rn.nodes.len()).choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            d.set_heuristic(l.landmark_heuristic( *dst));
            total_heuristic_calc_time += now.elapsed();
            let now = Instant::now();
            total_cost += d.compute_shortest_path(*src, Some(*dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("------  A-star Landmark (S) -------");
        println!("Landmark precompute time: {} s", landmarks_precompute_time.as_secs_f32());
        println!("Average heuristic calculation time: {} s", total_heuristic_calc_time.as_secs_f32() / 100.0);
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }

    #[test]
    fn saarland_arc_flags() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut arc_flags_precompute_time = Duration::ZERO;
        let mut a = ArcFlagsAlgorithm::new(&mut rn, 49.20, 49.25, 6.95, 7.05);
        let now = Instant::now();
        a.precompute_arc_flags();
        arc_flags_precompute_time += now.elapsed();
        let mut target_nodes = Vec::new();
        while target_nodes.len() < 100 {
            let node = (0..a.rn.nodes.len()).choose(&mut rng).unwrap();
            if a.contains(node) {
                target_nodes.push(node);
            }
        }
        for (src, dst) in (0..a.rn.nodes.len()).choose_multiple(&mut rng, 100).iter().zip(target_nodes) {
            let now = Instant::now();
            total_cost += a.compute_shortest_path(*src, dst).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += a.num_settled_nodes;
        }
        println!("--------  Arc Flags (S) ---------");
        println!("Arc flags precompute time: {} s", arc_flags_precompute_time.as_secs_f32());
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }

    #[test]
    fn saarland_node_contraction() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let num_contractions = 1000;
        let mut total_contraction_time = Duration::ZERO;
        let mut shortcuts_hist = [0; 5];
        let mut edge_difference_hist = [0; 5];
        let mut total_shortcuts = 0;
        let mut ch = ContractionHierarchies::new(&mut rn);
        ch.compute_random_node_ordering();
        ch.set_max_num_settled_nodes(20);
        for i in 0..num_contractions {
            let now = Instant::now();
            let (shortcuts, edge_differences) = ch.contract_node(i, false);
            total_contraction_time += now.elapsed();
            total_shortcuts += shortcuts;
            match shortcuts {
                0 => shortcuts_hist[0] += 1,
                1 => shortcuts_hist[1] += 1,
                2 => shortcuts_hist[2] += 1,
                3 => shortcuts_hist[3] += 1,
                4..=u32::MAX => shortcuts_hist[4] += 1
            }
            match edge_differences {
                i32::MIN..=-3 => edge_difference_hist[0] += 1,
                -2 => edge_difference_hist[1] += 1,
                -1..=1 => edge_difference_hist[2] += 1,
                2 => edge_difference_hist[3] += 1,
                3..=i32::MAX => edge_difference_hist[4] += 1
            }
        }
        println!("----- Node Contraction (S) -----");
        println!("Average contraction time:  {} µs", total_contraction_time.as_secs_f64() * 1_000_000.0 / num_contractions as f64);
        println!("Shortcuts histogram: {} / {} / {} / {} / {}", shortcuts_hist[0], shortcuts_hist[1], shortcuts_hist[2], shortcuts_hist[3], shortcuts_hist[4]);
        println!("Edge difference histogram: {} / {} / {} / {} / {}", edge_difference_hist[0], edge_difference_hist[1], edge_difference_hist[2], edge_difference_hist[3], edge_difference_hist[4]);
        println!("Shortcuts per contraction: {}", total_shortcuts as f64 / num_contractions as f64);
        assert_eq!(shortcuts_hist.iter().sum::<usize>(), num_contractions);
        assert_eq!(edge_difference_hist.iter().sum::<usize>(), num_contractions);
    }

    #[test]
    fn saarland_contraction_hierarchies() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut ch = ContractionHierarchies::new(&mut rn);
        ch.set_max_num_settled_nodes(20);
        let now = Instant::now();
        let shortcuts = ch.precompute();
        let precompute_time = now.elapsed();
        let mut rng = &mut rand::thread_rng();
        let mut total_cost = 0;
        let mut total_elapsed_time = Duration::ZERO;
        let mut total_settled_nodes = 0;
        for (src, dst) in (0..ch.rn.nodes.len()).choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            total_cost += ch.compute_shortest_path(*src, *dst).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += ch.num_settled_nodes;
        }
        println!("----- Contraction Hierarchies (S) -----");
        println!("Precompute time:  {} s", precompute_time.as_secs_f64());
        println!("Total shortcuts: {}", shortcuts);
        println!("Average query time: {} s", total_elapsed_time.as_secs_f64() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }

    #[test]
    fn baden_wuerttemberg_osm() {
        let now = Instant::now();
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/baden-wuerttemberg.osm.pbf").unwrap();
        let elapsed_time = now.elapsed();
        println!("Reading BW OSM took {} s", elapsed_time.as_secs_f32());
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 14_593_458);
        assert_eq!(rn.graph.iter().map(|v| v.len()).sum::<usize>() / 2, 2_642_949);
        let now = Instant::now();
        rn.reduce_to_largest_connected_component();
        let elapsed_time = now.elapsed();
        println!("BW LCC reduction took {} s", elapsed_time.as_secs_f32());
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 2_458_230);
        assert_eq!(rn.graph.iter().map(|v| v.len()).sum::<usize>() / 2, 2_613_338);
    }

    #[test]
    fn baden_wuerttemberg_dijkstra() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/baden-wuerttemberg.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time: Duration = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut d = DijkstrasAlgorithm::new(&rn);
        for (src, dst) in (0..rn.nodes.len()).choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            total_cost += d.compute_shortest_path(*src, Some(*dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("---------- Dijkstra (BW) -----------");
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }

    #[test]
    fn baden_wuerttemberg_a_star() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/baden-wuerttemberg.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut total_heuristic_calc_time = Duration::ZERO;
        let mut d = DijkstrasAlgorithm::new(&rn);
        for (src, dst) in (0..rn.nodes.len()).choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            d.set_heuristic(d.simple_heuristic(*dst));
            total_heuristic_calc_time += now.elapsed();
            let now = Instant::now();
            total_cost += d.compute_shortest_path(*src, Some(*dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("----------  A-star (BW) -----------");
        println!("Average heuristic calculation time: {} s", total_heuristic_calc_time.as_secs_f32() / 100.0);
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }

    #[test]
    fn baden_wuerttemberg_a_star_landmarks() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/baden-wuerttemberg.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut total_heuristic_calc_time = Duration::ZERO;
        let mut landmarks_precompute_time = Duration::ZERO;
        let mut d = DijkstrasAlgorithm::new(&rn);
        let mut l = LandmarkAlgorithm::new(&rn);
        l.select_landmarks(42);
        let now = Instant::now();
        l.precompute_landmark_distances();
        landmarks_precompute_time += now.elapsed();
        for (src, dst) in (0..rn.nodes.len()).choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            d.set_heuristic(l.landmark_heuristic( *dst));
            total_heuristic_calc_time += now.elapsed();
            let now = Instant::now();
            total_cost += d.compute_shortest_path(*src, Some(*dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("------  A-star Landmark (BW) -------");
        println!("Landmark precompute time: {} s", landmarks_precompute_time.as_secs_f32());
        println!("Average heuristic calculation time: {} s", total_heuristic_calc_time.as_secs_f32() / 100.0);
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }

    #[test]
    fn baden_wuerttemberg_arc_flags() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/baden-wuerttemberg.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut arc_flags_precompute_time = Duration::ZERO;
        let mut a = ArcFlagsAlgorithm::new(&mut rn, 47.95, 48.05, 7.75, 7.90);
        let now = Instant::now();
        a.precompute_arc_flags();
        arc_flags_precompute_time += now.elapsed();
        let mut target_nodes = Vec::new();
        while target_nodes.len() < 100 {
            let node = (0..a.rn.nodes.len()).choose(&mut rng).unwrap();
            if a.contains(node) {
                target_nodes.push(node);
            }
        }
        for (src, dst) in (0..a.rn.nodes.len()).choose_multiple(&mut rng, 100).iter().zip(target_nodes) {
            let now = Instant::now();
            total_cost += a.compute_shortest_path(*src, dst).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += a.num_settled_nodes;
        }
        println!("--------  Arc Flags (BW) ---------");
        println!("Arc flags precompute time: {} s", arc_flags_precompute_time.as_secs_f32());
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }

    #[test]
    fn baden_wuerttemberg_node_contraction() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/baden-wuerttemberg.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let num_contractions = 1000;
        let mut total_contraction_time = Duration::ZERO;
        let mut shortcuts_hist = [0; 5];
        let mut edge_difference_hist = [0; 5];
        let mut ch = ContractionHierarchies::new(&mut rn);
        ch.compute_random_node_ordering();
        ch.set_max_num_settled_nodes(20);
        for i in 0..num_contractions {
            let now = Instant::now();
            let (shortcuts, edge_differences) = ch.contract_node(i, false);
            total_contraction_time += now.elapsed();
            match shortcuts {
                0 => shortcuts_hist[0] += 1,
                1 => shortcuts_hist[1] += 1,
                2 => shortcuts_hist[2] += 1,
                3 => shortcuts_hist[3] += 1,
                4..=u32::MAX => shortcuts_hist[4] += 1
            }
            match edge_differences {
                i32::MIN..=-3 => edge_difference_hist[0] += 1,
                -2 => edge_difference_hist[1] += 1,
                -1..=1 => edge_difference_hist[2] += 1,
                2 => edge_difference_hist[3] += 1,
                3..=i32::MAX => edge_difference_hist[4] += 1
            }
        }
        println!("----- Node Contraction (BW) -----");
        println!("Average contraction time:  {} µs", total_contraction_time.as_secs_f64() * 1_000_000.0 / num_contractions as f64);
        println!("Shortcuts histogram: {} / {} / {} / {} / {}", shortcuts_hist[0], shortcuts_hist[1], shortcuts_hist[2], shortcuts_hist[3], shortcuts_hist[4]);
        println!("Edge difference histogram: {} / {} / {} / {} / {}", edge_difference_hist[0], edge_difference_hist[1], edge_difference_hist[2], edge_difference_hist[3], edge_difference_hist[4]);
        assert_eq!(shortcuts_hist.iter().sum::<usize>(), num_contractions);
        assert_eq!(edge_difference_hist.iter().sum::<usize>(), num_contractions);
    }
    
    #[test]
    fn baden_wuerttemberg_contraction_hierarchies() {
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/baden-wuerttemberg.osm.pbf").unwrap();
        rn.reduce_to_largest_connected_component();
        let mut ch = ContractionHierarchies::new(&mut rn);
        ch.set_max_num_settled_nodes(20);
        let now = Instant::now();
        let shortcuts = ch.precompute();
        let precompute_time = now.elapsed();
        let mut rng = &mut rand::thread_rng();
        let mut total_cost = 0;
        let mut total_elapsed_time = Duration::ZERO;
        let mut total_settled_nodes = 0;
        for (src, dst) in (0..ch.rn.nodes.len()).choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            total_cost += ch.compute_shortest_path(*src, *dst).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += ch.num_settled_nodes;
        }
        println!("----- Contraction Hierarchies (BW) -----");
        println!("Precompute time:  {} s", precompute_time.as_secs_f64());
        println!("Total shortcuts: {}", shortcuts);
        println!("Average query time: {} s", total_elapsed_time.as_secs_f64() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
    }
}
