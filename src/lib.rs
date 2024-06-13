use osmpbfreader::{NodeId, OsmObj};
use itertools::Itertools;
use geographiclib_rs::{Geodesic, InverseGeodesic};
use std::{cmp::Reverse, collections::{BinaryHeap, HashMap}, path::Path};
use rand::{rngs::ThreadRng, seq::IteratorRandom};

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
        let diff_lon = (n_1.lon - n_2.lat) * 71_695.0;
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
        let self.parents = vec![None; rn.nodes.len()];
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

#[cfg(test)]
mod tests {
    use crate::{DijkstrasAlgorithm, LandmarkAlgorithm, RoadNetwork};
    use std::time::{Duration, Instant};
    use rand::seq::IteratorRandom;
    use itertools::Itertools;

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
}
