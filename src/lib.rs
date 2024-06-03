use osmpbfreader::{NodeId, OsmObj};
use itertools::Itertools;
use geographiclib_rs::{Geodesic, InverseGeodesic};
use std::{cmp::Reverse, collections::{BinaryHeap, HashMap, HashSet}, path::Path};
use rand::{rngs::ThreadRng, seq::IteratorRandom};


pub struct RoadNetwork {
    graph: HashMap<NodeId, HashMap<NodeId, u32>>,
    nodes: HashMap<NodeId, (f64, f64)>
}

impl RoadNetwork {
    fn new() -> RoadNetwork {
        RoadNetwork{ graph: HashMap::new(), nodes: HashMap::new() }
    }

    pub fn read_from_osm_file<P: AsRef<Path>>(path: P) -> Result<RoadNetwork, std::io::Error> {
        let f = std::fs::File::open(path).unwrap();
        let mut pbf = osmpbfreader::OsmPbfReader::new(f);
        let mut rn = RoadNetwork::new();
        pbf.iter().map(Result::unwrap).for_each(|obj| {
            match obj {
                OsmObj::Node(node) => {
                    rn.nodes.insert(node.id, (node.lat(), node.lon()));
                },
                OsmObj::Way(way) => {
                    if let Some(v) = way.tags.into_inner().get("highway") {
                        if let Some(road_type) = RoadNetwork::classify_road(v) {
                            let speed = RoadNetwork::speed(road_type);
                            for (node_1, node_2) in way.nodes.iter().tuple_windows() {
                                let cost = (rn.approx_distance(node_1, node_2).unwrap() / speed) as u32;
                                rn.graph.entry(*node_1)
                                .and_modify(|edge_list| {edge_list.insert(*node_2, cost);})
                                .or_insert({
                                    let mut hm = HashMap::new();
                                    hm.insert(*node_2, cost);
                                    hm
                                });
                                rn.graph
                                .entry(*node_2)
                                .and_modify(|edge_list| {edge_list.insert(*node_1, cost);})
                                .or_insert({
                                    let mut hm = HashMap::new();
                                    hm.insert(*node_1, cost);
                                    hm
                                });
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
    pub fn distance(&self, id_1: &NodeId, id_2: &NodeId) -> Option<f64> {
        let g = Geodesic::wgs84();
        if let (Some(p_1), Some(p_2)) = (self.nodes.get(id_1), self.nodes.get(id_2)) {
            Some(g.inverse(p_1.0, p_1.1, p_2.0, p_2.1))
        } else {
            None
        }
    }

    // approximate distance relevant to germany
    fn approx_distance(&self, id_1: &NodeId, id_2: &NodeId) -> Option<f64> {
        if let (Some(p_1), Some(p_2)) = (self.nodes.get(id_1), self.nodes.get(id_2)) {
            let diff_lat = (p_1.0 - p_2.0) * 111_229.0;
            let diff_lon = (p_1.1 - p_2.1) * 71_695.0;
            Some(f64::sqrt(diff_lat * diff_lat + diff_lon * diff_lon))
        } else {
            None
        }
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
        let mut source_nodes: HashSet<NodeId> = self.graph.keys().copied().collect();
        while !source_nodes.is_empty() {
            let source_node = source_nodes.iter().next().unwrap();
            d.compute_shortest_path(*source_node, None, Some(source_node.0 as u64));
            let v = d.visited_nodes.keys().copied().collect();
            source_nodes = source_nodes.difference(&v).copied().collect();
            if v.len() > source_nodes.len() {
                break;
            }
        }
        let mut node_hist = HashMap::new();
        d.visited_nodes.values().for_each(|x| {node_hist.entry(x)
            .and_modify(|count| {*count += 1;})
            .or_insert(1);});
        let lcc_node = NodeId(**node_hist.iter().max_by(|a, b| a.1.cmp(b.1)).map(|(k, _v)| k).unwrap() as i64);
        d.visited_nodes.clear();
        d.compute_shortest_path(lcc_node, None, Some(lcc_node.0 as u64));
        let visited = d.visited_nodes;
        self.graph.retain(|k, _v| visited.contains_key(k));
        self.nodes.retain(|k, _v| visited.contains_key(k));
        self.graph.shrink_to_fit();
        self.nodes.shrink_to_fit();
    }
}

impl std::fmt::Display for RoadNetwork {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "# of Nodes: {}, # of Arcs: {}", self.nodes.len(), 
            self.graph.values().map(|v| v.len()).sum::<usize>() / 2)
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
    visited_nodes: HashMap<NodeId, u64>,
    num_settled_nodes: usize,
    heuristic: Option<HashMap<NodeId, u64>>,
}

impl DijkstrasAlgorithm<'_> {
    pub fn new(rn: &RoadNetwork) -> DijkstrasAlgorithm {
        DijkstrasAlgorithm{ rn, visited_nodes : HashMap::new(), num_settled_nodes: 0, heuristic: None}
    }

    // returns cost of shortest path to target if target exists.
    // marks visited nodes with marker if marker exists.
    pub fn compute_shortest_path(&mut self, source: NodeId, target: Option<NodeId>,
        marker: Option<u64>) -> Option<u64> {
        let mut settled_nodes = HashSet::new();
        let mut pq = BinaryHeap::new(); // defaults to max-heap
        let mut node_costs = HashMap::<NodeId, u64>::new();
        let mut h = 0;
        if self.heuristic.is_some() {
            h = *self.heuristic.as_ref().unwrap().get(&source).unwrap();
        }
        pq.push((Reverse(h), source));
        node_costs.insert(source, 0);
        if let Some(marker) = marker {
            self.visited_nodes.insert(source, marker);
        }
        while let Some((_, closest_node)) = pq.pop() {
            if settled_nodes.contains(&closest_node) {
                continue; // no point going back over a settled node
            }
            settled_nodes.insert(closest_node);
            if target.is_some_and(|target| closest_node == target) { // found target
                self.num_settled_nodes = settled_nodes.len();
                return Some(*node_costs.get(&closest_node).unwrap());
            }
            if let Some(edges) = self.rn.graph.get(&closest_node) {
                for (dest, cost) in edges {
                    if settled_nodes.contains(dest) {
                        continue; // no point touching a settled node
                    }
                    // marking visited nodes
                    if let Some(marker) = marker {
                        self.visited_nodes.insert(*dest, marker);
                    }
                    let cost_from_closest = node_costs.get(&closest_node).unwrap() + *cost as u64;
                    if let Some(current_best_cost) = node_costs.get(dest) {
                        if current_best_cost <= &cost_from_closest {
                            continue; // can't relax edge
                        }
                    }
                    if self.heuristic.is_some() {
                        h = *self.heuristic.as_ref().unwrap().get(dest).unwrap();
                    } else {
                        h = 0;
                    }
                    pq.push((Reverse(cost_from_closest + h), *dest));
                    node_costs.insert(*dest, cost_from_closest);
                }
            }
        }
        self.num_settled_nodes = settled_nodes.len();
        None
    }

    pub fn set_heuristic(&mut self, h: HashMap<NodeId, u64>) {
        self.heuristic = Some(h);
    }

    pub fn simple_heuristic(&self, target: NodeId) -> HashMap<NodeId, u64> {
        let mut h = HashMap::new();
        for id in self.rn.nodes.keys() {
            h.insert(*id, (self.rn.approx_distance(id, &target).unwrap() * 3600.0 / 110_000.0) as u64);
        }
        h
    }
}

pub struct LandmarkAlgorithm<'a> {
    rn: &'a RoadNetwork,
    landmarks: Vec<NodeId>,
    costs: HashMap<NodeId, Vec<u64>>,
    rng: ThreadRng,
}

impl LandmarkAlgorithm<'_> {
    pub fn new(rn: &RoadNetwork) -> LandmarkAlgorithm<'_> {
        LandmarkAlgorithm {
            rn,
            landmarks: Vec::new(),
            costs: HashMap::new(),
            rng: rand::thread_rng(),
        }
    }

    pub fn select_landmarks(&mut self, n: usize) {
        self.landmarks.clear();
        self.landmarks = self.rn.nodes.keys().choose_multiple(&mut self.rng, n).iter().map(|x| *(*x)).collect();
    }

    // simplified and slightly modified to give all distances from a single source
    fn dijkstra(&self, source: NodeId) -> HashMap<NodeId, u64> {
        let mut settled_nodes = HashSet::new();
        let mut pq: BinaryHeap<(Reverse<u64>, NodeId)> = BinaryHeap::new(); // defaults to max-heap
        let mut node_costs = HashMap::<NodeId, u64>::new();
        pq.push((Reverse(0), source));
        node_costs.insert(source, 0);
        while let Some((_, closest_node)) = pq.pop() {
            if settled_nodes.contains(&closest_node) {
                continue; // no point going back over a settled node
            }
            settled_nodes.insert(closest_node);
            if let Some(edges) = self.rn.graph.get(&closest_node) {
                for (dest, cost) in edges {
                    if settled_nodes.contains(dest) {
                        continue; // no point touching a settled node
                    }
                    let cost_from_closest = node_costs.get(&closest_node).unwrap() + *cost as u64;
                    if let Some(current_best_cost) = node_costs.get(dest) {
                        if current_best_cost <= &cost_from_closest {
                            continue; // can't relax edge
                        }
                    }
                    pq.push((Reverse(cost_from_closest), *dest));
                    node_costs.insert(*dest, cost_from_closest);
                }
            }
        }
        node_costs
    }

    pub fn precompute_landmark_distances(&mut self) {
        self.costs.clear();
        for l in self.landmarks.iter() {
            for (n, cost) in self.dijkstra(*l) {
                self.costs.entry(n)
                .and_modify(|x| {x.push(cost)})
                .or_insert(vec![cost]);
            }
            
        }
    }

    pub fn landmark_heuristic(&self, target: NodeId) -> HashMap<NodeId, u64> {
        let mut h = HashMap::<NodeId, u64>::new();
        let l_t_dists: &Vec<u64> = self.costs.get(&target).unwrap();
        for (node, _) in self.rn.nodes.iter() {
            h.insert(*node, self.costs.get(node).unwrap().iter().zip(l_t_dists.iter()).map(|(x, y)| x.abs_diff(*y)).max().unwrap());
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
    fn saarland() {
        let now = Instant::now();
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
        let elapsed_time = now.elapsed();
        println!("Reading Saarland OSM took {} s", elapsed_time.as_secs_f32());
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 1_119_289);
        //assert_eq!(rn.graph.iter().map(|(_, v)| v.len()).sum::<usize>() / 2, 227_826);
        let now = Instant::now();
        rn.reduce_to_largest_connected_component();
        let elapsed_time = now.elapsed();
        println!("Saarland LCC reduction took {} s", elapsed_time.as_secs_f32());
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 213_567);
        // assert_eq!(rn.graph.iter().map(|(_, v)| v.len()).sum::<usize>() / 2, 225_506);
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time: Duration = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut d = DijkstrasAlgorithm::new(&rn);
        for (src, dst) in rn.nodes.keys().collect::<Vec<_>>().iter().choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            total_cost += d.compute_shortest_path(***src, Some(***dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("---------- Dijkstra (S) -----------");
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
        total_elapsed_time = Duration::ZERO;
        total_cost = 0;
        total_settled_nodes = 0;
        let mut total_heuristic_calc_time = Duration::ZERO;
        for (src, dst) in rn.nodes.keys().collect::<Vec<_>>().iter().choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            d.set_heuristic(d.simple_heuristic(***dst));
            total_heuristic_calc_time += now.elapsed();
            let now = Instant::now();
            total_cost += d.compute_shortest_path(***src, Some(***dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("----------  A-star (S) -----------");
        println!("Average heuristic calculation time: {} s", total_heuristic_calc_time.as_secs_f32() / 100.0);
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
        total_elapsed_time = Duration::ZERO;
        total_cost = 0;
        total_settled_nodes = 0;
        total_heuristic_calc_time = Duration::ZERO;
        let mut landmarks_precompute_time = Duration::ZERO;
        let now = Instant::now();
        let mut l = LandmarkAlgorithm::new(&rn);
        l.select_landmarks(42);
        l.precompute_landmark_distances();
        landmarks_precompute_time += now.elapsed();
        for (src, dst) in rn.nodes.keys().collect::<Vec<_>>().iter().choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            d.set_heuristic(l.landmark_heuristic( ***dst));
            total_heuristic_calc_time += now.elapsed();
            let now = Instant::now();
            total_cost += d.compute_shortest_path(***src, Some(***dst), None).unwrap();
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
    fn baden_wuerttemberg() {
        let now = Instant::now();
        let mut rn = RoadNetwork::read_from_osm_file("rsrc/baden-wuerttemberg.osm.pbf").unwrap();
        let elapsed_time = now.elapsed();
        println!("Reading BW OSM took {} s", elapsed_time.as_secs_f32());
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 14_593_458);
        //assert_eq!(rn.graph.iter().map(|(_, v)| v.len()).sum::<usize>() / 2, 2_642_949);
        let now = Instant::now();
        rn.reduce_to_largest_connected_component();
        let elapsed_time = now.elapsed();
        println!("BW LCC reduction took {} s", elapsed_time.as_secs_f32());
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 2_458_230);
        //assert_eq!(rn.graph.iter().map(|(_, v)| v.len()).sum::<usize>() / 2, 2_613_338);
        let mut rng = &mut rand::thread_rng();
        let mut total_elapsed_time: Duration = Duration::ZERO;
        let mut total_cost = 0;
        let mut total_settled_nodes = 0;
        let mut d = DijkstrasAlgorithm::new(&rn);
        for (src, dst) in rn.nodes.keys().collect::<Vec<_>>().iter().choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            total_cost += d.compute_shortest_path(***src, Some(***dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("---------- Dijkstra (BW) -----------");
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
        total_elapsed_time = Duration::ZERO;
        total_cost = 0;
        total_settled_nodes = 0;
        let mut total_heuristic_calc_time = Duration::ZERO;
        for (src, dst) in rn.nodes.keys().collect::<Vec<_>>().iter().choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            d.set_heuristic(d.simple_heuristic(***dst));
            total_heuristic_calc_time += now.elapsed();
            let now = Instant::now();
            total_cost += d.compute_shortest_path(***src, Some(***dst), None).unwrap();
            total_elapsed_time += now.elapsed();
            total_settled_nodes += d.num_settled_nodes;
        }
        println!("----------  A-star (BW) -----------");
        println!("Average heuristic calculation time: {} s", total_heuristic_calc_time.as_secs_f32() / 100.0);
        println!("Average query time: {} s", total_elapsed_time.as_secs_f32() / 100.0);
        println!("Average cost: {}", total_cost / 100);
        println!("Average settled nodes: {}", total_settled_nodes / 100);
        total_elapsed_time = Duration::ZERO;
        total_cost = 0;
        total_settled_nodes = 0;
        total_heuristic_calc_time = Duration::ZERO;
        let mut landmarks_precompute_time = Duration::ZERO;
        let now = Instant::now();
        let mut l = LandmarkAlgorithm::new(&rn);
        l.select_landmarks(42);
        l.precompute_landmark_distances();
        landmarks_precompute_time += now.elapsed();
        for (src, dst) in rn.nodes.keys().collect::<Vec<_>>().iter().choose_multiple(&mut rng, 200).iter().tuples() {
            let now = Instant::now();
            d.set_heuristic(l.landmark_heuristic( ***dst));
            total_heuristic_calc_time += now.elapsed();
            let now = Instant::now();
            total_cost += d.compute_shortest_path(***src, Some(***dst), None).unwrap();
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
