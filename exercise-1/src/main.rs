use osmpbfreader::{OsmObj, NodeId};
use itertools::Itertools;
use geographiclib_rs::{Geodesic, InverseGeodesic};
use std::{collections::HashMap, path::Path};


struct RoadNetwork {
    graph: HashMap<NodeId, HashMap<NodeId, f64>>,
    nodes: HashMap<NodeId, (f64, f64)>
}

impl<'a> RoadNetwork {
    fn new() -> RoadNetwork {
        RoadNetwork{ graph: HashMap::new(), nodes: HashMap::new() }
    }

    fn read_from_osm_file<P: AsRef<Path>>(path: P) -> Result<RoadNetwork, std::io::Error> {
        let f = std::fs::File::open(path).unwrap();
        let mut pbf = osmpbfreader::OsmPbfReader::new(f);
        let mut rn = RoadNetwork::new();
        let _ = pbf.iter().map(Result::unwrap).for_each(|obj| {
            match obj {
                OsmObj::Node(node) => {
                    rn.graph.insert(node.id, HashMap::new());
                    rn.nodes.insert(node.id, (node.lat(), node.lon()));
                },
                OsmObj::Way(way) => {
                    if let Some(v) = way.tags.into_inner().get("highway") {
                        if let Some(road_type) = RoadNetwork::classify_road(v) {
                            let speed = RoadNetwork::speed(road_type);
                            for (node_1, node_2) in way.nodes.iter().tuple_windows() {
                                let cost = rn.distance(node_1, node_2).unwrap() / speed;
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
    fn distance(&self, id_1: &NodeId, id_2: &NodeId) -> Option<f64> {
        let g = Geodesic::wgs84();
        if let (Some(p_1), Some(p_2)) = (self.nodes.get(id_1), self.nodes.get(id_2)) {
            Some(g.inverse(p_1.0, p_1.1, p_2.0, p_2.1))
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
}

impl std::fmt::Display for RoadNetwork {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "# of Nodes: {}, # of Arcs: {}", self.nodes.len(), 
            self.graph.iter().map(|(_, v)| v.len()).sum::<usize>() / 2)
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

fn main() {
    let rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
    println!("{rn}");
}

#[cfg(test)]
mod tests {
    use crate::RoadNetwork;

    #[test]
    fn saarland_from_osm() {
        let rn = RoadNetwork::read_from_osm_file("rsrc/saarland.osm.pbf").unwrap();
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 1_119_289);
        //assert_eq!(rn.graph.iter().map(|(_, v)| v.len()).sum::<usize>() / 2, 227_826);
    }

    #[test]
    fn baden_wuerttemberg_from_osm() {
        let rn = RoadNetwork::read_from_osm_file("rsrc/baden-wuerttemberg.osm.pbf").unwrap();
        println!("{rn}");
        assert_eq!(rn.nodes.len(), 14_593_458);
        //assert_eq!(rn.graph.iter().map(|(_, v)| v.len()).sum::<usize>() / 2, 2_642_949);
    }
}
