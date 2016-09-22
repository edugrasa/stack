

//
// Link-state policy set for Routing
//
//    Eduard Grasa <eduard.grasa@i2cat.net>
//    Vincenzo Maffione <v.maffione@nextworks.it>
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301  USA
//

#include <assert.h>
#include <climits>
#include <set>
#include <sstream>
#include <string>

#define IPCP_MODULE "routing-ps-link-state"
#include "../../ipcp-logging.h"

#include <librina/timer.h>

#include "ipcp/components.h"
#include "routing-ps.h"
#include "common/encoders/FlowStateMessage.pb.h"
#include "common/encoders/FlowStateGroupMessage.pb.h"
#include "ipcp/ipc-process.h"

namespace rinad {

std::string LinkStateRoutingPs::LINK_STATE_POLICY = "LinkState";

LinkStateRoutingPs::LinkStateRoutingPs(IRoutingComponent * rc_) : rc(rc_)
{
	lsr_policy = 0;
}

void LinkStateRoutingPs::set_dif_configuration(const rina::DIFConfiguration& dif_configuration)
{
	lsr_policy = new LinkStateRoutingPolicy(rc->ipcp);
	lsr_policy->set_dif_configuration(dif_configuration);
}

int LinkStateRoutingPs::set_policy_set_param(const std::string& name,
                                            const std::string& value)
{
        LOG_IPCP_DBG("No policy-set-specific parameters to set (%s, %s)",
                        name.c_str(), value.c_str());
        return -1;
}

LinkStateRoutingPs::~LinkStateRoutingPs()
{
	delete lsr_policy;
}

extern "C" rina::IPolicySet *
createRoutingComponentPs(rina::ApplicationEntity * ctx)
{
		IRoutingComponent * rc = dynamic_cast<IRoutingComponent *>(ctx);

		if (!rc) {
			return NULL;
		}

		return new LinkStateRoutingPs(rc);
}

extern "C" void
destroyRoutingComponentPs(rina::IPolicySet * ps)
{
        if (ps) {
                delete ps;
        }
}

Edge::Edge(unsigned int address1, unsigned int address2, int weight)
{
	address1_ = address1;
	address2_ = address2;
	weight_ = weight;
}

bool Edge::isVertexIn(unsigned int address) const
{
	if (address == address1_) {
		return true;
	}

	if (address == address2_) {
		return true;
	}

	return false;
}

unsigned int Edge::getOtherEndpoint(unsigned int address)
{
	if (address == address1_) {
		return address2_;
	}

	if (address == address2_) {
		return address1_;
	}

	return 0;
}

std::list<unsigned int> Edge::getEndpoints()
{
	std::list<unsigned int> result;
	result.push_back(address1_);
	result.push_back(address2_);

	return result;
}

bool Edge::operator==(const Edge & other) const
{
	if (!isVertexIn(other.address1_)) {
		return false;
	}

	if (!isVertexIn(other.address2_)) {
		return false;
	}

	if (weight_ != other.weight_) {
		return false;
	}

	return true;
}

bool Edge::operator!=(const Edge & other) const
{
	return !(*this == other);
}

const std::string Edge::toString() const
{
	std::stringstream ss;

	ss << address1_ << " " << address2_;
	return ss.str();
}

Graph::Graph(const std::list<FlowStateObject>& flow_state_objects)
{
	set_flow_state_objects(flow_state_objects);
}

Graph::Graph()
{
}

void Graph::set_flow_state_objects(const std::list<FlowStateObject>& flow_state_objects)
{
	flow_state_objects_ = flow_state_objects;
	init_vertices();
	init_edges();
}

Graph::~Graph()
{
	std::list<CheckedVertex *>::iterator it;
	for (it = checked_vertices_.begin(); it != checked_vertices_.end(); ++it) {
		delete (*it);
	}

	std::list<Edge *>::iterator edgeIt;
	for (edgeIt = edges_.begin(); edgeIt != edges_.end(); ++edgeIt) {
		delete (*edgeIt);
	}
}

void Graph::init_vertices()
{
	std::list<FlowStateObject>::const_iterator it;
	for (it = flow_state_objects_.begin(); it != flow_state_objects_.end();
			++it) {
		if (!contains_vertex(it->address)) {
			vertices_.push_back(it->address);
		}

		if (!contains_vertex(it->neighbor_address)) {
			vertices_.push_back(it->neighbor_address);
		}
	}
}

bool Graph::contains_vertex(unsigned int address) const
{
	std::list<unsigned int>::const_iterator it;
	for (it = vertices_.begin(); it != vertices_.end(); ++it) {
		if ((*it) == address) {
			return true;
		}
	}

	return false;
}

bool Graph::contains_edge(unsigned int address1, unsigned int address2) const
{
	for(std::list<Edge *>::const_iterator eit = edges_.begin();
					eit != edges_.end(); ++eit) {
		if (*(*eit) == Edge(address1, address2, (*eit)->weight_)) {
			return true;
		}
	}

	return false;
}

void Graph::init_edges()
{
	std::list<unsigned int>::const_iterator it;
	std::list<FlowStateObject>::iterator flowIt;

	for (it = vertices_.begin(); it != vertices_.end(); ++it) {
		checked_vertices_.push_back(new CheckedVertex((*it)));
	}

	CheckedVertex * origin = 0;
	CheckedVertex * dest = 0;
	for (flowIt = flow_state_objects_.begin();
			flowIt != flow_state_objects_.end(); ++flowIt) {
		if (!flowIt->state) {
			continue;
		}

		LOG_IPCP_DBG("Processing flow state object: %s",
				flowIt->getObjectName().c_str());

		origin = get_checked_vertex(flowIt->address);
		if (origin == 0) {
			LOG_IPCP_WARN("Could not find checked vertex for address %ud",
					flowIt->address);
			continue;
		}

		dest = get_checked_vertex(flowIt->neighbor_address);
		if (dest == 0) {
			LOG_IPCP_WARN("Could not find checked vertex for address %ud",
					flowIt->neighbor_address);
			continue;
		}

		if (origin->connection_contains_address(dest->address_)
				&& dest->connection_contains_address(origin->address_)) {
			edges_.push_back(new Edge(origin->address_,
						  dest->address_,
						  1));
			origin->connections.remove(dest->address_);
			dest->connections.remove(origin->address_);
		} else {
			origin->connections.push_back(dest->address_);
			dest->connections.push_back(origin->address_);
		}
	}
}

Graph::CheckedVertex * Graph::get_checked_vertex(unsigned int address) const
{
	std::list<Graph::CheckedVertex *>::const_iterator it;
	for (it = checked_vertices_.begin(); it != checked_vertices_.end(); ++it) {
		if ((*it)->address_ == address) {
			return (*it);
		}
	}

	return 0;
}

void Graph::print() const
{
	LOG_IPCP_INFO("Graph edges:");

	for (std::list<Edge *>::const_iterator it = edges_.begin();
					it != edges_.end(); it++) {
		const Edge& e = **it;

		LOG_IPCP_INFO("    (%u --> %u, %d)", e.address1_,
			      e.address2_, e.weight_);
	}
}

PredecessorInfo::PredecessorInfo(unsigned int nPredecessor,
				 unsigned int cost)
{
	predecessor_ = nPredecessor;
	cost_ = cost;
}

DijkstraAlgorithm::DijkstraAlgorithm()
{
}

void DijkstraAlgorithm::clear()
{
	unsettled_nodes_.clear();
	settled_nodes_.clear();
	predecessors_.clear();
	distances_.clear();
}

void DijkstraAlgorithm::computeShortestDistances(const Graph& graph,
					unsigned int source_address,
					std::map<unsigned int, int>& distances)
{
	execute(graph, source_address);

	// Write back the result
	distances = distances_;

	clear();
}

void DijkstraAlgorithm::computeRoutingTable(const Graph& graph,
			 	 	    const std::list<FlowStateObject>& fsoList,
					    unsigned int source_address,
					    std::list<rina::RoutingTableEntry *>& rt)
{
	std::list<unsigned int>::const_iterator it;
	unsigned int nextHop;
	rina::RoutingTableEntry * entry;
	unsigned int cost = 0;

	execute(graph, source_address);

	for (it = graph.vertices_.begin(); it != graph.vertices_.end(); ++it) {
		if ((*it) != source_address) {
			nextHop = getNextHop((*it), source_address, cost);
			if (nextHop != 0) {
				entry = new rina::RoutingTableEntry();
				entry->address = (*it);
				entry->nextHopAddresses.push_back(rina::NHopAltList(nextHop));
				entry->qosId = 0;
				rt.push_back(entry);
				LOG_IPCP_DBG("Added entry to routing table: destination %u, next-hop %u",
						entry->address, nextHop);
			}
		}
	}

	clear();
}

void DijkstraAlgorithm::execute(const Graph& graph, unsigned int source)
{
	distances_[source] = 0;
	unsettled_nodes_.insert(source);

	unsigned int node;
	while (unsettled_nodes_.size() > 0) {
		node = getMinimum();
		settled_nodes_.insert(node);
		unsettled_nodes_.erase(node);
		findMinimalDistances(graph, node);
	}
}

unsigned int DijkstraAlgorithm::getMinimum() const
{
	unsigned int minimum = UINT_MAX;
	std::set<unsigned int>::iterator it;

	for (it = unsettled_nodes_.begin(); it != unsettled_nodes_.end(); ++it) {
		if (minimum == UINT_MAX) {
			minimum = (*it);
		} else {
			if (getShortestDistance((*it)) < getShortestDistance(minimum)) {
				minimum = (*it);
			}
		}
	}

	return minimum;
}

int DijkstraAlgorithm::getShortestDistance(unsigned int destination) const
{
	std::map<unsigned int, int>::const_iterator it;
	int distance = INT_MAX;

	it = distances_.find(destination);
	if (it != distances_.end()) {
		distance = it->second;
	}

	return distance;
}

void DijkstraAlgorithm::findMinimalDistances(const Graph& graph,
					     unsigned int node)
{
	std::list<unsigned int> adjacentNodes;
	std::list<Edge *>::const_iterator edgeIt;

	unsigned int target = 0;
	int shortestDistance;
	for (edgeIt = graph.edges_.begin(); edgeIt != graph.edges_.end();
			++edgeIt) {
		if (isNeighbor((*edgeIt), node)) {
			target = (*edgeIt)->getOtherEndpoint(node);
			shortestDistance = getShortestDistance(node) + (*edgeIt)->weight_;
			if (getShortestDistance(target) > shortestDistance) {
				distances_[target] = shortestDistance;
				predecessors_[target] = new PredecessorInfo(node,
							                    (*edgeIt)->weight_);
				unsettled_nodes_.insert(target);
			}
		}
	}
}

bool DijkstraAlgorithm::isNeighbor(Edge * edge, unsigned int node) const
{
	if (edge->isVertexIn(node)) {
		if (!isSettled(edge->getOtherEndpoint(node))) {
			return true;
		}
	}

	return false;
}

bool DijkstraAlgorithm::isSettled(unsigned int node) const
{
	std::set<unsigned int>::iterator it;

	for (it = settled_nodes_.begin(); it != settled_nodes_.end(); ++it) {
		if ((*it) == node) {
			return true;
		}
	}

	return false;
}

unsigned int DijkstraAlgorithm::getNextHop(unsigned int target,
					   unsigned int source,
					   unsigned int & cost)
{
	std::map<unsigned int, PredecessorInfo *>::iterator it;
	PredecessorInfo * step;
	unsigned int nextHop = target;

	it = predecessors_.find(target);
	if (it == predecessors_.end()) {
		return 0;
	} else {
		step = it->second;
	}

	it = predecessors_.find(step->predecessor_);
	cost = 0;
	while (it != predecessors_.end()) {
		nextHop = step->predecessor_;
		step = it->second;
		cost++;
		if (step->predecessor_ == source) {
			break;
		}

		it = predecessors_.find(step->predecessor_);
	}

	if (step->predecessor_ == target) {
		return 0;
	}

	return nextHop;
}

// ECMP Dijkstra algorithm
ECMPDijkstraAlgorithm::ECMPDijkstraAlgorithm()
{
	t = NULL;
}

void ECMPDijkstraAlgorithm::clear()
{
	unsettled_nodes_.clear();
	settled_nodes_.clear();
	predecessors_.clear();
	distances_.clear();
	delete t;
}

void ECMPDijkstraAlgorithm::computeShortestDistances(const Graph& graph,
						     unsigned int source_address,
						     std::map<unsigned int, int>& distances)
{
	execute(graph, source_address);

	// Write back the result
	distances = distances_;

	clear();
}

void ECMPDijkstraAlgorithm::computeRoutingTable(const Graph& graph,
			 	 	    	const std::list<FlowStateObject>& fsoList,
						unsigned int source_address,
						std::list<rina::RoutingTableEntry *>& rt)
{
	std::list<unsigned int>::const_iterator it;
	std::list<unsigned int>::const_iterator nextHopsIt;
	std::list<unsigned int> nextHops;
	rina::RoutingTableEntry * entry;

	(void)fsoList; // avoid compiler barfs

	execute(graph, source_address);

	for(std::set<TreeNode *>::iterator it = t->chl.begin(); it != t->chl.end(); it++){
		std::list<rina::RoutingTableEntry *>::iterator pos = findEntry(rt,
									       (*it)->addr);

		if(pos != rt.end()){
			(*pos)->nextHopAddresses.push_back((*it)->addr);
		}
		else{
			entry = new rina::RoutingTableEntry();
		        entry->address = (*it)->addr;
		        entry->qosId = 1;
		        entry->cost = (*it)->metric;
			entry->nextHopAddresses.push_back((*it)->addr);
			LOG_IPCP_DBG("Added entry to routing table: destination %u, next-hop %u",
                        entry->address, (*it)->addr);
			rt.push_back(entry);
		}
		addRecursive(rt, 1, (*it)->addr, *it);
	}
	clear();
}

void ECMPDijkstraAlgorithm::addRecursive(std::list<rina::RoutingTableEntry *> &table,
					 int qos,
					 unsigned int next,
					 TreeNode * node)
{
	for(std::set<TreeNode *>::iterator it = node->chl.begin(); it != node->chl.end(); it++){
		std::list<rina::RoutingTableEntry *>::iterator pos = findEntry(table,
									       (*it)->addr);

		if(pos != table.end()){
			(*pos)->nextHopAddresses.push_back(next);
		}
		else{
			rina::RoutingTableEntry * entry = new rina::RoutingTableEntry();
		        entry->address = (*it)->addr;
		        entry->qosId = 1;
		        entry->cost = (*it)->metric;
			entry->nextHopAddresses.push_back(next);
			LOG_IPCP_DBG("Added entry to routing table: destination %u, next-hop %u",
                        entry->address, next);
			table.push_back(entry);
		}
		addRecursive(table, 1, next, *it);
	}
}

std::list<rina::RoutingTableEntry *>::iterator ECMPDijkstraAlgorithm::findEntry(std::list<rina::RoutingTableEntry *> &table,
										unsigned int addr)
{
	std::list<rina::RoutingTableEntry *>::iterator it;
	for(it = table.begin(); it != table.end(); it++)
	{
		if((*it)->address == addr){
			return it;
		}
	}
	return it;
}

void ECMPDijkstraAlgorithm::execute(const Graph& graph,
				    unsigned int source)
{
	distances_[source] = 0;
	settled_nodes_.insert(source);
	t = new TreeNode(source, 0);

	std::list<Edge *>::const_iterator edgeIt;
	int cost;
	unsigned int target = 0;
	int shortestDistance;
	for (edgeIt = graph.edges_.begin();
			edgeIt != graph.edges_.end(); ++edgeIt) {
		if (isNeighbor((*edgeIt), source)) {
			target = (*edgeIt)->getOtherEndpoint(source);
			distances_[target]=(*edgeIt)->weight_;
			predecessors_[target].push_front(t);
			unsettled_nodes_.insert(target);
		}
	}

	while (unsettled_nodes_.size() > 0) {
		std::set<unsigned int>::iterator it;
		getMinimum();
		for(it = minimum_nodes_.begin(); it != minimum_nodes_.end(); ++it){
			settled_nodes_.insert(*it);

			TreeNode * nt = new TreeNode(*it, distances_.find(*it)->second);
			bool fPar = true;
			for(std::list<TreeNode *>::iterator par = predecessors_.find(*it)->second.begin();
					par != predecessors_.find(*it)->second.end(); ++par){
				if(fPar){
					(*par)->chldel.insert(nt);
					fPar = false;
				}
				(*par)->chl.insert(nt);
			}

			unsettled_nodes_.erase(*it);
			findMinimalDistances(graph, nt);
		}
	}
}

void ECMPDijkstraAlgorithm::getMinimum()
{
	unsigned int minimum = UINT_MAX;
	std::set<unsigned int>::iterator it;
	minimum_nodes_.clear();
	std::list<Edge *>::const_iterator edgeIt;
	for (it = unsettled_nodes_.begin(); it != unsettled_nodes_.end(); ++it) {
		if (minimum == UINT_MAX) {
			minimum_nodes_.insert(*it);
			minimum = (*it);
		} else {
			if (getShortestDistance((*it)) < getShortestDistance(minimum)) {
								
				minimum = (*it);
				minimum_nodes_.clear();
			}
			if (getShortestDistance((*it)) == getShortestDistance(minimum)) {
								
				minimum_nodes_.insert(*it);
			}
		}
	}
}

int ECMPDijkstraAlgorithm::getShortestDistance(unsigned int destination) const
{
	std::map<unsigned int, int>::const_iterator it;
	int distance = INT_MAX;

	it = distances_.find(destination);
	if (it != distances_.end()) {
		distance = it->second;
	}

	return distance;
}

void ECMPDijkstraAlgorithm::findMinimalDistances(const Graph& graph,
						 TreeNode * pred)
{
	std::list<unsigned int> adjacentNodes;
	std::list<Edge *>::const_iterator edgeIt;
	int cost;

	unsigned int target = 0;
	int shortestDistance;
	for (edgeIt = graph.edges_.begin(); edgeIt != graph.edges_.end();
			++edgeIt) {
		if (isNeighbor((*edgeIt), pred->addr)) {
			target = (*edgeIt)->getOtherEndpoint(pred->addr);
			cost = (*edgeIt)->weight_;
			shortestDistance = getShortestDistance(pred->addr) + cost;
			if (shortestDistance < getShortestDistance(target)) {
				distances_[target] = shortestDistance;
				predecessors_[target].clear();
			}
			if (shortestDistance == getShortestDistance(target)) {
				predecessors_[target].push_front(pred);
				unsettled_nodes_.insert(target);
			}
		}
	}
}

bool ECMPDijkstraAlgorithm::isNeighbor(Edge * edge,
				       unsigned int node) const
{
	if (edge->isVertexIn(node)) {
		if (!isSettled(edge->getOtherEndpoint(node))) {
			return true;
		}
	}

	return false;
}

bool ECMPDijkstraAlgorithm::isSettled(unsigned int node) const
{
	std::set<unsigned int>::iterator it;

	for (it = settled_nodes_.begin();
			it != settled_nodes_.end(); ++it) {
		if ((*it) == node) {
			return true;
		}
	}

	return false;
}

//Class IResiliencyAlgorithm
IResiliencyAlgorithm::IResiliencyAlgorithm(IRoutingAlgorithm& ra)
						: routing_algorithm(ra)
{
}

//Class LoopFreeAlternateAlgorithm
LoopFreeAlternateAlgorithm::LoopFreeAlternateAlgorithm(IRoutingAlgorithm& ra)
						: IResiliencyAlgorithm(ra)
{
}

void LoopFreeAlternateAlgorithm::extendRoutingTableEntry(
			std::list<rina::RoutingTableEntry *>& rt,
			unsigned int target_address, unsigned int nexthop)
{
	std::list<rina::RoutingTableEntry *>::iterator rit;
	bool found = false;

	// Find the involved routing table entry
	for (rit = rt.begin(); rit != rt.end(); rit++) {
		if ((*rit)->address == target_address) {
			break;
		}
	}

	if (rit == rt.end()) {
		LOG_WARN("LFA: Couldn't find routing table entry for "
			 "target address %u", target_address);
		return;
	}

	// Assume unicast and try to extend the routing table entry
	// with the new alternative 'nexthop'
	rina::NHopAltList& altlist = (*rit)->nextHopAddresses.front();

	for (std::list<unsigned int>::iterator
			hit = altlist.alts.begin();
				hit != altlist.alts.end(); hit++) {
		if (*hit == nexthop) {
			// The nexthop is already in the alternatives
			found = true;
			break;
		}
	}

	if (!found) {
		altlist.alts.push_back(nexthop);
		LOG_DBG("Node %u selected as LFA node towards the "
			 "destination node %u", nexthop, target_address);
	}
}

void LoopFreeAlternateAlgorithm::fortifyRoutingTable(const Graph& graph,
						unsigned int source_address,
						std::list<rina::RoutingTableEntry *>& rt)
{
	std::map<unsigned int, std::map< unsigned int, int > > neighbors_dist_trees;
	std::map<unsigned int, int> src_dist_tree;

	// TODO avoid this, can be computed when invoke computeRoutingTable()
	routing_algorithm.computeShortestDistances(graph, source_address, src_dist_tree);

	// Collect all the neighbors, and for each one use the routing algorithm to
	// compute the shortest distance map rooted at that neighbor
	for (std::list<unsigned int>::const_iterator it = graph.vertices_.begin();
						it != graph.vertices_.end(); ++it) {
		if ((*it) != source_address && graph.contains_edge(source_address, *it)) {
			neighbors_dist_trees[*it] = std::map<unsigned int, int>();
			routing_algorithm.computeShortestDistances(graph,
						*it, neighbors_dist_trees[*it]);
		}
	}

	// For each node X other than than the source node
	for (std::list<unsigned int>::const_iterator it = graph.vertices_.begin();
						it != graph.vertices_.end(); ++it) {
		if ((*it) == source_address) {
			continue;
		}

		// For each neighbor of the source node, excluding X
		for (std::map<unsigned int, std::map<unsigned int, int> >::iterator
			nit = neighbors_dist_trees.begin();
				nit != neighbors_dist_trees.end(); nit++) {
			// If this neighbor is a LFA node for the current
			// destination (*it) extend the routing table to take it
			// into account
			std::map< unsigned int, int>& neigh_dist_map = nit->second;
			unsigned int neigh = nit->first;

			if (neigh == *it) {
				continue;
			}

			// dist(neigh, target) < dist(neigh, source) + dist(source, target)
			if (neigh_dist_map[*it] < src_dist_tree[neigh] + src_dist_tree[*it]) {
				//LOG_DBG("Node %u is a possible LFA for destination %u",
				//	neigh, *it);
				extendRoutingTableEntry(rt, *it, neigh);
			}
		}
	}
}

// CLASS FlowStateObject
FlowStateObject::FlowStateObject()
{
	address = 0;
	neighbor_address = 0;
	cost = 0;
	state = false;
	sequence_number = 0;
	age = 0;
	modified = false;
	avoid_port = 0;
	deprecated = true;
}
FlowStateObject::FlowStateObject(unsigned int address_,
				 unsigned int neighbor_address_,
				 unsigned int cost_,
				 bool up_,
				 int sequence_number_,
				 unsigned int age_)
{
	address = address_;
	neighbor_address = neighbor_address_;
	cost = cost_;
	state = up_;
	sequence_number = sequence_number_;
	age = age_;
	modified = true;
	deprecated = false;
	avoid_port = 0;
}

FlowStateObject::~FlowStateObject()
{
}

const std::string FlowStateObject::toString()
{
	std::stringstream ss;
	ss << "Address: " << address << "; Neighbor address: " << neighbor_address
		<< "; cost: " << cost << std::endl;
	ss << "Up: " << state << "; Sequence number: " << sequence_number
		<< "; Age: " << age;

	return ss.str();
}

FlowStateObject& FlowStateObject::operator=(const FlowStateObject& other)
{
	address = other.address;
	neighbor_address = other.neighbor_address;
	cost = other.cost;
	state = other.state;
	sequence_number = other.sequence_number;
	age = other.age;
	modified = true;
	deprecated = false;

	return *this;
}


void FlowStateObject::deprecateObject(unsigned int max_age)
{
	LOG_IPCP_DBG("Object %s deprecated", getObjectName().c_str());
	state = false;
	age = max_age + 1;
	sequence_number ++;
	modified = true;
}

std::string FlowStateObject::getObjectName() const
{
	std::stringstream ss;
	ss << FlowStateRIBObject::object_name_prefix
	   << address << "-" << neighbor_address;

	return ss.str();
}

/*std::string FlowStateObject::getKey()
{
	std::stringstream ss;
	ss << address << "-" << neighbor_address;
	return ss.str();
}*/

// CLASS FlowStateRIBObject
const std::string FlowStateRIBObject::clazz_name = "FlowStateObject";
const std::string FlowStateRIBObject::object_name_prefix = "/resalloc/fsos/key=";


FlowStateRIBObject::FlowStateRIBObject(FlowStateObjects * fsos,
		   	   	       const std::string& fqn):
		rina::rib::RIBObj(clazz_name)
{
	fsos_ = fsos;
	fqn_ = fqn;
}

void FlowStateRIBObject::read(const rina::cdap_rib::con_handle_t &con, 
	const std::string& fqn, const std::string& clas, 
	const rina::cdap_rib::filt_info_t &filt, const int invoke_id,
	rina::ser_obj_t &obj_reply, rina::cdap_rib::res_info_t& res)
{
	FlowStateObjectEncoder encoder;
	FlowStateObject result;

	if (fsos_->getObjectCopy(fqn_, result)) {
		res.code_ = rina::cdap_rib::CDAP_SUCCESS;
		encoder.encode(result, obj_reply);
	} else {
		res.code_ = rina::cdap_rib::CDAP_ERROR;
	}
}

const std::string FlowStateRIBObject::get_displayable_value() const
{
	FlowStateObject result;

	if (fsos_->getObjectCopy(fqn_, result))
		return result.toString();
	else
		return "Error retrieving object";
}

// CLASS FlowStateObjects
FlowStateObjects::FlowStateObjects(LinkStateRoutingPolicy * lsr)
{
	modified_ = false;
	lsr_ = lsr;
	rina::rib::RIBObj *rib_objects = new FlowStateRIBObjects(lsr);
	IPCPRIBDaemon* rib_daemon = (IPCPRIBDaemon*)IPCPFactory::getIPCP()
		->get_rib_daemon();
	rib_daemon->addObjRIB(FlowStateRIBObjects::object_name, &rib_objects);
}

FlowStateObjects::~FlowStateObjects()
{
	for (std::map<std::string, FlowStateObject*>::iterator it = 
		objects.begin(); it != objects.end(); ++it)
	{
		delete it->second;
	}
	objects.clear();
	IPCPRIBDaemon* rib_daemon = (IPCPRIBDaemon*)IPCPFactory::getIPCP()
		->get_rib_daemon();
	rib_daemon->removeObjRIB(FlowStateRIBObjects::object_name);
}

bool FlowStateObjects::addObject(const FlowStateObject& object)
{
	rina::ScopedLock g(lock);

	if (objects.find(object.getObjectName()) == objects.end())
	{
		addCheckedObject(object);
		return true;
	}
	else
	{
		LOG_DBG("FlowStateObject with name %s already present in the database",
			object.getObjectName().c_str());
		return false;
	}
}

void FlowStateObjects::addCheckedObject(const FlowStateObject& object)
{
	FlowStateObject * fso = new FlowStateObject(object.address,
						    object.neighbor_address,
						    object.cost,
						    object.state,
						    object.sequence_number,
						    object.age);
	objects[object.getObjectName()] = fso;
	LOG_IPCP_INFO("Create new Flow State object with name %s and pointer %p",
			fso->getObjectName().c_str(),
			fso);
	rina::rib::RIBObj* rib_obj = new FlowStateRIBObject(this, object.getObjectName());
	IPCPRIBDaemon* rib_daemon = (IPCPRIBDaemon*)IPCPFactory::getIPCP()->get_rib_daemon();
	rib_daemon->addObjRIB(fso->getObjectName(), &rib_obj);
	modified_ = true;
}

void FlowStateObjects::deprecateObject(const std::string& fqn, 
				       unsigned int max_age)
{
	rina::ScopedLock g(lock);

	std::map<std::string, FlowStateObject*>::iterator it =
			objects.find(fqn);
	if(it != objects.end())
	{
		it->second->deprecateObject(max_age);
	}
}

void FlowStateObjects::deprecateObjects(unsigned int neigh_address,
                                        unsigned int address,
		     	     	        unsigned int max_age)
{
	rina::ScopedLock g(lock);

	std::map<std::string, FlowStateObject *>::iterator it;
	for (it = objects.begin(); it != objects.end();
			++it) {
		if (it->second->neighbor_address == neigh_address &&
				it->second->address == address) {
			it->second->deprecateObject(max_age);
			modified_ = true;
		}
	}
}

void FlowStateObjects::deprecateObjectsWithAddress(unsigned int address,
						   unsigned int max_age,
						   bool neighbor)
{
	rina::ScopedLock g(lock);

	std::map<std::string, FlowStateObject *>::iterator it;
	for (it = objects.begin(); it != objects.end();
			++it) {
		if (!neighbor && it->second->address == address) {
			it->second->deprecateObject(max_age);
			modified_ = true;
		} else if (neighbor && it->second->neighbor_address == address) {
			it->second->deprecateObject(max_age);
			modified_ = true;
		}
	}
}

void FlowStateObjects::removeObject(const std::string& fqn)
{
	rina::ScopedLock g(lock);

	LOG_IPCP_INFO("Trying to remove object %s", fqn.c_str());

	std::map<std::string, FlowStateObject*>::iterator it =
			objects.find(fqn);

	if (it == objects.end())
		return;

	LOG_IPCP_INFO("Object name before: %s", it->second->getObjectName().c_str());

	IPCPRIBDaemon* rib_daemon = (IPCPRIBDaemon*) IPCPFactory::getIPCP()->get_rib_daemon();
	rib_daemon->removeObjRIB(it->second->getObjectName());

	LOG_IPCP_INFO("Object name after: %s", it->second->getObjectName().c_str());

	objects.erase(it);
	LOG_IPCP_INFO("About to delete %p", it->second);
	//delete it->second;
}

bool FlowStateObjects::getObjectCopy(const std::string& fqn,
		   	   	     FlowStateObject& object)
{
	rina::ScopedLock g(lock);

	std::map<std::string, FlowStateObject*>::iterator it =
		objects.find(fqn);
	if (it == objects.end()) {
		LOG_IPCP_WARN("Could not find FSO with name %s", fqn.c_str());
		return false;
	}

	object.address = it->second->address;
	object.age = it->second->age;
	object.avoid_port = it->second->avoid_port;
	object.cost = it->second->cost;
	object.deprecated = it->second->deprecated;
	object.modified = it->second->modified;
	object.neighbor_address = it->second->neighbor_address;
	object.sequence_number = it->second->sequence_number;

	return true;
}

void FlowStateObjects::has_modified(bool modified)
{
	rina::ScopedLock g(lock);
	modified_ = modified;
}

void FlowStateObjects::modifyObject(const std::string& fqn,
		  	  	    FlowStateObject& object)
{
	std::map<std::string, FlowStateObject *>::iterator it;

	rina::ScopedLock g(lock);

	it = objects.find(fqn);
	if (it == objects.end()) {
		LOG_IPCP_WARN("Could not find FSO with name %s", fqn.c_str());
		return;
	}

	it->second->address = object.address;
	it->second->age = object.age;
	it->second->avoid_port = object.avoid_port;
	it->second->cost = object.cost;
	it->second->modified = true;
	it->second->neighbor_address = object.neighbor_address;
	it->second->state = object.state;
	it->second->sequence_number = object.sequence_number;
}

void FlowStateObjects::getModifiedFSOs(std::list<FlowStateObject >& result)
{
	rina::ScopedLock g(lock);

	for (std::map<std::string, FlowStateObject*>::iterator it
		= objects.begin(); it != objects.end();++it)
	{
		if (it->second->modified)
		{
			result.push_back(*(it->second));
		}
	}
}

void FlowStateObjects::getAllFSOs(std::list<FlowStateObject>& result)
{
	rina::ScopedLock g(lock);

	for (std::map<std::string, FlowStateObject*>::iterator it
			= objects.begin(); it != objects.end();++it)
	{
		result.push_back(*(it->second));
	}
}

void FlowStateObjects::incrementAge(unsigned int max_age,
			            unsigned long wait_until_remove,
				    rina::Timer* timer)
{
	rina::ScopedLock g(lock);

	for (std::map<std::string, FlowStateObject *>::iterator
		it = objects.begin(); it != objects.end(); ++it) 
	{
		if (it->second->age < UINT_MAX)
			it->second->age = it->second->age + 1;

		if (it->second->age >= max_age && !it->second->deprecated) {
			LOG_IPCP_INFO("Object %s will be removed %p. Age: %d",
					it->second->getObjectName().c_str(),
					it->second,
					it->second->age);
			it->second->deprecated = true;
			KillFlowStateObjectTimerTask* ksttask =
				new KillFlowStateObjectTimerTask(lsr_,
								 it->second->getObjectName());

			timer->scheduleTask(ksttask, wait_until_remove);
		}
	}
}

void FlowStateObjects::updateObject(const std::string& fqn, 
				    unsigned int avoid_port)
{
	rina::ScopedLock g(lock);
	std::map<std::string, FlowStateObject*>::iterator it =
		objects.find(fqn);
	if (it == objects.end())
		return;

	FlowStateObject* obj = it->second;
	obj->age = 0;
	obj->avoid_port = avoid_port;
	obj->deprecated = false;
	obj->state = true;
	obj->sequence_number = 1;
	obj->modified = true;
}

void FlowStateObjects::encodeAllFSOs(rina::ser_obj_t& obj)
{
	rina::ScopedLock g(lock);

	FlowStateObjectListEncoder encoder;
	std::list<FlowStateObject> result;
	if (!objects.empty())
	{
		for (std::map<std::string, FlowStateObject*>::iterator it
			= objects.begin(); it != objects.end();++it)
		{
			result.push_back(*(it->second));
		}
		encoder.encode(result, obj);
	}
	else
	{
		obj.size_ = 0;
		obj.message_ = 0;
	}
}

bool FlowStateObjects::is_modified()
{
	rina::ScopedLock g(lock);
	return modified_;
}

//Class FlowStateRIBObjects
const std::string FlowStateRIBObjects::clazz_name = "FlowStateObjects";
const std::string FlowStateRIBObjects::object_name= "/resalloc/fsos";

FlowStateRIBObjects::FlowStateRIBObjects(LinkStateRoutingPolicy * lsr) :
		rina::rib::RIBObj(clazz_name)
{
	lsr_ = lsr;
}

void FlowStateRIBObjects::read(const rina::cdap_rib::con_handle_t &con,
			       const std::string& fqn,
			       const std::string& clas,
			       const rina::cdap_rib::filt_info_t &filt,
			       const int invoke_id,
			       rina::ser_obj_t &obj_reply,
			       rina::cdap_rib::res_info_t& res)
{
	//TODO, implement following LSR policy specification
}

void FlowStateRIBObjects::write(const rina::cdap_rib::con_handle_t &con,
				const std::string& fqn,
				const std::string& clas,
				const rina::cdap_rib::filt_info_t &filt,
				const int invoke_id,
				const rina::ser_obj_t &obj_req,
				rina::ser_obj_t &obj_reply,
				rina::cdap_rib::res_info_t& res)
{
	FlowStateObjectListEncoder encoder;
	std::list<FlowStateObject> new_objects;
	encoder.decode(obj_req, new_objects);
	lsr_->updateObjects(new_objects,
			    con.port_id,
			    IPCPFactory::getIPCP()->get_address());
}

// ComputeRoutingTimerTask
ComputeRoutingTimerTask::ComputeRoutingTimerTask(
		LinkStateRoutingPolicy * lsr_policy, long delay)
{
	lsr_policy_ = lsr_policy;
	delay_ = delay;
}

void ComputeRoutingTimerTask::run()
{
	lsr_policy_->routingTableUpdate();

	//Re-schedule
	ComputeRoutingTimerTask * task = new ComputeRoutingTimerTask(
			lsr_policy_, delay_);

	lsr_policy_->timer_->scheduleTask(task, delay_);
}

KillFlowStateObjectTimerTask::KillFlowStateObjectTimerTask(LinkStateRoutingPolicy *lsr,
							   const std::string& fqn)
{
	lsr_ = lsr;
	fqn_ = fqn;
}

void KillFlowStateObjectTimerTask::run()
{
	lsr_->removeFlowStateObject(fqn_);
}

PropagateFSODBTimerTask::PropagateFSODBTimerTask(
		LinkStateRoutingPolicy * lsr_policy, long delay)
{
	lsr_policy_ = lsr_policy;
	delay_ = delay;
}

void PropagateFSODBTimerTask::run()
{
	lsr_policy_->propagateFSDB();

	//Re-schedule
	PropagateFSODBTimerTask * task = new PropagateFSODBTimerTask(
			lsr_policy_, delay_);
	lsr_policy_->timer_->scheduleTask(task, delay_);
}

UpdateAgeTimerTask::UpdateAgeTimerTask(
		LinkStateRoutingPolicy * lsr_policy, long delay)
{
	lsr_policy_ = lsr_policy;
	delay_ = delay;
}

void UpdateAgeTimerTask::run()
{
	lsr_policy_->updateAge();

	//Re-schedule
	UpdateAgeTimerTask * task = new UpdateAgeTimerTask(lsr_policy_,
			delay_);
	lsr_policy_->timer_->scheduleTask(task, delay_);
}

ExpireOldAddressTimerTask::ExpireOldAddressTimerTask(LinkStateRoutingPolicy * lsr_policy,
						     unsigned int addr,
						     bool neigh)
{
	lsr_policy_ = lsr_policy;
	address = addr;
	neighbor = neigh;
}

void ExpireOldAddressTimerTask::run()
{
	lsr_policy_->expireOldAddress(address, neighbor);
}

// CLASS LinkStateRoutingPolicy
const std::string LinkStateRoutingPolicy::OBJECT_MAXIMUM_AGE = "objectMaximumAge";
const std::string LinkStateRoutingPolicy::WAIT_UNTIL_READ_CDAP = "waitUntilReadCDAP";
const std::string LinkStateRoutingPolicy::WAIT_UNTIL_ERROR = "waitUntilError";
const std::string LinkStateRoutingPolicy::WAIT_UNTIL_PDUFT_COMPUTATION = "waitUntilPDUFTComputation";
const std::string LinkStateRoutingPolicy::WAIT_UNTIL_FSODB_PROPAGATION = "waitUntilFSODBPropagation";
const std::string LinkStateRoutingPolicy::WAIT_UNTIL_AGE_INCREMENT = "waitUntilAgeIncrement";
const std::string LinkStateRoutingPolicy::WAIT_UNTIL_REMOVE_OBJECT = "waitUntilRemoveObject";
const std::string LinkStateRoutingPolicy::WAIT_UNTIL_DEPRECATE_OLD_ADDRESS = "waitUntilDeprecateAddress";
const std::string LinkStateRoutingPolicy::ROUTING_ALGORITHM = "routingAlgorithm";
const int LinkStateRoutingPolicy::MAXIMUM_BUFFER_SIZE = 4096;
const std::string LinkStateRoutingPolicy::DIJKSTRA_ALG = "Dijkstra";
const std::string LinkStateRoutingPolicy::ECMP_DIJKSTRA_ALG = "ECMPDijkstra";

LinkStateRoutingPolicy::LinkStateRoutingPolicy(IPCProcess * ipcp)
{
	test_ = false;
	ipc_process_ = ipcp;
	rib_daemon_ = ipc_process_->rib_daemon_;
	routing_algorithm_ = 0;
	resiliency_algorithm_ = 0;
	wait_until_deprecate_address_ = 0;

	subscribeToEvents();
	timer_ = new rina::Timer();
	fsos = new FlowStateObjects(this);
	maximum_age = UINT_MAX;
	wait_until_remove_obj = 0;
}

LinkStateRoutingPolicy::~LinkStateRoutingPolicy()
{
	delete timer_;
	delete routing_algorithm_;
	delete resiliency_algorithm_;
	delete fsos;
}

void LinkStateRoutingPolicy::subscribeToEvents()
{
	ipc_process_->internal_event_manager_->
		subscribeToEvent(rina::InternalEvent::APP_N_MINUS_1_FLOW_DEALLOCATED, this);
	ipc_process_->internal_event_manager_->
		subscribeToEvent(rina::InternalEvent::APP_N_MINUS_1_FLOW_ALLOCATED, this);
	ipc_process_->internal_event_manager_->
		subscribeToEvent(rina::InternalEvent::APP_NEIGHBOR_ADDED, this);
	ipc_process_->internal_event_manager_->
		subscribeToEvent(rina::InternalEvent::APP_CONNECTIVITY_TO_NEIGHBOR_LOST, this);
	ipc_process_->internal_event_manager_->
		subscribeToEvent(rina::InternalEvent::ADDRESS_CHANGE, this);
	ipc_process_->internal_event_manager_->
		subscribeToEvent(rina::InternalEvent::NEIGHBOR_ADDRESS_CHANGE, this);
}

void LinkStateRoutingPolicy::set_dif_configuration(
		const rina::DIFConfiguration& dif_configuration)
{
	std::string routing_alg;
        rina::PolicyConfig psconf;
        long delay;

        psconf = dif_configuration.routing_configuration_.policy_set_;

        try {
        	routing_alg = psconf.get_param_value_as_string(ROUTING_ALGORITHM);
        } catch (rina::Exception &e) {
        	LOG_WARN("No routing algorithm specified, using Dijkstra");
        	routing_alg = DIJKSTRA_ALG;
        }

        if (routing_alg == DIJKSTRA_ALG) {
        	routing_algorithm_ = new DijkstraAlgorithm();
                LOG_IPCP_DBG("Using Dijkstra as routing algorithm");
        } else if (routing_alg == ECMP_DIJKSTRA_ALG)  {
                routing_algorithm_ = new ECMPDijkstraAlgorithm();
                LOG_IPCP_DBG("Using ECMP Dijkstra as routing algorithm");
        } else {
        	throw rina::Exception("Unsupported routing algorithm");
        }
#if 0
	resiliency_algorithm_ = new LoopFreeAlternateAlgorithm(*routing_algorithm_);
#endif


	if (!test_) {
		try {

			maximum_age = psconf.get_param_value_as_int(OBJECT_MAXIMUM_AGE);
		} catch (rina::Exception &e) {
			maximum_age = PULSES_UNTIL_FSO_EXPIRATION_DEFAULT;
		}

		try {
			wait_until_remove_obj = psconf.get_param_value_as_long(WAIT_UNTIL_REMOVE_OBJECT);
		} catch (rina::Exception &e) {
			wait_until_remove_obj = WAIT_UNTIL_REMOVE_OBJECT_DEFAULT;
		}

		try {
			wait_until_deprecate_address_ = psconf.get_param_value_as_long(WAIT_UNTIL_DEPRECATE_OLD_ADDRESS);
		} catch (rina::Exception &e) {
			wait_until_deprecate_address_ = WAIT_UNTIL_DEPRECATE_OLD_ADDRESS_DEFAULT;
		}

		// Task to compute PDUFT
		try {
			delay = psconf.get_param_value_as_long(WAIT_UNTIL_PDUFT_COMPUTATION);
		} catch (rina::Exception &e) {
			delay = WAIT_UNTIL_PDUFT_COMPUTATION_DEFAULT;
		}
		ComputeRoutingTimerTask * cttask = new ComputeRoutingTimerTask(this, delay);
		timer_->scheduleTask(cttask, delay);

		// Task to increment age
		try {
			delay = psconf.get_param_value_as_long(WAIT_UNTIL_AGE_INCREMENT);
		} catch (rina::Exception &e) {
			delay = WAIT_UNTIL_AGE_INCREMENT_DEFAULT;
		}
		UpdateAgeTimerTask * uattask = new UpdateAgeTimerTask(this, delay);
		timer_->scheduleTask(uattask, delay);

		// Task to propagate modified FSO
		try {
			delay = psconf.get_param_value_as_long(WAIT_UNTIL_FSODB_PROPAGATION);
		} catch (rina::Exception &e) {
			delay = WAIT_UNTIL_FSODB_PROPAGATION_DEFAULT;
		}
		PropagateFSODBTimerTask * pfttask = new PropagateFSODBTimerTask(this,
				delay);
		timer_->scheduleTask(pfttask, delay);
	}
}

const std::list<rina::FlowInformation>& LinkStateRoutingPolicy::get_allocated_flows() const
{
	return allocated_flows_;
}

void LinkStateRoutingPolicy::eventHappened(rina::InternalEvent * event)
{
	if (!event)
		return;

	if (event->type == rina::InternalEvent::APP_N_MINUS_1_FLOW_DEALLOCATED) {
		rina::NMinusOneFlowDeallocatedEvent * flowEvent =
				(rina::NMinusOneFlowDeallocatedEvent *) event;
		processFlowDeallocatedEvent(flowEvent);
	} else if (event->type == rina::InternalEvent::APP_N_MINUS_1_FLOW_ALLOCATED) {
		rina::NMinusOneFlowAllocatedEvent * flowEvent =
			(rina::NMinusOneFlowAllocatedEvent *) event;
		processFlowAllocatedEvent(flowEvent);
	} else if (event->type == rina::InternalEvent::APP_NEIGHBOR_ADDED) {
		rina::NeighborAddedEvent * neighEvent = (rina::NeighborAddedEvent *) event;
		processNeighborAddedEvent(neighEvent);
	} else if (event->type == rina::InternalEvent::APP_CONNECTIVITY_TO_NEIGHBOR_LOST) {
		rina::ConnectiviyToNeighborLostEvent * neighEvent =
			(rina::ConnectiviyToNeighborLostEvent *) event;
		processNeighborLostEvent(neighEvent);
	} else if (event->type == rina::InternalEvent::ADDRESS_CHANGE) {
		rina::AddressChangeEvent * addrEvent =
			(rina::AddressChangeEvent *) event;
		processAddressChangeEvent(addrEvent);
	} else if (event->type == rina::InternalEvent::NEIGHBOR_ADDRESS_CHANGE) {
		rina::NeighborAddressChangeEvent * addrEvent =
			(rina::NeighborAddressChangeEvent *) event;
		processNeighborAddressChangeEvent(addrEvent);
	}
}

void LinkStateRoutingPolicy::processAddressChangeEvent(rina::AddressChangeEvent * event)
{
	std::list<FlowStateObject> all_fsos;

	rina::ScopedLock g(lock_);

	fsos->getAllFSOs(all_fsos);

	//Add LSOs to reflect the routes to the new address
	for (std::list<FlowStateObject>::iterator it = all_fsos.begin();
			it != all_fsos.end(); ++it) {
		if (it->address == event->old_address) {
			FlowStateObject newObject(event->new_address,
						  it->neighbor_address,
						  1,
						  true,
						  1,
						  0);
			fsos->addObject(newObject);
		}
	}

	//Schedule task to remove all objects with old address
	ExpireOldAddressTimerTask * task = new ExpireOldAddressTimerTask(this,
									 event->old_address,
									 false);
	timer_->scheduleTask(task, event->deprecate_old_timeout);
}

void LinkStateRoutingPolicy::processNeighborAddressChangeEvent(rina::NeighborAddressChangeEvent * event)
{
	rina::ScopedLock g(lock_);

	std::list<FlowStateObject> all_fsos;

	fsos->getAllFSOs(all_fsos);

	LOG_IPCP_DBG("Neighbor address changed: old address %d, new address %d",
		     event->old_address,
		     event->new_address);
	//Add LSOs to reflect the routes to the new address
	for (std::list<FlowStateObject>::iterator it = all_fsos.begin();
			it != all_fsos.end(); ++it) {
		if (it->neighbor_address == event->old_address) {
			FlowStateObject newObject(ipc_process_->get_address(),
						  event->new_address,
						  1,
						  true,
						  1,
						  0);
			fsos->addObject(newObject);
		}
	}

	//Schedule task to remove all objects with old address
	ExpireOldAddressTimerTask * task = new ExpireOldAddressTimerTask(this,
									 event->old_address,
									 true);
	timer_->scheduleTask(task, wait_until_deprecate_address_);
}

void LinkStateRoutingPolicy::expireOldAddress(unsigned int address, bool neighbor)
{
	rina::ScopedLock g(lock_);
	fsos->deprecateObjectsWithAddress(address, maximum_age, neighbor);
}

void LinkStateRoutingPolicy::removeFlowStateObject(const std::string& fqn)
{
	rina::ScopedLock g(lock_);
	fsos->removeObject(fqn);
}

void LinkStateRoutingPolicy::updateObjects(const std::list<FlowStateObject>& newObjects,
		   	   	   	   unsigned int avoidPort,
					   unsigned int address)
{
	rina::ScopedLock g(lock_);
	LOG_IPCP_INFO("Received Flow State object WRITE message via N-1 port-id %d (%d objects)",
		       avoidPort,
		       newObjects.size());
	FlowStateObject aux_fso;

	for (std::list<FlowStateObject>::const_iterator
		newIt = newObjects.begin(); newIt != newObjects.end(); ++newIt)
	{
		//1 If the object exists update
		if (fsos->getObjectCopy(newIt->getObjectName(), aux_fso))
		{
			LOG_IPCP_DBG("Found the object in the DB. Object: %s",
				aux_fso.getObjectName().c_str());
			//1.1 If the object has a higher sequence number update
			if (newIt->sequence_number > aux_fso.sequence_number)
			{
				if (newIt->address == address)
				{
					LOG_IPCP_DBG("Object is self generated, updating the sequence number and age of %s to %d",
						      aux_fso.getObjectName().c_str(),
						      aux_fso.sequence_number);
					aux_fso.sequence_number = newIt->sequence_number + 1;
					aux_fso.avoid_port = NO_AVOID_PORT;
					aux_fso.age = 0;
				}
				else
				{
					LOG_IPCP_DBG("Update the object %s with seq num %d",
						      aux_fso.getObjectName().c_str(),
						      newIt->sequence_number);
					aux_fso.avoid_port = avoidPort;
					if (newIt->age >= maximum_age) {
						aux_fso.deprecateObject(maximum_age);
					} else {
						aux_fso.age = 0;
						aux_fso.sequence_number = newIt->sequence_number;
					}
				}

				aux_fso.modified = true;
				fsos->modifyObject(newIt->getObjectName(), aux_fso);
				fsos->has_modified(true);
			}
		}

		//2. If the object does not exist create
		else
		{
			if(newIt->address != address)
			{
				LOG_IPCP_DBG("New object added");
				FlowStateObject fso(*newIt);
				fso.avoid_port = avoidPort;
				fso.modified = true;
				fsos->addObject(fso);
				fsos->has_modified(true);
			}
		}
	}
}

void LinkStateRoutingPolicy::processFlowDeallocatedEvent(
		rina::NMinusOneFlowDeallocatedEvent * event)
{
	rina::ScopedLock g(lock_);

	for (std::list<rina::FlowInformation>::iterator it
		= allocated_flows_.begin(); it != allocated_flows_.end(); ++it) {
		if (it->portId == event->port_id_) {
			allocated_flows_.erase(it);
			return;
		}
	}

	LOG_IPCP_DBG("N-1 Flow with neighbor lost");
	//TODO update cost
}

void LinkStateRoutingPolicy::processNeighborLostEvent(
		rina::ConnectiviyToNeighborLostEvent* event)
{
	rina::ScopedLock g(lock_);
	fsos->deprecateObjects(event->neighbor_.address_,
			       ipc_process_->get_address(),
			       maximum_age);
}


void LinkStateRoutingPolicy::processFlowAllocatedEvent(
		rina::NMinusOneFlowAllocatedEvent * event)
{
	rina::ScopedLock g(lock_);

	if (ipc_process_->resource_allocator_->get_n_minus_one_flow_manager()->
			numberOfFlowsToNeighbour(event->flow_information_.remoteAppName.processName,
					event->flow_information_.remoteAppName.processInstance) > 1) {
		LOG_IPCP_DBG("Already had an N-1 flow with this neighbor IPCP");
		//TODO update the cost of the FlowStateObject
		return;
	}

	try 
	{
		FlowStateObject newObject(ipc_process_->get_address(),
				 	  ipc_process_->namespace_manager_->getAdressByname(event->flow_information_.remoteAppName),
					  1,
					  true,
					  1,
					  0);

		fsos->addObject(newObject);
	} catch (rina::Exception &e) 
	{
		LOG_IPCP_DBG("flow allocation waiting for enrollment");
		allocated_flows_.push_back(event->flow_information_);
	}
}

void LinkStateRoutingPolicy::processNeighborAddedEvent(
		rina::NeighborAddedEvent * event)
{
	rina::ScopedLock g(lock_);

	int portId = event->neighbor_.get_underlying_port_id();
	for (std::list<rina::FlowInformation>::iterator it = 
		allocated_flows_.begin(); it != allocated_flows_.end(); ++it) 
	{
		if (it->portId == portId)
		{
			LOG_IPCP_INFO("There was an allocation flow event waiting for enrollment, launching it");
			try {
				FlowStateObject newObject(ipc_process_->get_address(),
							  ipc_process_->namespace_manager_->getAdressByname(event->neighbor_.get_name()),
							  1,
							  true,
							  1,
							  0);

				fsos->addObject(newObject);
				allocated_flows_.erase(it);
				break;
			} catch (rina::Exception &e) {
				LOG_IPCP_ERR("Could not allocate the flow, no neighbor found");
			}
		}
	}

	try {
		rina::cdap_rib::obj_info_t obj;
		rina::cdap_rib::con_handle_t con;
		obj.class_ = FlowStateRIBObjects::clazz_name;
		obj.name_ = FlowStateRIBObjects::object_name;
		fsos->encodeAllFSOs(obj.value_);
		obj.inst_ = 0;
		rina::cdap_rib::flags_t flags;
		rina::cdap_rib::filt_info_t filt;
		con.port_id = portId;
		if (obj.value_.size_ != 0)
			rib_daemon_->getProxy()->remote_write(con,
							      obj,
							      flags,
							      filt,
							      0);
	} catch (rina::Exception &e) {
		LOG_IPCP_ERR("Problems encoding and sending CDAP message: %s", e.what());
	}
}

void LinkStateRoutingPolicy::propagateFSDB()
{
	rina::ScopedLock g(lock_);

	//1 Get the active flows
	std::list<rina::FlowInformation> nMinusOneFlows =
			ipc_process_->resource_allocator_->get_n_minus_one_flow_manager()->getAllNMinusOneFlowInformation();
	//2 Initilize the map
	std::map <int, std::list<FlowStateObject> > objectsToSend;
	for(std::list<rina::FlowInformation>::iterator it = nMinusOneFlows.begin();
		it != nMinusOneFlows.end(); ++it) 
	{
		objectsToSend[it->portId] = std::list<FlowStateObject>();
	}

	//3 Get the objects to send
	std::list<FlowStateObject> modifiedFSOs;
	fsos->getModifiedFSOs(modifiedFSOs);

	//4 add each modified object to its port list
	for (std::list<FlowStateObject>::iterator it = modifiedFSOs.begin();
			it != modifiedFSOs.end(); ++it)
	{
		LOG_DBG("Propagation: Check modified object %s with age %d and status %d",
			it->getObjectName().c_str(),
			it->age,
			it->state);

		for(std::map<int, std::list<FlowStateObject> >::iterator it2 =
				objectsToSend.begin(); it2 != objectsToSend.end(); ++it2)
		{
			if(it2->first != it->avoid_port)
			{
				it2->second.push_back(*it);
			}
		}

		it->modified = false;
		it->avoid_port = NO_AVOID_PORT;
		fsos->modifyObject(it->getObjectName(), *it);
	}

	if (objectsToSend.size() == 0) {
		return;
	}

	FlowStateObjectListEncoder encoder;
	rina::cdap_rib::con_handle_t con;
	for (std::map<int, std::list<FlowStateObject> >::iterator it = objectsToSend.begin();
		it != objectsToSend.end(); ++it)

	{
		if (it->second.size() == 0)
			continue;

		rina::cdap_rib::flags flags;
		rina::cdap_rib::filt_info_t filter;
		try
		{
			rina::cdap_rib::object_info obj;
			obj.class_ = FlowStateRIBObjects::clazz_name;
			obj.name_ = FlowStateRIBObjects::object_name;
			encoder.encode(it->second, obj.value_);
			con.port_id = it->first;
			rib_daemon_->getProxy()->remote_write(con,
							      obj,
							      flags,
							      filter,
							      0);
		}
		catch (rina::Exception &e) 
		{
			LOG_IPCP_ERR("Errors sending message: %s", e.what());
		}
	}
}

void LinkStateRoutingPolicy::updateAge()
{
	rina::ScopedLock g(lock_);
	fsos->incrementAge(maximum_age,
			   wait_until_remove_obj,
			   timer_);
}

void LinkStateRoutingPolicy::removeDuplicateEntries(std::list<rina::RoutingTableEntry *>& rt)
{
	std::list<rina::RoutingTableEntry *>::iterator it;
	std::list<rina::RoutingTableEntry *>::iterator jt;
	rina::RoutingTableEntry * current = 0;
	rina::RoutingTableEntry * candidate = 0;
	bool increment = true;

	it = rt.begin();
	while (it != rt.end()) {
		current = *it;
		increment = true;

		for (jt = rt.begin(); jt != rt.end(); ++jt) {
			if(it == jt)
				continue;

			candidate = *jt;

			//Detect multiple next hops to the same target, some with higher cost
			if (current->address == candidate->address &&
			    current->qosId == candidate->qosId &&
			    current->cost > candidate->cost) {
				it = rt.erase(it);
				delete current;
				increment = false;
				break;
			}

			//Duplicate detection
			if (current->address == candidate->address &&
			    current->qosId == candidate->qosId &&
			    current->cost == candidate->cost &&
			    current->nextHopAddresses.size() == candidate->nextHopAddresses.size() &&
			    current->nextHopAddresses.front().alts.front() == candidate->nextHopAddresses.front().alts.front()) {
				it = rt.erase(it);
				delete current;
				increment = false;
				break;
			}
		}

		if (increment)
			++it;
	}
}

void LinkStateRoutingPolicy::printNhopTable(std::list<rina::RoutingTableEntry *>& rt)
{
	std::list<rina::RoutingTableEntry *>::iterator it;
	std::list<rina::NHopAltList>::iterator jt;
	rina::RoutingTableEntry * current = 0;
	std::stringstream ss;

	for(it = rt.begin(); it != rt.end(); ++it) {
		ss.str(std::string());
		ss.clear();
		current = *it;
		ss << "Dest. address: " << current->address
		   << "; QoS-id: " << current->qosId
		   << "; Cost: " << current->cost
		   << "; Next hop addresses: ";
		for (jt = current->nextHopAddresses.begin();
				jt != current->nextHopAddresses.end(); ++jt) {
			ss << jt->alts.front() << "; ";
		}

		LOG_IPCP_INFO("%s", ss.str().c_str());
	}
}

void LinkStateRoutingPolicy::routingTableUpdate()
{
	std::list<FlowStateObject> all_fsos;
	std::list<FlowStateObject> g1_fsos;
	std::list<FlowStateObject> g2_fsos;
	std::list<rina::RoutingTableEntry *> rt;
	Graph g1;
	Graph g2;
	std::list<FlowStateObject>::iterator it;
	bool add = false;
	unsigned int address = 0;
	unsigned int old_address = 0;

	rina::ScopedLock g(lock_);

	if (!fsos->is_modified())
		return;
	fsos->has_modified(false);

	fsos->getAllFSOs(all_fsos);
	address = ipc_process_->get_address();
	old_address = ipc_process_->get_old_address();
	if (old_address != 0) {
		//Remove fsos with old_address to avoid routes to myself
		for (it=all_fsos.begin(); it != all_fsos.end(); ++it) {
			if (it->address != old_address && it->neighbor_address!= old_address)
				g1_fsos.push_back(*it);
		}

		g1.set_flow_state_objects(g1_fsos);
		routing_algorithm_->computeRoutingTable(g1,
							g1_fsos,
							address,
							rt);

		//Remove fsos with address to avoid routes to myself
		for (it=all_fsos.begin(); it != all_fsos.end(); ++it) {
			if (it->address != address && it->neighbor_address != address)
				g2_fsos.push_back(*it);
		}

		g2.set_flow_state_objects(g2_fsos);
		routing_algorithm_->computeRoutingTable(g2,
							g2_fsos,
							old_address,
							rt);
		removeDuplicateEntries(rt);
	} else {
		g1.set_flow_state_objects(all_fsos);
		routing_algorithm_->computeRoutingTable(g1,
							all_fsos,
							address,
							rt);
	}

	// Run the resiliency algorithm, if any, to extend the routing table
	if (resiliency_algorithm_) {
		resiliency_algorithm_->fortifyRoutingTable(g1,
							   ipc_process_->get_address(),
							   rt);
		if (old_address != 0) {
			resiliency_algorithm_->fortifyRoutingTable(g2,
								   ipc_process_->get_old_address(),
								   rt);
		}
	}

	LOG_IPCP_INFO("Computed new Next Hop and PDU Forwarding Tables");
	printNhopTable(rt);

	assert(ipc_process_->resource_allocator_->pduft_gen_ps);
	ipc_process_->resource_allocator_->pduft_gen_ps->routingTableUpdated(rt);
}

// CLASS FlowStateObjectEncoder
namespace fso_helpers{
void toGPB(	const FlowStateObject &fso, 
	rina::messages::flowStateObject_t &gpb_fso)
{
	gpb_fso.set_address(fso.address);
	gpb_fso.set_age(fso.age);
	gpb_fso.set_neighbor_address(fso.neighbor_address);
	gpb_fso.set_cost(fso.cost);
	gpb_fso.set_state(fso.state);
	gpb_fso.set_sequence_number(fso.sequence_number);
}

void toModel(
	const rina::messages::flowStateObject_t &gpb_fso, FlowStateObject &fso)
{
	fso.address = gpb_fso.address();
	fso.neighbor_address = gpb_fso.neighbor_address();
	fso.cost = gpb_fso.cost();
	fso.state = gpb_fso.state();
	fso.sequence_number = gpb_fso.sequence_number();
	fso.age = gpb_fso.age();
}
} //namespace fso_helpers

void FlowStateObjectEncoder::encode(const FlowStateObject &obj, 
				    rina::ser_obj_t &serobj)
{
	rina::messages::flowStateObject_t gpb;

	fso_helpers::toGPB(obj, gpb);

	serobj.size_ = gpb.ByteSize();
	serobj.message_ = new unsigned char[serobj.size_];
	gpb.SerializeToArray(serobj.message_, serobj.size_);
}

void FlowStateObjectEncoder::decode(const rina::ser_obj_t &serobj,
				    FlowStateObject &des_obj)
{
	rina::messages::flowStateObject_t gpb;
	gpb.ParseFromArray(serobj.message_, serobj.size_);

	fso_helpers::toModel(gpb, des_obj);
}

void FlowStateObjectListEncoder::encode(const std::list<FlowStateObject> &obj,
					rina::ser_obj_t& serobj)
{
	rina::messages::flowStateObjectGroup_t gpb;

	for (std::list<FlowStateObject>::const_iterator it= obj.begin();
		it != obj.end(); ++it) 
	{
		rina::messages::flowStateObject_t *gpb_fso;
		gpb_fso = gpb.add_flow_state_objects();
		fso_helpers::toGPB(*it, *gpb_fso);
	}

	serobj.size_ = gpb.ByteSize();
	serobj.message_ = new unsigned char[serobj.size_];
	gpb.SerializeToArray(serobj.message_, serobj.size_);
}

void FlowStateObjectListEncoder::decode(const rina::ser_obj_t &serobj,
					std::list<FlowStateObject> &des_obj)
{
	rina::messages::flowStateObjectGroup_t gpb;
	gpb.ParseFromArray(serobj.message_, serobj.size_);
	for(int i=0; i<gpb.flow_state_objects_size(); i++)
	{
		FlowStateObject fso;
		fso_helpers::toModel(gpb.flow_state_objects(i), fso);
		des_obj.push_back(fso);
	}
}

}// namespace rinad
