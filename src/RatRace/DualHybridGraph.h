#pragma once

#include <vector>
#include <array>
#include <map>
#include <experimental/filesystem>

#include <glog/logging.h>

#include "MeshHeader.h"
#include "PolySurfMesh.h"
#include <memory>

class DualHybridGraph
{
	
public:
	using NodeIdxT = long;
	using HalfedgeIdxT = long;

	class Halfedge
	{
		NodeIdxT from_ = -1;
		NodeIdxT to_ = -1;
		HalfedgeIdxT idx_ = -1;
		bool isValid_ = false;

		HalfedgeIdxT oppositeIdx_ = -1;
		HalfedgeIdxT nextIdx_ = -1;
		HalfedgeIdxT prevIdx_ = -1;

		PolyMesh::VertexHandle vertex_;

	public:
		Halfedge() {}
		Halfedge( HalfedgeIdxT idx ) : idx_{ idx }, isValid_{ true } {}
		
		auto& from() { return from_; }
		auto from() const { return from_; }
		auto& to() { return to_; }
		auto to() const { return to_; }
		auto idx() const { return idx_; }
		auto isValid() const { return isValid_; }
		auto& isValid() { return isValid_; }

		auto& oppositeIdx() { return oppositeIdx_; }
		auto& nextIdx() { return nextIdx_; }
		auto& prevIdx() { return prevIdx_; }

		auto& vertex() { return vertex_; }
		auto vertex() const { return vertex_; }
	};

	class : public std::vector<Halfedge>
	{
	public:
		Halfedge& add() {
			push_back( Halfedge( (HalfedgeIdxT)this->size() ) );
			return this->back();
		}

	} halfedges;

	class Node
	{
		NodeIdxT idx_ = -1;
		bool isValid_ = false;
		HalfedgeIdxT halfedgeIdx_ = -1;
		std::vector<PolyMesh::VertexHandle> vertices_;
	public:
		Node() {}
		Node( NodeIdxT idx ) : idx_{ idx }, isValid_{ true } {}
		auto idx() const { return idx_; }
		auto isValid() const { return isValid_; }
		auto& isValid() { return isValid_; }
		auto& halfedgeIdx() { return halfedgeIdx_; }
		auto halfedgeIdx() const { return halfedgeIdx_; }
		auto& vertices() { return vertices_; }
		auto vertices() const { return vertices_; }

	};

	class : public std::vector<Node>
	{
	public:
		Node& add() {
			push_back( Node((NodeIdxT)this->size()) );
			return this->back();
		}
	} nodes;

	PolyMesh mesh_;
	
	std::vector<bool> isTriangleVec_;

	DualHybridGraph( PolySurfMesh& mesh );
	
	// return ALL neighbors including boundaries, e.g. invalid node handles
	std::vector<NodeIdxT> neighborIdxs( const NodeIdxT& idx );
	// return neighbors without boundaries (no invalid node handles)
	std::vector<NodeIdxT> innerNeighborIdxs( const NodeIdxT& idx );

	/// <summary>
	/// Iterates around vertex and returns all outgoing halfedges
	/// </summary>
	/// <param name="idx">center node index</param>
	/// <returns>outgoing halfedges</returns>
	auto outgoingHalfedges( const NodeIdxT& idx );

	// return halfedge index
	HalfedgeIdxT findHalfedge( const NodeIdxT& nodeFrom, const NodeIdxT& nodeTo );

	PolyMesh getMesh();

	void printMesh( const std::experimental::filesystem::path& filename );

	/// <summary>
	/// Collapse two nodes in along halfedge. Multiple halfedges are also deleted.
	/// </summary>
	/// <param name="idx">halfedge pointing towards remaining node</param>
	/// <returns>Index of deleted node</returns>
	NodeIdxT halfedgeCollapse( const HalfedgeIdxT& idx );

	/// <summary>
	/// Split nodeFrom of halfedges[idx] s. th. it contains a triangle neighboring nodeTo.
	/// The split is performed s.th.the best quad is generated.
	/// The new halfedge points to nodeFrom, halfedges[idx] remains valid.
	/// </summary>
	/// <param name="idx"></param>
	/// <param name="invalidNodeIdx"></param>
	void split3to4( const HalfedgeIdxT& idx, const NodeIdxT& invalidNodeIdx );

	void garbageCollection();

	bool isTriangle( const NodeIdxT& idx );

	bool hasTriangles();

	int nTriangles();

	// Search for the first triangle in nodes and connect it with the nearest triangle
	void connectTriangles();

	void printGridInfo();
};


class OutgoingHalfedgeRange
{
	DualHybridGraph& dhg_;
	DualHybridGraph::NodeIdxT nodeIdx_;
	DualHybridGraph::HalfedgeIdxT first_;
public:
	OutgoingHalfedgeRange( DualHybridGraph::NodeIdxT nodeIdx, DualHybridGraph& dhg ) 
		: dhg_{ dhg }, nodeIdx_{ nodeIdx }, first_{ dhg.nodes[nodeIdx].halfedgeIdx() } {}

	class OutgoingHalfedgeIterator
	{
		DualHybridGraph& dhg_;
		int halfedgeIdx_;
		bool isFirst_ = false;
	public:
		OutgoingHalfedgeIterator( DualHybridGraph::HalfedgeIdxT halfedgeIdx, DualHybridGraph& dhg, bool isFirst = false ) 
			: dhg_{ dhg }, halfedgeIdx_{ halfedgeIdx }, isFirst_{ isFirst } {}
		OutgoingHalfedgeIterator operator++() {
			auto gheh = dhg_.halfedges[halfedgeIdx_];
			gheh = dhg_.halfedges[gheh.prevIdx()];
			gheh = dhg_.halfedges[gheh.oppositeIdx()];
			halfedgeIdx_ = gheh.idx();
			isFirst_ = false;
			return *this;
		}
		bool operator!=( const OutgoingHalfedgeIterator& other ) const { 
			return ( halfedgeIdx_ != other.halfedgeIdx_ ) ^ ( isFirst_ != other.isFirst_ );
		}
		DualHybridGraph::Halfedge& operator*() const { return dhg_.halfedges[halfedgeIdx_]; }
	};
public:
	OutgoingHalfedgeIterator begin() const { return OutgoingHalfedgeIterator( first_, dhg_, true ); }
	OutgoingHalfedgeIterator end() const { return OutgoingHalfedgeIterator( first_, dhg_ ); }
};
