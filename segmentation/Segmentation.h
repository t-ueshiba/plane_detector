/*!
  \file		Segmentation.h
  \author	Toshio UESHIBA
  \brief	クラス TU::Segmentation の定義と実装
*/
#pragma once

#include <list>
#include <vector>
#include <map>
#include <cstdint>
#include <string>
#include <limits>

namespace plane_detector
{
/************************************************************************
*  class Segmentation<R>						*
************************************************************************/
template <class R>
class Segmentation
{
  public:
    class Vertex;
    class Region;

    using riterator	  = typename std::list<Region>::iterator;
    using const_riterator = typename std::list<Region>::const_iterator;
    using viterator	  = typename std::list<Vertex>::iterator;
    using const_riterator = typename std::list<Vertex>::const_iterator;

    class Edge;

    class Vertex
    {
	friend class	Edge;

      public:
			Vertex(riterator r[], viterator vend)
			{
			    for (size_t e = 0; e < 4; ++e)
				_r[e] = r[e];
			}

	size_t		valence(viterator vend) const
			{
			    size_t	n = 0;
			    for (auto vt : _v)
				if (vt != vend)
				    ++n;
			    return n;
			}
	int		x()				const	{ return _x; }
	int		y()				const	{ return _y; }

      private:
	viterator	self(viterator vend) const
			{
			    for (auto vt : _v)
				if (vt != vend)
				{
				    for (auto v : vt->_v)
					if (v != vend && &(*v) == this)
					    return v;

				    throw std::runtime_error("Segmentation<R>::Vertex::self(): Internal error!");
				}
			    return vend;
			}

      private:
	riterator	_r[4];
	viterator	_v[4];
	int		_x;
	int		_y;

	friend std::ostream&
			operator <<(std::ostream& out, const Vertex& v)
			{
			    out << "Vertex[" << &v << "]:" << std::endl;
			    for (_size_t e = 0; e < 4; ++e)
				if (v._v[e] != vend)
				    out << "  v[" << &(*v._v[e])
					<< "], r[" << &(*v._r[e])
					<< "]:" << *v._r[e]
					<< std::endl;
			    return out;
			}
    };

    class Edge
    {
      public:
			Edge(const Vertex& v, viterator vend)
			    :Edge(vertex.self(vend), vend)		{}

	bool		operator ==(const Edge& edge) const
			{
			    return (_e == edge._e) && (_v == edge._v);
			}
	bool		operator !=(const Edge& edge) const
			{
			    return !(*this == edge);
			}
	size_t		valence() const
			{
			    return _v->valence(_vend);
			}
	Edge&		operator ++()
			{
			    do
			    {
				_e = (_e == 3 ? 0 : _e + 1);
			    } while (vt() == _vend);

			    return *this;
			}
	Edge&		operator --()
			{
			    do
			    {
				_e = (_e == 0 ? 3 : _e - 1);
			    } while (vt() == _vend);

			    return *this;
			}
	Edge&		operator ~()
			{
			    const auto	v = _v;
			    _v = vt();
			    for (_e = 0; _e < 4; ++_e)
				if (vt() == v)
				    return *this;

			    throw std::runtime_error("Segmentation<R>::Edge::operator ~(): Internal error!");
			    return *this;
			}
	Edge		next()		const	{ return ++Edge(*this); }
	Edge		prev()		const	{ return --Edge(*this); }
	Edge		conj()		const	{ return  ~Edge(*this); }

	friend class	Segmentation;

      private:
			Edge(viterator v, viterator vend)
			    :_v(v), _e(0), _vend(vend)
			{
			    while (vt() == _vend)
				++_e;
			}

	riterator	r()		const	{ return _v->_r[_e]; }
	viterator	vs()		const	{ return _v; }
	viterator&	vt()		const	{ return _v->_v[_e]; }
	void		pair(const Edge& edge) const
			{
			    vt() = edge._v;
			    edge.vt() = _v;
			}
	void		replaceRegion(riterator r, const Edge& edgeE) const
			{
			    auto	edge = *this;
			    do
			    {
				edge._v->_r[edge._e] = r;
			    } while (--(~edge) != edgeE);

			    resetRegion();
			}
	void		resetRegion() const
			{
			    r()->setEdge(*this);
			}

      private:
	viterator	_v;		//!< 親の頂点を指す反復子
	size_t		_e;		//!< 辺の番号
	const viterator	_vend;
    };

    class Region : public R
    {
      public:
	Region(const R& region, const Edge& edge)
	    :R(region), _edge(edge)			{}

	Edge	edge()				const	{ return _edge; }
	void	setEdge(const Edge& edge)		{ _edge = edge; }

      private:
	Edge	_edge;	// parnet edge
    };

  public:
    Edge		initialize(const R region[])			;
    void		clear()						;

    Edge		split(Edge& edge, const Vertex& vertex)		;
    bool		merge(Edge& edge)				;
    Edge		kill(Edge& edge)				;
    Edge		make(const Edge& edge0,
			     const Edge& edge1, const R& r)		;

    riterator		begin()			{ return _regions.begin(); }
    const_riterator	begin()		const	{ return _regions.begin(); }
    riterator		end()			{ return _regions.end(); }
    const_riterator	end()		const	{ return _regions.end(); };

  private:
    viterator	vend()			const	{ return _vertices.vend() }
    riterator	newRegion(const R& region)	{ _regions.push_front(region);
						  return _regions.begin(); }
    viterator	newVertex(const Vertex& vertex)	{ _vertices.push_front(vertex);
						  return _vertices.begin(); }
    void	deleteRegion(riterator r)	{ _regions.erase(r); }
    void	deleteVertex(viterator v)	{ _vertices.erase(v); }

  private:
    std::list<Region>	_regions;			//!< 領域のリスト
    std::list<Vertex>	_vertices;			//!< 頂点のリスト
};

template <class R> typename Segmentation<R>::Edge
Segmentation<R>::initialize(const R region[4])
{
  // 表の頂点を生成する．
    riterator	r[4];
    for (size_t e = 0; e < 4; ++e)
	r[e] = newRegion(region[e]);
    viterator	f = newVertex(V(r, vend()));

  // 裏の頂点を生成する．
    riterator	rC[4];
    for (size_t e = 0; e < 4; ++e)
	rC[e] = r[4-1-e];
    viterator	vC = newVertex(V(rC, vend()));

  // 表と裏を貼り合わせる．
    Edge	edge0(f), edge(edge0), edgeC(fC);
    --edgeC;
    do
    {
	edge.pair(--edgeC);
    } while (++edge != edge0);

    return edge0;
}

template <class R> inline void
Segmentation<R>::clear()
{
    _regions.clear();
    _vertices.clear();
}

template <class R> typename Segmentation<R>::Edge
Segmentation<R>::split(Edge& edge, const Vertex& vertex)
{
    const auto	v = newVertex(vertex);
}

template <class R> bool
Segmentation<R>::merge(Edge& edge)
{
    if (edge.valence() == 2)
    {
	const auto	edgeC = edge.conj();
	const auto	v     = edge.vs();
	~(++edge);
	edge.pair(edgeC);
	deleteVertex(v);

	edge.resetRegion();
	edgeC.resetRegion();

	return true;
    }
    else
	return false;
}

template <class R> typename Segmentation<R>::Edge
Segmentation<R>::kill(Edge& edge)
{
    const auto	edge1 = edge.conj();

    if (edge.valence() < 3 || edge1.valence() < 3)
    {
	std::cerr << "TU::Segmentation<R, V, 3u>::kill(): Too small valence!"
		  << std::endl;
	return edge;
    }

  // Remove a link between source and target vertices of edge.
    edge.vt()  = vend();
    edge1.vt() = vend();

  // Replace regions
    const auto	r = edge.r();
    --edge;
    --edge1;
    edge1.replaceRegion(edge.r(), edge);
    deleteRegion(r);

    reduce(edge);
    reduce(edge1);

    return edge1;
}

template <class R> typename Segmentation<R>::Edge
Segmentation<R>::make(const Edge& edge0, const Edge& edge1, const R& r)
{
  // Check whether edge0 and edge1 share their parent region.
    if (edge0.r() != edge1.r())
	throw std::domain_error("Segmentation<R>::make(): Given two edges share a region!");

  // Check whether edge0 and edge1 are different.
    if (edge0 == edge1)
	throw std::domain_error("Segmentation<R>::make(): Given two edges are identical!");

  // Forward edge0 and edge1 to the next link.
    auto	edge = edge0;
    edge._e = (edge._e == 3 ? 0 : edge._e + 1);
    if (edge.vt() != vend())
	throw std::domain_error("Segmentation<R>::make(): edge._v[" +
				std::to_string(edge._e) +
				"] is already occupied!");
    auto	edgeC = edge1;
    edgeC._e = (edgeC._e == 3 ? 0 : edgeC._e + 1);
    if (edgeC.vt() != vend())
	throw std::domain_error("Segmentation<R>::make(): edgeC._v[" +
				std::to_string(edgeC._e) +
				"] is already occupied!");

  // Make a link between source vertices of edge0 and edge1.
    edge.pair(edgeC);

  // Create a new region and make it a parent of vertices surrounding it.
    const auto	rnew = newRegion(r);
    edge.replaceRegion(rnew, edge);

    return edge;
}

}	// namespace plane_detector
