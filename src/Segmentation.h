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

    using riterator	  = typename std::list<R>::iterator;
    using const_riterator = typename std::list<R>::const_iterator;
    using viterator	  = typename std::list<Vertex>::iterator;
    using const_riterator = typename std::list<Vertex>::const_iterator;

    class Edge;

    class Vertex
    {
	friend class	Edge;

      public:
#ifndef TU_MESH_DEBUG
	Vertex(riterator r[], viterator vend)				;
#else
	Vertex(riterator r[], viterator vend, size_t vn)		;
#endif

      private:
	viterator	self(viterator vend)			const	;

      private:
	riterator	_r[4];
	viterator	_v[4];

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

#ifdef TU_MESH_DEBUG
      public:
	const size_t	vnum;
#endif
    };

    class Edge
    {
      public:
	Edge(const Vertex& v, viterator vend)				;

	bool		operator ==(const Edge& edge)		const	;
	bool		operator !=(const Edge& edge)		const	;
	bool		commonRegion(const Edge& edge)		const	;
	size_t		valence()				const	;
	Edge&		operator ++()					;
	Edge&		operator --()					;
	Edge&		operator ~()					;
	Edge		next()					const	;
	Edge		prev()					const	;
	Edge		conj()					const	;

	friend class	Segmentation;

      private:
	Edge(viterator v, viterator vend)				;

	riterator	r()					const	;
	viterator	vt()					const	;
	viterator&	vt()						;
	void		pair(const Edge& edge)			const	;
	void		replaceRegion(riterator r,
				      const Edge& edgeE)	const	;
	
      private:
	viterator	_v;		//!< 親の頂点を指す反復子
	size_t		_e;		//!< 辺の番号
	const viterator	_vend;
    };

  public:
    Edge		initialize(const R region[])			;
    void		clear()						;

    bool		reduce(Edge& edge)				;
    Edge		kill(Edge& edge)				;
    Edge		make(const Edge& edge0,
			     const Edge& edge1, const R& r)		;

    riterator		begin()						;
    const_riterator	begin()					const	;
    riterator		end()						;
    const_riterator	end()					const	;
#ifdef TU_MESH_DEBUG
    std::ostream&	showTopology(std::ostream& out)		const	;
#endif

  private:
    viterator		vend()					const	;
    riterator		newRegion(const R& r)				;
    viterator		newVertex(const V& v)				;
    void		deleteRegion(riterator r)			;
    void		deleteVertex(viterator f)			;

  private:
    std::list<R>	_regions;			//!< 領域のリスト
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

template <class R> bool
Segmentation<R>::reduce(Edge& edge)
{
    if (edge.valence() == 2)
    {
	const auto	edgeC = edge.conj();
	const auto	v     = edge.v();
	~(++edge);
	edge.pair(edgeC);
	deleteVertex(v);

	return true;
    }
    else
	return false;
}

template <class R> typename Segmentation<R>::Edge
Segmentation<R>::kill(Edge& edge)
{
    const auto	edgeC	 = edge.conj();
    const auto	valence	 = edge.valence();
    const auto	valenceC = edgeC.valence();

    if (valence < 3 || valenceC < 3)
    {
	std::cerr << "TU::Segmentation<R, V, 3u>::kill(): Too small valence!"
		  << std::endl;
	return edge;
    }

    edge.vt()  = vend();
    edgeC.vt() = vend();

  // Replace regions
    const auto	edgeP = edge.prev();
    const auto	rP    = edgeP.r();
    edge = edgeC.prev();
    edgeP.replaceRegion(edge.r(), edge);
    deleteRegion(rP);
    
    reduce(edgeP);
    reduce(edge);
    
    return edgeP;
}

template <class R> typename Segmentation<R>::Edge
Segmentation<R>::make(const Edge& edge0, const Edge& edge1, const R& v)
{
    using namespace	std;

  // edge0とedge1が始点を共有しているかチェック．
    if (!edge0.commonRegion(edge1))
	throw domain_error("TU::Segmentation<R, V, 3u>::make(): Given two edges have no common region!");

  // edge0とedge1が同一でないかチェック．
    if (edge0 == edge1)
	throw domain_error("TU::Segmentation<R, V, 3u>::make(): Given two edges are identical!");

  // 新しい領域を作る．
    riterator	rnew = newRegion(v);

  // 新しい頂点を2つ作る．
    riterator	rp[3];
    rp[0] = rnew;
    rp[1] = edge0.r();
    rp[2] = edge0.next().r();
    viterator	v = newVertex(V(rp, vend()));
    rp[0] = edge1.r();
    rp[1] = rnew;
    rp[2] = edge1.next().r();
    viterator	vC = newVertex(V(rp, vend()));

  // edge0の始点を置き換える前にedge0とedge1の裏を保持しておく．
    Edge	edge0C(edge0.conj()), edge1C(edge1.conj());

  // [edge0, edge1)の範囲の辺の始点を新しい領域に置き換える.
    edge0.replaceRegion(rnew, edge1);

  // winged-edge構造を作る．
    Edge	edge(v), edgeC(vC);
    edge.pair(edgeC);
    (--edge ).pair(edge0);
    (--edge ).pair(edge0C);
    (--edgeC).pair(edge1);
    (--edgeC).pair(edge1C);

    return --edge;
}

template <class R> inline typename Segmentation<R>::riterator
Segmentation<R>::begin()
{
    return _regions.begin();
}

template <class P> inline typename Segmentation<R>::const_riterator
Segmentation<R>::begin() const
{
    return _regions.begin();
}

template <class R> inline typename Segmentation<R>::riterator
Segmentation<R>::end()
{
    return _regions.end();
}

template <class R> inline typename Segmentation<R>::const_riterator
Segmentation<R>::end() const
{
    return _regions.end();
}

#ifdef TU_MESH_DEBUG
template <class R> std::ostream&
Segmentation<R>::showTopology(std::ostream& out) const
{
    for (const auto& v : _vertices)
    {
	out << "Vertex[" << v->vnum << "]:";
	for (size_t e = 0; e < 4; ++e)
	    out << ' ' << v.v(e).vnum;
	out << std::endl;
    }

    return out;
}
#endif

template <class R> inline typename Segmentation<R>::viterator
Segmentation<R>::vend() const
{
    _vertices.end();
}

template <class R> inline typename Segmentation<R>::riterator
Segmentation<R>::newRegion(const R& r)
{
    _regions.push_front(r);
    return _regions.begin();
}

template <class R> inline typename Segmentation<R>::viterator
Segmentation<R>::newVertex(const V& v)
{
    _vertices.push_front(v);
    return _vertices.begin();
}

template <class R> inline void
Segmentation<R>::deleteRegion(riterator r)
{
    _regions.erase(r);
}

template <class R> inline void
Segmentation<R>::deleteVertex(viterator v)
{
    _vertices.erase(v);
}

/************************************************************************
*  class Segmentation<R>::Vertex					*
************************************************************************/
template <class R> inline
#ifndef TU_MESH_DEBUG
Segmentation<R>::Vertex::Vertex(riterator r[])
#else
Segmentation<R>::Vertex::Vertex(riterator r[], size_t fn)
    :vnum(fn)
#endif
{
    for (size_t e = 0; e < 4; ++e)
	_r[e] = r[e];
}

template <class R> typename Segmentation<R>::viterator
Segmentation<R>::Vertex::self(viterator vend) const
{
    for (auto vt : _v)
	if (vt != vend)
	{
	    for (auto v : vt->_v)		// vcのe番目の辺を介して
		if (v != vend && &(*v) == this)	// この頂点を指していたら
		    return v;			// vがこの頂点への反復子．
    
	    throw std::runtime_error("Segmentation<R>::Vertex::self(): Internal error!");
	}
    
    return vend;
}

/************************************************************************
*  class Segmentation<R>::Edge						*
************************************************************************/
template <class R> inline
Segmentation<R>::Edge::Edge(const Vertex& vertex, viterator vend)
    :Edge(vertex.self(vend), vend)
{
}

template <class R> inline bool
Segmentation<R>::Edge::operator ==(const Edge& edge) const
{
    return (_e == edge._e) && (_v == edge._v);
}

template <class R> inline bool
Segmentation<R>::Edge::operator !=(const Edge& edge) const
{
    return !(*this == edge);
}

template <class R> bool
Segmentation<R>::Edge::commonRegion(const Edge& edge) const
{
    Edge	tmp(*this);
    do
    {
	if (tmp == edge)
	    return true;
    } while (~(--tmp) != *this);

    return false;
}

template <class R> size_t
Segmentation<R>::Edge::valence() const
{
    size_t	n = 0;
    for (size_t e = 0; e < 4; ++e)
	if (_v->_v[e] != _vend)
	    ++n;
    return n;
}
    
template <class R> inline typename Segmentation<R>::Edge&
Segmentation<R>::Edge::operator ++()
{
    do
    {
	if (_e == 3)
	    _e = 0;
	else
	    ++_e;
    } while (vt() == _vend);

    return *this;
}

template <class R> inline typename Segmentation<R>::Edge&
Segmentation<R>::Edge::operator --()
{
    do
    {
	if (_e == 0)
	    _e = 3;
	else
	    --_e;
    } while (vt() == _vend);

    return *this;
}

template <class R> inline typename Segmentation<R>::Edge&
Segmentation<R>::Edge::operator ~()
{
    const auto	v = _v;
    _v = vt();
    for (_e = 0; _e < 4; ++_e)
	if (vt() == v)
	    return *this;

    throw std::runtime_error("Segmentation<R>::Edge::operator ~(): Internal error!");
    return *this;
}

template <class R> inline typename Segmentation<R>::Edge
Segmentation<R>::Edge::next() const
{
    auto	edge(*this);
    return ++edge;
}

template <class R> inline typename Segmentation<R>::Edge
Segmentation<R>::Edge::prev() const
{
    auto	edge(*this);
    return --edge;
}

template <class R> typename Segmentation<R>::Edge
Segmentation<R>::Edge::conj() const
{
    auto	edge(*this);
    return ~edge;
}

template <class R> inline
Segmentation<R>::Edge::Edge(viterator v, viterator vend)
    :_v(v), _e(0), _vend(vend)
{
    while (vt() == _vend)
	++_e;
}

template <class R> typename Segmentation<R>::riterator
Segmentation<R>::Edge::r() const
{
    return _v->_r[_e];
}

template <class R> typename Segmentation<R>::viterator
Segmentation<R>::Edge::vt() const
{
    return _v->_v[_e];
}

template <class R> typename Segmentation<R>::viterator&
Segmentation<R>::Edge::vt()
{
    return _v->_v[_e];
}

template <class R> inline void
Segmentation<R>::Edge::pair(const Edge& edge) const
{
    vt() = edge._v;	// Set source of this to target of 'edge'
    edge.vt() = _v;	// Set target of 'edge' to source of this
}

template <class R> void
Segmentation<R>::Edge::replaceRegion(riterator r, const Edge& edgeE) const
{
    auto	edge = *this;
    do
    {
	edge._v->_r[edge._e] = r;
    } while (--(~edge) != edgeE)
}

}
