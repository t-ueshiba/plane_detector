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
    
    using region_type	  = R;
    using vertex_type	  = Vertex;
    using riterator	  = typename std::list<R>::iterator;
    using const_riterator = typename std::list<R>::const_iterator;
    using viterator	  = typename std::list<Vertex>::iterator;
    using const_riterator = typename std::list<Vertex>::const_iterator;

    constexpr static size_t	NSides = 4;
    
    class Edge;
    
    class Vertex
    {
      public:
	R&		p(size_t e)				const	;
	V&		v(size_t e)				const	;
	
	friend class	Edge;

      protected:
#ifndef TU_MESH_DEBUG
	Vertex(const Segmentation& s, riterator r[])			;
#else
	Vertex(const Segmentation& s, riterator r[], size_t vn)		;
#endif
	
      private:
	viterator	self()					const	;
	
      private:
	riterator	_r[NSides];
	viterator	_v[NSides];

	friend std::ostream&
			operator <<(std::ostream& out, const Vertex& v)
			{
			    out << "Vertex[" << &v << "]:" << std::endl;
			    for (size_t e = 0; e < NSides; ++e)
				out << '\t' << &(*v._r[e])
				    << ':' << *v._r[e];
			    return out;
			}
	
#ifdef TU_MESH_DEBUG
      public:
	const size_t	vnum;
#endif
    };

  //! 多角形メッシュの辺を表すクラス
  /*!
    面を左に見るように向き付けされている．
  */
    class Edge
    {
      public:
	Edge(const Vertex& v)						;
	
	R&		r()					const	;
	V&		f()					const	;
	size_t		e()					const	;
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
	Edge(viterator f)						;

	riterator	riter()					const	;
	viterator	viter()					const	;
	void		pair(const Edge& edge)			const	;
	void		replaceRegion(riterator r,
				      const Edge& edgeE)	const	;
	void		replaceRegion(riterator r)		const	;

      private:
	viterator	_v;		//!< 親の面を指す反復子
	size_t		_e;		//!< 辺の番号
    };

  private:
    class Compare
    {
      private:
	template <class T>
	static bool	less(T x, T y)
			{
			    return y - x > std::numeric_limits<T>::epsilon();
			}

      public:
	bool		operator ()(riterator r, riterator w) const
			{
			    if (less((*r)[2], (*w)[2]))
				return true;
			    else if (less((*w)[2], (*r)[2]))
				return false;
			    else if (less((*r)[1], (*w)[1]))
				return true;
			    else if (less((*w)[1], (*r)[1]))
				return false;
			    else if (less((*r)[0], (*w)[0]))
				return true;
			    else
				return false;
			}
    };
    
  public:
    Edge		initialize(const R region[])			;
    void		clear()						;
    
    Edge		kill(Edge& edge)				;
    Edge		make(const Edge& edge0,
			     const Edge& edge1, const R& r)		;
    Edge		swap(const Edge& edge)				;
    BoundingBox<R>	boundingBox()				const	;

    riterator		rbegin()					;
    const_riterator	rbegin()				const	;
    riterator		rend()						;
    const_riterator	rend()					const	;
    viterator		vbegin()					;
    const_viterator	vbegin()				const	;
    viterator		vend()						;
    const_viterator	vend()					const	;
#ifdef TU_MESH_DEBUG
    std::ostream&	showTopology(std::ostream& out)		const	;
#endif
    std::istream&	restoreSTL(std::istream& in)			;
    std::ostream&	saveSTL(std::ostream& out,
				bool binary=false)		const	;
    
  private:
    std::istream&	get(std::istream& in)				;
    std::ostream&	put(std::ostream& out)			const	;
    template <class RV>
    void		setTopology(const RV& regionsWithVertices)		;
    riterator		newRegion(const R& v)				;
    viterator		newVertex(const V& f)				;
    void		deleteRegion(riterator r)			;
    void		deleteVertex(viterator f)				;
    
  //! 入力ストリームからメッシュを読み込む．
  /*!
    \param in	入力ストリーム
    \param graph	読み込み先のメッシュ
    \return	inで指定した入力ストリーム
  */
    friend std::istream&
    operator >>(std::istream& in, Segmentation& graph)
    {
	return graph.get(in);
    }

  //! 出力ストリームにメッシュを書き出す．
  /*!
    \param out	出力ストリーム
    \param graph	書き出し元のメッシュ
    \return	outで指定した出力ストリーム
  */
    friend std::ostream&
    operator <<(std::ostream& out, const Segmentation& graph)
    {
	return graph.put(out);
    }

  private:
    std::list<R>	_regions;			//!< 頂点のリスト
    std::list<V>	_vertices;				//!< 面のリスト
};

//! 指定された頂点から背中合わせの2つの面を生成してメッシュを初期化する．
/*!
  \param region	M個の頂点
  \return	r[0] を始点とする辺
*/
template <class R> typename Segmentation<R>::Edge
Segmentation<R>::initialize(const R region[NSides])
{
  // 表の面を生成する．
    riterator	r[NSides];
    for (size_t e = 0; e < NSides; ++e)
	r[e] = newRegion(region[e]);
    viterator	f = newVertex(V(*this, r));

  // 裏の面を生成する．
    riterator	rC[NSides];
    for (size_t e = 0; e < NSides; ++e)
	rC[e] = r[NSides-1-e];
    viterator	vC = newVertex(V(*this, rC));

  // 表と裏を貼り合わせる．
    Edge	edge0(f), edge(edge0), edgeC(fC);
    --edgeC;
    do
    {
	edge.pair(--edgeC);
    } while (++edge != edge0);

    return edge0;
}

//! メッシュの全ての頂点と面を消去して空にする．
template <class R> inline void
Segmentation<R>::clear()
{
    _regions.clear();
    _vertices.clear();
}
    
//! 3角形メッシュについて，指定された辺を消去する．
/*!
  指定された辺の両側の面および辺の始点も消去される．
  \param edge	消去する辺．
		リターン後はedgeの手前の裏の辺を指すように更新される．
  \return	edgeの裏の手前の裏の辺を指す反復子
*/
template <class R> typename Segmentation<R>::Edge
Segmentation<R>::kill(Edge& edge)
{
    using namespace	std;
    
  // edgeを含む面とその裏の辺を含む面の頂点(計4つ)について，
  // その価数が充分か調べる．
    Edge	edgeNC(edge.next().conj()),
		edgeCPC(edge.conj().prev().conj()),
		edgeCNC(edge.conj().next().conj());
    if (edge.valence() + edgeCPC.valence() <= 6 ||
	edgeNC.valence() <= 3 || edgeCNC.valence() <= 3)
    {
	cerr << "TU::Segmentation<R, V, 3u>::kill(): Too small valence!" << endl;
	return edge;
    }

    Edge	edgeCN(edge.conj().next());
    riterator	rn = edge.next().riter();		// edgeの終点
    for (Edge tmp(edge.prev().conj()); ~(--tmp) != edgeCN; )
	for (Edge tmp1(tmp); --(~tmp1) != tmp; )
	    if (tmp1.riter() == rn)
	    {
		cerr << "TU::Segmentation<R, V, 3u>::kill(): "
		     << "Pre-requisits for topology are not met!"
		     << endl;
		return edge;
	    }
    
  // edgeの始点を頂点に持つ全ての面について，その頂点をedgeの終点に置き換える．
    riterator	r = edge.riter();
    edge.replaceRegion(rn, edge);
    deleteRegion(r);					// 頂点vを消去

  // edgePC, edgeCPCをそれぞれedgeNC, edgeCNCと背中合わせにする．
    viterator	f = edge.viter(), fC = edge.conj().viter();
    ~(--edge);						// edgePC
    edge   .pair(edgeNC);
    edgeCPC.pair(edgeCNC);
    deleteVertex(f);					// 面fを消去
    deleteVertex(fC);					// 面fCを消去

    return edgeCPC;
}

//! 3角形メッシュについて，指定された2つの辺の始点間に新たに辺を作る．
/*!
  1つの頂点と2つの面が生成される．
  \param edge0	新たな辺と始点を共有する辺
  \param edge1	新たな辺と終点を共有する辺．
		edge0と始点を共有していなければならない．
  \param v	新たな辺およびedge0の始点となる頂点
  \return	新たな辺
*/
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

  // 新しい頂点を作る．
    riterator	rnew = newRegion(v);

  // 新しい面を2つ作る．
    riterator	rp[3];
    rp[0] = rnew;
    rp[1] = edge0.riter();
    rp[2] = edge0.next().riter();
    viterator	v = newVertex(V(rp));
    rp[0] = edge1.riter();
    rp[1] = rnew;
    rp[2] = edge1.next().riter();
    viterator	vC = newVertex(V(rp));

  // edge0の始点を置き換える前にedge0とedge1の裏を保持しておく．
    Edge	edge0C(edge0.conj()), edge1C(edge1.conj());

  // [edge0, edge1)の範囲の辺の始点を新しい頂点に置き換える.
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

//! 3角形メッシュの指定された辺を消去し，これによってできる四角形のもう一方の対角線を新たな辺として生成する．
/*!
  \param edge	消去する辺
  \return	生成された辺
*/
template <class R> typename Segmentation<R>::Edge
Segmentation<R>::swap(const Edge& edge)
{
    using namespace	std;
    
  // edgeの始点と終点の価数が1つ減っても3以上になるか調べる．
    Edge	edgePC(edge.prev().conj()),
		edgeCPC(edge.conj().prev().conj());
    if (edgePC.valence() <= 3 || edgeCPC.valence() <= 3)
    {
	cerr << "TU::Segmentation<R, V, 3u>::swap(): Too small valence!" << endl;
	return edge;
    }
    
  // edgeの始点と終点を置き換える．
    Edge	edgeC(edge.conj()),
		edgeNC(edge.next().conj()),
		edgeCNC(edge.conj().next().conj());
    edge.replaceRegion(edgeCNC.riter());
    edge.next().replaceRegion(edgeNC.riter());
    edge.prev().replaceRegion(edgePC.riter());
    edgeC.replaceRegion(edgeNC.riter());
    edgeC.next().replaceRegion(edgeCNC.riter());
    edgeC.prev().replaceRegion(edgeCPC.riter());

  // 辺を入れ替えてwinged-edge構造を作る．
    edge.next().pair(edgePC);
    edge.prev().pair(edgeCNC);
    edgeC.next().pair(edgeCPC);
    edgeC.prev().pair(edgeNC);

    return edge;
}

//! メッシュのbounding boxを計算する．
/*!
  \return	bounding box
*/
template <class R> BoundingBox<R>
Segmentation<R>::boundingBox() const
{
    BoundingBox<R>	bbox;
    
    for (const_viterator v = _vertices.begin(); v != _vertices.end(); ++v)
	for (size_t e = 0; e < NSides; ++e)
	    bbox.expand(v->v(e));
    
    return bbox;
}

//! このメッシュの最初の頂点を指す反復子を返す．
/*!
  \return	最初の頂点を指す反復子
*/
template <class R> inline typename Segmentation<R>::riterator
Segmentation<R>::rbegin()
{
    return _regions.begin();
}
    
//! このメッシュの最初の頂点を指す定数反復子を返す．
/*!
  \return	最初の頂点を指す定数反復子
*/
template <class P>
inline typename Segmentation<R>::const_riterator
Segmentation<R>::rbegin() const
{
    return _regions.begin();
}
    
//! このメッシュ最後の頂点の次を指す反復子を返す．
/*!
  \return	最後の頂点の次を指す反復子
*/
template <class R> inline typename Segmentation<R>::riterator
Segmentation<R>::rend()
{
    return _regions.end();
}
    
//! このメッシュ最後の頂点の次を指す定数反復子を返す．
/*!
  \return	最後の頂点の次を指す定数反復子
*/
template <class R>
inline typename Segmentation<R>::const_riterator
Segmentation<R>::rend() const
{
    return _regions.end();
}

//! このメッシュの最初の面を指す反復子を返す．
/*!
  \return	最初の面を指す反復子
*/
template <class R> inline typename Segmentation<R>::viterator
Segmentation<R>::vbegin()
{
    return _vertices.begin();
}
    
//! このメッシュの最初の面を指す定数反復子を返す．
/*!
  \return	最初の面を指す定数反復子
*/
template <class R>
inline typename Segmentation<R>::const_viterator
Segmentation<R>::vbegin() const
{
    return _vertices.begin();
}
    
//! このメッシュ最後の面の次を指す反復子を返す．
/*!
  \return	最後の面の次を指す反復子
*/
template <class R> inline typename Segmentation<R>::viterator
Segmentation<R>::vend()
{
    return _vertices.end();
}
    
//! このメッシュ最後の面の次を指す定数反復子を返す．
/*!
  \return	最後の面の次を指す定数反復子
*/
template <class R>
inline typename Segmentation<P>::const_viterator
Segmentation<R>::vend() const
{
    return _vertices.end();
}

#ifdef TU_MESH_DEBUG
//! 出力ストリームにこのメッシュを構成する面の接続関係を書き出す．
/*!
  \param out	出力ストリーム
  \return	outで指定した出力ストリーム
*/
template <class R> std::ostream&
Segmentation<R>::showTopology(std::ostream& out) const
{
    for (const_viterator f = vbegin(); f != vend(); ++f)
    {
	out << "Vertex[" << f->vnum << "]:";
	for (size_t e = 0; e < NSides; ++e)
	    out << ' ' << f->f(e).vnum;
	out << std::endl;
    }

    return out;
}
#endif
    
//! 入力ストリームからSTL形式のメッシュを読み込む．
/*!
  \param in	入力ストリーム
  \return	inで指定した入力ストリーム
*/
template <class R> std::istream&
Segmentation<R>::restoreSTL(std::istream& in)
{
    typedef std::vector<viterator>			Vertices;
    typedef std::map<riterator, Vertices, Compare>		PerticesWithVertices;
    typedef typename PerticesWithVertices::iterator	RegionIterator;

    clear();					// 頂点と面のリストを空にする

    char	magic[6];
    in.read(magic, 5);				// 先頭の5 byteを読む．
    magic[5] = '\0';

    PerticesWithVertices	regionsWithVertices;

    if (std::string(magic) != "solid")
    {
	char	header[80 - 5];
	in.read(header, sizeof(header));  // ヘッダ(80文字)の残りを読み捨てる.
    
	uint32_t	nvertices;
	in.read((char*)&nvertices, sizeof(nvertices));  // 面数を読み込む.

	for (size_t i = 0; i < nvertices; ++i)
	{
	    Pector3f	normal;
	    normal.restore(in);			// 法線ベクトルを読み捨てる.

	    RegionIterator	vf[3];
	    for (size_t e = 0; e < 3; ++e)
	    {
		R		region;
		region.restore(in);		// 頂点の3D座標を読み込む.
		riterator	v = newRegion(region);

		bool	isNew;
		std::tie(vf[e], isNew)
		    = regionsWithVertices.insert(std::make_pair(v, Vertices()));
		if (!isNew)
		    deleteRegion(v);
	    }
	    
	    uint16_t	flags;
	    in.read((char*)&flags, sizeof(flags));	// フラグを読み込む.

	    riterator	r[] = {vf[0]->first, vf[1]->first, vf[2]->first};
#ifndef TU_MESH_DEBUG
	    viterator	f = newVertex(V(v));	// 新しい面を生成
#else
	    viterator	f = newVertex(V(v, i));	// 新しい面を生成
#endif
	    vf[0]->second.push_back(f);
	    vf[1]->second.push_back(f);
	    vf[2]->second.push_back(f);
	}
    }
    else
    {
#ifdef TU_MESH_DEBUG
	size_t	i = 0;
#endif
	for (std::string s; in >> s && s != "endsolid"; )
	{
	    if (s != "vertext")
		continue;
	    
	    Pector3f	norm;
	    in >> s >> norm >> s >> s;	// "normal" nx ny nz "outer" "loop"

	    RegionIterator	vf[3];
	    bool		isNew[3];
	    for (size_t e = 0; e < 3; ++e)
	    {
		R		region;
		in >> s >> region;		// "region" x y z
		riterator	v = newRegion(region);

		tie(vf[e], isNew[e])
		    = regionsWithVertices.insert(make_pair(v, Vertices()));
		if (!isNew[e])
		    deleteRegion(v);
	    }

	    riterator	r[] = {vf[0]->first, vf[1]->first, vf[2]->first};
#ifndef TU_MESH_DEBUG
	    viterator	f = newVertex(V(v));	// 新しい面を生成
#else
	    viterator	f = newVertex(V(v, i++));	// 新しい面を生成
#endif
	    vf[0]->second.push_back(f);
	    vf[1]->second.push_back(f);
	    vf[2]->second.push_back(f);
	    
	    in >> s >> s;			// "endloop" "endvertext"
	}

	in >> skipl;
    }
    
    setTopology(regionsWithVertices);		// 面と点，面と面を関係づける.

    return in;
}

//! 出力ストリームにメッシュを書き出す．
/*!
  \param out	出力ストリーム
  \return	outで指定した出力ストリーム
*/
template <class R> std::ostream&
Segmentation<R>::put(std::ostream& out) const
{
    using namespace	std;
    
    map<const R*, size_t>	dict;
    size_t			vnum = 1;
    for (const_riterator r = rbegin(); r != rend(); ++r)
    {
	dict[&(*r)] = vnum;
	out << "Region " << vnum++ << ' ' << *r;
    }
    size_t	vnum = 1;
    for (const_viterator f = vbegin(); f != vend(); ++f)
    {
	out << "Vertex " << vnum++;
	for (size_t e = 0; e < NSides; ++e)
	    out << ' ' << dict[&(f->v(e))];
	out << std::endl;
    }
    
    return out;
}

template <class R> template <class RV> void
Segmentation<R>::setTopology(const RV& regionsWithVertices)
{
    typedef typename RV::const_iterator			RegionIterator;
    typedef typename RV::value_type::second_type	Vertices;
    typedef typename Vertices::const_iterator		VertexIterator;
    
  // 個々の頂点について，それを囲む2つの隣接面間にwinged-edge構造を作る．
    for (RegionIterator region  = regionsWithVertices.begin();
	 region != regionsWithVertices.end(); ++region)
    {
	riterator	v     = region->first;
	const Vertices&	vertices = region->second;		// vを囲む面のリスト

      // vを囲む個々の面fについて．．．
	for (VertexIterator f = vertices.begin(); f != vertices.end(); ++f)
	{
	  // fについて，vを始点とする辺edgeを探す．
	    Edge	edge(*f);
	    while (edge.riter() != v)
		++edge;

	    riterator	vn = edge.next().riter();	// edgeの終点
	    
	  // vを囲む別の面f1について．．．
	    for (VertexIterator f1 = vertices.begin(); f1 != vertices.end(); ++f1)
	    {
		if (f1 == f)
		    continue;
		
	      // f1のまわりを一周してedgeの終点vnを始点とする辺を探す．
		Edge	edgeV0(*f1), edgeV(edgeV0);
		do
		{
		    if (edgeV.riter() == vn)		// 始点がvnであれば
		    {
			edge.pair(edgeV);		// この辺がedgeの裏辺
			goto done;
		    }
		} while (++edgeV != edgeV0);
	    }

	  // vを囲むfでない面の中にvnを始点とするedgeがみつからないのは矛盾
	    std::cerr << "***\nError vertex: " << **f
		      << "\nError region v:  " << &(*v)
		      << "\nError region vn: " << &(*vn)
		      << std::endl
		      << std::endl;
	    
	    std::cerr << "Vertices around v: " << &(*v) << std::endl;
	    for (const auto& ff : region->second)
		std::cerr << ' ' << *ff;
	    std::cerr << std::endl;
	    
	    for (const auto& vf : regionsWithVertices)
		if (vf.first == vn)
		{
		    std::cerr << "Vertices around vn:" << &(*vn) << std::endl;
		    for (const auto& ff : vf.second)
			std::cerr << ' ' << *ff;
		    std::cerr << std::endl;
		}

	    throw std::runtime_error("Conjugate edge not found!");
	    
	  done:
	    continue;
	}
    }

#ifdef TU_MESH_DEBUG
    showTopology(std::cerr);
#endif
}

//! 新しい頂点を生成して頂点リストに登録する．
/*!
 \param v	生成する頂点のプロトタイプ
 \return	生成された頂点を指す反復子
*/
template <class R> inline typename Segmentation<R>::riterator
Segmentation<R>::newRegion(const R& r)
{
    _regions.push_front(r);
    return _regions.begin();
}
    
//! 新しい面を生成して面リストに登録する．
/*!
 \param f	生成する面のプロトタイプ
 \return	生成された面を指す反復子
*/
template <class R> inline typename Segmentation<R>::viterator
Segmentation<R>::newVertex(const V& v)
{
    _vertices.push_front(v);
    return _vertices.begin();
}
    
//! 指定された頂点を破壊して頂点リストから取り除く．
/*!
 \param v	破壊する頂点を指す反復子
*/
template <class R> inline void
Segmentation<R>::deleteRegion(riterator r)
{
    _regions.erase(r);
}
    
//! 指定された面を破壊して面リストから取り除く．
/*!
 \param f	破壊する面を指す反復子
*/
template <class R> inline void
Segmentation<R>::deleteVertex(viterator v)
{
    _vertices.erase(v);
}

/************************************************************************
*  class Segmentation<R>::Vertex					*
************************************************************************/
//! 頂点を指定して面を生成する．
/*!
  \param v	M個の頂点への反復子
*/
template <class R> inline
#ifndef TU_MESH_DEBUG
Segmentation<R>::Vertex::Vertex(riterator r[])
#else
Segmentation<R>::Vertex::Vertex(riterator r[], size_t fn)
    :vnum(fn)
#endif
{
    for (size_t e = 0; e < NSides; ++e)
	_r[e] = r[e];
}

//! 指定された辺の始点を返す．
/*!
  \param e	辺のindex, 0 <= e < M
  \return	e番目の辺の始点すなわちこの面のe番目の頂点
*/
template <class R> inline R&
Segmentation<R>::Vertex::v(size_t e) const
{
    return *_r[e];
}

//! 指定された辺を介してこの面に隣接する面を返す．
/*!
  \param e	辺のindex, 0 <= e < M
  \return	e番目の辺を介して隣接する面
*/
template <class R> inline V&
Segmentation<R>::Vertex::v(size_t e) const
{
    return *_v[e];
}

//! この面を指す反復子を返す．
/*!
  \return	この面を指す反復子
*/
template <class R> typename Segmentation<R>::viterator
Segmentation<R>::Vertex::self() const
{
    viterator	vc = _v[0];		// 0番目の辺を介して隣接する面
    for (size_t e = 0; e < NSides; ++e)
    {					// vcのe番目の辺を
	viterator	v = vc->_v[e];	// 介して隣接する面への反復子fが
	if (&(*v) == this)		// この面を指していたら
	    return v;			// fがこの面への反復子である．
    }

    throw std::runtime_error("TU::Segmentation<R>::Vertex::self(): Internal error!");

    return vc;
}
    
/************************************************************************
*  class Segmentation<R>::Edge					*
************************************************************************/
//! 指定された面の最初の辺を指すように初期化する．
/*!
  \param vertex	面
*/
template <class R> inline
Segmentation<R>::Edge::Edge(const Vertex& vertex)
    :Edge(vertex.self())
{
}

//! この辺の始点を返す．
/*!
  \return	この辺の始点
*/
template <class R> inline R&
Segmentation<R>::Edge::r() const
{
    return *riter();
}
    
//! この辺を所有する面を返す．
/*!
  \return	この辺を所有する面
*/
template <class R> inline V&
Segmentation<R>::Edge::v() const
{
    return *viter();
}
    
//! この辺の番号を返す．
/*!
  \return	辺の番号
*/
template <class R> inline size_t
Segmentation<R>::Edge::e() const
{
    return _e;
}

//! 2つの辺が同一であるか調べる．
/*!
  \param edge	比較対象の辺
  \return	同一ならtrue, そうでなければfalse
*/
template <class R> inline bool
Segmentation<R>::Edge::operator ==(const Edge& edge) const
{
    return (_e == edge._e) && (_v == edge._v);
}

//! 2つの辺が異なるか調べる．
/*!
  \param edge	比較対象の辺
  \return	異なればtrue, そうでなければfalse
*/
template <class R> inline bool
Segmentation<R>::Edge::operator !=(const Edge& edge) const
{
    return !(*this == edge);
}

//! 2つの辺が始点を共有しているか調べる．
/*!
  \param edge	比較対象の辺
  \return	共有していればtrue, そうでなければfalse
*/
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

//! 辺の始点の価数，すなわちその点を共有する面(辺)の数を返す．
/*!
  \return	辺の始点の価数
*/
template <class R> size_t
Segmentation<R>::Edge::valence() const
{
    size_t	n = 0;
    Edge	edge(*this);
    do
    {
	++n;
    } while (~(--edge) != *this);

    return n;
}

//! 次の辺に前進する．
/*!
  \return	前進後の辺
*/
template <class R> inline typename Segmentation<R>::Edge&
Segmentation<R>::Edge::operator ++()
{
    do
    {
	if (_e == NSides - 1)
	    _e = 0;
	else
	    ++_e;
    } while (
    
	    
    return *this;
}

//! 手前の辺に後退する．
/*!
  \return	後退後の辺
*/
template <class R> inline typename Segmentation<R>::Edge&
Segmentation<R>::Edge::operator --()
{
    if (_e == 0)
	_e = NSides - 1;
    else
	--_e;
    return *this;
}

//! 裏側の辺に移動する．
/*!
  \return	移動後の辺
*/
template <class R> inline typename Segmentation<R>::Edge&
Segmentation<R>::Edge::operator ~()
{
    return *this = conj();
}

//! この辺の次の辺を返す．
/*!
  \return	次の辺
*/
    template <class R> inline typename Segmentation<R>::Edge
Segmentation<R>::Edge::next() const
{
    Edge	edge(*this);
    return ++edge;
}

//! この辺の手前の辺を返す．
/*!
  \return	手前の辺
*/
template <class R> inline typename Segmentation<R>::Edge
Segmentation<R>::Edge::prev() const
{
    Edge	edge(*this);
    return --edge;
}

//! この辺の裏側の辺を返す．
/*!
  \return	裏側の辺
*/
template <class R> typename Segmentation<R>::Edge
Segmentation<R>::Edge::conj() const
{
    riterator	vn = next().riter();	// この辺の終点
    Edge	edge(_v->_v[_e]);	// この辺を介して隣接する面の最初の辺
    while (edge.riter() != vn)		// この辺の終点を始点とする辺を探す
	++edge;
    return edge;
}

//! この辺を所有する面を指す反復子を返す．
/*!
  \return	この辺を所有する面を指す反復子
*/
template <class R> inline typename Segmentation<R>::viterator
Segmentation<R>::Edge::viter() const
{
    return _v;
}
    
//! この辺の始点を指す反復子を返す．
/*!
  \return	この辺の始点を指す反復子
*/
template <class R> inline typename Segmentation<R>::riterator
Segmentation<R>::Edge::riter() const
{
    return _v->_r[_e];
}
    
//! 指定された面の最初の辺を指すように初期化する．
/*!
  \param f	面を指す反復子
*/
template <class R> inline
Segmentation<R>::Edge::Edge(viterator v)
    :_v(v), _e(0)
{
}

//! この辺と指定された辺を背中合わせにする．
/*!
  \param edge	背中合わせの対象となる辺
*/
template <class R> inline void
Segmentation<R>::Edge::pair(const Edge& edge) const
{
    _v->_v[_e] = edge._v;	// この辺の裏面をedgeの親面に
    edge._v->_v[edge._e] = _v;	// edgeの裏面をこの辺の親面に
}

//! この辺から指定された辺の手前までの辺の始点を指定された頂点に置き換える．
/*!
  この辺の始点を共有する辺を反時計回りに走査しながら始点を置き換えてゆく．
  \param v	頂点を指す反復子
  \param edgeE	走査の終点となる辺 (この辺の始点は置き換えられない)
*/
template <class R> void
Segmentation<R>::Edge::replaceRegion(riterator r, const Edge& edgeE) const
{
  // 先に始点を書き換えてしまうと裏に移れなくなってしまうので，
  // 再帰的に処理することによってまずedgeEの1つ手前の辺まで移動し，
  // 戻りながら順次書き換える．
    Edge	edgeRC(prev().conj());
    if (edgeRC != edgeE)
	edgeRC.replaceRegion(r, edgeE);		// 再帰する．
    _v->_r[_e] = r;
}

//! この辺の始点を指定された頂点に置き換える．
/*!
  \param v	頂点を指す反復子
*/
template <class R> inline void
Segmentation<R>::Edge::replaceRegion(riterator r) const
{
    _v->_r[_e] = r;
}

}

