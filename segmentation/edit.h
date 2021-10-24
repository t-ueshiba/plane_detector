/*
 *  $Id: edit.h,v 1.3 2011-01-30 23:39:46 ueshiba Exp $
 */
#ifndef __edit_h
#define __edit_h

#include "TU/Segmentation.h"
#include "TU/Manip.h"
#include "TU/v/OglDC.h"

namespace plane_detector
{
/************************************************************************
*  class TriFace							*
************************************************************************/
class Region : public Segmentation<TriVertex, TriFace>::
{
  private:
    typedef Mesh<TriVertex, TriFace>	mesh_t;
    typedef mesh_t::Face		super;

  public:
    typedef mesh_t::viterator		viterator;
    typedef mesh_t::Edge		Edge;

  public:
#ifndef TU_MESH_DEBUG
    TriFace(viterator v[])
	:super(v), mark(0)					{}
#else
    TriFace(viterator v[], int fn)
	:super(v, fn), mark(0)					{}

    void		print(std::ostream& out)	const	;
#endif
    Vector3f		normal()			const	;
    Vector3f		G()				const	;

  public:
    mutable u_int	mark;
};

#ifdef TU_MESH_DEBUG
inline void
TriFace::print(std::ostream& out) const
{
    out << "  face " << std::setw(4) << fnum << ':' << G();
}
#endif

inline Vector3f
TriFace::normal() const
{
    Vector3f	normal = (v(1) - v(0)) ^ (v(2) - v(0));
    return normalize(normal);
}

inline Vector3f
TriFace::G() const
{
    return (v(0) + v(1) + v(2))/3;
}

using MySegmentation	= Segmentation<>;

}	// namespace plane_detector


namesace TU
{
namespace v
{
/************************************************************************
*  drawing functions							*
************************************************************************/
OglDC&	operator <<(OglDC& dc,
		    const plane_detector::MySegmentation& segmentation)	;
OglDC&	operator <<(OglDC& dc, const TriFace& face)			;
OglDC&	operator <<(OglDC& dc, const TriMesh::Edge& edge)		;
OglDC&	operator <<(OglDC& dc, const TriVertex& v)			;

OglDC&	drawColoredMeshFaces(OglDC& dc, const TriMesh& mesh)		;

OglDC&	draw(OglDC& dc, const TriMesh::Edge& edge)			;
OglDC&	erace(OglDC& dc, const TriMesh::Edge& edge)			;
}
}

#endif	/* !__edit_h	*/
