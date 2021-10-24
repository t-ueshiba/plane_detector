/*
 *  $Id: draw.cc,v 1.2 2011-01-30 23:39:46 ueshiba Exp $
 */
#include "edit.h"

namespace TU
{
namespace v
{
/************************************************************************
*  global functions							*
************************************************************************/
OglDC&
operator <<(OglDC& dc, const TriMesh& mesh)
{
    typedef TriMesh::const_fiterator	const_fiterator;
    
    for (const_fiterator face = mesh.fbegin(); face != mesh.fend(); ++face)
	dc << *face;
    
    return dc;
}

OglDC&
operator <<(OglDC& dc, const TriFace& face)
{
    glBegin(GL_TRIANGLES);
      glNormal3fv(face.normal().data());
      glVertex3fv(face.v(0).data());
      glVertex3fv(face.v(1).data());
      glVertex3fv(face.v(2).data());
    glEnd();

    return dc;
}

OglDC&
operator <<(OglDC& dc, const TriMesh::Edge& edge)
{
    glBegin(GL_LINES);
      glVertex3fv(edge.v().data());
      glVertex3fv(edge.next().v().data());
    glEnd();

    return dc;
}

OglDC&
operator <<(OglDC& dc, const TriVertex& v)
{
    glBegin(GL_POINTS);
      glVertex3fv(v.data());
    glEnd();

    return dc;
}

OglDC&
drawColoredMeshFaces(OglDC& dc, const TriMesh& mesh)
{
    typedef TriMesh::const_fiterator	const_fiterator;
    
    glPushAttrib(GL_CURRENT_BIT);

    struct Color
    {
	float	r, g, b;
	bool	used;
    };
    Color	color[] = {{0.0, 0.0, 1.0, false}, {0.0, 1.0, 1.0, false},
			   {1.0, 0.0, 1.0, false}, {1.0, 1.0, 0.0, false}};

    for (const_fiterator face = mesh.fbegin(); face != mesh.fend(); ++face)
	face->mark = 4;
    for (const_fiterator face = mesh.fbegin(); face != mesh.fend(); ++face)
    {
	color[0].used = color[1].used = color[2].used = color[3].used = false;
	for (int e = 0; e < 3; ++e)
	    if (face->f(e).mark != 4)
		color[face->f(e).mark].used = true;
	int	index = 0;
	while (color[index].used)
	    ++index;
	glColor3f(color[index].r, color[index].g, color[index].b);
    	dc << *face;
	face->mark = index;
    }
    for (const_fiterator face = mesh.fbegin(); face != mesh.fend(); ++face)
	face->mark = 0;

    glPopAttrib();
    
    return dc;
}

OglDC&
draw(OglDC& dc, const TriMesh::Edge& edge)
{
    glPushAttrib(GL_CURRENT_BIT);
    glPushAttrib(GL_POINT_BIT);

    glColor3f(1.0, 0.0, 0.0);
    dc << edge;
    glColor3f(0.0, 1.0, 1.0);	// light blue
    glPointSize(3.0);
    dc << edge.v();

    glPopAttrib();
    glPopAttrib();
    
    return dc;
}

OglDC&
erace(OglDC& dc, const TriMesh::Edge& edge)
{
    glPushAttrib(GL_CURRENT_BIT);
    glPushAttrib(GL_POINT_BIT);
    glColor3f(1.0, 1.0, 1.0);
    glPointSize(3.0);
    dc << edge << edge.v();
    glPopAttrib();
    glPopAttrib();
    
    return dc;
}
 
}
}
