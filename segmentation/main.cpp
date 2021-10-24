/*
 *  $Id: main.cc,v 1.3 2011-01-30 23:39:47 ueshiba Exp $
 */
#include <fstream>
#include <iomanip>
#include "TU/v/App.h"
#include "TU/v/CmdWindow.h"
#include "TU/v/CmdPane.h"
#include "TU/v/OglDC.h"
#include "edit.h"

namespace TU
{
namespace v
{
static const u_int	WinWidth = 640, WinHeight = 480;
    
/************************************************************************
*  menus and commands							*
************************************************************************/
enum	{c_Prev=100, c_Conj, c_Next, c_Kill, c_Make, c_Swap,
	 c_FaceNum, c_Edge, c_Face, c_DrawMode};

static MenuDef FileMenu[] =
{
    {"New",  M_New,  false, noSub},
    {"Open", M_Open, false, noSub},
    {"-",    M_Line, false, noSub},
    {"Quit", M_Exit, false, noSub},
    EndOfMenu
};

static CmdDef MainMenu[] =
{
    {C_MenuButton, M_File, 0, "File", FileMenu, CA_None, 0, 0, 1, 1, 0},
    EndOfCmds
};

static CmdDef RadioButtonCmds[] =
{
    {C_RadioButton, c_Edge,    0, "Edge",	    noProp, CA_None,
     0, 0, 1, 1},
    {C_RadioButton, c_Face,    0, "Face",	    noProp, CA_None,
     1, 0, 1, 1},
    EndOfCmds
};
	
static CmdDef Cmds[] =
{
    {C_Button, c_Prev,   0, "prev",   noProp, CA_None, 0, 0, 1, 1, 0},
    {C_Button, c_Conj,   0, "conj",   noProp, CA_None, 1, 0, 1, 1, 0},
    {C_Button, c_Next,   0, "next",   noProp, CA_None, 2, 0, 1, 1, 0},
  //    {C_Text,   c_FaceNum,0, "face#",  noProp, CA_None, 3, 0, 1, 1, 0},
    {C_Button, c_Kill,   0, "kill",   noProp, CA_None, 3, 0, 1, 1, 0},
    {C_Button, c_Make,   0, "make",   noProp, CA_None, 4, 0, 1, 1, 0},
    {C_Button, c_Swap,   0, "swap",   noProp, CA_None, 5, 0, 1, 1, 0},
    {C_ChoiceFrame, c_DrawMode, c_Edge, "draw mode", RadioButtonCmds, CA_None,
     6, 0, 1, 1, 0},
    EndOfCmds
};

/************************************************************************
*  class MyCanvasPane							*
************************************************************************/
class MyCanvasPane : public CanvasPane
{
  public:
    MyCanvasPane(Window& parentWin, TriMesh& mesh)
	:CanvasPane(parentWin, WinWidth, WinHeight),
	 _dc(*this),
	 _mesh(mesh),
	 _edge(*_mesh.fbegin()),
	 _edge0(_edge),
	 _edge1(_edge),
	 _v(_edge.v())							{}

    void		prev()						;
    void		next()						;
    void		conj()						;
    void		kill()						;
    void		make()						;
    void		swap()						;
    void		edgeMode()					;
    void		faceMode()					;
    
    OglDC&				dc()			{return _dc;}
    const TriMesh::Edge&	edge()		const	{return _edge;}
    
    virtual void	repaintUnderlay()				;

  protected:
    virtual void	initializeGraphics()				;
    
  private:
    enum		{MyMesh = 1};

    OglDC			_dc;
    TriMesh&		_mesh;
    TriMesh::Edge	_edge, _edge0, _edge1;
    TriVertex			_v;
};

void
MyCanvasPane::repaintUnderlay()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glCallList(MyMesh);
    glDisable(GL_DEPTH_TEST);
    draw(_dc, _edge);
    glEnable(GL_DEPTH_TEST);

    glFlush();
}

void
MyCanvasPane::initializeGraphics()
{
    using namespace	std;
    
    BoundingBox<TriVertex>	bbox = _mesh.boundingBox();
    cerr << "Mesh BoundingBox:   [" << bbox.min(0) << ", " << bbox.max(0)
	 << "]x[" << bbox.min(1) << ", " << bbox.max(1)
	 << "]x[" << bbox.min(2) << ", " << bbox.max(2) << ']'
	 << endl;

    _dc.setInternal(WinWidth/2, WinHeight/2, 1000.0, 1000.0, 0.01, 1000.0);

    Vector3d	c = (bbox.min() + bbox.max()) / 2.0, t = c;
    t[2] += 2.0 * bbox.width();
    Matrix33d	Rt;
    Rt[0][0] =  1.0;
    Rt[1][1] = -1.0;
    Rt[2][2] = -1.0;
    _dc.setExternal(t, Rt) << distance(length(t-c));

    GLfloat     position[] = {1.0, 1.0, 1.0, 0.0};
    glLightfv(GL_LIGHT0, GL_POSITION, position);
    glShadeModel(GL_FLAT);

    glEnable(GL_DEPTH_TEST);		// Normally, depth-test is enabled.
    glDisable(GL_AUTO_NORMAL);

    edgeMode();				// Make a display list for the mesh.
}

void
MyCanvasPane::prev()
{
    glDisable(GL_DEPTH_TEST);
    erace(_dc, _edge);
    draw(_dc, --_edge);
    glEnable(GL_DEPTH_TEST);
    glFlush();
}

void
MyCanvasPane::next()
{
    glDisable(GL_DEPTH_TEST);
    erace(_dc, _edge);
    draw(_dc, ++_edge);
    glEnable(GL_DEPTH_TEST);
    glFlush();
}

void
MyCanvasPane::conj()
{
    glDisable(GL_DEPTH_TEST);
    erace(_dc, _edge);
    draw(_dc, ~_edge);
    glEnable(GL_DEPTH_TEST);
    glFlush();
}

void
MyCanvasPane::kill()
{
    std::cerr << "kill: edge.f() = " << std::hex << &_edge.f() << std::endl;
    _v = _edge.v();
    _edge1 = _mesh.kill(_edge);
    _edge0 = _edge;
}

void
MyCanvasPane::make()
{
    if (_edge0 != _edge1)
	_edge = _edge0 = _edge1 = _mesh.make(_edge0, _edge1, _v);
    else
	std::cerr << "Cannot make an edge!" << std::endl;
}

void
MyCanvasPane::swap()
{
    glDisable(GL_DEPTH_TEST);
    erace(_dc, _edge);
    _mesh.swap(_edge);
    draw(_dc, _edge);
    glEnable(GL_DEPTH_TEST);
    glFlush();
}

void
MyCanvasPane::edgeMode()
{
  /* Light settings. */
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_COLOR_MATERIAL);
    
    glDeleteLists(MyMesh, 1);
    glNewList(MyMesh, GL_COMPILE);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glColor3f(1.0, 1.0, 1.0);
      _dc << _mesh;
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.0, 1.0);
      glColor3f(0.0, 0.0, 0.0);
      _dc << _mesh;
      glDisable(GL_POLYGON_OFFSET_FILL);
    glEndList();
}

void
MyCanvasPane::faceMode()
{
  /* Light settings. */
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
    glDeleteLists(MyMesh, 1);
    glNewList(MyMesh, GL_COMPILE);
      drawColoredMeshFaces(_dc, _mesh);
    glEndList();
}

/************************************************************************
*  class MyCmdWindow							*
************************************************************************/
class MyCmdWindow : public CmdWindow
{
  public:
    MyCmdWindow(App&				parentApp,
		const char*			name,
		const XVisualInfo*		vinfo,
		TriMesh&		mesh)		;

    virtual void	callback(CmdId, CmdVal)			;

  private:
    CmdPane		_menu;
    CmdPane		_cmd;
    MyCanvasPane	_canvas;
};

MyCmdWindow::MyCmdWindow(App& parentApp, const char* name,
			 const XVisualInfo* vinfo,
			 TriMesh& mesh)
    :CmdWindow(parentApp, name, vinfo, Colormap::RGBColor, 16, 0, 0),
     _menu(*this, MainMenu), _cmd(*this, Cmds), _canvas(*this, mesh)
{
    _menu.place(0, 0, 1, 1);
    _cmd.place(0, 1, 1, 1);
    _canvas.place(0, 2, 1, 1);
    
    show();

    _canvas.repaintUnderlay();
}

void
MyCmdWindow::callback(CmdId id, CmdVal val)
{
    bool	remake = false;
    
    switch (id)
    {
      case M_Exit:
	app().exit();
	break;

      case c_Prev:
	_canvas.prev();
	break;

      case c_Conj:
	_canvas.conj();
	break;
	
      case c_Next:
	_canvas.next();
	break;
	
      case c_Kill:
	_canvas.kill();
	remake = true;
	break;

      case c_Make:
	_canvas.make();
	remake = true;
	break;

      case c_Swap:
	_canvas.swap();
	remake = true;
	break;

      case c_DrawMode:
	remake = true;
	break;
    }

    if (remake)
    {
	switch (_cmd.getValue(c_DrawMode))
	{
	  case c_Edge:
	    _canvas.edgeMode();
	    break;
	
	  case c_Face:
	    _canvas.faceMode();
	    break;
	}
	_canvas.dc() << repaint;
    }
  /*
    ostrstream	s;
    s << setw(4) << _canvas.edge().f().fnum << ends;
    _cmd.setString(c_FaceNum, s.str());*/
}
 
}
}
/************************************************************************
*  global functions							*
************************************************************************/
int
main(int argc, char* argv[])
{
    using namespace	std;
    using namespace	TU;
    
  /* Process command line. */
    char		*infile = 0;
    extern char*	optarg;
    for (int c; (c = getopt(argc, argv, "i:S")) != EOF; )
	switch (c)
	{
	  case 'i':
	    infile = optarg;
	    break;
	}

  /* Read mesh file. */    
    TriMesh	mesh;
    ifstream	fin(infile);
    istream&	in = (infile == 0 ? (istream&)cin : (istream&)fin);
    if (!in)
    {
	cerr << "Cannot open the mesh file (" << infile << ")!" << endl;
	return 1;
    }
    cerr << "Reading meshes...  ";
    in >> mesh;
    cerr << "Done!" << endl;
    
  /* Select visual for OpenGL. */
    v::App		vapp(argc, argv);
    int			attrs[] = {GLX_RGBA,
				   GLX_RED_SIZE,	1,
				   GLX_GREEN_SIZE,	1,
				   GLX_BLUE_SIZE,	1,
				   GLX_DEPTH_SIZE,	1,
				   None};
    XVisualInfo*	vinfo = glXChooseVisual(vapp.colormap().display(),
						vapp.colormap().vinfo().screen,
						attrs);
    if (vinfo == 0)
    {
	cerr << "No appropriate visual!!" << endl;
	return 1;
    }

  /* Main loop. */
    v::MyCmdWindow	myWin0(vapp, "mesh viewer", vinfo, mesh);
    vapp.run();
    XFree(vinfo);
    
    cout << mesh;
  //mesh.saveSTL(cout);
    
    return 0;
}
