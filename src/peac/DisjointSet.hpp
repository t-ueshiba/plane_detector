//
// Copyright 2014 Mitsubishi Electric Research Laboratories All
// Rights Reserved.
//
// Permission to use, copy and modify this software and its
// documentation without fee for educational, research and non-profit
// purposes, is hereby granted, provided that the above copyright
// notice, this paragraph, and the following three paragraphs appear
// in all copies.
//
// To request permission to incorporate this software into commercial
// products contact: Director; Mitsubishi Electric Research
// Laboratories (MERL); 201 Broadway; Cambridge, MA 02139.
//
// IN NO EVENT SHALL MERL BE LIABLE TO ANY PARTY FOR DIRECT,
// INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
// LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
// DOCUMENTATION, EVEN IF MERL HAS BEEN ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGES.
//
// MERL SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN
// "AS IS" BASIS, AND MERL HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE,
// SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
//
#pragma once

#include <vector>

class DisjointSet
{
  public:
		DisjointSet(int n)
		    :_parent(n), _size(n)
		{
		    for (int i = 0; i < n; ++i)
		    {
			_parent[i] = i;
			_size[i]   = 1;
		    }
		}

		~DisjointSet()			{}

    void	remove(int x)
		{
		    if(_parent[x] != x)
		    {
			--_size[Find(x)];
			_parent[x] = x;
			_size[x]   = 1;
		    }
		}

    int		getSetSize(int x)		{ return _size[Find(x)]; }

    int		Union(int x, int y)
		{
		    const auto	xRoot = Find(x);
		    const auto	yRoot = Find(y);

		    if (xRoot == yRoot)
			return xRoot;

		    const auto	xRootSize = _size[xRoot];
		    const auto	yRootSize = _size[yRoot];

		    if (xRootSize < yRootSize)
		    {
			_parent[xRoot] = yRoot;
			_size[yRoot]  += _size[xRoot];
			return yRoot;
		    }
		    else
		    {
			_parent[yRoot] = xRoot;
			_size[xRoot]  += _size[yRoot];
			return xRoot;
		    }
		}

    int		Find(int x)
		{
		    if (_parent[x] != x)
			_parent[x] = Find(_parent[x]);

		    return _parent[x];
		}

			DisjointSet()				= delete;
			DisjointSet(const DisjointSet& rhs)	= delete;
    DisjointSet&	operator=(const DisjointSet& rhs)	= delete;

  private:
    std::vector<int>	_parent;
    std::vector<int>	_size;
};
