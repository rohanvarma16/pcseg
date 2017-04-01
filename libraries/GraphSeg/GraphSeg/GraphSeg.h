//Composed by Su Dongcai on 2009/11/15
//If you have any suggestions, questions, and bug reports etc please feel free
//to contact me (suntree4152@gmail.com)

//Copyright (c) 2009, Su Dongcai
//All rights reserved.

//Redistribution and use in source and binary forms, with or without 
//modification, are permitted provided that the following conditions are 
//met:

//    * Redistributions of source code must retain the above copyright 
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright 
//      notice, this list of conditions and the following disclaimer in 
//      the documentation and/or other materials provided with the distribution

//GraphSeg.h
#ifndef GRAPHSEG_H
#define GRAPHSEG_H

#include "binaryHeap.h"
#include <fstream>
#include <iostream>
using namespace std;

//const bool DEBUG_STATE=0;

#define GRAPHSEG_TEMPLATE template<class ImgType, class NodeType>
#define GRAPHSEG_WITH_TEMPLATE GraphSeg<ImgType, NodeType>
#define EDGE_TEMPLATE template<class ImgType, class NodeType>
#define EDGE_WITH_TEMPLATE Edge<ImgType, NodeType>
#define RF_TEMPLATE template<class NodeType>
#define RF_WITH_TEMPLATE RegionForest<NodeType>

//******************************************The Edge class***************************************
EDGE_TEMPLATE
class Edge
{
public:
    Edge()
    {
    }
public:
    //members:
    //start poNodeType:
    NodeType a;
    NodeType b;
    ImgType w;
};
//******************************************End of the Edge class********************************

//******************************************The RegionForest class*******************************
RF_TEMPLATE
struct RF_elt
{
    NodeType rank;
    NodeType p;
    NodeType size;
};

RF_TEMPLATE
class RegionForest
{
public:
    //methods:
    RegionForest(NodeType numOelements);
    ~RegionForest();
    NodeType find(NodeType x);
    void join(NodeType x, NodeType y);
    NodeType size(NodeType x);
    NodeType num_sets();
private:
    RF_elt<NodeType> *m_elts;
    NodeType m_numOregions;
};
RF_TEMPLATE
NodeType RF_WITH_TEMPLATE::size(NodeType x)
{
    return m_elts[x].size;
}

RF_TEMPLATE
NodeType RF_WITH_TEMPLATE::num_sets()
{
    return m_numOregions;
}

RF_TEMPLATE
RF_WITH_TEMPLATE::RegionForest(NodeType numOelements)
{
    m_elts = new RF_elt<NodeType>[numOelements];
    m_numOregions = numOelements;
    NodeType i;
    for(i=0; i<numOelements; i++)
    {
        m_elts[i].rank = 0;
        m_elts[i].size = 1;
        m_elts[i].p = i;
    }
}

RF_TEMPLATE
RF_WITH_TEMPLATE::~RegionForest()
{
    delete [] m_elts;
}

RF_TEMPLATE
NodeType RF_WITH_TEMPLATE::find(NodeType x)
{
    NodeType RegionId = x;
    NodeType next;
    while( RegionId!= m_elts[RegionId].p )
    {
        RegionId = m_elts[RegionId].p;
    }
    
    while( m_elts[x].p!=RegionId)
    {
        next = m_elts[x].p;
        m_elts[x].p=RegionId;
        x=next;
    }
    
    return RegionId;
}

RF_TEMPLATE
void RF_WITH_TEMPLATE::join(NodeType x, NodeType y)
{
    m_elts[y].p = x;
    m_elts[x].size += m_elts[y].size;

    m_numOregions--;
}
//******************************************End of the RegionForest class************************

//******************************************The GraphSeg class***********************************
GRAPHSEG_TEMPLATE
class GraphSeg
{
public:
    //methods:
    GraphSeg(ImgType* in_Img, NodeType* labeledImg, NodeType nImgHeight, NodeType nImgWidth,NodeType nRadius, double threshold, NodeType min_size);
    GraphSeg(ImgType* in_Img, NodeType* labeledImg, NodeType nImgHeight, NodeType nImgWidth,NodeType numOneighbors, double threshold, NodeType min_size, double* KNNG, double* KNNG_DIST);
    ~GraphSeg();
    void construct_graph();
    void construct_graph(NodeType numOnns);
    void segment_graph();
    void label();
    void HeapSort_edges();
	ImgType pixelDistance(NodeType p, NodeType q);
	bool InRange(NodeType x, NodeType y);
	double cmp_Threshold(NodeType N, double c);
public:
    //members:
    //the input image
    ImgType* m_Image_in;
    //the output image
    NodeType* m_Image_seg;
    //the edges:
    EDGE_WITH_TEMPLATE* m_edges;
    //the Regions:
    RF_WITH_TEMPLATE* m_Regions;
    //the thresholds for regions:
    double* m_RegionThresh;
    //the number of vertex and edges:
    NodeType m_numOvertex;
    NodeType m_numOedges;
    NodeType m_numOnns;
    //the minimum size of a component:
    NodeType m_minSize;
    //the conditional parameters:
    NodeType m_neighbor_radius;
    NodeType m_ImgHeight;
    NodeType m_ImgWidth;
    double m_threshold;
    //the k nearest neighbourhood
    double *m_knng;
    double *m_knng_dist;
};

GRAPHSEG_TEMPLATE
GRAPHSEG_WITH_TEMPLATE::GraphSeg(ImgType* in_Img, NodeType* labeledImg, NodeType nImgHeight, NodeType nImgWidth, NodeType nRadius, double threshold, NodeType min_size)
{
    int i;
    m_Image_in = in_Img;
    m_Image_seg = labeledImg;


    m_ImgHeight = nImgHeight;
    m_ImgWidth = nImgWidth;
    m_numOvertex = m_ImgHeight*m_ImgWidth;
    m_numOedges = 0;

    m_neighbor_radius = nRadius;
    m_threshold = threshold;
    m_minSize = min_size;

    //build graph
    construct_graph();
    //segmentation
    segment_graph();
    //labelling
    label();
    
    
}

GRAPHSEG_TEMPLATE
GRAPHSEG_WITH_TEMPLATE::GraphSeg(ImgType* in_Img, NodeType* labeledImg, NodeType nImgHeight, NodeType nImgWidth, NodeType numOneighbors, double threshold, NodeType min_size, double* knng, double* knng_dist)
{
    int i;
    m_Image_in = in_Img;
    m_Image_seg = labeledImg;
    m_knng=knng;
    m_knng_dist = knng_dist;


    m_ImgHeight = nImgHeight;
    m_ImgWidth = nImgWidth;
    m_numOvertex = m_ImgHeight*m_ImgWidth;
    m_numOnns = numOneighbors;
    m_numOedges = 0;

    m_neighbor_radius = 0;
    m_threshold = threshold;
    m_minSize = min_size;

    //build graph
    construct_graph(m_numOnns);
    //segmentation
    segment_graph();
    //labelling
    label();  
    
}

//free memories
GRAPHSEG_TEMPLATE
GRAPHSEG_WITH_TEMPLATE::~GraphSeg()
{
    delete [] m_edges;
    delete m_Regions;
    delete [] m_RegionThresh;
}

//Building graph:
GRAPHSEG_TEMPLATE
void GRAPHSEG_WITH_TEMPLATE::construct_graph()
{
    NodeType numOedges_max = ( (2*m_neighbor_radius+1)*(2*m_neighbor_radius+1)-1 )*m_numOvertex;
    m_edges = new EDGE_WITH_TEMPLATE[numOedges_max];
    NodeType x, y, nx_idx, ny_idx, nx, ny, p, q;
    //for each pixel p in the image
    for(x=0; x<m_ImgWidth; x++)
    {
        for(y=0; y<m_ImgHeight; y++)
        {
            p=y+m_ImgHeight*x;
            //for each neighbouring q of p
            for(nx_idx=-m_neighbor_radius; nx_idx<=m_neighbor_radius; nx_idx++)
            {
                for(ny_idx=-m_neighbor_radius; ny_idx<=m_neighbor_radius; ny_idx++)
                {
                    nx=x+nx_idx;
                    ny=y+ny_idx;
                    q=ny+m_ImgHeight*nx;
                    if( InRange(nx, ny) )
                    {
                        if(q==p) continue;
                        m_edges[m_numOedges].a = p;
                        m_edges[m_numOedges].b = q;
                        m_edges[m_numOedges].w = pixelDistance(p, q);
                        m_numOedges++;
                    }
                }
            }
        }
    }
    HeapSort_edges();
}

GRAPHSEG_TEMPLATE
void GRAPHSEG_WITH_TEMPLATE::construct_graph(NodeType numOneighbors)
{
    int i;
    NodeType x, y, p, q, q_idx;
    double w_pq;
    m_edges = new EDGE_WITH_TEMPLATE[numOneighbors*m_numOvertex];
    for(x=0; x<m_numOnns; x++)
    {
        for(y=0; y<m_numOvertex; y++)
        {
            p=y;
            q_idx = y+x*m_numOvertex;
            q = (NodeType)m_knng[q_idx];
            if(p!=q)
            {
                w_pq = m_knng_dist[q_idx];
                m_edges[m_numOedges].a = p;
                m_edges[m_numOedges].b = q;
                m_edges[m_numOedges].w = w_pq;
                m_numOedges++;
            }
        }
    }
    HeapSort_edges();
}

GRAPHSEG_TEMPLATE
void GRAPHSEG_WITH_TEMPLATE::segment_graph()
{
    //initiate a disjoNodeType-set forest:
    m_Regions = new RF_WITH_TEMPLATE(m_numOvertex);
    //initiate the thresholds:
    m_RegionThresh = new double[m_numOvertex];
    NodeType i, a, b;
    //initiate threshold
    for(i=0; i<m_numOvertex; i++)
    {
        m_RegionThresh[i] = cmp_Threshold(1, m_threshold);
    }
    //for each edge, in non-decreasing weight order...
    for(i=0; i<m_numOedges; i++)
    {
        EDGE_WITH_TEMPLATE* pEdge = &m_edges[i];

        //components connected by this edge:
        a = m_Regions->find(pEdge->a);
        b = m_Regions->find(pEdge->b);
        
        if(a!=b)
        {
            if( (pEdge->w <= m_RegionThresh[a]) && (pEdge->w <= m_RegionThresh[b]))
            {
                m_Regions->join(a, b);
                a = m_Regions->find(a);
                m_RegionThresh[a] = pEdge->w + cmp_Threshold(m_Regions->size(a), m_threshold);
            }
        }
    }
}

GRAPHSEG_TEMPLATE
void GRAPHSEG_WITH_TEMPLATE::label()
{
    //post process small components
    NodeType i, a, b;
    NodeType x, y, p;
    for(i=0; i<m_numOedges; i++)
    {
        a = m_Regions->find(m_edges[i].a);
        b = m_Regions->find(m_edges[i].b);
        if( a!=b )
        {
			if( (m_Regions->size(a)<m_minSize) && (m_Regions->size(b)<m_minSize) )
			{
				m_Regions->join(a, b);
			}
        }
    }
    //labelling:
    for(x=0; x<m_ImgWidth; x++)
    {
        for(y=0; y<m_ImgHeight; y++)
        {
            p=y+m_ImgHeight*x;
            m_Image_seg[p] = m_Regions->find(p);
        }
    }

}

GRAPHSEG_TEMPLATE
void GRAPHSEG_WITH_TEMPLATE::HeapSort_edges()
{
    CBinaryHeap<EDGE_WITH_TEMPLATE, ImgType> priorityQueue(m_numOedges);
    HeapNode<EDGE_WITH_TEMPLATE, ImgType> pair;
    NodeType i;
    for(i=0; i<m_numOedges; i++)
    {
        priorityQueue.Insert(m_edges[i], m_edges[i].w);
    }

    for(i=0; i<m_numOedges; i++)
    {
        priorityQueue.Extract(&pair);
        m_edges[i] = pair.content;
    }
}


GRAPHSEG_TEMPLATE
ImgType GRAPHSEG_WITH_TEMPLATE::pixelDistance(NodeType p, NodeType q)
{
	ImgType p_intensity, q_intensity, intensity_Distance;
	p_intensity = m_Image_in[p];
	q_intensity = m_Image_in[q];
    intensity_Distance = p_intensity-q_intensity;
    if(intensity_Distance<0)
    {
        intensity_Distance = -intensity_Distance;
    }
	return intensity_Distance;
}

GRAPHSEG_TEMPLATE
bool GRAPHSEG_WITH_TEMPLATE::InRange(NodeType x, NodeType y)
{
	bool flag_x, flag_y;
	flag_x = (x>=0 && x<m_ImgWidth);
	flag_y = (y>=0 && y<m_ImgHeight);
	return (flag_x && flag_y);
}

GRAPHSEG_TEMPLATE
double GRAPHSEG_WITH_TEMPLATE::cmp_Threshold(NodeType N, double c)
{
	return c/N;
}
//**********************************************End of the GrapSeg class**************************
#endif 