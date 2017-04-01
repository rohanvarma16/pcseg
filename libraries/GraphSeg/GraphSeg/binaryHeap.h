//Copyright (c) 2009, Yan Zhang
//All rights reserved.

//Redistribution and use in source and binary forms, with or without 
//modification, are permitted provided that the following conditions are 
//met:
//    * Redistributions of source code must retain the above copyright 
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright 
//      notice, this list of conditions and the following disclaimer in 
//      the documentation and/or other materials provided with the distribution
//Note:
//This script was originally composed by Yan Zhang and modified by Su Dongcai on 2009/11/15

//BinaryHeap.h
#ifndef __ZY_CBinaryHeap_H__
#define __ZY_CBinaryHeap_H__
#define HEAPTEMPLATE template<class ContentType, class WeightType>
#include <iostream>
#include <exception>
using namespace std;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//     CBinaryHeap
//   BinaryHeap implementation 
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
HEAPTEMPLATE
struct HeapNode
{
    ContentType content; 
    WeightType weight; 
};

HEAPTEMPLATE
class CBinaryHeap
{    
private:
    HeapNode<ContentType, WeightType>* m_heap; 
    int m_heap_idx; // indicates the heap position for next insertion.
    int m_heap_size; 
public:
    CBinaryHeap(int number_of_points);
    CBinaryHeap(ContentType* contenArray, WeightType* weightArray, int number_of_points);
    ~CBinaryHeap();
public:
    void Insert(ContentType content, WeightType weight);

    bool Extract(HeapNode<ContentType, WeightType>* heap_node);

    bool IsEmpty();

    void ReplaceData(ContentType content, WeightType weight);

    void Reset();
}; 

HEAPTEMPLATE
CBinaryHeap<ContentType, WeightType>::CBinaryHeap(int number_of_points)
{
    m_heap = new HeapNode<ContentType, WeightType>[number_of_points];             
    m_heap_idx = 0; 
    m_heap_size = number_of_points; 
}

HEAPTEMPLATE
CBinaryHeap<ContentType, WeightType>::CBinaryHeap(ContentType* contenArray, WeightType* weightArray, int number_of_points)
{
    m_heap = new HeapNode<ContentType, WeightType>[number_of_points];
    m_heap_idx = 0;
    m_heap_size = number_of_points;

    int i;
    
    for(i=0; i<number_of_points; i++)
    {
        Insert(contenArray[i], weightArray[i]);
    }
}

HEAPTEMPLATE
CBinaryHeap<ContentType, WeightType>::~CBinaryHeap()
{
    if(m_heap) delete []m_heap; 
}

HEAPTEMPLATE
void CBinaryHeap<ContentType, WeightType>::Insert(ContentType content, WeightType weight)
{
    m_heap[m_heap_idx].content = content; 
    m_heap[m_heap_idx].weight = weight; 

    int idx = m_heap_idx; 
    m_heap_idx++; 
    if( m_heap_idx > m_heap_size ) throw std::exception("CFastMarching: Heap overflow!"); 

    while(idx != 0)
    {
        int parent_idx = int((idx-1)/2); 
        if( m_heap[idx].weight < m_heap[parent_idx].weight )
        {
            HeapNode<ContentType, WeightType> tmp = m_heap[idx]; 
            m_heap[idx] = m_heap[parent_idx]; 
            m_heap[parent_idx] = tmp; 

            idx = parent_idx; 
        }
        else break; 
    }
}

HEAPTEMPLATE
bool CBinaryHeap<ContentType, WeightType>::Extract(HeapNode<ContentType, WeightType>* heap_node)
{
    if( this->IsEmpty() ) return false; 

    heap_node->content = m_heap[0].content; 
    heap_node->weight = m_heap[0].weight; 

    m_heap_idx --; 
    int idx = m_heap_idx; 
    if( idx > 0) // if still exists node in the heap
    {
        m_heap[0] = m_heap[idx]; // replace the first node by the last node
        idx = 0; 
        while(1)
        {
            int child_idx0 = 2*idx + 1; 
            int child_idx1 = child_idx0 +1; 
            if( child_idx0 >= m_heap_idx )  break; // no children for idx, i.e. idx has reached a leave, break. 
            else
            {
                int child_idx = child_idx0; 
                if( child_idx1 < m_heap_idx ) // if idx has two children.
                    if( m_heap[child_idx].weight > m_heap[child_idx1].weight ) // find the child with larger value. 
                        child_idx = child_idx1; 
                    
                // if idx has smaller value than the child found, swap.
                if( m_heap[idx].weight > m_heap[child_idx].weight ) 
                {
                    HeapNode<ContentType, WeightType> tmp = m_heap[idx]; 
                    m_heap[idx] = m_heap[child_idx]; 
                    m_heap[child_idx] = tmp; 
                    idx = child_idx; 
                }
                else break; // otherwise idx is in the right place, break; 
            }
        } 
    }

    return true; 
}

HEAPTEMPLATE
bool CBinaryHeap<ContentType, WeightType>::IsEmpty()
{
    return m_heap_idx == 0; 
}

HEAPTEMPLATE
void CBinaryHeap<ContentType, WeightType>::ReplaceData(ContentType content, WeightType weight)
{
    Insert(content, weight); 
}

HEAPTEMPLATE
void CBinaryHeap<ContentType, WeightType>::Reset()
{
    m_heap_idx = 0; 
}
#endif // __ZY_CBinaryHeap_H__ 
