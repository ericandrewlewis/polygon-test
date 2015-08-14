// ===================================================================
//simple_Polygon.cpp - Check if a polygon is simple
// Written by Dan Sunday (2001) (http://geomalgorithms.com/a09-_intersect-3.html)
// Modified by Glenn Burkhardt (2014) to integrate it with the AVL balanced tree

// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

// http://geomalgorithms.com/a09-_intersect-3.html#Simple-Polygons

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Avl.h"
#include "simple_polygon.h"


// Assume that classes are already given for the objects:
//    Point with 2D coordinates {float x, y;}
//    Polygon with n vertices {int n; Point *V;} with V[n]=V[0]
//    Tnode is a node element structure for a BBT
//    BBT is a class for a Balanced Binary Tree
//        such as an AVL, a 2-3, or a red-black tree
//===================================================================

/*
 * An enumerator type to define whether an "event" (line segment endpoint)
 * is a "left" endpoint (lower x coordinate) or a "right" endpoint
 * (higher x coordinate).
 */
enum SEG_SIDE { LEFT, RIGHT };

/*
 * Determine the xy lexicographical order of two points.
 *
 * @return 1 if p1 > p2, -1 if p1 < p2, and 0 if equal.
 */
int xyorder( Point* p1, Point* p2 )
{
    if (p1->x > p2->x) return 1;
    if (p1->x < p2->x) return (-1);
    if (p1->y > p2->y) return 1;
    if (p1->y < p2->y) return (-1);
    return 0;
}

/*
 * Check if a point is left, on, or right of a line.
 *(see the January 2001 Algorithm on Area of Triangles)
 *
 * @param Point P0 The point to test.
 * @param Point P1 One endpoint of the line.
 * @param Point P2 The other endpoint of the line.
 * @return >0 if the point is left of the line, 0 if it is on, <0 if right.
 */
inline double isLeft( Point P0, Point P1, Point P2 )
{
    return (P1.x - P0.x)*(P2.y - P0.y) - (P2.x - P0.x)*(P1.y - P0.y);
}

class SLseg;

// EventQueue Class

/*
 * An Event.
 */
typedef struct _event Event;
struct _event {
    int           edge;     // polygon edge i is V[i] to V[i+1]
    enum SEG_SIDE type;     // event type: LEFT or RIGHT vertex
    Point*        eV;       // event vertex
    SLseg*        seg;      // segment in tree
    Event*        otherEnd; // segment is [this.eV, otherEnd.Ev]
};

int E_compare( const void* v1, const void* v2 ) // qsort compare two events
{
    Event**    pe1 = (Event**)v1;
    Event**    pe2 = (Event**)v2;

    int r = xyorder( (*pe1)->eV, (*pe2)->eV );
    if (r == 0) {
        if ((*pe1)->type == (*pe2)->type) return 0;
        if ((*pe1)->type == LEFT) return -1;
        else return 1;
    } else
        return r;
}

// the EventQueue is a presorted array (no insertions needed)
class EventQueue {
    int      ne;               // total number of events in array
    int      ix;               // index of next event on queue
    Event*   Edata;            // array of all events
    Event**  Eq;               // sorted list of event pointers
public:
    EventQueue(Polygon &P);    // constructor
    ~EventQueue(void)          // destructor
    {
        delete[] Eq;
        delete[] Edata;
    }

    Event*   next();                    // next event on queue
};

// EventQueue Routines
EventQueue::EventQueue( Polygon &P )
{
    ix = 0;
    ne = 2 * P.n;          // 2 vertex events for each edge
    Edata = (Event*)new Event[ne];
    Eq = (Event**)new Event*[ne];
    for (int i=0; i < ne; i++)          // init Eq array pointers
        Eq[i] = &Edata[i];

    // Initialize event queue with edge segment endpoints
    for (int i=0; i < P.n; i++) {       // init data for edge i
        Eq[2*i]->edge = i;
        Eq[2*i+1]->edge = i;
        Eq[2*i]->eV   = &(P.V[i]);
        Eq[2*i]->otherEnd = Eq[2*i+1];
        Eq[2*i+1]->otherEnd = Eq[2*i];
        Eq[2*i]->seg = Eq[2*i+1]->seg = 0;

        Point *pi1 = (i+1 < P.n) ? &(P.V[i+1]):&(P.V[0]);
        Eq[2*i+1]->eV = pi1;
        if (xyorder( &P.V[i], pi1) < 0) { // determine type
            Eq[2*i]->type   = LEFT;
            Eq[2*i+1]->type = RIGHT;
        }
        else {
            Eq[2*i]->type   = RIGHT;
            Eq[2*i+1]->type = LEFT;
        }
    }

    // Sort Eq[] by increasing x and y
    ::qsort( Eq, ne, sizeof(Event*), E_compare );
}

Event* EventQueue::next()
{
    if (ix >= ne)
        return (Event*)0;
    else
        return Eq[ix++];
}

/**
 * SweepLine segment data struct
 */
class SLseg : public Comparable<SLseg*> {
public:
    int      edge;         // polygon edge i is V[i] to V[i+1]
    Point    lP;           // leftmost vertex point
    Point    rP;           // rightmost vertex point
    Point*   lPp;          // pointer to leftmost vertex in poly V
    SLseg*   above;        // segment above this one
    SLseg*   below;        // segment below this one

    SLseg() : Comparable<SLseg*>(this) {}
    ~SLseg() {}

    // return true if 'this' is below 'a'
    bool operator< (const SLseg& a)
    {
        // First check if these two segments share a left vertex
        if (this->lPp == a.lPp) {
            // Same point - the two segments share a vertex.
            // use y coord of right end points
            if (this->rP.y < a.rP.y)
                return true;
            else
                return false;
        }

        return isLeft(this->lP, this->rP, a.lP) > 0;
    }

    bool operator== (const SLseg& a)
    {
        return this->edge == a.edge;
    }

    cmp_t Compare(SLseg* key) const
    {
        return (*key == *this) ? EQ_CMP
            : ((*key < *this) ? MIN_CMP : MAX_CMP);
    }
};

typedef AvlNode<SLseg*> Tnode;

// the Sweep Line itself
class SweepLine {
    int      nv;           // number of vertices in polygon
    Polygon* Pn;           // initial Polygon
    AvlTree<SLseg*> Tree;  // balanced binary tree
public:
    SweepLine(Polygon &P)          // constructor
    { nv = P.n; Pn = &P; }

    ~SweepLine(void)               // destructor
    {
        cleanTree(Tree.myRoot);
    }

    void cleanTree(Tnode *p)
    {
        if (!p) return;
        delete p->Data();
        cleanTree(p->Subtree(AvlNode<SLseg*>::LEFT));
        cleanTree(p->Subtree(AvlNode<SLseg*>::RIGHT));
    }

    SLseg*   add( Event* );
    SLseg*   find( Event* );
    bool     intersect( SLseg*, SLseg* );
    void     remove( SLseg* );
};

SLseg* SweepLine::add( Event* E )
{
    // fill in SLseg element data
    SLseg* s = new SLseg;
    s->edge  = E->edge;
    E->seg = s;

    // if it is being added, then it must be a LEFT edge event
    // but need to determine which endpoint is the left one
    Point* v1 = &(Pn->V[s->edge]);
    Point* eN = (s->edge+1 < Pn->n ? &(Pn->V[s->edge+1]):&(Pn->V[0]));
    Point* v2 = eN;
    if (xyorder( v1, v2) < 0) { // determine which is leftmost
        s->lPp = v1;
        s->lP = *v1;
        s->rP = *v2;
    }
    else {
        s->rP = *v1;
        s->lP = *v2;
        s->lPp = v2;
    }
    s->above = (SLseg*)0;
    s->below = (SLseg*)0;

    // add a node to the balanced binary tree
    Tnode* nd = Tree.Insert(s);
    Tnode* nx = Tree.Next(nd);
    Tnode* np = Tree.Prev(nd);

    if (nx != (Tnode*)0) {
        s->above = (SLseg*)nx->Data();
        s->above->below = s;
    }
    if (np != (Tnode*)0) {
        s->below = (SLseg*)np->Data();
        s->below->above = s;
    }
    return s;
}

void SweepLine::remove( SLseg* s )
{
    // remove the node from the balanced binary tree
    Tnode* nd = Tree.Search(s);
    if (nd == (Tnode*)0)
        return;      // not there !

    // get the above and below segments pointing to each other
    Tnode* nx = Tree.Next(nd);
    if (nx != (Tnode*)0) {
        SLseg* sx = (SLseg*)(nx->Data());
        sx->below = s->below;
    }
    Tnode* np = Tree.Prev(nd);
    if (np != (Tnode*)0) {
        SLseg* sp = (SLseg*)(np->Data());
        sp->above = s->above;
    }
    Tree.Delete(nd->Key());       // now can safely remove it
    delete s;                     // note:  s == nd->Data()
}

// test intersect of 2 segments and return: 0=none, 1=intersect
bool SweepLine::intersect( SLseg* s1, SLseg* s2)
{
    if (s1 == (SLseg*)0 || s2 == (SLseg*)0)
        return false;      // no intersect if either segment doesn't exist

    // check for consecutive edges in polygon
    int e1 = s1->edge;
    int e2 = s2->edge;
    if (((e1+1)%nv == e2) || (e1 == (e2+1)%nv))
        return false;      // no non-simple intersect since consecutive

    // test for existence of an intersect point
    double lsign, rsign;
    lsign = isLeft(s1->lP, s1->rP, s2->lP);    // s2 left point sign
    rsign = isLeft(s1->lP, s1->rP, s2->rP);    // s2 right point sign
    if (lsign * rsign > 0) // s2 endpoints have same sign relative to s1
        return false;      // => on same side => no intersect is possible
    lsign = isLeft(s2->lP, s2->rP, s1->lP);    // s1 left point sign
    rsign = isLeft(s2->lP, s2->rP, s1->rP);    // s1 right point sign
    if (lsign * rsign > 0) // s1 endpoints have same sign relative to s2
        return false;      // => on same side => no intersect is possible
    // the segments s1 and s2 straddle each other
    return true;           // => an intersect exists
}
//===================================================================


// simple_Polygon(): test if a Polygon P is simple or not
//     Input:  Pn = a polygon with n vertices V[]
//     Return: FALSE(0) = is NOT simple
//             TRUE(1)  = IS simple

bool simple_Polygon( Polygon &Pn )
{
    EventQueue  Eq(Pn);
    SweepLine   SL(Pn);
    Event*      e;                 // the current event
    SLseg*      s;                 // the current SL segment

    // This loop processes all events in the sorted queue
    // Events are only left or right vertices since
    // No new events will be added (an intersect => Done)
    while ((e = Eq.next())) {      // while there are events
        if (e->type == LEFT) {     // process a left vertex
            s = SL.add(e);         // add it to the sweep line
            if (SL.intersect( s, s->above))
                return false;      // Pn is NOT simple
            if (SL.intersect( s, s->below))
                return false;      // Pn is NOT simple
        }
        else {                     // process a right vertex
            s = e->otherEnd->seg;
            if (SL.intersect( s->above, s->below))
                return false;      // Pn is NOT simple
            SL.remove(s);          // remove it from the sweep line
        }
    }
    return true;      // Pn is simple
}
//===================================================================