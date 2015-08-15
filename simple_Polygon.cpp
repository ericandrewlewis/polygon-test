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

/*
 * An enumerator type to define whether an "event" (line segment endpoint)
 * is a "left" endpoint (lower x coordinate) or a "right" endpoint
 * (higher x coordinate).
 */
enum SEG_SIDE { LEFT, RIGHT };

/*
 * Determine the xy order of two points.
 *
 * Ordered first by ascending x coordinate, and then ascending y coordinate.
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

// This is oddly separated from the class definition so _event is aware of it.
class SLseg;

// Make `Event` an alias to the `_event` data structure. Weird.
typedef struct _event Event;

/**
 * An event, which is a point
 */
struct _event {
    /*
     * polygon edge i is V[i] to V[i+1]
     */
    int edge;

    /*
     * Whether the event is the LEFT or RIGHT vertex.
     */
    enum SEG_SIDE type;
    
    /*
     * The event vertex?
     */
    Point* eV;

    /*
     * The segment in tree.
     */
    SLseg* seg;
    
    /*
     * The segment is [this.eV, otherEnd.Ev].
     */
    Event* otherEnd;
};

/*
 * qsort compare two events.
 */
int E_compare( const void* v1, const void* v2 )
{
    Event** pe1 = (Event**)v1;
    Event** pe2 = (Event**)v2;

    int r = xyorder( (*pe1)->eV, (*pe2)->eV );
    if (r == 0) {
        if ((*pe1)->type == (*pe2)->type) return 0;
        if ((*pe1)->type == LEFT) return -1;
        else return 1;
    } else
        return r;
}

/*
 * The EventQueue is a presorted array (no insertions needed).
 */
class EventQueue {
    /*
     * Total number of events in array.
     */
    int numberOfEvents;

    /*
     * The index of the next event on the queue.
     */
    int nextEventIndex;

    /*
     * The array of all events.
     */
    Event* Edata;

    /*
     * Sorted list of event pointers.
     */
    Event** Eq;
public:
    // constructor
    EventQueue(Polygon &P);
    // destructor
    ~EventQueue(void)
    {
        delete[] Eq;
        delete[] Edata;
    }

    // next event on queue
    Event* next();
};

/*
 * @constructor
 * @param Polygon The polygon.
 */
EventQueue::EventQueue( Polygon &P )
{
    nextEventIndex = 0;
    // 2 vertex events for each edge.
    numberOfEvents = 2 * P.n;
    Edata = new Event[numberOfEvents];
    Eq = new Event*[numberOfEvents];

    // Initialize Eq array pointers
    for ( int i = 0; i < numberOfEvents; i++ )
        Eq[i] = &Edata[i];

    // Initialize the event queue with edge segment endpoints.
    for ( int i = 0; i < P.n; i++ ) {
        // init data for edge i
        Eq[2*i]->edge = i;
        Eq[2*i+1]->edge = i;
        Eq[2*i]->eV   = &(P.V[i]);
        Eq[2*i]->otherEnd = Eq[2*i+1];
        Eq[2*i+1]->otherEnd = Eq[2*i];
        Eq[2*i]->seg = Eq[2*i+1]->seg = 0;

        Point *pi1 = ( i + 1 < P.n ) ? &( P.V[i+1] ) : &( P.V[0] );
        Eq[2*i+1]->eV = pi1;

        // Set the event as either left or right bound.
        if ( xyorder( &P.V[i], pi1) < 0 ) {
            Eq[2*i]->type   = LEFT;
            Eq[2*i+1]->type = RIGHT;
        }
        else {
            Eq[2*i]->type   = RIGHT;
            Eq[2*i+1]->type = LEFT;
        }
    }

    // Sort the queue array by ascending x and then ascending y coordinates.
    ::qsort( Eq, numberOfEvents, sizeof(Event*), E_compare );
}

/*
 * Get a the next event on the queue.
 *
 * @return Event*
 */
Event* EventQueue::next()
{
    if (nextEventIndex >= numberOfEvents)
        return (Event*)0;
    else
        return Eq[nextEventIndex++];
}

/*
 * SweepLine segment data struct
 */
class SLseg : public Comparable<SLseg*> {
public:
    /*
     * polygon edge i is V[i] to V[i+1]
     */
    int      edge;
    
    /*
     * leftmost vertex point.
     */
    Point    lP;
    
    /*
     * rightmost vertex point.
     */
    Point    rP;
    
    /*
     * The pointer to leftmost vertex in poly V.
     */
    Point*   lPp;
    
    /*
     * The segment above this one.
     */
    SLseg*   above;
    
    /*
     * The segment below this one.
     */
    SLseg*   below;

    SLseg() : Comparable<SLseg*>(this) {}
    ~SLseg() {}

    /*
     * Define a custom comparator for whether
     * 
     * @param SLseg a
     * @return true if 'this' is below 'a'
     */
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

    /*
     * Whether the instance is equal to another SLseg
     *
     * @param SLseg segment to compare.
     */
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

// Make an AvlNode type specific for SLsegments avaible under the alias Tnode.
typedef AvlNode<SLseg*> Tnode;

/*
 * The SweepLine.
 */
class SweepLine {
    /*
     * Number of vertices in polygon.
     */
    int      nv;
    /*
     * Initial Polygon.
     */
    Polygon* Pn;
    /*
     * Balanced binary tree.
     */
    AvlTree<SLseg*> Tree;
public:
    // constructor
    SweepLine(Polygon &P)
    { nv = P.n; Pn = &P; }

    // destructor
    ~SweepLine(void)
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

/*
 *
 *
 * @param Event* E
 */
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

/*
 *
 *
 * @param SLseg* s
 */
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

/*
 * Check whether 2 segments intersect.
 *
 * @param SLseg* s1 First segment.
 * @param SLseg* s2 Other segment.
 * @return True if segments intersect, false if not.
 */
bool SweepLine::intersect( SLseg* s1, SLseg* s2)
{
    // Bail early if either segment doesn't exist.
    if (s1 == (SLseg*)0 || s2 == (SLseg*)0)
        return false;

    // check for consecutive edges in polygon.
    int e1 = s1->edge;
    int e2 = s2->edge;

    // no non-simple intersect since consecutive
    if ( ( (e1+1)%nv == e2 ) || ( e1 == (e2+1)%nv ) )
        return false;

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

/*
 * Check whether a polygon is simple (none of its lines intersect) or not.
 *
 * @param Polygon
 * @return True if the polygon is simple, false if not.
 */
bool simple_Polygon( Polygon &Pn )
{
    EventQueue  Eq(Pn);
    SweepLine   SL(Pn);
    // the current event
    Event*      e;
    // the current SL segment
    SLseg*      s;

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