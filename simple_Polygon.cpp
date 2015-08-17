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

// Make `Event` an alias to the `_event` data structure. Weird.
typedef struct _event Event;

class SweepLineSegment;

/*
 * Determine the xy order of two points.
 *
 * Ordered first by ascending x coordinate, and then ascending y coordinate.
 *
 * @param Point* p1
 * @param Point* p2
 * @return 1 if p1 > p2, -1 if p1 < p2, and 0 if equal.
 */
int xyorder( Point p1, Point p2 )
{
    if (p1.x > p2.x) return 1;
    if (p1.x < p2.x) return (-1);
    if (p1.y > p2.y) return 1;
    if (p1.y < p2.y) return (-1);
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

/**
 * An event, which relates to a line segment endpoint.
 */
struct _event {
    /*
     * The numerical representation of what "side" of the shape the event is on.
     */
    int edge;
    
    /*
     * Whether the event is the LEFT or RIGHT vertex.
     */
    enum SEG_SIDE type;
    
    /*
     * The x,y coordinate of the line segment's endpoint.
     */
    Point* point;
    
    /*
     * The segment in tree.
     */
    SweepLineSegment* seg;
    
    /*
     * The segment is [this.point, otherEnd.Ev].
     */
    Event* otherEnd;
};

/*
 * qsort compare two events.
 *
 * @param Event* v1
 * @param Event* v2
 */
int E_compare( const void* v1, const void* v2 )
{
    Event** pe1 = (Event**)v1;
    Event** pe2 = (Event**)v2;
    
    int r = xyorder( *(*pe1)->point, *(*pe2)->point );
    if (r == 0) {
        if ((*pe1)->type == (*pe2)->type) return 0;
        if ((*pe1)->type == LEFT) return -1;
        else return 1;
    } else
        return r;
}

/*
 * The event queue.
 *
 * Given a polygon, for every line segment endpoint that makes up the shape
 * an event is created and sorted for a sweepline to pass through.
 */
class EventQueue {
    /*
     * Total number of events in the queue.
     */
    int numberOfEvents;
    
    /*
     * The index of the next event on the queue.
     */
    int nextEventIndex;
    
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
    Eq = new Event*[numberOfEvents];
    
    // Populate each event's data from the Polygon's points.
    for ( int i = 0; i < P.n; i++ ) {
        Eq[2*i] = new Event;
        Eq[2*i+1] = new Event;

        Eq[2*i]->edge = i;
        Eq[2*i]->point = &(P.V[i]);
        Eq[2*i]->otherEnd = Eq[2*i+1];
        Eq[2*i]->seg = 0;
        
        Eq[2*i+1]->edge = i;
        Eq[2*i+1]->otherEnd = Eq[2*i];
        Eq[2*i+1]->seg = 0;
        
        Point *pi1 = ( i + 1 < P.n ) ? &( P.V[i+1] ) : &( P.V[0] );
        Eq[2*i+1]->point = pi1;
        
        // Set the event as either left or right bound.
        if ( xyorder( P.V[i], *pi1 ) < 0 ) {
            Eq[2*i]->type = LEFT;
            Eq[2*i+1]->type = RIGHT;
        }
        else {
            Eq[2*i]->type = RIGHT;
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

class SweepLineSegment : public Comparable<SweepLineSegment*> {
public:
    /*
     * polygon edge i is V[i] to V[i+1]
     */
    int edge;
    
    /*
     * leftmost vertex point.
     */
    Point* leftPoint;
    
    /*
     * rightmost vertex point.
     */
    Point* rightPoint;
    
    /*
     * The segment above this one.
     */
    SweepLineSegment*   above;
    
    /*
     * The segment below this one.
     */
    SweepLineSegment*   below;
    
    SweepLineSegment() : Comparable<SweepLineSegment*>(this) {}
    ~SweepLineSegment() {}
    
    /*
     * Compare two sweep line segments by which is more "left".
     *
     * @param SweepLineSegment a The other sweep line segment to compare to.
     * @return True if the class instance is below
     */
    bool operator< (const SweepLineSegment& otherLineSegment)
    {
        // If these two segments share a left vertex, compare the right points.
        if ( this->leftPoint == otherLineSegment.leftPoint ) {
            // Same point - the two segments share a vertex.
            // use y coord of right end points
            if ( this->rightPoint->y < otherLineSegment.rightPoint->y )
                return true;
            else
                return false;
        }
        
        return isLeft(*this->leftPoint, *this->rightPoint, *otherLineSegment.leftPoint) > 0;
    }
    
    /*
     * Whether the instance is equal to another SweepLineSegment
     *
     * @param SweepLineSegment segment to compare.
     */
    bool operator== (const SweepLineSegment& a)
    {
        return this->edge == a.edge;
    }
    
    cmp_t Compare(SweepLineSegment* key) const
    {
        return (*key == *this) ? EQ_CMP
        : ((*key < *this) ? MIN_CMP : MAX_CMP);
    }
};

// Make an AvlNode type specific for SweepLineSegmentments avaible under the alias Tnode.
typedef AvlNode<SweepLineSegment*> Tnode;

/*
 * The sweep line keeps track of line segments that are intersecting
 * at the currently processing line segment endpoint.
 */
class SweepLine {
    /*
     * Number of vertices in polygon.
     */
    int nv;
    
    /*
     * Initial Polygon.
     */
    Polygon* polygon;
    
    /*
     * Balanced binary tree.
     */
    AvlTree<SweepLineSegment*> Tree;
public:
    // constructor
    SweepLine(Polygon &P)
    {
        nv = P.n; polygon = &P;
    }
    
    // destructor
    ~SweepLine(void)
    {
        cleanTree(Tree.myRoot);
    }
    
    void cleanTree(Tnode *p)
    {
        if (!p) return;
        delete p->Data();
        cleanTree(p->Subtree(AvlNode<SweepLineSegment*>::LEFT));
        cleanTree(p->Subtree(AvlNode<SweepLineSegment*>::RIGHT));
    }
    
    SweepLineSegment* add( Event* );
    SweepLineSegment* find( Event* );
    bool intersect( SweepLineSegment*, SweepLineSegment* );
    void remove( SweepLineSegment* );
};

/*
 * Add an event (line segment endpoint) to the sweep line.
 *
 * @param Event* E Event to add to the sweep line.
 * @return SweepLineSegment* s The new line segment.
 */
SweepLineSegment* SweepLine::add( Event* event )
{
    // Create a line segment from the event.
    SweepLineSegment* lineSegment = new SweepLineSegment;
    lineSegment->edge = event->edge;
    event->seg = lineSegment;
    
    // If it is being added, then it must be a LEFT edge event
    // but need to determine which endpoint is the left one.
    Point* endpoint1 = &( polygon->V[lineSegment->edge] );
    Point* eN = (lineSegment->edge+1 < polygon->n ? &(polygon->V[lineSegment->edge+1]) : &(polygon->V[0]));
    Point* endpoint2 = eN;
    
    // Determine which is endpoint is leftmost.
    if ( xyorder( *endpoint1, *endpoint2 ) < 0 ) {
        lineSegment->leftPoint = endpoint1;
        lineSegment->rightPoint = endpoint2;
    }
    else {
        lineSegment->rightPoint = endpoint1;
        lineSegment->leftPoint = endpoint2;
    }
    lineSegment->above = (SweepLineSegment*)0;
    lineSegment->below = (SweepLineSegment*)0;
    
    // Add a node to the balanced binary tree.
    Tnode* node = Tree.Insert(lineSegment);
    Tnode* nextNode = Tree.Next(node);
    Tnode* previousNode = Tree.Prev(node);
    
    if ( nextNode != (Tnode*)0 ) {
        lineSegment->above = (SweepLineSegment*)nextNode->Data();
        lineSegment->above->below = lineSegment;
    }
    if ( previousNode != (Tnode*)0 ) {
        lineSegment->below = (SweepLineSegment*)previousNode->Data();
        lineSegment->below->above = lineSegment;
    }
    return lineSegment;
}

/*
 * Remove a line segment from the sweep line.
 *
 * @param SweepLineSegment* s The line segment to be removed.
 */
void SweepLine::remove( SweepLineSegment* lineSegment )
{
    // remove the node from the balanced binary tree
    Tnode* node = Tree.Search(lineSegment);
    
    // If the node can't be found, bail.
    if ( node == (Tnode*)0 )
        return;
    
    // Get the above and below segments pointing to each other.
    Tnode* nextNode = Tree.Next(node);
    if ( nextNode != (Tnode*)0 ) {
        SweepLineSegment* sx = (SweepLineSegment*)(nextNode->Data());
        sx->below = lineSegment->below;
    }
    Tnode* previousNode = Tree.Prev(node);
    if ( previousNode != (Tnode*)0 ) {
        SweepLineSegment* sp = (SweepLineSegment*)(previousNode->Data());
        sp->above = lineSegment->above;
    }
    // Now can safely remove it.
    Tree.Delete( node->Key() );
    // note:  lineSegment == nd->Data()
    delete lineSegment;
}

/*
 * Check whether 2 segments intersect.
 *
 * @param SweepLineSegment* s1 First segment.
 * @param SweepLineSegment* s2 Other segment.
 * @return True if segments intersect, false if not.
 */
bool SweepLine::intersect( SweepLineSegment* s1, SweepLineSegment* s2 )
{
    // Bail early if either segment doesn't exist.
    if (s1 == (SweepLineSegment*)0 || s2 == (SweepLineSegment*)0)
        return false;
    
    // check for consecutive edges in polygon.
    int e1 = s1->edge;
    int e2 = s2->edge;
    
    // no non-simple intersect since consecutive
    if ( ( (e1+1)%nv == e2 ) || ( e1 == (e2+1)%nv ) )
        return false;
    
    // test for existence of an intersect point
    double lsign, rsign;
    lsign = isLeft(*s1->leftPoint, *s1->rightPoint, *s2->leftPoint);    // s2 left point sign
    rsign = isLeft(*s1->leftPoint, *s1->rightPoint, *s2->rightPoint);    // s2 right point sign
    if (lsign * rsign > 0) // s2 endpoints have same sign relative to s1
        return false;      // => on same side => no intersect is possible
    lsign = isLeft(*s2->leftPoint, *s2->rightPoint, *s1->leftPoint);    // s1 left point sign
    rsign = isLeft(*s2->leftPoint, *s2->rightPoint, *s1->rightPoint);    // s1 right point sign
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
bool simple_Polygon( Polygon &polygon )
{
    EventQueue eventQueue( polygon );
    SweepLine sweepline( polygon );
    Event* currentEvent;
    SweepLineSegment* currentSegment;
    
    // Loop through all sorted events in the queue.
    // Events are only left or right vertices since no new events will be added (an intersect => Done).
    while ( (currentEvent = eventQueue.next()) ) {
        // Process a left endpoint.
        if ( currentEvent->type == LEFT ) {
            // Add it to the sweep line.
            currentSegment = sweepline.add(currentEvent);
            // If the current segment intersects with its sweep line neighbor, the polygon is not simple.
            if ( sweepline.intersect( currentSegment, currentSegment->above ) )
                return false;
            if ( sweepline.intersect( currentSegment, currentSegment->below ) )
                return false;
        }
        // process a right endpoint.
        else {
            currentSegment = currentEvent->otherEnd->seg;
            // If the new neighbor segments intersect, the polygon is not simple.
            if ( sweepline.intersect( currentSegment->above, currentSegment->below ) )
                return false;
            // Remove it from the sweep line.
            sweepline.remove( currentSegment );
        }
    }
    // If no intersections were found after the sweepline has passed all events, the polygon is simple.
    return true;
}