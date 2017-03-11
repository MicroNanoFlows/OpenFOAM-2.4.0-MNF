/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "visibilityGraphAndDjikstra.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(visibilityGraphAndDjikstra, 0);

addToRunTimeSelectionTable(scheduleModel, visibilityGraphAndDjikstra, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



void visibilityGraphAndDjikstra::setBoundBoxes()
{
    PtrList<entry> boxList(propsDict_.lookup("boxes"));

    boxes_.setSize(boxList.size());
    destinations_.setSize(boxList.size());
    nScheduledDestinations_ = destinations_.size();
    
    if (propsDict_.found("timeAtDestination"))
    {
        timeAtDestination_ = scalarList(propsDict_.lookup("timeAtDestination"));
    }
    else
    {
        timeAtDestination_.setSize(nScheduledDestinations_, 0.0);
    }    
    
    forAll(boxList, b)
    {
        const entry& boxI = boxList[b];
        const dictionary& dict = boxI.dict();

        vector startPoint = dict.lookup("startPoint");
        vector endPoint = dict.lookup("endPoint");
        boxes_[b].resetBoundedBox(startPoint, endPoint);
        destinations_[b] = boxes_[b].midpoint();        
    }
    
    if( (nScheduledDestinations_ <= 0) || 
        (nScheduledDestinations_ != timeAtDestination_.size())
    )
    {
        FatalError
            << "Something went wrong with the localGoalsSchedule " << endl
            << " lists must match in size and need to have sizes greater than zero"
            << endl << "destinations = " << destinations_
            << endl << "timeAtDestination = " << timeAtDestination_
            << nl << abort(FatalError);  
        
    }     
}

// Construct from components
visibilityGraphAndDjikstra::visibilityGraphAndDjikstra
(
    Time& time,
    agentCloud& cloud,
    const dictionary& dict
)
:
    scheduleModel(time, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    timeIndex_(0.0),
    counter_(-1),
    deltaT_(time_.deltaT().value())
{
    // ADDED - read borders list
    borderList_ = List<vectorList>(propsDict_.lookup("bordersList"));

    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();  
    

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

visibilityGraphAndDjikstra::~visibilityGraphAndDjikstra()
{
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void visibilityGraphAndDjikstra::initialConfiguration()
{
    setBoundBoxes();
    
    label nextDest = 0;
    
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {    
        if(findIndex(agentIds_, mol().id()) != -1)
        {
            mol().d()= destinations_[nextDest];
        }
    } 
    
    
    // Generate a list of all the points to construct the visibility graph
    // Points need to be offset from the corners of borders in the right direction
    // Need to get the normals to the edges of borders (like in borders code)
    label nBorders = borderList_.size();
    normalVectors_.setSize(nBorders);
    midPoints_.setSize(nBorders);
    
    label nSegments = 4*nBorders;
    
    //offsetPointsPtr = new vectorList(nSegments);
    //vectorList offsetPoints_ = *offsetPointsPtr;
    
    offsetPoints_.setSize(nSegments);
    
    {
	SquareMatrix<scalar> graph(nSegments+2);
    
	graph_ = graph;
    }
    
//     offsetPoints_.setSize(nSegments);
    
    // adjacency matrix
    //graphPtr = new scalarSquareMatrix(nSegments+2); // +2 because of agent position and goal position
    //scalarSquareMatrix graph_ = *graphPtr;
    
    scalar offset = 0.4;
    
    forAll(borderList_, i)
    {
// 	const vectorList& borderI = borderList_[i];
	
	label sz=borderList_[i].size()-1;
        normalVectors_[i].setSize(sz);
        midPoints_[i].setSize(sz);
        vector centrePoint = vector::zero;
        
        // set midpoints
        for (int p = 0; p < sz; p++)
        {
            const vector& vI=borderList_[i][p];
            centrePoint += vI;
        }
        
        centrePoint /= scalar(sz);
        
//         Info << "centrePoint = " << centrePoint << endl;
        
        // set normals
        for (int p = 0; p < sz; p++)
        {
            const vector& vI=borderList_[i][p];
            const vector& vJ=borderList_[i][p+1];
            vector vIJ=vJ-vI;
            vector midPoint = vI + vIJ*0.5;
            midPoints_[i][p]=midPoint;
            
            vector n = midPoint-centrePoint; // pointing outwards
            n /= mag(n);
            
            // deal with skewness
            scalar theta = acos(vIJ & n / mag(vIJ));

            n *= sin(theta);
            n /= mag(n);
            
            normalVectors_[i][p] = n;
	    
        }
        
        // Offsetting the points by a specified value "offset"
        for (int p = 0; p < sz; p++)
        {
	    if(p == 0)
	    {
		vector offsetDirection = normalVectors_[i][p] + normalVectors_[i][p+3];
		offsetDirection /= mag(offsetDirection);
		offsetPoints_[4*i+p] = borderList_[i][p] + offset*offsetDirection;
		
// 		borderList_[i][p] = borderList_[i][p] + (0.99*offset)*offsetDirection;
// 		// 0.99 multiplied to avoid a colinear scenario later
		
// 		borderList_[i][4] = borderList_[i][p];
	    }
	    else
	    {
		vector offsetDirection = normalVectors_[i][p] + normalVectors_[i][p-1];
		offsetDirection /= mag(offsetDirection);
		offsetPoints_[4*i+p] = borderList_[i][p] + offset*offsetDirection;
		
// 		borderList_[i][p] = borderList_[i][p] + (0.99*offset)*offsetDirection;
		
	    }
	}
	     
        
    }
    
//     Info << "offsetPoints = " << offsetPoints_ << endl;
    
    
    
    // Create the adjacency matrix to construct the graph
    // (no agent position and no goal position yet)
    
    for (int i = 0; i < graph_.n(); i++)
    {
	graph_[i][0] = 0;
	graph_[0][i] = 0;
	graph_[i][graph_.m()-1] = 0;
	graph_[graph_.n()-1][i] = 0;
	
// 	Info << ">>>>>>>>>>>> i = " << i << endl;
    }
    
    forAll(offsetPoints_, i)
    {
	forAll(offsetPoints_, j)
	{
	    bool intersectionYes1 = false;
	    
	    vector p1 = offsetPoints_[i];
	    vector q1 = offsetPoints_[j];
	    
	    if(i == j)
	    {
		graph_[i+1][j+1] = 0; // +1 as index 0 will be the agent's position
		graph_[j+1][i+1] = 0;
	    }
	    else
	    {
	    
		forAll(borderList_, k)
		{
		    label sz=borderList_[k].size()-1;
		    for(int p = 0; p < sz; p++)
		    {
			vector p2 = borderList_[k][p];
			vector q2 = borderList_[k][p+1];
			
// 			if (i==0 && j==1 && p==0)
// 			{
// 			    Info << "i = " << i << endl;
// 			    Info << "j = " << j << endl;
// 			    Info << "p1 = " << p1 << endl;
// 			    Info << "q1 = " << q1 << endl;
// 			    Info << "p2 = " << p2 << endl;
// 			    Info << "q2 = " << q2 << endl;
// 			    Info << "doIntersect = " << doIntersect(p1, q1, p2, q2) << endl;
// 			}
			
			if(doIntersect(p1, q1, p2, q2))
			{
			    graph_[i+1][j+1] = 0;
			    graph_[j+1][i+1] = 0;
			    intersectionYes1 = true;
			    break;
			}

		    }
		    
		    if(intersectionYes1)
		    {
			break;
		    }
		    
		}
		
		if(!intersectionYes1)
		{
		    graph_[i+1][j+1] = mag(q1-p1);
		    graph_[j+1][i+1] = mag(q1-p1);
		}
		
	    }
	}
    }
    
//     Info << "adjacencyMatrix = " << graph_ << endl;
    
}









// MEMBER FUNCTIONS USED TO CHECK IF TWO SEGMENTS INTERSECT

// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool visibilityGraphAndDjikstra::onSegment(vector p, vector q, vector r)
{
    if (q.x() <= max(p.x(), r.x()) && q.x() >= min(p.x(), r.x()) &&
        q.y() <= max(p.y(), r.y()) && q.y() >= min(p.y(), r.y()))
    {
       return true;
    }
    else
    {
       return false;
    }
}
 
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int visibilityGraphAndDjikstra::orientation(vector p, vector q, vector r)
{
    // See http://www.geeksforgeeks.org/orientation-3-ordered-points/
    // for details of below formula.
    scalar val = (q.y() - p.y()) * (r.x() - q.x()) -
              (q.x() - p.x()) * (r.y() - q.y());
	      
 
    if (val == 0)
    {
	return 0;  // colinear
    }
    else
    { 
        return (val > 0)? 1: 2; // clock or counterclock wise
    }
}
 
// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool visibilityGraphAndDjikstra::doIntersect(vector p1, vector q1, vector p2, vector q2)
{
    // Find the four orientations needed for general and// 	Info << "i = " << i << endl;
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
    
 
    // General case
    if (o1 != o2 && o3 != o4)
        return true;
 
    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
 
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
 
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
 
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;
 
    return false; // Doesn't fall in any of the above cases
}










// is called every time-step to set new destinations of agents
void visibilityGraphAndDjikstra::setSchedule()
{
    
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {    
        if(findIndex(agentIds_, mol().id()) != -1)
        {
            if(mol().t() > 0 )
            {
		mol().v() = vector::zero;
                mol().t() -= deltaT_;
                
                if(mol().t() < 0)
                {
		    
                    mol().t() = 0.0;
                }
            }
            else
            {
            
//                 label index = -1;

                vector finalGoal = destinations_[0];
// 		Info << "finalGoal = " << finalGoal << endl;
		
		vector agentPosition = mol().position();
// 		Info << "agentPosition = " << agentPosition << endl;
		
		// Include agent position and finalGoal in the visibility graph
// 		scalarSquareMatrix graph_ = *graphPtr;
// 		vectorList offsetPoints_ = *offsetPointsPtr;
		
// 		Info << "graph_ = " << graph_ << endl;
// 		Info << "offsetPoints_= " << offsetPoints_ << endl;
// 		Info << "graph_ = " << graph_ << endl;
// 		Info << "offsetPoints_ = " << offsetPoints_ << endl;

		graph_[0][0] = 0;
		graph_[graph_.n()-1][graph_.m()-1] = 0;
		

// 		Info << "offsetPoints_= " << offsetPoints_ << endl;
		
		
		
		// Check the visibility between the agent and his destination first
		vector p1 = agentPosition;
		vector q1 = finalGoal;
		bool intersectionYes1 = false;
		forAll(borderList_, k)
		{
		    label sz=borderList_[k].size()-1;
		    for(int p = 0; p < sz; p++)
		    {
			vector p2 = borderList_[k][p];
			vector q2 = borderList_[k][p+1];
			
			if(doIntersect(p1, q1, p2, q2))
			{
			    graph_[0][graph_.m()-1] = 0;
			    graph_[graph_.n()-1][0] = 0;
			    intersectionYes1 = true;
			}
			
			if(intersectionYes1)
			{
			    break;
			}
		    }
		    
		    if(intersectionYes1)
		    {
			break;
		    }
		}
		
		if(!intersectionYes1)
		{
		    graph_[0][graph_.m()-1] = mag(q1-p1);
		    graph_[graph_.n()-1][0] = mag(q1-p1);
		    mol().d() = finalGoal;
		}
		else
		{

		    // Check visibility with other points and run dijkstra
		    forAll(offsetPoints_, j)
		    {
			bool intersectionYes1 = false;
			bool intersectionYes2 = false;
			
			vector p1 = agentPosition;
			vector q1 = offsetPoints_[j];
			
			vector p1F = finalGoal;
			
			forAll(borderList_, k)
			{
			    label sz=borderList_[k].size()-1;
			    for(int p = 0; p < sz; p++)
			    {
				vector p2 = borderList_[k][p];
				vector q2 = borderList_[k][p+1];
				
				if(doIntersect(p1, q1, p2, q2))
				{
				    graph_[0][j+1] = 0;
				    graph_[j+1][0] = 0;
				    intersectionYes1 = true;
				}
				
				if(doIntersect(p1F, q1, p2, q2))
				{
				    graph_[graph_.n()-1][j+1] = 0;
				    graph_[j+1][graph_.m()-1] = 0;
				    intersectionYes2 = true;
				}
				
				if(intersectionYes1 && intersectionYes2)
				{
				    break;
				}

			    }
			    
			    if(intersectionYes1 && intersectionYes2)
			    {
				break;
			    }
			    
			}
			
			if (j==0)
			{
// 			    Info << "intersectionYes1 = " << intersectionYes1 << endl;
// 			    Info << "intersectionYes2 = " << intersectionYes2 << endl;
			}
			
			if(!intersectionYes1)
			{
			    graph_[0][j+1] = mag(q1-p1);
			    graph_[j+1][0] = mag(q1-p1);
			}
			
			if(!intersectionYes2)
			{
			    graph_[graph_.n()-1][j+1] = mag(q1-p1F);
			    graph_[j+1][graph_.m()-1] = mag(q1-p1F);
			}
			
		    }
		    
// 		    Info << "graph_ = " << graph_ << endl;
		    
		    scalar nextDestinationIndex = dijkstra(graph_, 0);
		    
// 		    Info << "nextDestinationIndex = " << nextDestinationIndex << endl;
// 		    Info << "nextDestination = " << offsetPoints_[nextDestinationIndex-1] << endl;
		    
		    mol().d() = offsetPoints_[nextDestinationIndex-1];

		}
		
            }
  
        }
    }

}










// DJIKSTRA ALGORITHM + SUPPORTING MEMBER FUNCTIONS

// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
int visibilityGraphAndDjikstra::minDistance(scalar dist[], bool sptSet[], int V)
{
   // Initialize min value
   scalar min = GREAT;
   int min_index = -1;
  
   for (int v = 0; v < V; v++)
   {
     if (sptSet[v] == false && dist[v] <= min)
     {
         min = dist[v];
         min_index = v;
     }
   }
  
   return min_index;
}
  
// A utility function to print the constructed distance array
// int visibilityGraphAndDjikstra::printSolution(int dist[], int n)
// {
//    printf("Vertex   Distance from Source\n");
//    for (int i = 0; i < V; i++)
//       printf("%d \t\t %d\n", i, dist[i]);
// }
  
// Funtion that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
int visibilityGraphAndDjikstra::dijkstra(scalarSquareMatrix graph, int src)
{
     int V = graph.n();  
    
     scalar dist[V];     // The output array.  dist[i] will hold the shortest
                      // distance from src to i
                      
     scalar prev[V]; // prev[i] holds the id of the graph point previous to graph
		     // point i in the shortest path
  
     bool sptSet[V]; // sptSet[i] will true if vertex i is included in shortest
                     // path tree or shortest distance from src to i is finalized
  
     // Initialize all distances as INFINITE and stpSet[] as false
     for (int i = 0; i < V; i++)
     {
        dist[i] = GREAT;
        sptSet[i] = false;
	prev[i] = -1;
     }
  
     // Distance of source vertex from itself is always 0
     dist[src] = 0;
     prev[src] = src;
  
     // Find shortest path for all vertices
     for (int count = 0; count < V-1; count++)
     {
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in first iteration.
       int u = minDistance(dist, sptSet, V);
  
       // Mark the picked vertex as processed
       sptSet[u] = true;
  
       // Update dist value of the adjacent vertices of the picked vertex.
       for (int v = 0; v < V; v++)
       {
  
         // Update dist[v] only if is not in sptSet, there is an edge from 
         // u to v, and total weight of path from src to  v through u is 
         // smaller than current value of dist[v]
         if (!sptSet[v] && graph[u][v] != 0 && dist[u] != GREAT 
                                       && dist[u]+graph[u][v] < dist[v])
	 {
            dist[v] = dist[u] + graph[u][v];
	    prev[v] = u;
	    
// 	    Info << "v = " << v << endl;
// 	    Info << "dist = " << dist[v] << endl;
// 	    Info << "prev = " << prev[v] << endl;
	 }
       }
     }
  
     // print the constructed distance array
//      printSolution(dist, V);
     
     int a = V-1;
     int b = a;
     while(a > 0)
     {
	 b = a;
	 a = prev[a];
     }
 
//      for (int i = 0; i < V; i++)
//      {
// 	 Info << "dist = " << dist[i] << endl;
// 	 Info << "prev = " << prev[i] << endl;
//      }
     
//      Info << "b = " << b << endl;
     
     return b; // returns the index of the next agent's destination
}









} // End namespace Foam

// ************************************************************************* //
