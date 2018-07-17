/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::boundsBox

Description
    A bounding box defined in terms of the points at its extremities.

\*---------------------------------------------------------------------------*/

#ifndef boundsBox_H
#define boundsBox_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class boundsBox;
template<class T> class tmp;

Istream& operator>>(Istream&, boundsBox&);
Ostream& operator<<(Ostream&, const boundsBox&);


/*---------------------------------------------------------------------------*\
                          Class boundsBox Declaration
\*---------------------------------------------------------------------------*/

class boundsBox
{
    // Private data

        //- Minimum and maximum describing the bounding box
        point min_, max_;

    // Private Member Functions

        //- Calculate the bounding box from the given points.
        //  Does parallel communication (doReduce = true)
        void calculate(const UList<point>&, const bool doReduce = true);

public:

    // Static data members

        //- The great value used for greatBox and invertedBox
        static const scalar great;

        //- A very large boundsBox: min/max == -/+ VGREAT
        static const boundsBox greatBox;

        //- A very large inverted boundsBox: min/max == +/- VGREAT
        static const boundsBox invertedBox;


    // Constructors

        //- Construct null, setting points to zero
        inline boundsBox();

        //- Construct from components
        inline boundsBox(const point& min, const point& max);

        //- Construct as the bounding box of the given points
        //  Does parallel communication (doReduce = true)
        boundsBox(const UList<point>&, const bool doReduce = true);

        //- Construct as the bounding box of the given temporary pointField.
        //  Does parallel communication (doReduce = true)
        boundsBox(const tmp<pointField>&, const bool doReduce = true);

        //- Construct bounding box as subset of the pointField.
        //  The indices could be from cell/face etc.
        //  Does parallel communication (doReduce = true)
        boundsBox
        (
            const UList<point>&,
            const labelUList& indices,
            const bool doReduce = true
        );

        
        //- Construct bounding box as subset of the pointField.
        //  The indices could be from edge/triFace etc.
        //  Does parallel communication (doReduce = true)
        template<unsigned Size>
        boundsBox
        (
            const UList<point>&,
            const FixedList<label, Size>& indices,
            const bool doReduce = true
        );

        //- Construct from Istream
        inline boundsBox(Istream&);


    // Member functions

        // Access

            //- Minimum describing the bounding box
            inline const point& min() const;

            //- Maximum describing the bounding box
            inline const point& max() const;

            //- Minimum describing the bounding box, non-const access
            inline point& min();

            //- Maximum describing the bounding box, non-const access
            inline point& max();

            //- The midpoint of the bounding box
            inline point midpoint() const;

            //- The bounding box span (from minimum to maximum)
            inline vector span() const;

            //- The magnitude of the bounding box span
            inline scalar mag() const;

            //- The volume of the bound box
            inline scalar volume() const;

            //- Smallest length/height/width dimension
            inline scalar minDim() const;

            //- Largest length/height/width dimension
            inline scalar maxDim() const;

            //- Average length/height/width dimension
            inline scalar avgDim() const;

            //- Return corner points in an order corresponding to a 'hex' cell
            tmp<pointField> points() const;


        // Manipulate

            //- Inflate box by factor*mag(span) in all dimensions
            void inflate(const scalar s);

            // reset bounded box from two points 
            // (no need for them to be ordered)
            void resetBoundedBox  //***
            (
                const vector& startPoint,
                const vector& endPoint
            );  
            
            //- expands box only in direction of min and max by offset s
            void expand(const scalar s);
            
            void expandII(const scalar s);
            void contractII(const scalar s);
            
        // Query

            //- Overlaps/touches boundingBox?
            inline bool overlaps(const boundsBox&) const;
            
            //- Just overlaps boundingBox? (does not include cells with just edges) 
            inline bool justOverlaps(const boundsBox&) const; //***
            
            //- Contains point? (inside or on edge)
            inline bool contains(const point&) const;

            //- Fully contains other boundingBox?
            inline bool contains(const boundsBox&) const;

            //- Contains point? (inside only)
            inline bool containsInside(const point&) const;

            //- Contains all of the points? (inside or on edge)
            bool contains(const UList<point>&) const;

            //- Contains all of the points? (inside or on edge)
            bool contains
            (
                const UList<point>&,
                const labelUList& indices
            ) const;

            //- Contains all of the points? (inside or on edge)
            template<unsigned Size>
            bool contains
            (
                const UList<point>&,
                const FixedList<label, Size>& indices
            ) const;


            //- Contains any of the points? (inside or on edge)
            bool containsAny(const UList<point>&) const;

            //- Contains any of the points? (inside or on edge)
            bool containsAny
            (
                const UList<point>&,
                const labelUList& indices
            ) const;

            //- Contains any of the points? (inside or on edge)
            template<unsigned Size>
            bool containsAny
            (
                const UList<point>&,
                const FixedList<label, Size>& indices
            ) const;


    // Friend Operators

        inline friend bool operator==(const boundsBox&, const boundsBox&);
        inline friend bool operator!=(const boundsBox&, const boundsBox&);


    // IOstream operator

        friend Istream& operator>>(Istream&, boundsBox&);
        friend Ostream& operator<<(Ostream&, const boundsBox&);
};


//- Data associated with boundsBox type are contiguous
template<>
inline bool contiguous<boundsBox>() {return contiguous<point>();}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "boundsBoxI.H"

#ifdef NoRepository
#   include "boundsBoxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
