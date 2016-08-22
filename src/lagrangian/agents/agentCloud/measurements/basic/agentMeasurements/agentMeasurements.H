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

Class
    agentMeasurements

Description

\*---------------------------------------------------------------------------*/

#ifndef agentMeasurements_H
#define agentMeasurements_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "agentMeasurement.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentMeasurements Declaration
\*---------------------------------------------------------------------------*/

class agentMeasurements
{

    // Private data

        Time& time_;

        IOdictionary measurementsDict_;

        PtrList<entry> fieldList_;
        List<word> fieldNames_;
        List<label> fieldIds_;
        List< autoPtr<agentMeasurement> > fields_;

        labelList measurementsDuringForceComp_;
        labelList measurementsDuringForceCompSite_;

public:

    // Constructors

        //- Null Constructor (mdInitialise)
        agentMeasurements
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFoam
        agentMeasurements
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud
        );

        ~agentMeasurements();

    // Member Functions


        void createFields();

        void calculateFields();

        void writeFields();

        void measurementsDuringForceComputation
        (
            agent* molI,
            agent* molJ
        );

        void measurementsDuringForceComputationSite
        (
            agent* molI,
            agent* molJ,
            label sI,
            label sJ
        );


        inline List< autoPtr<agentMeasurement> >& fields();
        inline const List< autoPtr<agentMeasurement> >& fields() const;
        inline const List<word>& fieldNames() const;
        inline const List<label>& fieldIds() const;

    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const agentMeasurements&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "agentMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
