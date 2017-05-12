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
    utilMeasurements

Description

	Stores all the information for the utilMeasurements models (state and fluxes),
    defined within the utilMeasurementsDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef utilMeasurements_H
#define utilMeasurements_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "utilField.H"
#include "reducedUnits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class utilMeasurements Declaration
\*---------------------------------------------------------------------------*/

class utilMeasurements
{

    // Private data

        //- The entire field-properties dictionary 
        //  (it contains multiple subDictionaries)
        Time& time_;

        IOdictionary dict_;

        PtrList<entry> fieldList_;
        List<word> fieldNames_;
        List<label> fieldIds_;
        List< autoPtr<utilField> > fields_;

public:

    // Constructors

        //- Null Constructor (mdInitialise)
        utilMeasurements
        (
            Time& t,
            const polyMesh& mesh,
            const reducedUnits& rU 
        );

        ~utilMeasurements();

    // Member Functions
        
        void createFields();
        
        void calculateFields();

        void writeFields();

        inline List< autoPtr<utilField> >& fields();
        inline const List< autoPtr<utilField> >& fields() const;
        inline const List<word>& fieldNames() const;
        inline const List<label>& fieldIds() const;

    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const utilMeasurements&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
