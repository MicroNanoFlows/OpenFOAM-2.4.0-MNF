/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "agent.H"
#include "IOstreams.H"
#include "agentCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::agent::agent
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    v_(vector::zero),
    d_(vector::zero),
    f_(vector::zero),    
    specialPosition_(vector::zero),
    mass_(0.0),
    potentialEnergy_(0.0),
    R_(GREAT),
    frac_(1.0),
    special_(0),
    id_(0),
    trackingNumber_(-1)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> v_;
            is  >> d_;
            is  >> f_;
            is  >> specialPosition_;
            is >> mass_;
            is >> potentialEnergy_;
            is >> R_;
            is >> frac_;
            special_ = readLabel(is);
            id_ = readLabel(is);
            trackingNumber_ = readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&v_),
                sizeof(v_)
              + sizeof(d_)
              + sizeof(f_)
              + sizeof(specialPosition_)
              + sizeof(mass_)
              + sizeof(potentialEnergy_)
              + sizeof(R_)
              + sizeof(frac_)              
              + sizeof(special_)
              + sizeof(id_)
              + sizeof(trackingNumber_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "Foam::agent::agent"
        "(const Cloud<agent>& cloud, Foam::Istream&), bool"
    );
}


void Foam::agent::readFields(Cloud<agent>& mC)
{
    if (!mC.size())
    {
        return;
    }

    particle::readFields(mC);
    
    Info << "Reading fields" << endl;



    // MB: I have removed reading/writing these fields because no one ever uses them
    //    and they by doing this we are saving a lot on hard disk space
    
    IOField<vector> v(mC.fieldIOobject("v", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, v);

//     IOField<vector> f(mC.fieldIOobject("f", IOobject::MUST_READ));
//     mC.checkFieldIOobject(mC, f);

    IOField<vector> specialPosition
    (
        mC.fieldIOobject("specialPosition", IOobject::MUST_READ)
    );
    
    mC.checkFieldIOobject(mC, specialPosition);

    IOField<scalar> mass(mC.fieldIOobject("mass", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, mass);    
    
    IOField<label> special(mC.fieldIOobject("special", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, special);
    
    IOField<label> id(mC.fieldIOobject("id", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, id);
    
    IOField<label> trackingNumber(mC.fieldIOobject("trackingNumber", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, trackingNumber);    

    label i = 0;
    forAllIter(agentCloud, mC, iter)
    {
        agent& mol = iter();

        mol.v_ = v[i];
//         mol.a_ = a[i];
        mol.specialPosition_ = specialPosition[i];
        mol.mass_ = mass[i];
        mol.special_ = special[i];
        mol.id_ = id[i];
        mol.trackingNumber_ = trackingNumber[i];        
        i++;
    }
}


void Foam::agent::writeFields(const Cloud<agent>& mC)
{
    particle::writeFields(mC);
    
    label np = mC.size();

//     IOField<tensor> Q(mC.fieldIOobject("Q", IOobject::NO_READ), np);
//     IOField<tensor> rf(mC.fieldIOobject("rf", IOobject::NO_READ), np);
    IOField<vector> v(mC.fieldIOobject("v", IOobject::NO_READ), np);
//     IOField<vector> a(mC.fieldIOobject("a", IOobject::NO_READ), np);
//     IOField<vector> pi(mC.fieldIOobject("pi", IOobject::NO_READ), np);
//     IOField<vector> tau(mC.fieldIOobject("tau", IOobject::NO_READ), np);
    IOField<vector> specialPosition
    (
        mC.fieldIOobject("specialPosition", IOobject::NO_READ),
        np
    );
    IOField<label> special(mC.fieldIOobject("special", IOobject::NO_READ), np);
    IOField<label> id(mC.fieldIOobject("id", IOobject::NO_READ), np);
    IOField<scalar> mass(mC.fieldIOobject("mass", IOobject::NO_READ), np);    
    IOField<label> trackingNumber(mC.fieldIOobject("trackingNumber", IOobject::NO_READ), np);
    
    // Post processing fields

    label i = 0;
    forAllConstIter(agentCloud, mC, iter)
    {
        const agent& mol = iter();

        v[i] = mol.v_;
        mass[i] = mol.mass_;
        specialPosition[i] = mol.specialPosition_;
        special[i] = mol.special_;
        id[i] = mol.id_;
        trackingNumber[i] = mol.trackingNumber_;
        i++;
    }

    v.write();
    mass.write();
    specialPosition.write();
    special.write();
    id.write();
    trackingNumber.write();

    Info<< "writeFields " << mC.name() << endl;
    
    if (isA<agentCloud>(mC))
    {
        const agentCloud& m = dynamic_cast<const agentCloud&>(mC);

        m.writeXYZ
        (
            m.mesh().time().timePath()/cloud::prefix/"agentCloud.xmol"
        );
    }
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const agent& mol)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << token::SPACE << static_cast<const particle&>(mol)
            << token::SPACE << mol.face()
            << token::SPACE << mol.stepFraction()
            << token::SPACE << mol.v_
            << token::SPACE << mol.d_
            << token::SPACE << mol.f_
            << token::SPACE << mol.specialPosition_
            << token::SPACE << mol.mass_
            << token::SPACE << mol.potentialEnergy_
            << token::SPACE << mol.R_
            << token::SPACE << mol.frac_
            << token::SPACE << mol.special_
            << token::SPACE << mol.id_
            << token::SPACE << mol.trackingNumber_;
    }
    else
    {
        os  << static_cast<const particle&>(mol);
        os.write
        (
            reinterpret_cast<const char*>(&mol.v_),
            sizeof(mol.v_)
          + sizeof(mol.d_)
          + sizeof(mol.f_)
          + sizeof(mol.specialPosition_)
          + sizeof(mol.mass_)
          + sizeof(mol.potentialEnergy_)
          + sizeof(mol.R_)
          + sizeof(mol.frac_)
          + sizeof(mol.special_)
		  + sizeof(mol.id_)
          + sizeof(mol.trackingNumber_)
        );
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::agent&)"
    );

    return os;
}


// ************************************************************************* //
