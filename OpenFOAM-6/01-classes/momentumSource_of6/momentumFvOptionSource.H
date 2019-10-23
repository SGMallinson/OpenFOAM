/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

Description

    momentumSource
    {
        type            momentumSourceSource;

        active          true;
        selectionMode   all;

        momentumSourceCoeffs
        {
            fieldNames      (U);
            nLoc = 25.0; // number of injection locations
            freq = 1.0e4; // injection frequency, Hz
            Ainj = 1.0e-7; // injector area, m^2
            velocity = 10.0; // injection velocity, m/s
            nuAir = 1.5e-5; // air kinematic viscosity, m/s^2
            dp = 20e-6; // injected particle diameter, m
        }
    }

SourceFiles
    momentumFvOptionSource.C

\*---------------------------------------------------------------------------*/

#ifndef momentumFvOptionSource_H
#define momentumFvOptionSource_H

#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"
#include "fvMatrix.H"
#include "fvMatrices.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                         momentumFvOptionSource
\*---------------------------------------------------------------------------*/

class momentumFvOptionSource
:
    public cellSetOption
{
protected:

    // Protected data

        // Number of injection locations
        scalar nLoc_;
 
        // Injection frequency, Hz
        dimensionedScalar freq_;

        // Injector area, m^2
        dimensionedScalar Ainj_;
 
        // Injection velocity, m/s
        dimensionedVector velocity_;

        // Air kinematic viscosity, m/s^2
        dimensionedScalar nuAir_;
 
        // Injected particle diameter, m
        dimensionedScalar dp_; 


    // Protected functions


public:

    //- Runtime type information
    TypeName("momentumSource");


    // Constructors

        //- Construct from patch and internal field
        momentumFvOptionSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~momentumFvOptionSource();


    // Member functions

            //- Calculate forces
            template<class RhoFieldType>
            void calculate
            (
                const RhoFieldType& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            ) const;

            //- Explicit and implicit matrix contributions
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add explicit contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv


template<class RhoFieldType>
void Foam::fv::momentumFvOptionSource::calculate
(
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
) const
{
    const scalar& Pi = Foam::constant::mathematical::pi;

    DimensionedField<vector, volMesh> UdL
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "UdL",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        velocity_,
        false
    );

    DimensionedField<vector, volMesh> Udiff
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "UDiff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<vector>("zero", dimVelocity, vector::zero),
        false
    );

    Udiff =  - UdL; // minimc vectorSemiImplicitSource

    DimensionedField<scalar, volMesh> Red
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Red",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>("zero", dimless, 0.0),
        false
    );

    Red = mag(Udiff) * dp_ / nuAir_; // local particle Reynolds number

    DimensionedField<scalar, volMesh> phi
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>("zero", dimless, 0.0),
        false
    );

    phi =  Red/(4.0 * (1.0 + Foam::sqrt(Red))) + Red/60.0;

    dimensionedScalar nde = nLoc_ * freq_ / (Ainj_ * mag(velocity_));

    DimensionedField<vector, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Su",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<vector>("zero", eqn.dimensions()/dimVolume, vector::zero),
        false
    );

    DimensionedField<scalar, volMesh> Sp
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Sp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>("zero", Su.dimensions()/eqn.psi().dimensions(), 0.0),
        false
    );

    DimensionedField<vector, volMesh> f0
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "f0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<vector>("zero", Su.dimensions(), vector::zero),
        false
    );

    f0 = 3.0*Pi*nde*nuAir_*dp_*velocity_*(1.0 + phi)*rho; 

    DimensionedField<scalar, volMesh> f1
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "f1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>("zero", Sp.dimensions(), 0.0),
        false
    );

    f1 = -3.0 * Pi * nde * nuAir_ * dp_ * (1.0 + phi)*rho;

    forAll(cells_, i)
    {
         Su[cells_[i]] = f0[cells_[i]];
         Sp[cells_[i]] = f1[cells_[i]];
    }

    eqn += Su + fvm::SuSp(Sp, eqn.psi());

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

