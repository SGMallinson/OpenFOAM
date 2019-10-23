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

\*---------------------------------------------------------------------------*/

#include "momentumFvOptionSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "DimensionedField.H"
#include "geometricOneField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(momentumFvOptionSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    momentumFvOptionSource,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

momentumFvOptionSource::
momentumFvOptionSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    nLoc_(readScalar(coeffs_.lookup("nLoc"))),
    freq_(coeffs_.lookup("frequency")),
    Ainj_(coeffs_.lookup("area")),
    velocity_(coeffs_.lookup("velocity")),
    nuAir_(coeffs_.lookup("nuAir")),
    dp_(coeffs_.lookup("particleDiameter")) 
{
    coeffs_.lookup("frequency") >> freq_;
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

momentumFvOptionSource::~momentumFvOptionSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void momentumFvOptionSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    Info<< "momentumSource::addSup(eqn, fieldI) for source " << name_ << endl;

    this->calculate(geometricOneField(), eqn, fieldI);
}

void momentumFvOptionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{

    Info<< "momentumSource::addSup(rho, eqn, fieldI) for source " << name_ << endl;

    tmp<volScalarField> trhoRef
    (
        new volScalarField
        (
            IOobject
            (
                "rhoRef",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("rho", dimDensity, scalar(1))
        )
    );

    this->calculate(trhoRef(), eqn, fieldI);

}

bool momentumFvOptionSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        //nLoc_ = readScalar(coeffs_.lookup("nLoc"));
        //freq_ = coeffs_.lookup("frequency");
        //Ainj_ = coeffs_.lookup("area");
        //velocity_ = coeffs_.lookup("velocity");
        //nuAir_ = coeffs_.lookup("nuAir");
        dp_ = coeffs_.lookup("particleDiameter"); 

        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv
// ************************************************************************* //

