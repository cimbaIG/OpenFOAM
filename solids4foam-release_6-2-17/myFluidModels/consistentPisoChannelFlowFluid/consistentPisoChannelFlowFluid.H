/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    consistentPisoChannelFlowFluid

Description
    Consistent PISO solver for LES channel flow.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    consistentPisoChannelFlowFluid.C

\*---------------------------------------------------------------------------*/

#ifndef consistentPisoChannelFlowFluid_H
#define consistentPisoChannelFlowFluid_H

#include "fluidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "singlePhaseTransportModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class consistentPisoChannelFlowFluid Declaration
\*---------------------------------------------------------------------------*/

class consistentPisoChannelFlowFluid
:
    public fluidModel
{
    // Private data

        //- Velocity field
        volVectorField U_;

        //- Pressure field
        volScalarField p_;

        //- Pressure field
        volVectorField gradp_;

        //- Pressure gradient force term
        dimensionedScalar gradPft_;

        //- Flux field
        surfaceScalarField phi_;

        //- Divergence free face velocity field
        surfaceVectorField Uf_;

        //- Transport model
        singlePhaseTransportModel laminarTransport_;

        //- Turbulence model
        autoPtr<incompressible::LESModel> sgsModel_;

        //- Sgs viscosity field
        volScalarField nuSgs_;

        //- Density
        dimensionedScalar rho_;

        //- Ubar value
        dimensionedVector Ubar_;

        //- Consistency in time
        const Switch extrapolateFlux_;
  
    // Private Member Functions

        //- Disallow default bitwise copy construct
        consistentPisoChannelFlowFluid(const consistentPisoChannelFlowFluid&);

        //- Disallow default bitwise assignment
        void operator=(const consistentPisoChannelFlowFluid&);

public:

    //- Runtime type information
    TypeName("consistentPisoChannelFlowFluid");

    // Constructors

        //- Construct from components
        consistentPisoChannelFlowFluid(const fvMesh& mesh);


    // Destructor

        virtual ~consistentPisoChannelFlowFluid()
        {}


    // Member Functions

        // Access

            //- Return velocity field
            virtual const volVectorField& U() const;

            //- Return velocity field
            volVectorField& U()
            {
                return U_;
            }

            //- Return pressure field
            virtual const volScalarField& p() const;

            //- Return pressure field
            volScalarField& p()
            {
                return p_;
            }

            //- Return pressure gradient
            volVectorField& gradp()
            {
                return gradp_;
            }

            //- Return pressure gradient force term
            const dimensionedScalar& gradPft()
            {
                return gradPft_;
            }

            //- Return const reference to flux field
            virtual const surfaceScalarField& phi() const
            {
                return phi_;
            }

            //- Return reference to flux field
            virtual surfaceScalarField& phi()
            {
                return phi_;
            }

            //- Return const reference to face velocity field
            const surfaceVectorField& Uf() const
            {
                return Uf_;
            }

            //- Return reference to face velocity field
            surfaceVectorField& Uf()
            {
                return Uf_;
            }

            //- Return effective viscosity field
            volScalarField& nuSgs()
            {
                return nuSgs_;
            }
  
            //- Density
            const dimensionedScalar& rho()
            {
                return rho_;
            }

            //- Ubar
            const dimensionedVector& Ubar()
            {
                return Ubar_;
            }

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> faceZoneViscousForce
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> faceZonePressureForce
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff
            (
                const label zoneID,
                const label patchID
            ) const;

        // Edit

            //- Evolve the fluid model
            virtual void evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
