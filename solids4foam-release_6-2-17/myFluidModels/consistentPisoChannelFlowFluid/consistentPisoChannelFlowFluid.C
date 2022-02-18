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

\*---------------------------------------------------------------------------*/

#include "consistentPisoChannelFlowFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
#include "steadyStateDdtScheme.H"

#include "inletOutletFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "skewCorrectionVectors.H"
#include "linear.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(consistentPisoChannelFlowFluid, 0);
addToRunTimeSelectionTable(fluidModel, consistentPisoChannelFlowFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

consistentPisoChannelFlowFluid::consistentPisoChannelFlowFluid(const fvMesh& mesh)
:
    fluidModel(this->typeName, mesh),
    U_
    (
        IOobject
        (
            "U",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    gradp_(fvc::grad(p_)),
    gradPft_
    (
            "gradPft",
            dimensionSet(0, 1, -2, 0, 0),
            0.0
    ),
    phi_
    (
        IOobject
        (
            "phi",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_) & mesh.Sf()
    ),
    Uf_
    (
        IOobject
        (
            "Uf",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_)
    ),
    laminarTransport_(U_, phi_),
    sgsModel_
    (
        incompressible::LESModel::New
        (
            U_, phi_, laminarTransport_
        )
    ),
    rho_
    (
        IOdictionary
        (
            IOobject
            (
                "transportProperties",
                runTime().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("rho")
    ),
    Ubar_
    (
        IOdictionary
        (
            IOobject
            (
                "channelFlowProperties",
                runTime().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("Ubar")
    ),
    extrapolateFlux_
    (
        fluidProperties().lookupOrDefault<Switch>("extrapolateFlux", false)
    )
{
  // Info << "Uf: " << Uf_ << endl;
  // Info << "Uf.oldTime: " << Uf_.oldTime() << endl;
    phi_.oldTime().oldTime();
    Uf_.oldTime().oldTime();
    phi_.oldTime().oldTime();
    Uf_.oldTime().oldTime();
    phi_.oldTime().oldTime();
    Uf_.oldTime().oldTime();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const volVectorField& consistentPisoChannelFlowFluid::U() const
{
    return U_;
}


const volScalarField& consistentPisoChannelFlowFluid::p() const
{
    return p_;
}


tmp<vectorField> consistentPisoChannelFlowFluid::patchViscousForce
(
    const label patchID
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & sgsModel_->devBeff()().boundaryField()[patchID]
        );

//     vectorField n = mesh().boundary()[patchID].nf();
//     tvF() -= n*(n&tvF());

    return tvF;
}


tmp<scalarField> consistentPisoChannelFlowFluid::patchPressureForce
(
    const label patchID
) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<vectorField> consistentPisoChannelFlowFluid::faceZoneViscousForce
(
    const label zoneID,
    const label patchID
) const
{
    vectorField pVF = patchViscousForce(patchID);

    tmp<vectorField> tvF
    (
        new vectorField(mesh().faceZones()[zoneID].size(), vector::zero)
    );
    vectorField& vF = tvF();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pVF, i)
    {
        vF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            pVF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(vF, sumOp<vectorField>());


    return tvF;
}


tmp<scalarField> consistentPisoChannelFlowFluid::faceZonePressureForce
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pPF = patchPressureForce(patchID);

    tmp<scalarField> tpF
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& pF = tpF();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pPF, i)
    {
        pF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            pPF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(pF, sumOp<scalarField>());

    return tpF;
}


tmp<scalarField> consistentPisoChannelFlowFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pMuEff =
        rho_.value()*sgsModel_->nuEff()().boundaryField()[patchID];

    tmp<scalarField> tMuEff
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& muEff = tMuEff();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pMuEff, i)
    {
        muEff[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            pMuEff[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(muEff, sumOp<scalarField>());

    return tMuEff;
}


void consistentPisoChannelFlowFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    const fvMesh& mesh = fluidModel::mesh();

    const int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    const int nCorr(readInt(fluidProperties().lookup("nCorrectors")));
    
    int nOuterCorr =
        readInt(fluidProperties().lookup("nOuterCorrectors"));
    if (nOuterCorr < 1)
    {
        nOuterCorr = 1;
    }

    scalar convergenceCriterion = 0;
    fluidProperties().readIfPresent("convergence", convergenceCriterion);

    // Read channel flow properties and calculate flow direction
    dimensionedScalar magUbar = mag(Ubar_);
    vector flowDirection = (Ubar_/magUbar).value();

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p_, fluidProperties(), pRefCell, pRefValue);

    phi_.oldTime().oldTime();
    Uf_.oldTime().oldTime();

    // Info << runTime().timeIndex() << endl;
    
    if (extrapolateFlux_)
    {
        p_ = 2*p_.oldTime() - p_.oldTime().oldTime();
        p_.correctBoundaryConditions();
#       include "calcPressureGradient.H"

        U_ = 2*U_.oldTime() - U_.oldTime().oldTime();
        U_.correctBoundaryConditions();
        
        phi_ = 2*phi_.oldTime() - phi_.oldTime().oldTime();
        // Correct phi at the boundary
        forAll(phi_.boundaryField(), patchI)
        {
            if (U_.boundaryField()[patchI].fixesValue())
            {
                phi_.boundaryField()[patchI] =
                (
                    U_.boundaryField()[patchI]
                  & mesh.Sf().boundaryField()[patchI]
                );
            }
        }
        
        
        if (nOuterCorr == 1 && runTime().timeIndex() < 3)
        {
            nOuterCorr = 100;
        }
    }

    // if (runTime().timeIndex() == 1)
    // {
    //     nOuterCorr = 10;
    // }

    sgsModel_->correct();

    for (int oCorr = 0; oCorr < nOuterCorr; oCorr++)
    {
        scalar eqnResidual = 1, maxResidual = 0;

        if (mesh.moving())
        {
            // Make the fluxes relative
            phi_ -= fvc::meshPhi(U_);
        }
        
        // Calculate CourantNo
        {
            scalar CoNum = 0.0;
            scalar meanCoNum = 0.0;
            scalar velMag = 0.0;

            if (mesh.nInternalFaces())
            {
                surfaceScalarField SfUfbyDelta =
                    mesh.surfaceInterpolation::deltaCoeffs()*mag(phi_);

                CoNum =
                    max(SfUfbyDelta/mesh.magSf()).value()
                   *runTime().deltaT().value();

                meanCoNum =
                    (sum(SfUfbyDelta)/sum(mesh.magSf())).value()
                   *runTime().deltaT().value();

                velMag = max(mag(phi_)/mesh.magSf()).value();
            }

            Info<< "Courant Number mean: " << meanCoNum
                << " max: " << CoNum
                << " velocity magnitude: " << velMag << endl;
        }
        
        // Time derivative matrix
        fvVectorMatrix UEqnT(fvm::ddt(U_));

        // Construct momentum equation
        // Convection-diffusion matrix
        fvVectorMatrix UEqnH
        (
            fvm::div(phi_, U_)
          + sgsModel_->divDevBeff(U_)
          ==
            flowDirection*gradPft_
        );

        fvVectorMatrix ueqn = UEqnT + UEqnH;

        tmp<fvVectorMatrix> UEqn(UEqnT + UEqnH);

        // Solve momentum equation
        eqnResidual =
            solve(UEqn() == -gradp_).initialResidual();
        // solve(UEqn() == -gradp_);
        maxResidual = max(eqnResidual, maxResidual);

        volScalarField raU = 1/ueqn.A();

        UEqn.clear();

        for (int corr = 0; corr < nCorr; corr++)
        {
            volScalarField rAU = 1/UEqnH.A();
            surfaceScalarField rAUf = fvc::interpolate(rAU);

            U_ = UEqnH.H()*rAU;
#           include "calcConsistentPhiPISO.H"

            adjustPhi(phi_, U_, p_);

            for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
            {
                // Construct pressure equation
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAUf, p_) == fvc::div(phi_)
                );

                // Solve pressure equation

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi_ -= pEqn.flux();
                }
            }

            // Calculate Continuity error
            {
                volScalarField contErr = fvc::div(phi_);

                scalar sumLocalContErr =
                    runTime().deltaT().value()
                   *mag(contErr)().weightedAverage(mesh.V()).value();

                scalar globalContErr =
                    runTime().deltaT().value()
                   *contErr.weightedAverage(mesh.V()).value();

                Info<< "time step continuity errors : sum local = "
                    << sumLocalContErr << ", global = "
                    << globalContErr << endl;
            }

#           include "calcPressureGradient.H"
            // gradp_ = fvc::grad(p_);

            U_ = 1.0/(1/rAU + UEqnT.A())*
                (
                    U_/rAU + UEqnT.H() - gradp_
                );

            U_.correctBoundaryConditions();
        }

        // Correct driving force for a constant mass flow rate

        // Extract the velocity in the flow direction
        dimensionedScalar magUbarStar =
            (flowDirection & U_)().weightedAverage(mesh.V());

        // Calculate the pressure gradient increment needed to
        // adjust the average flow-rate to the correct value
        dimensionedScalar gragPplus =
            (magUbar - magUbarStar)/(raU).weightedAverage(mesh.V());

        U_ += flowDirection*(raU)*gragPplus;

        gradPft_ += gragPplus;

        Info << "Uncorrected Ubar = " << magUbarStar.value() << tab
            << "Pressure gradient = " << gradPft_.value() << endl;

        if (maxResidual < convergenceCriterion)
        {
            Info<< "reached convergence criterion: "
                << convergenceCriterion << endl;
            Info<< "Number of iterations: " << oCorr << endl;
            break;
        }
    }

    // Update divergence free face velocity field
    Uf_ = fvc::interpolate(U_);
    Uf_ -= (mesh.Sf() & Uf_)*mesh.Sf()/magSqr(mesh.Sf());
    Uf_ += phi_*mesh.Sf()/magSqr(mesh.Sf());    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
