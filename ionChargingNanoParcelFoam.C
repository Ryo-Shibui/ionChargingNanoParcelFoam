/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "fvCFD.H"
#include "rhoReactionThermo.H"
#include "SLGThermo.H"
#include "basicReactingMultiphaseCloud.H"
#include "IonCharging.H"

int main(int argc, char *argv[])
{
    using namespace Foam;
    argList::addNote("One-way aerosol tracking with stochastic bipolar ion attachment.");
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word fieldTime(runTime.controlDict().getOrDefault<word>("carrierFieldTime", "0"));
    const volVectorField U
    (
        IOobject("U", fieldTime, mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh
    );
    const volVectorField E
    (
        IOobject("E", fieldTime, mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh
    );
    const volScalarField nP
    (
        IOobject("nP", fieldTime, mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh
    );
    const volScalarField nN
    (
        IOobject("nN", fieldTime, mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh
    );
    if
    (
        nP.dimensions() != dimless/dimVolume || nN.dimensions() != dimless/dimVolume
     || E.dimensions() != dimensionSet(1, 1, -3, 0, 0, -1, 0)
    )
    {
        FatalErrorInFunction
            << "Use number densities in m^-3 and an electric field in V/m."
            << exit(FatalError);
    }
    // In m^-3. Suppress only physically negligible CFD undershoot in the rate.
    const scalar negativeRoundoff = 1e-6;
    for (const volScalarField* n : {&nP, &nN})
    {
        for (const scalar value : n->primitiveField())
        {
            if (!std::isfinite(value) || value < -negativeRoundoff)
            {
                FatalErrorInFunction
                    << n->name() << " has a non-finite or negative internal value."
                    << exit(FatalError);
            }
        }
    }

    const IOdictionary physicalProperties
    (
        IOobject("physicalProperties", runTime.constant(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE)
    );

    autoPtr<rhoReactionThermo> thermo(rhoReactionThermo::New(mesh));
    thermo->validate(args.executable(), "h", "e");
    SLGThermo slgThermo(mesh, thermo());
    const volScalarField rho
    (
        IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        thermo->rho()
    );
    #include "readGravitationalAcceleration.H"

    basicReactingMultiphaseCloud parcels("reactingCloud1", rho, U, g, slgThermo);
    const dictionary& cloud = parcels.particleProperties();
    const dictionary& constants = cloud.subDict("constantProperties");
    const dictionary& models = cloud.subDict("subModels");
    if (parcels.solution().coupled() || !parcels.solution().transient())
    {
        FatalIOErrorInFunction(cloud)
            << "Use a transient cloud with coupled false." << exit(FatalIOError);
    }
    if (!constants.getOrDefault<Switch>("diameterChange", false))
    {
        for (const word& key : wordList
            ({"phaseChangeModel", "devolatilisationModel", "surfaceReactionModel"}))
        {
            if (models.getOrDefault<word>(key, "none") != "none")
            {
                FatalIOErrorInFunction(models)
                    << key << " must be none when diameterChange is false."
                    << exit(FatalIOError);
            }
        }
    }

    using KinematicCloudType = basicReactingMultiphaseCloud::kinematicCloudType;
    auto& functions = parcels.functions();
    functions.append
    (
        new IonCharging<KinematicCloudType>(static_cast<KinematicCloudType&>(parcels), nP, nN, physicalProperties)
    );

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        parcels.evolve();
        runTime.write();
    }
    Info<< "End" << endl;
    return 0;
}
