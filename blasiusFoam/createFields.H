/*
#include "fvCFD.H"

Info << "Creating fields\n" << endl;

// Read the x-velocity field `u`
volScalarField u
(
    IOobject
    (
        "u",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,  // Read the field from the time directory
        IOobject::AUTO_WRITE  // Write the field during the simulation
    ),
    mesh
);

// Read the y-velocity field `v`, initialized to zero
volScalarField v
(
    IOobject
    (
        "v",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,     // Do not read from a file, initialize inside the solver
        IOobject::AUTO_WRITE   // Write the field during the simulation
    ),
    mesh,
    dimensionedScalar("v", dimensionSet(0, 1, -1, 0, 0), 0.0)  // Initialize with zero
);

// Create the velocity vector field `U` from `u` and `v`
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,     // Do not read from a file, calculated from `u` and `v`
        IOobject::AUTO_WRITE   // Write during the simulation
    ),
    mesh,
    //dimensionedVector(U.dimensions(), Zero)  // Initialize `U` to zero vector field
    //vector(Zero, Zero, Zero)  // Initialize U to a zero vector (optional initialization)
    dimensionedVector("U", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0, 0, 0))  // Initialize U as a velocity field (m/s)

);

// Set up the flux field `phi` based on the velocity field `U`
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,     // Flux is calculated in the solver
        IOobject::AUTO_WRITE   // Write the flux field during the simulation
    ),
    fvc::flux(U)  // Initialize phi as the flux of the velocity field `U`
);

// Read or construct the temperature field `T`
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,   // Read the temperature field from the case directory
        IOobject::AUTO_WRITE   // Write the field during the simulation
    ),
    mesh
);

// Attempt to read from transportProperties
IOdictionary transportProperties
(
 IOobject
 (
     "transportProperties",
     runTime.constant(),
     mesh,
     IOobject::MUST_READ,
     IOobject::NO_WRITE
 )
);

dimensionedScalar Pr
(
     "Pr",
     dimless,
     transportProperties
);

dimensionedScalar nu
(
     "nu",
     pow(dimLength,2)/dimTime,
     transportProperties
);

// Construct a pressure field
// If it is available read it otherwise construct from the velocity BCs
// converting fixed-value BCs to zero-gradient and vice versa.

// Allow override from command-line -pName option
const word pName = args.getOrDefault<word>("pName", "p");

// Infer the pressure BCs from the velocity
wordList pBCTypes
(
    U.boundaryField().size(),
    fixedValueFvPatchScalarField::typeName
);

forAll(U.boundaryField(), patchi)
{
    if (U.boundaryField()[patchi].fixesValue())
    {
        pBCTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
}

Info<< "Constructing pressure field " << pName << nl << endl;
volScalarField p
(
    IOobject
    (
        pName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(sqr(dimVelocity), Zero),
    pBCTypes
);


// Infer the velocity potential BCs from the pressure
wordList PhiBCTypes
(
    p.boundaryField().size(),
    zeroGradientFvPatchScalarField::typeName
);

forAll(p.boundaryField(), patchi)
{
    if (p.boundaryField()[patchi].fixesValue())
    {
        PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}

// Construct the velocity potential field `Phi`
// Boundary conditions for Phi inferred from the pressure field if it exists
Info << "Constructing velocity potential field Phi\n" << endl;

// Construct boundary conditions for Phi from pressure `p`
volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, Zero),
    PhiBCTypes
);

label PhiRefCell = 0;
scalar PhiRefValue = 0;


// Ensure that the mesh requires the flux for Phi
mesh.setFluxRequired(Phi.name());

*/

#include "fvCFD.H"

Info << "Creating fields\n" << endl;

// Read the x-velocity field `u`
volScalarField u
(
    IOobject
    (
        "u",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,  // Read the field from the time directory
        IOobject::AUTO_WRITE  // Write the field during the simulation
    ),
    mesh
);

// Read the y-velocity field `v`, initialized to zero
volScalarField v
(
    IOobject
    (
        "v",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,     // Do not read from a file, initialize inside the solver
        IOobject::AUTO_WRITE   // Write the field during the simulation
    ),
    mesh,
    dimensionedScalar("v", dimensionSet(0, 1, -1, 0, 0), 0.0)  // Initialize with zero
);

// Create the velocity vector field `U` from `u` and `v`
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,     // Do not read from a file, calculated from `u` and `v`
        IOobject::AUTO_WRITE   // Write during the simulation
    ),
    mesh,
    dimensionedVector("U", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0, 0, 0))  // Initialize U as a velocity field (m/s)
);

// Set up the flux field `phi` based on the velocity field `U`
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,     // Flux is calculated in the solver
        IOobject::AUTO_WRITE   // Write the flux field during the simulation
    ),
    fvc::flux(U)  // Initialize phi as the flux of the velocity field `U`
);

// Read or construct the temperature field `T`
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,   // Read the temperature field from the case directory
        IOobject::AUTO_WRITE   // Write the field during the simulation
    ),
    mesh
);

// Read transport properties (nu, Pr) from transportProperties
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

/*

dimensionedScalar nu
(
    transportProperties.lookup("nu")
);

dimensionedScalar Pr
(
    transportProperties.lookup("Pr")
);

*/

dimensionedScalar Pr
(
     "Pr",
     dimless,
     transportProperties
);

dimensionedScalar nu
(
     "nu",
     pow(dimLength,2)/dimTime,
     transportProperties
);


/*
// Construct the velocity potential field `Phi`
Info << "Constructing velocity potential field Phi\n" << endl;

wordList PhiBCTypes
(
    U.boundaryField().size(),
    zeroGradientFvPatchScalarField::typeName   // Default BC for Phi: zeroGradient
);

forAll(U.boundaryField(), patchi)
{
    if (U.boundaryField()[patchi].fixesValue())
    {
        PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}

volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength * dimVelocity, Zero),  // Initialize Phi with zero
    PhiBCTypes   // Assign boundary conditions to Phi
);

label PhiRefCell = 0;
scalar PhiRefValue = 0;

// Ensure that the mesh requires the flux for Phi
mesh.setFluxRequired(Phi.name());
*/

// Construct a pressure field
// If it is available read it otherwise construct from the velocity BCs
// converting fixed-value BCs to zero-gradient and vice versa.

// Allow override from command-line -pName option
const word pName = args.getOrDefault<word>("pName", "p");

// Infer the pressure BCs from the velocity
wordList pBCTypes
(
    U.boundaryField().size(),
    fixedValueFvPatchScalarField::typeName
);

forAll(U.boundaryField(), patchi)
{
    if (U.boundaryField()[patchi].fixesValue())
    {
        pBCTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
}

Info<< "Constructing pressure field " << pName << nl << endl;
volScalarField p
(
    IOobject
    (
        pName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(sqr(dimVelocity), Zero),
    pBCTypes
);


// Infer the velocity potential BCs from the pressure
wordList PhiBCTypes
(
    p.boundaryField().size(),
    zeroGradientFvPatchScalarField::typeName
);

forAll(p.boundaryField(), patchi)
{
    if (p.boundaryField()[patchi].fixesValue())
    {
        PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}

// Construct the velocity potential field `Phi`
// Boundary conditions for Phi inferred from the pressure field if it exists
Info << "Constructing velocity potential field Phi\n" << endl;

// Construct boundary conditions for Phi from pressure `p`
volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, Zero),
    PhiBCTypes
);

label PhiRefCell = 0;
scalar PhiRefValue = 0;


// Ensure that the mesh requires the flux for Phi
mesh.setFluxRequired(Phi.name());
