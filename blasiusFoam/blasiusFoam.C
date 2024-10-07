/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    simpleFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flows.

    \heading Solver details
    The solver uses the SIMPLE algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    simpleControl simple(mesh);

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << endl;

            // Insert Info statements here to print out the dimensions of phi, u, and nu
		/*
        Info << "phi dimensions: " << phi.dimensions() << endl;
        Info << "u dimensions: " << u.dimensions() << endl;
        Info << "nu dimensions: " << nu.dimensions() << endl;
        Info << "U dimensions: " << U.dimensions() << endl;
		*/


        // Solve the x-momentum equation for u
        tmp<fvScalarMatrix> tuEqn
        (
            fvm::div(phi, u) == fvm::laplacian(nu, u) // or fvm::ddt(u) for transient problems
        );

		fvScalarMatrix& uEqn = tuEqn.ref();

        uEqn.relax();
        //uEqn.solve();
        solve(uEqn);
        //u.correctBoundaryConditions();


        // Correct boundary conditions for U
        //U.correctBoundaryConditions();
		

        
        // Correct the velocity flux (phi)
        //surfaceScalarField phi = fvc::flux(U);  // Define phi


		

        phi = fvc::flux(U);  // Recalculate phi from U
		

        // Solve the continuity equation using the Laplace equation for Phi
        //fvScalarMatrix PhiEqn(fvm::laplacian(Phi));

		
        tmp<fvScalarMatrix> tPhiEqn
        (
            fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
         ==
            fvc::div(phi)
        );
		/*

        fvScalarMatrix PhiEqn
        (
            //fvm::laplacian(Phi) == fvc::div(phi)  // Solve Poisson equation

            fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
         ==
            fvc::div(phi)

        );
		*/

		fvScalarMatrix& PhiEqn = tPhiEqn.ref();
        //PhiEqn.setReference(PhiRefCell, PhiRefValue);
        PhiEqn.relax();
        PhiEqn.solve();
        //Phi.correctBoundaryConditions();

        phi -= PhiEqn.flux();

        // Recalculate velocity flux (phi) and correct U
        //phi = fvc::flux(U);  // Correct the face flux

        volVectorField UIntermediate(fvc::reconstruct(phi));
        v = UIntermediate.component(vector::Y);
        v.correctBoundaryConditions();  // Correct the boundary conditions for v
        // Update U vector field from the new x-component u
	
	
        U = vector(1, 0, 0) * u + vector(0, 1, 0) * v;
		/*

        // Solve the temperature equation

        //fvScalarMatrix TEqn
        //(
        //    fvm::div(phi, T) == fvm::laplacian(nu / Pr, T)
        //);
        //TEqn.relax();
        //TEqn.solve();
        //T.correctBoundaryConditions();

		*/
        tmp<fvScalarMatrix> tTEqn
        ( 
            fvm::div(phi, T)
            == 
            fvm::laplacian(nu / Pr, T)
        ); 

        fvScalarMatrix& TEqn = tTEqn.ref();

        TEqn.relax();

        solve(TEqn);

		T.correctBoundaryConditions();
		
        

        // Write results for this time step
        runTime.write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
