/**
 * A code to read a mesh and write the number of cells to the terminal
 *
 * @author Maalik, ali@tensorfields.com
 */

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IOobject.H"

using namespace Foam;

int main(int argc, char* argv[])
{
    argList args(argc, argv);

    Time runTime("controlDict", args);

	fvMesh mesh
	(
	    IOobject
		(
		    fvMesh::defaultRegion,
			runTime.name(),
			runTime,
			IOobject::MUST_READ
		)
	);

    Info << "nCells = " << mesh.nCells() << endl;
}
