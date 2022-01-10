/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 2021 AUTHOR,AFFILIATION
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
    BarycentricTurbulence

Description
    Utility for calculating barycentric coordinates of volSymmTensorField

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "SymmTensor.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "tensor.H"
#include "Identity.H"
#include "argList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcRGB
(
    const fvMesh& mesh,
    const Time& runTime,
    const volSymmTensorField& field,
    const word fieldName
)
{
          volVectorField RGB
          (
              IOobject
              (
                  fieldName + ".barycentricTurbulence",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              mesh,
              dimensionedVector(dimless, Zero)
          );


	forAll( RGB,i)
        {
        
		auto b = field[i]/(tr(field[i]))-I/3; // Anysotropy tensor 
               	auto ev = eigenValues(b);             // Vector of eigenvalues in ascending order

		// Eigenvalues sorted in descending order
		scalar l1 = ev.z(); 
		scalar l2 = ev.y(); 
		scalar l3 = ev.x(); 

        	// Calculating barycentric coordinates
        	RGB[i].x() = l1 - l2;
        	RGB[i].y() = 2*(l2 - l3);
        	RGB[i].z() = 3*l3 + 1;
	}
          
        Info << "Writting " << RGB.name() <<" to " << runTime.timeName() << endl;
	
	RGB.write();

}





int main(int argc, char *argv[])
{
   


    Foam::argList::addOption
    (
        "fields",
        "<wordList>",
        "specify names of symmTensor fields for which the barycentric coordinates will be calculated "
        "example: -fields (UPrime2Mean RMean)"
    );

    
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
 
   
    List<word> selectedFields;
    args.lookup("fields")() >> selectedFields;
   
   	
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        
	Info<< "Time = " << runTime.timeName() << endl;

	forAll(selectedFields,i)
	{
		word fieldName = selectedFields[i];
		
		IOobject fieldHeader
        	(
            		fieldName,
            		runTime.timeName(),
            		mesh,
            		IOobject::MUST_READ,
            		IOobject::NO_WRITE
       		);


		autoPtr<volSymmTensorField> fieldPtr_(nullptr);

		bool foundField = fieldHeader.typeHeaderOk<volSymmTensorField>();
	
		if(!foundField)
		{  
			Info <<	"Not found field " << fieldName << " at time: " << runTime.timeName() << endl;
		}
		else
        	{
			Info << "Reading field " << fieldName << endl;

                	fieldPtr_.reset
                	(
                     		new volSymmTensorField
                     		(
                        		fieldHeader,
                        		mesh
                     		)
                	);

                	calcRGB(mesh,runTime, fieldPtr_.ref(),fieldName);

       		}


 
	}


    }

    Info<< "End\n" << endl;

    return 0;
        
    
}


// ************************************************************************* //
