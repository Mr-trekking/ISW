/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    setFields

Group
    grpPreProcessingUtilities

Description
    Set values on a selected set of cells/patch-faces via a dictionary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "faceSet.H"
#include "volFields.H"

using namespace Foam;

template<class Type>
bool setCellFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedCells,
    Istream& fieldValueStream
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);

    // Check the current time directory
    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check the "constant" directory
    if (!fieldHeader.typeHeaderOk<fieldType>(true))
    {
        fieldHeader = IOobject
        (
            fieldName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        );
    }

    // Check field exists
    if (fieldHeader.typeHeaderOk<fieldType>(true))
    {
        Info<< "    Setting internal values of "
            << fieldHeader.headerClassName()
            << " " << fieldName << endl;

        fieldType field(fieldHeader, mesh, false);

        const Type& value = pTraits<Type>(fieldValueStream);

        // Parameter input
        const scalar h1 = 0.1; 

        const scalar h2 = 0.4;

        const scalar rho1 = 996;

        const scalar rho2 = 1030;

        // the internal solitary wave amplitude
        const scalar a = 0.065; 

        // initial position
        const scalar x0 = 12.501;

        const scalar g = 9.81;
        		
        // Velocity field calculation
        scalar c02 = g * h1 * h2 * (rho2 - rho1) / (rho1 * h2 + rho2 * h1);
        const scalar c0 = Foam::sqrt(c02);
    
        const scalar c1 = -1.5 * c0 * (rho1 * h2 * h2 - rho2 * h1 * h1) / ( rho1 * h1 * h2 * h2 + rho2 * h1 * h1 * h2);

        const scalar c2 = c0/6.0 * (rho1 * h1 * h1 * h2 + rho2 * h1 * h2 * h2) / ( rho1  * h2 + rho2 * h1 );

        const scalar term1 = 3 * c0 /( h1 * h1 * h2 * h2);

        const scalar term2 = 0.875 * pow((rho1 * h2 * h2 - rho2 * h1 * h1)/(rho1 * h2 + rho2 * h1),2);

        const scalar term3 = (rho1 * h2 * h2 * h2 + rho2 * h1 * h1 * h1)/(rho1 * h2 + rho2 * h1);

        const scalar c3 =  term1 * (term2 - term3);
        
        const scalar H = -a;

        const scalar B = -1 * H * c3 / (2 * c1 + H * c3);
       
        const scalar lambda2 = H * (c1 + 0.5 * c3 * H) / (12 * c2);

        const scalar lambda = Foam::sqrt(lambda2);
    
        const scalar ce = c0 + H / 3.0 * (c1 + 0.5 * c3 * H);

        forAll(field,cellI)
        {
            scalar pos = lambda*(mesh.C()[cellI].x()-x0);
            const scalar eta = H/(B + (1 - B)*sqr(Foam::cosh(pos)));


            //R TO L
            if (mesh.C()[cellI].z() - h2 > eta)
            {
                field[cellI] = 1 * ce * eta/(h1 - eta) * value;
            }
	    else
            {
                field[cellI] = -1 * ce * eta/(h2 + eta) * value;
            }
        }
/*
            //L TO R
            if (mesh.C()[cellI].z() - h2 > eta)
            {
                field[cellI] = -1 * ce * eta/(h1 - eta) * value;
            }
            else
            {
                field[cellI] = 1 * ce * eta/(h2 + eta) * value;
            }
        }
*/
        typename GeometricField<Type, fvPatchField, volMesh>::
            Boundary& fieldBf = field.boundaryFieldRef();

        forAll(field.boundaryField(), patchi)
        {
            fieldBf[patchi] = fieldBf[patchi].patchInternalField();
        }

        if (!field.write())
        {
            FatalErrorInFunction
              << "Failed writing field " << fieldName << endl;
        }
    }
    else
    {
        WarningInFunction
          << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    return true;
}


class setCellField
{

public:

    setCellField()
    {}

    autoPtr<setCellField> clone() const
    {
        return autoPtr<setCellField>::New();
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedCells_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedCells)
        :
            mesh_(mesh),
            selectedCells_(selectedCells)
        {}

        iNew(const fvMesh& mesh, labelList&& selectedCells)
        :
            mesh_(mesh),
            selectedCells_(std::move(selectedCells))
        {}

        autoPtr<setCellField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setCellFieldType<scalar>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<vector>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<symmTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<tensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                )
            )
            {
                WarningInFunction
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setCellField>::New();
        }
    };
};


template<class Type>
bool setFaceFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedFaces,
    Istream& fieldValueStream
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);

    // Check the current time directory
    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check the "constant" directory
    if (!fieldHeader.typeHeaderOk<fieldType>(true))
    {
        fieldHeader = IOobject
        (
            fieldName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        );
    }

    // Check field exists
    if (fieldHeader.typeHeaderOk<fieldType>(true))
    {
        Info<< "    Setting patchField values of "
            << fieldHeader.headerClassName()
            << " " << fieldName << endl;

        fieldType field(fieldHeader, mesh);

        const Type& value = pTraits<Type>(fieldValueStream);

        // Create flat list of selected faces and their value.
        Field<Type> allBoundaryValues(mesh.nBoundaryFaces());
        forAll(field.boundaryField(), patchi)
        {
            SubField<Type>
            (
                allBoundaryValues,
                field.boundaryField()[patchi].size(),
                field.boundaryField()[patchi].patch().start()
              - mesh.nInternalFaces()
            ) = field.boundaryField()[patchi];
        }

        // Override
        bool hasWarned = false;
        labelList nChanged
        (
            returnReduce(field.boundaryField().size(), maxOp<label>()),
            0
        );
        forAll(selectedFaces, i)
        {
            label facei = selectedFaces[i];
            if (mesh.isInternalFace(facei))
            {
                if (!hasWarned)
                {
                    hasWarned = true;
                    WarningInFunction
                        << "Ignoring internal face " << facei
                        << ". Suppressing further warnings." << endl;
                }
            }
            else
            {
                label bFacei = facei-mesh.nInternalFaces();
                allBoundaryValues[bFacei] = value;
                nChanged[mesh.boundaryMesh().patchID()[bFacei]]++;
            }
        }

        Pstream::listCombineGather(nChanged, plusEqOp<label>());
        Pstream::listCombineScatter(nChanged);

        typename GeometricField<Type, fvPatchField, volMesh>::
            Boundary& fieldBf = field.boundaryFieldRef();

        // Reassign.
        forAll(field.boundaryField(), patchi)
        {
            if (nChanged[patchi] > 0)
            {
                Info<< "    On patch "
                    << field.boundaryField()[patchi].patch().name()
                    << " set " << nChanged[patchi] << " values" << endl;
                fieldBf[patchi] == SubField<Type>
                (
                    allBoundaryValues,
                    fieldBf[patchi].size(),
                    fieldBf[patchi].patch().start()
                  - mesh.nInternalFaces()
                );
            }
        }

        if (!field.write())
        {
            FatalErrorInFunction
                << "Failed writing field " << field.name() << exit(FatalError);
        }
    }
    else
    {
        WarningInFunction
          << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    return true;
}


class setFaceField
{

public:

    setFaceField()
    {}

    autoPtr<setFaceField> clone() const
    {
        return autoPtr<setFaceField>::New();
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedFaces_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedFaces)
        :
            mesh_(mesh),
            selectedFaces_(selectedFaces)
        {}

        iNew(const fvMesh& mesh, labelList&& selectedFaces)
        :
            mesh_(mesh),
            selectedFaces_(std::move(selectedFaces))
        {}

        autoPtr<setFaceField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setFaceFieldType<scalar>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<vector>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<symmTensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<tensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                )
            )
            {
                WarningInFunction
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setFaceField>::New();
        }
    };
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Set values on a selected set of cells/patch-faces via a dictionary"
    );

    argList::addOption("dict", "file", "Use alternative setFieldsDict");

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("setUFieldsDict");
    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictName << "\n" << endl;

    IOdictionary setFieldsDict(dictIO);

    if (setFieldsDict.found("defaultFieldValues"))
    {
        Info<< "Setting field default values" << endl;
        PtrList<setCellField> defaultFieldValues
        (
            setFieldsDict.lookup("defaultFieldValues"),
            setCellField::iNew(mesh, labelList(mesh.nCells()))
        );
        Info<< endl;
    }


    Info<< "Setting field region values" << endl;

    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    forAll(regions, regionI)
    {
        const entry& region = regions[regionI];

        autoPtr<topoSetSource> source =
            topoSetSource::New(region.keyword(), mesh, region.dict());

        if (source().setType() == topoSetSource::CELLSETSOURCE)
        {
            cellSet selectedCellSet
            (
                mesh,
                "cellSet",
                mesh.nCells()/10+1  // Reasonable size estimate.
            );

            source->applyToSet
            (
                topoSetSource::NEW,
                selectedCellSet
            );

            PtrList<setCellField> fieldValues
            (
                region.dict().lookup("fieldValues"),
                setCellField::iNew(mesh, selectedCellSet.sortedToc())
            );
        }
        else if (source().setType() == topoSetSource::FACESETSOURCE)
        {
            faceSet selectedFaceSet
            (
                mesh,
                "faceSet",
                mesh.nBoundaryFaces()/10+1
            );

            source->applyToSet
            (
                topoSetSource::NEW,
                selectedFaceSet
            );

            PtrList<setFaceField> fieldValues
            (
                region.dict().lookup("fieldValues"),
                setFaceField::iNew(mesh, selectedFaceSet.sortedToc())
            );
        }
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
