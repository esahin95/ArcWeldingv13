const label NC = mesh.nCells();
simpleMatrix<scalar> A(NC);

scalarField& src = A.source();
for(label i=0; i<NC; i++)
{
    src[i] = 0.0;
    for(label j=0; j<NC; j++)
    {
        A[i][j] = 0.0;
    }
}

const scalarField& TDiag = TEqn.diag();
for(label i=0; i<NC; i++)
{
    A[i][i] = TDiag[i];
    src[i] = TEqn.source()[i];
}

const lduAddressing& ldu = TEqn.lduAddr();
for(label i=0; i<ldu.lowerAddr().size(); i++)
{
    const label u = ldu.upperAddr()[i];
    const label l = ldu.lowerAddr()[i];
    A[l][u] = TEqn.upper()[i];
    A[u][l] = TEqn.lower()[i];
}

const volScalarField::Boundary& TBoundary = T.boundaryField();
forAll(TBoundary, patchI)
{
    const fvPatch& TPatch = TBoundary[patchI].patch();
    forAll(TPatch, faceI)
    {
        label cellI = TPatch.faceCells()[faceI];
        A[cellI][cellI] += TEqn.internalCoeffs()[patchI][faceI];
        src[cellI] += TEqn.boundaryCoeffs()[patchI][faceI];
    }
}

/*
// Print Indexing
const scalarField& upper = TEqn.upper();
forAll(upper, faceI)
{
    Info<<upper[faceI]<< " ";
}
Info<<endl;

const scalarField& lower = TEqn.lower();
forAll(lower, faceI)
{
    Info<<lower[faceI]<< " ";
}
Info<<endl;

const labelUList& lowerAddr = ldu.lowerAddr();
forAll(lowerAddr, faceI)
{
    Info<<lowerAddr[faceI]<< " ";
}
Info<<endl;

const labelUList& upperAddr = ldu.upperAddr();
forAll(upperAddr, faceI)
{
    Info<<upperAddr[faceI]<< " ";
}
Info<<endl;
*/

// Print source vector
for(label i=0; i<NC; i++)
{
    Info<<src[i] << " ";
}
Info<<endl;

// Print matrix coefficients
for(label i=0; i<NC; i++)
{
    for(label j=0; j<NC; j++)
    {
        Info<<A[i][j]<<" ";
    }
    Info<<endl;
}