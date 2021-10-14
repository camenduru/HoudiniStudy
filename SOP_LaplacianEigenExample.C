/*
 * Copyright (c) 2019
 *	Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 * DetailAttrib SOP
 */

#include "SOP_LaplacianEigenExample.h"

#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_DSOVersion.h>

#include <Eigen/Sparse>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <iostream>

using namespace Eigen;
using namespace Spectra;
using namespace std;
void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "my_laplacianeigensolver",
        "LaplacianEigensolverExample",
        SOP_LaplacianEigenExample::myConstructor,
        SOP_LaplacianEigenExample::myTemplateList,
        1,
        1));
}


static PRM_Name sop_names[] = {
    PRM_Name("vecattribname", "Eigenvalue Name"),
    PRM_Name("valattribname", "Eigenvector Name"),
    PRM_Name("returncount", "Return Count"),
};

static PRM_Default prm_defaults[] = {
    PRM_Default(0, "eigenVectors"),
    PRM_Default(0, "eigenValues"),
    PRM_Default(4, "")
} ;

static PRM_Range sop_valueRange(PRM_RANGE_RESTRICTED, 2, PRM_RANGE_UI, 100);

PRM_Template SOP_LaplacianEigenExample::myTemplateList[] =
{
        PRM_Template(PRM_STRING, 1, &sop_names[0], &prm_defaults[0]),
        PRM_Template(PRM_STRING, 1, &sop_names[1], &prm_defaults[1]),
        PRM_Template(PRM_INT, 1, &sop_names[2], &prm_defaults[2], 0,
                        &sop_valueRange),
        PRM_Template()
};


OP_Node *
SOP_LaplacianEigenExample::myConstructor(OP_Network *net,const char *name,OP_Operator *entry)
{
    return new SOP_LaplacianEigenExample(net, name, entry);
}


SOP_LaplacianEigenExample::SOP_LaplacianEigenExample(OP_Network *net, const char *name, OP_Operator *entry)
    : SOP_Node(net, name, entry)
{
    // This indicates that this SOP manually manages its data IDs,
    // so that Houdini can identify what attributes may have changed,
    // e.g. to reduce work for the viewport, or other SOPs that
    // check whether data IDs have changed.
    // By default, (i.e. if this line weren't here), all data IDs
    // would be bumped after the SOP cook, to indicate that
    // everything might have changed.
    // If some data IDs don't get bumped properly, the viewport
    // may not update, or SOPs that check data IDs
    // may not cook correctly, so be *very* careful!
    mySopFlags.setManagesDataIDs(true);
}

SOP_LaplacianEigenExample::~SOP_LaplacianEigenExample()
{
}


OP_ERROR
SOP_LaplacianEigenExample::cookMySop(OP_Context &context)
{
    fpreal t = context.getTime();

    // We must lock our inputs before we try to access their geometry.
    // OP_AutoLockInputs will automatically unlock our inputs when we return.
    // NOTE: Don't call unlockInputs yourself when using this!
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    duplicateSource(0, context);

    UT_String eVectName;
    evalString(eVectName, 0, 0, t);
    UT_String eValName;
    evalString(eValName, 1, 0, t);
    exint eigencount;
    eigencount = evalInt(2, 0, t);

    UT_String neighboursname = UT_String("neighbours");
    UT_String weightsumname = UT_String("weightsum");
    UT_String weightsname = UT_String("weights");

    GA_RWHandleIA nAttr(gdp->findIntArray(GA_ATTRIB_POINT, neighboursname));
    GA_ROHandleF wsAttr(gdp->findFloatTuple(GA_ATTRIB_POINT, weightsumname));
    GA_ROHandleFA wAttr(gdp->findFloatArray(GA_ATTRIB_POINT, weightsname));

    int numPoints = (int)gdp->getNumPoints();
    SparseMatrix<double> mat = SparseMatrix<double, RowMajor>(numPoints, numPoints);
    int numEigenVs = Index((numPoints - 1 < (int)eigencount ? numPoints - 1 : (int)eigencount));

    if(nAttr.isValid() && wsAttr.isValid() && wAttr.isValid() && numPoints >= 2) {
        UT_Int32Array nArray;
        UT_FloatArray wArray;
        
        GA_Range ptRange = gdp->getPointRange();
        int counter = 0;

        for (GA_Iterator it(ptRange); !it.atEnd(); ++it)
        {
            GA_Offset offset = *it;
            nAttr.get(offset, nArray);
            wAttr.get(offset, wArray);
            fpreal32 wsVal = wsAttr.get(offset); 

            bool addingLower = false;
            for (int i = 0; i < nArray.size(); i++)
            {   
                if(nArray[i] > counter && !addingLower) {
                    mat.insert(counter, counter) = (double)wsVal;
                    addingLower = true;
                } 
                mat.insert(counter, nArray[i]) = (double)(wArray[i]);
            }
            if (!addingLower) mat.insert(counter, counter) = (double)wsVal;  
            counter++;
        }

        Spectra::SparseSymShiftSolve<double> op(mat);
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(op, numEigenVs, numEigenVs + 1, 0.0);
        eigs.init();
        eigs.compute();

        if (eigs.info() == Spectra::CompInfo::Successful)
        {
            VectorXd evalues = eigs.eigenvalues();
            MatrixXd evectors = eigs.eigenvectors();

            int eVectOutSize = (int)evectors.outerSize();
            int eValInSize = (int)evalues.innerSize();

            auto eVectAttr = gdp->addFloatArray(GA_ATTRIB_POINT, eVectName, eVectOutSize);
            GA_RWHandleFA eVectAttrHandle = GA_RWHandleFA(eVectAttr);

            auto eValAttr = gdp->addFloatArray(GA_ATTRIB_DETAIL, eValName, eValInSize);
            GA_RWHandleFA eValAttrHandle = GA_RWHandleFA(eValAttr);
            UT_FloatArray eValArr = UT_FloatArray();           

            if (eVectAttrHandle.isValid())
            {
                UT_FloatArray evArr = UT_FloatArray(exint(eVectOutSize));
                evArr.setSize(exint(eVectOutSize));
                counter = 0;
                for (GA_Iterator it(ptRange); !it.atEnd(); ++it)
                {
                    GA_Offset offset = *it;
                    for (int j = 0; j < eVectOutSize; j++)
                    {
                        evArr[j] = evectors(counter, j);
                    } 
                    eVectAttrHandle.set(GA_Offset(*it), evArr);
                    counter++;
                }
                eVectAttrHandle.bumpDataId();
            }

            if(eValAttrHandle.isValid()) 
            {
                eValArr.setSize(eValInSize);
                for(int i = 0; i < eValInSize; i++) 
                {
                    eValArr[i] = evalues[i];
                } 
                eValAttrHandle.set(0, eValArr);                 
                eValAttrHandle.bumpDataId();  
            }
        } 
    }

    return error();
}

