//============================================================================
//  The contents of this file are covered by the Viskores license. See
//  LICENSE.txt for details.
//
//  By contributing to this file, all contributors agree to the Developer
//  Certificate of Origin Version 1.1 (DCO 1.1) as stated in DCO.txt.
//============================================================================

//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#include <iostream>
#include <viskores/cont/Initialize.h>
#include <viskores/cont/DataSet.h>
#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/Invoker.h>
#include <viskores/worklet/WorkletPointNeighborhood.h>
#include <viskores/io/VTKDataSetReader.h>
#include <viskores/io/VTKDataSetWriter.h>

// Worklet: computes the mean and variance of the divergence of a 2D vector field by operating on a neighborhood around each point to compute finite differences.
struct ComputeDivergence : public viskores::worklet::WorkletPointNeighborhood
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInNeighborhood meanX,
                                    FieldInNeighborhood varX,
                                    FieldInNeighborhood meanY,
                                    FieldInNeighborhood varY,
                                    FieldOut div_mean,
                                    FieldOut div_variance);
    using ExecutionSignature = void(Boundary, _2, _3, _4, _5, _6, _7);
    using InputDomain = _1;

    template<typename BoundaryType, typename NeighborhoodType, typename OutType>
    VISKORES_EXEC void operator()(const BoundaryType& boundary,
                                    const NeighborhoodType& uMean,
                                    const NeighborhoodType& uVar,
                                    const NeighborhoodType& vMean,
                                    const NeighborhoodType& vVar,
                                    OutType& divergenceMean,
                                    OutType& divergenceVariance) const
    {
        OutType div_u_mean, var_squared_u;
        OutType div_v_mean, var_squared_v;

        if (boundary.MinNeighborIndices(1)[0] == 0)
        {
            div_u_mean = uMean.Get(1, 0, 0) - uMean.Get(0, 0, 0);
            var_squared_u = uVar.Get(1, 0, 0) + uVar.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[0] == 0)
        {
            div_u_mean = uMean.Get(0, 0, 0) - uMean.Get(-1, 0, 0);
            var_squared_u = uVar.Get(0, 0, 0) + uVar.Get(-1, 0, 0);
        }
        else
        {
            div_u_mean = (uMean.Get(1, 0, 0) - uMean.Get(-1, 0, 0)) / 2.0f;
            var_squared_u = (uVar.Get(1, 0, 0) + uVar.Get(-1, 0, 0)) / 4.0f;
        }
        
        if (boundary.MinNeighborIndices(1)[1] == 0)
        {
            div_v_mean = vMean.Get(0, 1, 0) - vMean.Get(0, 0, 0);
            var_squared_v = vVar.Get(0, 1, 0) + vVar.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[1] == 0) {
            div_v_mean = vMean.Get(0, 0, 0) - vMean.Get(0, -1, 0);
            var_squared_v = vVar.Get(0, 0, 0) + vVar.Get(0, -1, 0);
        }
        else {
            div_v_mean = (vMean.Get(0, 1, 0) - vMean.Get(0, -1, 0)) / 2.0f;
            var_squared_v = (vVar.Get(0, 1, 0) + vVar.Get(0, -1, 0)) / 4.0f;;
        }

        divergenceMean = div_u_mean + div_v_mean;
        divergenceVariance = viskores::Sqrt(var_squared_u + var_squared_v);
    }
};

// Worklet: computes the crossing probability of the divergence field at a given isovalue given the mean and variance of the divergence at each point in the cell.
/* struct CrossingProbability : viskores::worklet::WorkletVisitCellsWithPoints
{
    using ControlSignature = void(CellSetIn,
                                    FieldInPoint div_mean,
                                    FieldInPoint div_variance,
                                    FieldOutCell crossing_prob);
    using ExecutionSignature = void(_2, _3, _4);
    using InputDomain = _1;

    VISKORES_CONT
    CrossingProbability(viskores::FloatDefault isovalue) : Isovalue(isovalue) {}

    template <typename DivMeanVec, typename DivVarianceVec>
    VISKORES_EXEC void operator()(const DivMeanVec& div_means,
                                  const DivVarianceVec& div_variances,
                                  viskores::FloatDefault& crossing_prob) const
    {
        // 4 corners of the square
        constexpr viskores::IdComponent num_corners = 4;

        viskores::Vec<viskores::FloatDefault, num_corners> p_pos;
        viskores::Vec<viskores::FloatDefault, num_corners> p_neg;

        // part of CDF formula
        const viskores::FloatDefault inv_sqrt_2 = viskores::RSqrt(2.0f);

        for (viskores::IdComponent i = 0; i < num_corners; ++i) {
            viskores::FloatDefault mu = divMeans[i];
            viskores::FloatDefault sigma = divVariances[i];
            viskores::FloatDefault p = 0.5 * (1.0 + viskores::Erf((this->Isovalue - mu) * inv_sqrt_2 / sigma));
            p_pos[i] = p;
            p_neg[i] = 1.0 - p;
        }

        viskores::FloatDefault product_p_pos = 1.0f;
        viskores::FloatDefault product_p_neg = 1.0f;
        for (viskores::IdComponent i = 0; i < num_corners; ++i) {
            product_p_pos *= p_pos[i];
            product_p_neg *= p_neg[i];
        }

        crossing_prob = 1.0f - product_p_pos - product_p_neg;
    }
}; */

// Main Function: initializes Viskores, reads a .vtk dataset, computes divergence statistics and displays crossing probabilities.
int main(int argc, char* argv[])
{
    viskores::cont::Initialize(argc, argv);

    viskores::io::VTKDataSetReader reader("data/uncertainVectorField.vtk");
    viskores::cont::DataSet ds = reader.ReadDataSet();

    viskores::cont::Invoker invoker;
    
    using ArrayType = viskores::cont::ArrayHandle<viskores::FloatDefault>;
    ArrayType meanX_handle, varX_handle, meanY_handle, varY_handle;
    ds.GetField("meanX").GetData().AsArrayHandle(meanX_handle);
    ds.GetField("varX").GetData().AsArrayHandle(meanY_handle);
    ds.GetField("meanY").GetData().AsArrayHandle(varX_handle);
    ds.GetField("varY").GetData().AsArrayHandle(varY_handle);
    ArrayType div_mean_handle, div_variance_handle;

    invoker(ComputeDivergence{},
            ds.GetCellSet(),
            meanX_handle,
            varX_handle,
            meanY_handle,
            varY_handle,
            div_mean_handle,
            div_variance_handle
    );
    ds.AddPointField("divergence_mean", div_mean_handle);
    ds.AddPointField("divergence_variance", div_variance_handle);

    viskores::io::VTKDataSetWriter writer("out_vf_divergence.vtk");
    writer.WriteDataSet(ds);

    return 0;
}
