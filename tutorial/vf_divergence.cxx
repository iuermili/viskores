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
#include <viskores/Math.h>
#include <viskores/Types.h>
#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/DataSet.h>
#include <viskores/cont/Field.h>
#include <viskores/cont/Initialize.h>
#include <viskores/cont/Invoker.h>
#include <viskores/worklet/WorkletPointNeighborhood.h>
#include <viskores/worklet/WorkletMapTopology.h>
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

        if (boundary.MinNeighborIndices(1)[1] == 0)
        {
            div_u_mean = uMean.Get(0, 1, 0) - uMean.Get(0, 0, 0);
            var_squared_u = uVar.Get(0, 1, 0) + uVar.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[1] == 0) {
            div_u_mean = uMean.Get(0, 0, 0) - uMean.Get(0, -1, 0);
            var_squared_u = uVar.Get(0, 0, 0) + uVar.Get(0, -1, 0);
        }
        else {
            div_u_mean = (uMean.Get(0, 1, 0) - uMean.Get(0, -1, 0)) / 2.0;
            var_squared_u = uVar.Get(0, 1, 0) + uVar.Get(0, -1, 0);
        }

        if (boundary.MinNeighborIndices(1)[0] == 0)
        {
            div_v_mean = vMean.Get(1, 0, 0) - vMean.Get(0, 0, 0);
            var_squared_v = vVar.Get(1, 0, 0) + vVar.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[0] == 0)
        {
            div_v_mean = vMean.Get(0, 0, 0) - vMean.Get(-1, 0, 0);
            var_squared_v = vVar.Get(0, 0, 0) + vVar.Get(-1, 0, 0);
        }
        else
        {
            div_v_mean = (vMean.Get(1, 0, 0) - vMean.Get(-1, 0, 0)) / 2.0;
            var_squared_v = vVar.Get(1, 0, 0) + vVar.Get(-1, 0, 0);
        }

        divergenceMean = div_u_mean + div_v_mean;
        OutType total_variance = var_squared_u + var_squared_v;
        divergenceVariance = viskores::Sqrt(viskores::Max(total_variance, static_cast<OutType>(0.0)));
    }
};

// Helper Function: computes Gaussian CDF for a given value x, mean mu, and standard deviation sigma.
template <typename T>
VISKORES_EXEC T GaussianCDF(T x, T mu, T sigma)
{
    return static_cast<T>(0.5) * (static_cast<T>(1.0) + erf((x - mu) / (sigma * viskores::Sqrt(2.0))));
}

// Worklet: computes the crossing probability a 2D vector field at each point based on the mean and variance of the divergence.
struct CrossingProbability : public viskores::worklet::WorkletVisitCellsWithPoints
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInPoint div_mean,
                                    FieldInPoint div_variance,
                                    FieldOutCell crossing_prob);
    using ExecutionSignature = _4(_2, _3);
    using InputDomain = _1;

    viskores::Float64 Isovalue;

    template <typename MeanVecType, typename VarianceVecType>
    VISKORES_EXEC viskores::Float64 operator()(const MeanVecType& divMeanAtCorners,
                                    const VarianceVecType& divVarianceAtCorners) const
    {
        viskores::Float64 prob_pos = 1.0;
        viskores::Float64 prob_neg = 1.0;

        for (viskores::IdComponent i = 0; i < divMeanAtCorners.GetNumberOfComponents(); ++i)
        {
            auto mean = static_cast<viskores::Float64>(divMeanAtCorners[i]);
            auto variance = static_cast<viskores::Float64>(divVarianceAtCorners[i]);
            
            auto cdf_value = viskores::Sqrt(viskores::Max(variance, static_cast<viskores::Float64>(0.0)));

            auto prob_below = GaussianCDF(this->Isovalue, mean, cdf_value);
            auto prob_above = 1.0 - prob_below;

            prob_neg *= prob_below;
            prob_pos *= prob_above;
        }

        return 1.0f - prob_pos - prob_neg;
    }
};

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
    ds.GetField("varX").GetData().AsArrayHandle(varX_handle);
    ds.GetField("meanY").GetData().AsArrayHandle(meanY_handle);
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

    CrossingProbability crossingprobWorklet;
    crossingprobWorklet.Isovalue = -5.05;
    viskores::cont::ArrayHandle<viskores::Float64> crossing_prob_handle;
    invoker(crossingprobWorklet,
                ds.GetCellSet(),
                div_mean_handle,
                div_variance_handle,
                crossing_prob_handle
    );

    ds.AddPointField("divergence_mean", div_mean_handle);
    ds.AddPointField("divergence_variance", div_variance_handle);
    ds.AddCellField("crossing_probability", crossing_prob_handle);

    viskores::io::VTKDataSetWriter writer("out_vf_divergence.vtk");
    writer.WriteDataSet(ds);

    return 0;
}
