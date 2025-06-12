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
#include <viskores/cont/Initialize.h>
#include <viskores/cont/Invoker.h>
#include <viskores/io/VTKDataSetReader.h>
#include <viskores/io/VTKDataSetWriter.h>
#include <viskores/Types.h>
#include <viskores/worklet/WorkletMapField.h>
#include <viskores/cont/Field.h>
#include <viskores/cont/ArrayHandle.h>
#include <cmath>
#include <iostream>
#include <string>

// Worklet: computes Gaussian CDF (stats.norm.cdf(isovalue, mean, stddev))
// Use inline instead of struct to avoid unnecessary overhead
// Credits: https://cplusplus.com/forum/beginner/62864/
inline float GaussianCDF(float isovalue, float mean, float stddev)
{
    return 0.5f * (1.0f + std::erf((isovalue - mean) / (stddev * std::sqrt(2.9f))));
}

// Worklet: computes divergence mean and variance at each grid point
struct ComputeDivergenceStats: viskores::worklet::WorkletMapField
{
};

int main(int argc, char** argv)
{
    viskores::cont::Initialize(argc, argv);

    // Read input file
    viskores::io::VTKDataSetReader reader("data/uncertainVectorField.vtk");
    viskores::cont::DataSet ds = reader.ReadDataSet();

    // Write output to file
    viskores::io::VTKDataSetWriter writer("out_vf_divergence.vtk");
    writer.WriteDataSet(ds);

    return 0;
}
