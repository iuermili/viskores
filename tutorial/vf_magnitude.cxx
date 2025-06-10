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
#include <iostream>
#include <string>

// Worklet: combines two scalar fields into a Vec2f field
struct CombineToVec2: viskores::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn meanX, FieldIn meanY, FieldOut outputVectors);

    VISKORES_EXEC void operator()(viskores::FloatDefault x,
                                  viskores::FloatDefault y,
                                  viskores::Vec2f& outVec) const
    {
        // Combine x and y 
        outVec = viskores::Vec2f(static_cast<float>(x), static_cast<float>(y));
    }
};

// Worklet: copmutes magnitude of a Vec2f field
struct ComputeMagnitude: viskores::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn inputVec, FieldOut outputMagnitude);

    VISKORES_EXEC void operator()(const viskores::Vec2f& inVec, viskores::FloatDefault& outMagnitude) const
    {
        // Compute magnitude
        outMagnitude = viskores::Magnitude(inVec);
    }
};

int main(int argc, char** argv)
{
    viskores::cont::Initialize(argc, argv);

    // Read input file
    viskores::io::VTKDataSetReader reader("data/uncertainVectorField.vtk");
    viskores::cont::DataSet ds = reader.ReadDataSet();

    // Retrieve meanX and meanY as ArrayHandles
    auto meanXArray = ds.GetField("meanX").GetData().AsArrayHandle<viskores::cont::ArrayHandle<viskores::FloatDefault>>();
    auto meanYArray = ds.GetField("meanY").GetData().AsArrayHandle<viskores::cont::ArrayHandle<viskores::FloatDefault>>();
    
    // Combine meanX and meanY into a Vec2f vector field 
    viskores::cont::ArrayHandle<viskores::Vec2f> vectorField;
    viskores::cont::Invoker invoker;
    invoker(CombineToVec2{}, meanXArray, meanYArray, vectorField);

    // Add new vector field to dataset
    viskores::cont::Field vecField("meanXY", ds.GetField("meanX").GetAssociation(), vectorField);
    ds.AddField(vecField);

    // Compute vector magnitude
    viskores::cont::ArrayHandle<viskores::FloatDefault> magnitudeValues;
    invoker(ComputeMagnitude{}, vectorField, magnitudeValues);

    // Add magnitude field to dataset
    viskores::cont::Field magnitudeField("magnitude", ds.GetField("meanX").GetAssociation(), magnitudeValues);
    ds.AddField(magnitudeField);

    // Write output to file
    viskores::io::VTKDataSetWriter writer("out_vf_magnitude.vtk");
    writer.WriteDataSet(ds);

    return 0;
}
