//
//  Ramps.swift
//  Cusp
//
//  Created by June Russell on 28/02/2025.
//

import Foundation
import Plinth

extension Matrix where Scalar == Double {
    
    // Create an x-axis ramp matrix for fractional transforms
    internal static func frftXRamp(shape: Shape) -> Matrix {
        let width = shape.columns / 2
        let range = Scalar(-width)...Scalar(width - 1)
        return Matrix.xRamp(shape: shape, range: range)
    }
    
    // Create a y-axis ramp matrix for fractional transforms
    internal static func frftYRamp(shape: Shape) -> Matrix {
        let height = shape.rows / 2
        let range = Scalar(-height)...Scalar(height - 1)
        return Matrix.yRamp(shape: shape, range: range)
    }
    
}