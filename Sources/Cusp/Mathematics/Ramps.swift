//
//  Ramps.swift
//  Cusp
//
//  Created by June Russell on 28/02/2025.
//

import Foundation
import Plinth

extension Matrix where Scalar == Double {
    
    public static func centeredXRamp(shape: Shape) -> Matrix {
        let width = shape.columns / 2
        let range = Scalar(-width)...Scalar(width - 1)
        return Matrix.xRamp(shape: shape, range: range)
    }
    
    public static func centeredYRamp(shape: Shape) -> Matrix {
        let height = shape.rows / 2
        let range = Scalar(-height)...Scalar(height - 1)
        return Matrix.yRamp(shape: shape, range: range)
    }
    
}
