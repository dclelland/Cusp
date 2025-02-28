//
//  File.swift
//  Cusp
//
//  Created by June Russell on 01/03/2025.
//

import Foundation
import Plinth

extension ComplexMatrix where Scalar == Double {
    
    public static func frftInputChirp(shape: Shape, a: Scalar) -> ComplexMatrix {
        let alpha = a * .pi / 2.0
        
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi / Scalar.tan(alpha) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public static func frftOutputChirp(shape: Shape, a: Scalar) -> ComplexMatrix {
        let alpha = a * .pi / 2.0
        
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi / Scalar.sin(alpha) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public static func lctInputChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = (matrix.a / (2.0 * matrix.b)) * (2.0 * .pi) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public static func lctOutputChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = (matrix.d / (2.0 * matrix.b)) * (2.0 * .pi) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}
