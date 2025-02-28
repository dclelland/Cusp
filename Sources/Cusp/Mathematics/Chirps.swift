//
//  File.swift
//  Cusp
//
//  Created by June Russell on 01/03/2025.
//

import Foundation
import Plinth

extension ComplexMatrix where Scalar == Double {
    
    public static func frftInputChirp(shape: Shape, order: Scalar) -> ComplexMatrix {
        let alpha = order * .pi / 2.0
        
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi / tan(alpha) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public static func frftOutputChirp(shape: Shape, order: Scalar) -> ComplexMatrix {
        let alpha = order * .pi / 2.0
        
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi / sin(alpha) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public static func lctInputChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.a / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public static func lctOutputChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.d / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}
