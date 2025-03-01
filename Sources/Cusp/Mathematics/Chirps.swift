//
//  File.swift
//  Cusp
//
//  Created by June Russell on 01/03/2025.
//

import Foundation
import Plinth

extension ComplexMatrix where Scalar == Float {
    
    public static func lctPreChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.a / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public static func lctPostChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.d / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public static func lctPreChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.a / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public static func lctPostChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.d / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}
