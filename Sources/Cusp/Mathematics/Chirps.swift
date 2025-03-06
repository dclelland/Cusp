//
//  File.swift
//  Cusp
//
//  Created by June Russell on 01/03/2025.
//

import Foundation
import Numerics
import Plinth

extension ComplexMatrix where Scalar == Float {
    
    public static func lctPreChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        if matrix.b.isApproximatelyEqual(to: 0.0) {
            return .ones(shape: shape)
        }
        
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.a / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public static func lctPostChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        if matrix.b.isApproximatelyEqual(to: 0.0) {
            return .ones(shape: shape)
        }
        
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.d / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public static func frftPreChirp(shape: Shape, order: Scalar) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let angle = order * .pi / 2.0
        let factor = -.pi * tan(angle / 2.0) / Scalar(shape.length)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public static func frftPostChirp(shape: Shape, order: Scalar) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: .row(length: shape.length * 2))
        let angle = order * .pi / 2.0
        let factor = .pi / sin(angle) / Scalar(shape.length)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public static func lctPreChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        if matrix.b.isApproximatelyEqual(to: 0.0) {
            return .ones(shape: shape)
        }
        
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.a / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public static func lctPostChirp(shape: Shape, matrix: LCTMatrix<Scalar>) -> ComplexMatrix {
        if matrix.b.isApproximatelyEqual(to: 0.0) {
            return .ones(shape: shape)
        }
        
        let ramp = Matrix.centeredXRamp(shape: shape)
        let factor = .pi * matrix.d / matrix.b / Scalar(shape.columns)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}
