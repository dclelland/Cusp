//
//  Chirps.swift
//  Cusp
//
//  Created by Daniel Clelland on 01/03/2025.
//

import Foundation
import Plinth

extension ComplexMatrix where Scalar == Float {
    
    internal static func frftPreChirp(shape: Shape, order: Scalar) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let alpha = order * .pi / 2.0
        let factor = -.pi * tan(alpha / 2.0) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    internal static func frftPostChirp(shape: Shape, order: Scalar, scale: Int = 1) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: .init(rows: shape.rows, columns: shape.columns * scale))
        let alpha = order * .pi / 2.0
        let factor = .pi / sin(alpha) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    internal static func frftPreChirp(shape: Shape, order: Scalar) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: shape)
        let alpha = order * .pi / 2.0
        let factor = -.pi * tan(alpha / 2.0) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    internal static func frftPostChirp(shape: Shape, order: Scalar, scale: Int = 1) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: .init(rows: shape.rows, columns: shape.columns * scale))
        let alpha = order * .pi / 2.0
        let factor = .pi / sin(alpha) / Scalar(shape.columns)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}
