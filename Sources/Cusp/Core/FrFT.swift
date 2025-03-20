//
//  FrFT.swift
//  Cusp
//
//  Created by Daniel Clelland on 11/03/2025.
//

import Foundation
import Numerics
import Plinth

public struct FrFT<Scalar> where Scalar: Real {
    
    public var length: Int
    public var scale: Int
    public var order: Scalar
    
    public init(length: Int, scale: Int = 1, order: Scalar) {
        self.length = length
        self.order = order
        self.scale = scale
    }
    
}

extension FrFT where Scalar == Float {
    
    public func preChirp() -> ComplexMatrix<Scalar> {
        let ramp = Matrix<Scalar>.centeredXRamp(shape: .row(length: length))
        let alpha = order * .pi / 2.0
        let factor = -.pi * tan(alpha / 2.0) / Scalar(length)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public func postChirp() -> ComplexMatrix<Scalar> {
        let ramp = Matrix<Scalar>.centeredXRamp(shape: .row(length: length * scale))
        let alpha = order * .pi / 2.0
        let factor = .pi / sin(alpha) / Scalar(length)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}

extension FrFT where Scalar == Double {
    
    public func preChirp() -> ComplexMatrix<Scalar> {
        let ramp = Matrix<Scalar>.centeredXRamp(shape: .row(length: length))
        let alpha = order * .pi / 2.0
        let factor = -.pi * tan(alpha / 2.0) / Scalar(length)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    public func postChirp() -> ComplexMatrix<Scalar> {
        let ramp = Matrix<Scalar>.centeredXRamp(shape: .row(length: length * scale))
        let alpha = order * .pi / 2.0
        let factor = .pi / sin(alpha) / Scalar(length)
        let phase = ramp.square() * factor
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}
