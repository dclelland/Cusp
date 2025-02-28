//
//  LCT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 27/02/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Double {
    
    public func lct1D(matrix: ComplexMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        precondition(matrix.shape == .square(length: 2))
        return ComplexMatrix(real: self).lct1D(matrix: matrix, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func lct1D(matrix: ComplexMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        precondition(matrix.shape == .square(length: 2))
        
        let a = matrix[0, 0].real
        let b = matrix[0, 1].real
        let c = matrix[1, 0].real
        let d = matrix[1, 1].real
        
        let determinant = a * d - b * c
        precondition(abs(determinant - 1.0) < 1e-10, "LCT matrix must be symplectic (ad - bc = 1)")
        
        let inputChirp = ComplexMatrix.lct1DInputChirp(shape: shape, a: a, b: b, c: c, d: d)
        let outputChirp = ComplexMatrix.lct1DOutputChirp(shape: shape, a: a, b: b, c: c, d: d)
        
        let multiplied = self * inputChirp
        let transformed = multiplied.fft1D(setup: setup).fftShifted()
        let result = transformed * outputChirp
        
        return result / Scalar.sqrt(Scalar(shape.count))
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    fileprivate static func lct1DInputChirp(shape: Shape, a: Scalar, b: Scalar, c: Scalar, d: Scalar) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: .row(length: shape.count))
        let factor = (a / (2.0 * b)) * (2.0 * .pi) / Scalar(shape.count)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
    fileprivate static func lct1DOutputChirp(shape: Shape, a: Scalar, b: Scalar, c: Scalar, d: Scalar) -> ComplexMatrix {
        let ramp = Matrix.centeredXRamp(shape: .row(length: shape.count))
        let factor = (d / (2.0 * b)) * (2.0 * .pi) / Scalar(shape.count)
        let phase = ramp.square() * factor
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}
