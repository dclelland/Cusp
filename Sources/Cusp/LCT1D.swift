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
        
//        // Handle special cases - identity transform
//        if matrix.isApproximately(LCTMatrix<Scalar>.identityTransform) {
//            return self
//        }
        
        // Extract matrix components
        let a = matrix[0, 0].real
        let b = matrix[0, 1].real
        let c = matrix[1, 0].real
        let d = matrix[1, 1].real
        
        // Check if the matrix is symplectic (ad - bc = 1)
        let determinant = a * d - b * c
        precondition(abs(determinant - 1.0) < 1e-10, "LCT matrix must be symplectic (ad - bc = 1)")
        
        // Special cases
        
//        // Pure scaling (b = 0)
//        if abs(b) < 1e-10 {
//            return handlePureScaling(a: a, c: c, d: d)
//        }
        
        // General case - chirp method
        let inputChirp = ComplexMatrix.lct1DInputChirp(shape: shape, a: a, b: b, c: c, d: d)
        let outputChirp = ComplexMatrix.lct1DOutputChirp(shape: shape, a: a, b: b, c: c, d: d)
        
        let multiplied = self * inputChirp
        let transformed = multiplied.fft1D(setup: setup).fftShifted()
        let result = transformed * outputChirp
        
        // Scale factor based on the b parameter
        return result * (abs(b) != 1.0 ? Complex(1.0 / Scalar.sqrt(abs(b)), 0) : Complex(1, 0))
    }
    
//    // Handle the case when b = 0 (pure scaling)
//    private func handlePureScaling(a: Scalar, c: Scalar, d: Scalar) -> ComplexMatrix<Scalar> {
//        // For b = 0, the transform is a scaling and chirp multiplication
//        // The scaling factor is 1/sqrt(|a|)
//        let scaleFactor = 1.0 / Scalar.sqrt(abs(a))
//        
//        // Calculate the scaling phase factors if needed
//        let chirpFactor = ComplexMatrix(shape: shape) { i, j in
//            let x = Matrix.frftXRamp(shape: .row(length: shape.count))[0, i]
//            let phase = (c / (2 * a)) * x * x
//            return Complex(cos(phase), sin(phase))
//        }
//        
//        // Scale the signal and apply the chirp
//        // Note: actual scaling would require interpolation which is beyond the scope here
//        return self * chirpFactor * Complex(scaleFactor, 0)
//    }
}

extension ComplexMatrix where Scalar == Double {
    
    fileprivate static func lct1DInputChirp(shape: Shape, a: Scalar, b: Scalar, c: Scalar, d: Scalar) -> ComplexMatrix {
        let xRamp = Matrix.frftXRamp(shape: .row(length: shape.count))
        let chirpFactor = (.pi / b) * (a / 2.0)
        let chirpPhase = xRamp.square() * chirpFactor
        
        return ComplexMatrix(real: chirpPhase.cos(), imaginary: chirpPhase.sin())
    }
    
    fileprivate static func lct1DOutputChirp(shape: Shape, a: Scalar, b: Scalar, c: Scalar, d: Scalar) -> ComplexMatrix {
        let xRamp = Matrix.frftXRamp(shape: .row(length: shape.count))
        let chirpFactor = (.pi / b) * (d / 2.0)
        let chirpPhase = xRamp.square() * chirpFactor
        
        return ComplexMatrix(real: chirpPhase.cos(), imaginary: chirpPhase.sin())
    }
    
}

//// Add a helper method to check if two complex matrices are approximately equal
//extension ComplexMatrix where Scalar == Double {
//    
//    fileprivate func isApproximately(_ other: ComplexMatrix<Scalar>, tolerance: Scalar = 1e-10) -> Bool {
//        guard shape == other.shape else { return false }
//        
//        for i in 0..<shape.rows {
//            for j in 0..<shape.columns {
//                let diff = self[i, j] - other[i, j]
//                if diff.real * diff.real + diff.imaginary * diff.imaginary > tolerance * tolerance {
//                    return false
//                }
//            }
//        }
//        
//        return true
//    }
//    
//}
