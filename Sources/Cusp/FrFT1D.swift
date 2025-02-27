//
//  FrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 27/02/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Double {
    
    public func frft1D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let epsilon = 1e-8
        let aModulo = a.truncatingRemainder(dividingBy: 4)
        if abs(aModulo) < epsilon {
            return ComplexMatrix(real: self)
        }
        if abs(aModulo - 1) < epsilon {
            return fft1D(setup: setup)
        }
        if abs(aModulo - 2) < epsilon {
            return ComplexMatrix(real: self.reversed())
        }
        if abs(aModulo - 3) < epsilon {
            return ifft1D(setup: setup)
        }
        
        // General case: noninteger order.
        let theta = a * .pi / 2
        let cot = 1 / Scalar.tan(theta)
        let csc = 1 / Scalar.sin(theta)
        
        // Create index vector n = [0, 1, ..., N-1]
//        let n = (0..<shape.count).map { Scalar($0) - Scalar(shape.count - 1) / 2.0 }
        let n: Matrix<Scalar> = .fftXRamp(shape: shape)//.fftShifted()
        
        // Build the kernel matrix.
        // Each element is computed as:
        //   exp(-iπ*(n_i^2 + n_j^2)*cot/N + i2π*n_i*n_j*csc/N)
        // We assume that a 2D array supports pointwise arithmetic and that the `<*>` operator
        // performs matrix multiplication with a vector.
        var kernel = ComplexMatrix<Scalar>(shape: .square(length: shape.count)) { i, j in
            let phase = -Scalar.pi * (n[i] * n[i] + n[j] * n[j]) * cot / Scalar(shape.count) + 2.0 * .pi * n[i] * n[j] * csc / Scalar(shape.count)
            return Complex<Scalar>(length: 1.0, phase: phase)
        }
        
        // Compute the normalization constant:
        //   A = exp(-i*(π/4)*sgn(sinθ) - iθ/2) / sqrt(N * |sinθ|)
        let sgn: Scalar = Scalar.sin(theta)
        let length: Scalar = 1.0 / Scalar.sqrt(Scalar(shape.count) * abs(Scalar.sin(theta)))
        let phase: Scalar = -.pi / 4.0 * sgn - theta / 2.0
        let A = Complex<Scalar>(length: length, phase: phase)
        
        // Compute the transform by matrix–vector multiplication.
        // Here we assume that the operator `<*>` multiplies our 2D array (matrix) by the vector.
        let transformed = A * (kernel <*> self.asColumn())
        return transformed.asRow()
    }
    
//    public func frft1D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
//        let epsilon = 1e-8
//        let aModulo = a.truncatingRemainder(dividingBy: 4)
//        if abs(aModulo) < epsilon {
//            return ComplexMatrix(real: self)
//        }
//        if abs(aModulo - 1) < epsilon {
//            return fft1D(setup: setup)
//        }
//        if abs(aModulo - 2) < epsilon {
//            return ComplexMatrix(real: self.reversed())
//        }
//        if abs(aModulo - 3) < epsilon {
//            return ifft1D(setup: setup)
//        }
//        
//        let alpha = a * .pi / 2.0
//        let cotAlpha = 1.0 / Scalar.tan(alpha)
//        let sinAlpha = Scalar.sin(alpha)
//          
//        // Create an index vector centered around zero.
//        let n: Matrix<Scalar> = .fftXRamp(shape: shape)//.fftShifted()
//
//        // Pre-chirp multiplication: multiply input by a quadratic phase factor.
//        let preChirpPhase = Scalar.pi * cotAlpha * n.square() / Scalar(shape.count)
//        let preChirp = ComplexMatrix<Scalar>(real: preChirpPhase.cos(), imaginary: preChirpPhase.sin())
//        let xPre = self * preChirp
//        
//        // Compute FFT of the pre-chirped signal.
//        let X = xPre.fft1D(setup: setup)//.fftShifted()
//        
//        // Multiply in the Fourier domain by the chirp kernel.
//        let kernelPhase = Scalar.pi * n.square() / (Scalar(shape.count) * sinAlpha)
//        let kernel = ComplexMatrix<Scalar>(real: kernelPhase.cos(), imaginary: kernelPhase.sin())
//        let XKernel = X * kernel
//        
//        // Inverse FFT to complete the convolution.
//        let xIfft = XKernel.ifft1D(setup: setup)//.fftShifted()
//        
//        // Post-chirp multiplication.
//        let postChirpPhase = Scalar.pi * cotAlpha * n.square() / Scalar(shape.count)
//        let postChirp = ComplexMatrix<Scalar>(real: postChirpPhase.cos(), imaginary: postChirpPhase.sin())
//        let y = xIfft * postChirp
//        
//        // Overall scaling factor.
//        let phaseA = -(.pi / 4.0 - alpha / 2.0)
//        let A = Complex<Scalar>(Scalar.cos(phaseA), Scalar.sin(phaseA)) / Complex(Scalar.sqrt(abs(sinAlpha)))
//        
//        return y * A
//    }
    
}
