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
        
        let theta = a * .pi / 2.0
        
        let xRamp: Matrix<Scalar> = .fftXRamp(shape: .square(length: shape.count))
        let yRamp: Matrix<Scalar> = .fftYRamp(shape: .square(length: shape.count))
        let quadraticPhase: Matrix<Scalar> = (-Scalar.pi * (xRamp.square() + yRamp.square())) / Scalar.tan(theta)
        let crossCoupling: Matrix<Scalar> = (2.0 * .pi * xRamp * yRamp) / Scalar.sin(theta)
        let phase: Matrix<Scalar> = (quadraticPhase + crossCoupling) / Scalar(shape.count)
        let kernel: ComplexMatrix<Scalar> = .init(real: phase.cos(), imaginary: phase.sin())
        
        let sgn = Scalar.sin(theta) >= 0 ? 1.0 : -1.0
        let normLength: Scalar = 1.0 / Scalar.sqrt(Scalar(shape.count) * abs(Scalar.sin(theta)))
        let normPhase: Scalar = -.pi / 4.0 * sgn - theta / 2.0
        let norm = Complex<Scalar>(length: normLength, phase: normPhase)
        
        let transformed = (kernel <*> self.asColumn()) * norm
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
