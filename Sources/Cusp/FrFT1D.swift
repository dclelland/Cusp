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
        
        let alpha = a * .pi / 2.0
        let cotAlpha = 1.0 / Scalar.tan(alpha)
        let sinAlpha = Scalar.sin(alpha)
          
        // Create an index vector centered around zero.
        let n = (0..<shape.count).map { Scalar($0) - Scalar(shape.count / 2) }

        // Pre-chirp multiplication: multiply input by a quadratic phase factor.
        let preChirp: [Complex<Scalar>] = n.map { nVal in
            let phase = Scalar.pi * cotAlpha * (nVal * nVal) / Scalar(shape.count)
            return Complex<Scalar>(Scalar.cos(phase), Scalar.sin(phase))
        }
        let xPre = ComplexMatrix<Scalar>(shape: shape, elements: zip(self.elements, preChirp).map { Complex($0) * $1 })
        
        // Compute FFT of the pre-chirped signal.
        let X = xPre.fft1D(setup: setup)
        
        // Multiply in the Fourier domain by the chirp kernel.
        let kernel: [Complex<Scalar>] = n.map { nVal in
            let phase = Scalar.pi * (nVal * nVal) / (Scalar(shape.count) * sinAlpha)
            return Complex<Scalar>(Scalar.cos(phase), Scalar.sin(phase))
        }
        let XKernel = ComplexMatrix<Scalar>(shape: shape, elements: zip(X, kernel).map { $0 * $1 })
        
        // Inverse FFT to complete the convolution.
        let xIfft = XKernel.ifft1D(setup: setup)
        
        // Post-chirp multiplication.
        let postChirp: [Complex<Scalar>] = n.map { nVal in
            let phase = Scalar.pi * cotAlpha * (nVal * nVal) / Scalar(shape.count)
            return Complex<Scalar>(Scalar.cos(phase), Scalar.sin(phase))
        }
        let y = ComplexMatrix<Scalar>(shape: shape, elements: zip(xIfft, postChirp).map { $0 * $1 })
        
        // Overall scaling factor.
        let phaseA = -(.pi / 4.0 - alpha / 2.0)
        let A = Complex<Scalar>(Scalar.cos(phaseA), Scalar.sin(phaseA)) / Complex(Scalar.sqrt(abs(sinAlpha)))
        
        return y * A
    }
    
}
