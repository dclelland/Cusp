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
        return ComplexMatrix(real: self).frft1D(a: a, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func frft1DSlop(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let epsilon = 1e-6
        let aModulo = a.mod(4)
        if abs(aModulo) < epsilon {
            return self
        }
        if abs(aModulo - 1) < epsilon {
            return fft1D(setup: setup)
        }
        if abs(aModulo - 2) < epsilon {
            return self.reversed()
        }
        if abs(aModulo - 3) < epsilon {
            return ifft1D(setup: setup)
        }
        
        // Compute the rotation angle theta (in radians)
        let N = Scalar(shape.count)  // signal length
        let theta = a * .pi / 2.0
        
        // Make sure sin(theta) is not too small.
        guard abs(Scalar.sin(theta)) >= epsilon else {
            fatalError("Angle leads to a degenerate transform.")
        }
        
        // Create a centered index ramp: values from -N/2 to N/2 - 1.
        // (Assume .centeredIndexRamp(length:) is provided.)
        let n: Matrix = .fftXRamp(shape: .row(length: shape.count)).fftShifted()
        
        // Pre-Chirp: Compute C1(n) = exp(-i π (n²/N) cot(θ))
        let cotTheta = cos(theta) / sin(theta)
        let nSquared = n.square()  // elementwise square of the ramp
        let phasePre = (.pi * (nSquared / N)) * cotTheta
        // We want exp(-i * phasePre): cos(phasePre) - i*sin(phasePre)
        let preChirp = ComplexMatrix(real: phasePre.cos(), imaginary: -phasePre.sin())
        
        // Multiply input elementwise by the pre-chirp:
        let y = self * preChirp
        
        // Chirp Convolution Kernel: h(n) = exp(i π (n²/N) csc(θ))
        let cscTheta = 1.0 / sin(theta)
        let phaseKernel = (.pi * (nSquared / N)) * cscTheta
        let kernel = ComplexMatrix(real: phaseKernel.cos(), imaginary: phaseKernel.sin())
        
        // Convolution via FFT:
        // Compute the FFTs (assumed unitary) of y and the kernel.
        let Y = y.fft1D(setup: setup)
        let H = kernel.fft1D(setup: setup)
        // Elementwise multiplication in the Fourier domain:
        let convFFT = Y * H
        // Inverse FFT to perform the convolution:
        let z = convFFT.ifft1D(setup: setup)
        
        // Post-Chirp: Multiply by the same pre-chirp factor (elementwise)
        let zChirped = z * preChirp
        
        // Normalization factor:
        let sgn = sin(theta) >= 0 ? 1.0 : -1.0
        let normLength = 1.0 / Scalar.sqrt(N * abs(sin(theta)))
        let normPhase = -(theta / 2.0) - (.pi / 4.0) * sgn
        let A = Complex(length: normLength, phase: normPhase)
        
        // Final output:
        let X = A * zChirped
        return X
        
    }
    
    public func frft1D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let epsilon = 1e-6
        let aModulo = a.mod(4)
        if abs(aModulo) < epsilon {
            return self
        }
        if abs(aModulo - 1) < epsilon {
            return fft1D(setup: setup)
        }
        if abs(aModulo - 2) < epsilon {
            return self.reversed()
        }
        if abs(aModulo - 3) < epsilon {
            return ifft1D(setup: setup)
        }
        
        let theta = a * .pi / 2.0
        let cotTheta = 1.0 / Scalar.tan(theta)
        let cscTheta = 1.0 / Scalar.sin(theta)
        
        let xRamp: Matrix = .fftXRamp(shape: .square(length: shape.count))
        let yRamp: Matrix = .fftYRamp(shape: .square(length: shape.count))
        
        let quadratic: Matrix = (xRamp.square() + yRamp.square()) * cotTheta
        let cross = (2.0 * xRamp * yRamp) * cscTheta
        let kernelExponent = (Scalar.pi / Scalar(shape.count)) * (quadratic - cross)
        let kernel: ComplexMatrix = .init(real: kernelExponent.cos(), imaginary: kernelExponent.sin())
        
        let sgn = Scalar.sin(theta) >= 0 ? 1.0 : -1.0
        let aThetaLength: Scalar = 1.0 / Scalar.sqrt(Scalar(shape.count) * abs(Scalar.sin(theta)))
        let aThetaExponent: Scalar = -.pi / 4.0 * sgn - theta / 2.0
        let aTheta = Complex(length: aThetaLength, phase: aThetaExponent)
        
        let transformed = (kernel <*> self.asColumn()) * aTheta
        return transformed.asRow()
    }
    
    
    public func frft1DChirp(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let epsilon = 1e-8
        let aModulo = a.truncatingRemainder(dividingBy: 4)
        if abs(aModulo) < epsilon {
            return self
        }
        if abs(aModulo - 1) < epsilon {
            return fft1D(setup: setup)
        }
        if abs(aModulo - 2) < epsilon {
            return self.reversed()
        }
        if abs(aModulo - 3) < epsilon {
            return ifft1D(setup: setup)
        }

        let alpha = a * .pi / 2.0
        let cotAlpha = 1.0 / Scalar.tan(alpha)
        let sinAlpha = Scalar.sin(alpha)

        // Create an index vector centered around zero.
        let n: Matrix = .fftXRamp(shape: shape)//.fftShifted()

        // Pre-chirp multiplication: multiply input by a quadratic phase factor.
        let preChirpPhase = Scalar.pi * cotAlpha * n.square() / Scalar(shape.count)
        let preChirp = ComplexMatrix(real: preChirpPhase.cos(), imaginary: preChirpPhase.sin())
        let xPre = self * preChirp

        // Compute FFT of the pre-chirped signal.
        let X = xPre.fft1D(setup: setup)//.fftShifted()

        // Multiply in the Fourier domain by the chirp kernel.
        let kernelPhase = Scalar.pi * n.square() / (Scalar(shape.count) * sinAlpha)
        let kernel = ComplexMatrix(real: kernelPhase.cos(), imaginary: kernelPhase.sin())
        let XKernel = X * kernel

        // Inverse FFT to complete the convolution.
        let xIfft = XKernel.ifft1D(setup: setup)//.fftShifted()

        // Post-chirp multiplication.
        let postChirpPhase = Scalar.pi * cotAlpha * n.square() / Scalar(shape.count)
        let postChirp = ComplexMatrix(real: postChirpPhase.cos(), imaginary: postChirpPhase.sin())
        let y = xIfft * postChirp

        // Overall scaling factor.
        let phaseA = -(.pi / 4.0 - alpha / 2.0)
        let A = Complex(Scalar.cos(phaseA), Scalar.sin(phaseA)) / Complex(Scalar.sqrt(abs(sinAlpha)))

        return y * A
    }
    
}

private extension FloatingPoint {
    
    fileprivate func mod(_ modulus: Self) -> Self {
        let remainder = self.truncatingRemainder(dividingBy: modulus)
        return remainder < 0 ? remainder + modulus : remainder
    }
    
}
