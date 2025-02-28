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
    
    public func frft1DMatrix(a: Scalar) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft1DMatrix(a: a)
    }
    
    public func frft1D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft1D(a: a, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func frft1DMatrix(a: Scalar) -> ComplexMatrix {
        guard a != 0 else {
            return self
        }
        
        let kernel = ComplexMatrix.frft1DKernel(shape: shape, a: a)
        let transformed = kernel <*> self.asColumn()
        return transformed.asRow() / Scalar.sqrt(Scalar(shape.count))
    }
    
    public func frft1D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        // Prepare parameters for the chirp method
        let phi = a * .pi / 2
        let cotPhi = 1 / Scalar.tan(phi)
        let cscPhi = 1 / Scalar.sin(phi)
        
        let N = Scalar(shape.count)
        let normFactor = Scalar.sqrt(N * Scalar.sin(phi).magnitude)
        
        // Create chirp vectors
        let n = Array(0..<shape.count).map { Scalar($0 - shape.count / 2) }
        let chirp1 = n.map { Complex(Scalar.cos(.pi * cotPhi * $0 * $0 / N), Scalar.sin(.pi * cotPhi * $0 * $0 / N)) }
        let chirp2 = n.map { Complex(Scalar.cos(.pi * cscPhi * $0 * $0 / N), Scalar.sin(.pi * cscPhi * $0 * $0 / N)) }
        
        // Create ComplexMatrix objects from the chirp vectors
        let chirpMatrix1 = ComplexMatrix(shape: .row(length: shape.count), elements: chirp1)
        let chirpMatrix2 = ComplexMatrix(shape: .row(length: shape.count), elements: chirp2)
        
        // Step 1: Multiply input by first chirp
        let multiplied = self * chirpMatrix1
        
        // Step 2: Compute FFT
        let transformed = multiplied.fft1D(setup: setup)
        
        // Step 3: Multiply by second chirp
        let result = transformed * (chirpMatrix2)
        
        // Apply normalization
        return result / Scalar.sqrt(Scalar(shape.count))
    }
    
//    public func frft1D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
//        let N = Scalar(shape.count)
//        let alpha = a * .pi / 2.0
//        
//        // Create a centered index ramp: n = -N/2, ..., N/2 - 1.
//        let n: Matrix = .frftXRamp(shape: .row(length: shape.count)) // .fftShifted()
//        
//        // Step 1: Pre-Chirp Multiplication
//        // Multiply by: exp(-i π n² (cot(alpha))/N)
//        let cotAlpha: Scalar = 1.0 / Scalar.tan(alpha)
//        let prePhase = (.pi / N) * (n.square() * cotAlpha)
//        let preChirp = ComplexMatrix(real: prePhase.cos(), imaginary: -prePhase.sin())
//        
//        // Step 2: Convolution via FFT with the Chirp Kernel
//        // Define the chirp kernel: h(n) = exp(i π n² (csc(alpha))/N)
//        let cscAlpha: Scalar = 1.0 / Scalar.sin(alpha)
//        let hPhase = (.pi / N) * (n.square() * cscAlpha)
//        let h = ComplexMatrix(real: hPhase.cos(), imaginary: hPhase.sin())
//        
//        // Convolve y with h via FFT:
//        // Note: These FFT routines should be unitary.
//        let y = self * preChirp
//        let Y = y.fft1D(setup: setup)
//        let H = h.fft1D(setup: setup)
//        let convFFT = Y * H
//        let z = convFFT.ifft1D(setup: setup)
//        
//        // Step 3: Post-Chirp Multiplication
//        // Multiply the result by the same chirp factor as in the pre-step.
//        return z * preChirp / Scalar.sqrt(N)
//    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    fileprivate static func frft1DKernel(shape: Shape, a: Scalar) -> ComplexMatrix {
        let alpha = a * .pi / 2.0
        let cotAlpha = 1.0 / Scalar.tan(alpha)
        let cscAlpha = 1.0 / Scalar.sin(alpha)
        
        let xRamp: Matrix = .frftXRamp(shape: .square(length: shape.count))
        let yRamp: Matrix = .frftYRamp(shape: .square(length: shape.count))
        
        let quadratic: Matrix = (xRamp.square() + yRamp.square()) * cotAlpha
        let cross = (2.0 * xRamp * yRamp) * cscAlpha
        let phase = (.pi / Scalar(shape.count)) * (quadratic - cross)
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}

extension Matrix where Scalar == Double {
        
    fileprivate static func frftXRamp(shape: Shape) -> Matrix {
        let width = shape.columns / 2
        let range = Scalar(-width)...Scalar(width - 1)
        return Matrix.xRamp(shape: shape, range: range)
    }
    
    fileprivate static func frftYRamp(shape: Shape) -> Matrix {
        let height = shape.rows / 2
        let range = Scalar(-height)...Scalar(height - 1)
        return Matrix.yRamp(shape: shape, range: range)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
//    public func frft1DSlop(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
//        let epsilon = 1e-6
//        let aModulo = a.mod(4)
//        if abs(aModulo) < epsilon {
//            return self
//        }
//        if abs(aModulo - 1) < epsilon {
//            return fft1D(setup: setup)
//        }
//        if abs(aModulo - 2) < epsilon {
//            return self.reversed()
//        }
//        if abs(aModulo - 3) < epsilon {
//            return ifft1D(setup: setup)
//        }
//        
//        // Compute the rotation angle theta (in radians)
//        let N = Scalar(shape.count)  // signal length
//        let theta = a * .pi / 2.0
//        
//        // Make sure sin(theta) is not too small.
//        guard abs(Scalar.sin(theta)) >= epsilon else {
//            fatalError("Angle leads to a degenerate transform.")
//        }
//        
//        // Create a centered index ramp: values from -N/2 to N/2 - 1.
//        // (Assume .centeredIndexRamp(length:) is provided.)
//        let n: Matrix = .fftXRamp(shape: .row(length: shape.count)).fftShifted()
//        
//        // Pre-Chirp: Compute C1(n) = exp(-i π (n²/N) cot(θ))
//        let cotTheta = cos(theta) / sin(theta)
//        let nSquared = n.square()  // elementwise square of the ramp
//        let phasePre = (.pi * (nSquared / N)) * cotTheta
//        // We want exp(-i * phasePre): cos(phasePre) - i*sin(phasePre)
//        let preChirp = ComplexMatrix(real: phasePre.cos(), imaginary: -phasePre.sin())
//        
//        // Multiply input elementwise by the pre-chirp:
//        let y = self * preChirp
//        
//        // Chirp Convolution Kernel: h(n) = exp(i π (n²/N) csc(θ))
//        let cscTheta = 1.0 / sin(theta)
//        let phaseKernel = (.pi * (nSquared / N)) * cscTheta
//        let kernel = ComplexMatrix(real: phaseKernel.cos(), imaginary: phaseKernel.sin())
//        
//        // Convolution via FFT:
//        // Compute the FFTs (assumed unitary) of y and the kernel.
//        let Y = y.fft1D(setup: setup)
//        let H = kernel.fft1D(setup: setup)
//        // Elementwise multiplication in the Fourier domain:
//        let convFFT = Y * H
//        // Inverse FFT to perform the convolution:
//        let z = convFFT.ifft1D(setup: setup)
//        
//        // Post-Chirp: Multiply by the same pre-chirp factor (elementwise)
//        let zChirped = z * preChirp
//        
//        // Normalization factor:
//        let sgn = sin(theta) >= 0 ? 1.0 : -1.0
//        let normLength = 1.0 / Scalar.sqrt(N * abs(sin(theta)))
//        let normPhase = -(theta / 2.0) - (.pi / 4.0) * sgn
//        let A = Complex(length: normLength, phase: normPhase)
//        
//        // Final output:
//        let X = A * zChirped
//        return X
//        
//    }
    
//    public func frft1DChirp(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
//        let epsilon = 1e-8
//        let aModulo = a.truncatingRemainder(dividingBy: 4)
//        if abs(aModulo) < epsilon {
//            return self
//        }
//        if abs(aModulo - 1) < epsilon {
//            return fft1D(setup: setup)
//        }
//        if abs(aModulo - 2) < epsilon {
//            return self.reversed()
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
//        let n: Matrix = .fftXRamp(shape: shape)//.fftShifted()
//
//        // Pre-chirp multiplication: multiply input by a quadratic phase factor.
//        let preChirpPhase = Scalar.pi * cotAlpha * n.square() / Scalar(shape.count)
//        let preChirp = ComplexMatrix(real: preChirpPhase.cos(), imaginary: preChirpPhase.sin())
//        let xPre = self * preChirp
//
//        // Compute FFT of the pre-chirped signal.
//        let X = xPre.fft1D(setup: setup)//.fftShifted()
//
//        // Multiply in the Fourier domain by the chirp kernel.
//        let kernelPhase = Scalar.pi * n.square() / (Scalar(shape.count) * sinAlpha)
//        let kernel = ComplexMatrix(real: kernelPhase.cos(), imaginary: kernelPhase.sin())
//        let XKernel = X * kernel
//
//        // Inverse FFT to complete the convolution.
//        let xIfft = XKernel.ifft1D(setup: setup)//.fftShifted()
//
//        // Post-chirp multiplication.
//        let postChirpPhase = Scalar.pi * cotAlpha * n.square() / Scalar(shape.count)
//        let postChirp = ComplexMatrix(real: postChirpPhase.cos(), imaginary: postChirpPhase.sin())
//        let y = xIfft * postChirp
//
//        // Overall scaling factor.
//        let phaseA = -(.pi / 4.0 - alpha / 2.0)
//        let A = Complex(Scalar.cos(phaseA), Scalar.sin(phaseA)) / Complex(Scalar.sqrt(abs(sinAlpha)))
//
//        return y * A
//    }
    
}

private extension FloatingPoint {
    
    fileprivate func mod(_ modulus: Self) -> Self {
        let remainder = self.truncatingRemainder(dividingBy: modulus)
        return remainder < 0 ? remainder + modulus : remainder
    }
    
}
