//
//  Scalar1D.swift
//  Cusp
//
//  Created by June Russell on 02/03/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Double {
    
    public func frft1D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft1D(order: order, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    /// Performs a fractional Fourier transform
    public func frft1D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        // Check if the signal size is even
        precondition(shape.length % 2 == 0, "Signal size must be even")
        
        // Normalize order to [-2, 2] range
        var a = order.truncatingRemainder(dividingBy: 4.0)
        if a > 2.0 {
            a -= 4.0
        } else if a < -2.0 {
            a += 4.0
        }
        
        // Handle special cases
        if a.isApproximatelyEqual(to: 0.0) {
            return self
        } else if a.isApproximatelyEqual(to: 2.0) || a.isApproximatelyEqual(to: -2.0) {
            // Reverse the signal (reflection)
            return self.reversed()
        }
        
        // Perform band-limited interpolation
        let biz = bizInterpolate()
        let N = shape.length
        
        // Create zeros matrices
        let zeros = ComplexMatrix.zeros(shape: biz.shape)
        
        // Concatenate [zeros, biz, zeros]
        var concatenated = ComplexMatrix.zeros(shape: .row(length: 3 * N))
        for i in 0..<N {
            concatenated[0, i] = Complex.zero
        }
        for i in 0..<N {
            concatenated[0, i + N] = biz[i]
        }
        for i in 0..<N {
            concatenated[0, i + 2 * N] = Complex.zero
        }
        
        // Handle problematic ranges with decomposition approach
        var result = concatenated
        
        if (0.0 < a && a < 0.5) || (1.5 < a && a < 2.0) {
            // First apply a full Fourier transform, then adjust the order
            result = coreFrFT(signal: concatenated, order: 1.0, setup: setup)
            a -= 1.0
        }
        
        if (-0.5 < a && a < 0.0) || (-2.0 < a && a < -1.5) {
            // First apply an inverse Fourier transform, then adjust the order
            result = coreFrFT(signal: concatenated, order: -1.0, setup: setup)
            a += 1.0
        }
        
        // Apply the core FrFT for the remaining order
        result = coreFrFT(signal: result, order: a, setup: setup)
        
        // Extract the middle part
        var extracted = ComplexMatrix<Scalar>.zeros(shape: .row(length: N))
        for i in 0..<N {
            extracted[0, i] = result[i + N]
        }
        
        // Band-limited decimation
        let decimated = bizDecimate(extracted)
        
        // Double the first entry
        var final = decimated
        final[0, 0] = decimated[0, 0] * Complex(2.0)
        
        return final
    }
    
    // Core fractional Fourier transform implementation
    private func coreFrFT(signal: ComplexMatrix<Scalar>, order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let N = signal.shape.length
        let Nend = N / 2
        let Nstart = -(N % 2 + Nend)
        let deltax = Scalar.sqrt(Scalar(N))
        
        let phi = order * .pi / 2.0
        let alpha = Complex(0, -Scalar.pi * Scalar.tan(phi / 2.0))
        let beta = Complex(0, Scalar.pi / Scalar.sin(phi))
        
        let aphiNum = Complex.exp(-Complex(0, .pi * (phi < 0 ? -1.0 : 1.0) / 4.0 - phi / 2.0))
        let aphiDenom = Scalar.sqrt(Swift.abs(Scalar.sin(phi)))
        let aphi = aphiNum / Complex(aphiDenom)
        
        // Chirp Multiplication
        var chirp = ComplexMatrix.zeros(shape: signal.shape)
        for i in 0..<N {
            let x = Scalar(Nstart + i) / deltax
            chirp[0, i] = Complex.exp(alpha * Complex(x * x))
        }
        
        let multiplied = signal * chirp
        
        // Chirp Convolution
        var hlptc = ComplexMatrix.zeros(shape: .row(length: 2 * N - 1))
        for i in 0..<(2 * N - 1) {
            let t = Scalar(-N + 1 + i) / deltax
            hlptc[0, i] = Complex.exp(beta * Complex(t * t))
        }
        
        // Find next power of two for FFT
        let N2 = hlptc.shape.length
        let nextPowerTwo = Int(pow(2.0, ceil(log2(Double(N2 + N - 1)))))
        
        // Perform FFT-based convolution
        let multipFFT = multiplied.padded(to: .row(length: nextPowerTwo)).fft1D(setup: setup)
        let hlptcFFT = hlptc.padded(to: .row(length: nextPowerTwo)).fft1D(setup: setup)
        let convResult = (multipFFT * hlptcFFT).ifft1D(setup: setup)
        
        // Extract the relevant part of the convolution result
        var Hc = ComplexMatrix.zeros(shape: .row(length: N))
        for i in 0..<N {
            Hc[0, i] = convResult[i + N - 1]
        }
        
        // Final chirp multiplication
        let result = (Hc * chirp * aphi) / deltax
        
        // Apply adjustment for odd N (not needed since we enforce even N)
        return result
    }
    
    // Band-limited interpolation
    private func bizInterpolate() -> ComplexMatrix {
        // Handle complex matrices by processing real and imaginary parts separately
        let realPart = bizInterpolateReal(self.real)
        let imagPart = bizInterpolateReal(self.imaginary)
        
        // Combine the results
        return ComplexMatrix(real: realPart, imaginary: imagPart)
    }
    
    // Band-limited interpolation for real matrices
    private func bizInterpolateReal(_ x: Matrix) -> Matrix {
        let N = x.shape.length
        let N1 = N / 2 + (N % 2)
        let N2 = 2 * N - (N / 2)
        
        // Upsample by factor of 2 (insert zeros)
        let upsampled = upsample2(x)
        
        // FFT of upsampled signal
        let xf = ComplexMatrix(real: upsampled).fft1D()
        
        // Zero out high frequencies
        var xfFiltered = xf
        for i in N1..<N2 {
            xfFiltered[0, i] = Complex.zero
        }
        
        // IFFT to get interpolated signal
        let result = xfFiltered.ifft1D().real * 2.0
        
        return result
    }
    
    // Upsample by factor of 2 (insert zeros)
    private func upsample2(_ x: Matrix) -> Matrix {
        let N = x.shape.length
        var upsampled = Matrix.zeros(shape: .row(length: 2 * N))
        
        for i in 0..<N {
            upsampled[0, 2 * i] = x[0, i]
            upsampled[0, 2 * i + 1] = 0.0
        }
        
        return upsampled
    }
    
    // Band-limited decimation (take every other sample)
    private func bizDecimate(_ x: ComplexMatrix) -> ComplexMatrix {
        let N = x.shape.length
        var decimated = ComplexMatrix.zeros(shape: .row(length: N / 2))
        
        for i in 0..<(N / 2) {
            decimated[0, i] = x[0, 2 * i]
        }
        
        return decimated
    }
    
}
