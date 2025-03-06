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
        let biz = self.interpolated1D(setup: setup)
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
            result = concatenated._frft1D(order: 1.0, setup: setup)
            a -= 1.0
        }
        
        if (-0.5 < a && a < 0.0) || (-2.0 < a && a < -1.5) {
            // First apply an inverse Fourier transform, then adjust the order
            result = concatenated._frft1D(order: -1.0, setup: setup)
            a += 1.0
        }
        
        // Apply the core FrFT for the remaining order
        result = result._frft1D(order: a, setup: setup)
        
        // Extract the middle part
        var extracted = ComplexMatrix<Scalar>.zeros(shape: .row(length: N))
        for i in 0..<N {
            extracted[0, i] = result[0, i + N]
        }
        
        // NOTE: We're skipping the decimation step to maintain the original size
        // This is a key difference from the PyTorch implementation
        
        // Double the first entry (scaling factor for correct normalization)
        extracted[0, 0] = extracted[0, 0] * Complex(2.0)
        
        return extracted
        
        /* Previous implementation; doesn't work, incorrect shapes
        // Band-limited decimation
        let decimated = bizDecimate(extracted)
        
        // Double the first entry
        var final = decimated
        final[0, 0] = decimated[0, 0] * Complex(2.0)
        
        return final
         */
    }
    
    // Core fractional Fourier transform implementation
    private func _frft1D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let N = shape.length
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
        var chirp = ComplexMatrix.zeros(shape: shape)
        for i in 0..<N {
            let x = Scalar(Nstart + i) / deltax
            chirp[0, i] = Complex.exp(alpha * Complex(x * x))
        }
        
        let multiplied = self * chirp
        
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
            Hc[0, i] = convResult[0, i + N - 1]
        }
        
        // Final chirp multiplication
        let result = (Hc * chirp * aphi) / deltax
        
        // Apply adjustment for odd N (not needed since we enforce even N)
        return result
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    // Band-limited interpolation
    fileprivate func interpolated1D(setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix {
        return ComplexMatrix(
            real: real.interpolated1D(setup: setup),
            imaginary: imaginary.interpolated1D(setup: setup)
        )
    }
    
}

extension Matrix where Scalar == Double {
    
    // Band-limited interpolation for real matrices
    fileprivate func interpolated1D(setup: FFT<Scalar>.Setup? = nil) -> Matrix {
        // FFT of upsampled signal
        var fft = upsampled(columns: 2).fft1D(setup: setup)
        
        // Zero out high frequencies
        let n = shape.length
        let n1 = n / 2 + (n % 2)
        let n2 = 2 * n - (n / 2)
        for i in n1..<n2 {
            fft[0, i] = .zero
        }
        
        // IFFT to get interpolated signal
        return fft.ifft1D(setup: setup).real * 2.0
    }
    
}
    
extension ComplexMatrix where Scalar == Double {
    
    fileprivate func upsampled(rows: Int = 1, columns: Int = 1) -> ComplexMatrix {
        return ComplexMatrix(
            real: real.upsampled(rows: rows, columns: columns),
            imaginary: imaginary.upsampled(rows: rows, columns: columns)
        )
    }
    
    fileprivate func decimated(rows: Int = 1, columns: Int = 1) -> ComplexMatrix {
        return ComplexMatrix(
            real: real.decimated(rows: rows, columns: columns),
            imaginary: imaginary.decimated(rows: rows, columns: columns)
        )
    }
    
}

extension Matrix where Scalar == Double {
    
    fileprivate func upsampled(rows: Int = 1, columns: Int = 1) -> Matrix {
        return Matrix(shape: .init(rows: shape.rows * rows, columns: shape.columns * columns)) { row, column in
            guard row % rows == 0 && column % columns == 0 else {
                return 0.0
            }
            
            return self[row / rows, column / columns]
        }
    }
    
    fileprivate func decimated(rows: Int = 1, columns: Int = 1) -> Matrix {
        return Matrix(shape: .init(rows: shape.rows / rows, columns: shape.columns / columns)) { row, column in
            return self[row * rows, column * columns]
        }
    }
    
}
