//
//  FrFT1D.swift
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
        // Check if signal size is even
        let N = shape.length
        precondition(N % 2 == 0, "Signal size must be even")
        
        // Normalize order to [-2, 2] range
        var a = order.truncatingRemainder(dividingBy: 4.0)
        if a > 2.0 {
            a -= 4.0
        } else if a < -2.0 {
            a += 4.0
        }
        
        // Handle special cases
        if a.isApproximatelyEqual(to: 0.0) {
            // For order=0, return the input signal unchanged (identity)
            return self
        } else if a.isApproximatelyEqual(to: 2.0) || a.isApproximatelyEqual(to: -2.0) {
            // For order=±2, return the signal flipped (reflection)
            return self.reversed()
        }
        
        // Perform band-limited interpolation
        let biz = interpolated1D(setup: setup)
        
        // Create zeros matrices of the same shape as biz
        let zeros = ComplexMatrix.zeros(shape: shape)
        
        // Concatenate [zeros, biz, zeros]
        let concatenated: ComplexMatrix = [zeros, biz, zeros].concatenatedColumns()
        
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
        
        // Extract the middle part of the signal
        let extracted = result.cropped(to: biz.shape)
        
        // Decimate the result (take every other sample)
        let decimated = extracted.decimated(columns: 2)
        
        // Double the first entry (scaling factor for correct normalization)
        var final = decimated
        final[0, 0] = decimated[0, 0] * Complex(2.0)
        
        return final
    }
    
    private func _frft1D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        // Constants
        let N = shape.length
        let Nend = N / 2
        let Nstart = -(N % 2 + N / 2)
        let deltax = Scalar.sqrt(Scalar(N))
        
        // Calculate parameters
        let phi = order * .pi / 2.0
        let sinPhi = Scalar.sin(phi)
        let alpha = Complex(0.0, -Scalar.pi * Scalar.tan(phi / 2.0))
        let beta = Complex(0.0, Scalar.pi / sinPhi)
        
        // Calculate amplitude scale factor
        let aphiNum = Complex.exp(Complex(0.0, -.pi * (sinPhi < 0.0 ? -1.0 : 1.0) / 4.0 - phi / 2.0))
        let aphiDenom = Scalar.sqrt(Swift.abs(sinPhi))
        let aphi = aphiNum / Complex(aphiDenom)
        
        // First chirp multiplication
        let x = Matrix.centeredXRamp(shape: shape) / deltax
        let chirp = (ComplexMatrix(real: x.square()) * alpha).exp()
        let multiplied = self * chirp
        
        // Chirp for convolution
        let t = Matrix.centeredXRamp(shape: .row(length: shape.length * 2)) / deltax
        let hlptc = (ComplexMatrix(real: t.square()) * beta).exp()
        
        // Find next power of two for FFT
        let N2 = hlptc.shape.length
        let nextPowerTwo = Int(pow(2.0, ceil(log2(Double(N2 + N - 1)))))
        
        // Perform FFT-based convolution
        /* I think this padding is incorrect. */
        let multipFFT = multiplied.padded(right: nextPowerTwo - multiplied.shape.columns).fft1D(setup: setup)
        let hlptcFFT = hlptc.padded(right: nextPowerTwo - hlptc.shape.columns).fft1D(setup: setup)
        let convResult = (multipFFT * hlptcFFT).ifft1D(setup: setup)
        
        // Extract the relevant part of the convolution result
//        print(order, N, N - 1, 2 * N - 2, (N - 1)...(2 * N - 2))
//        let Hc = convResult[columns: (N - 1)...(2 * N - 1)]
        var Hc = ComplexMatrix<Scalar>.zeros(shape: .row(length: N))
        for i in 0..<N {
            Hc[0, i] = convResult[0, i + N - 1]
        }
        
        // Final chirp multiplication
        let result = (Hc * chirp * aphi) / deltax
        
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
