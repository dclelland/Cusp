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
    
    public func frft1D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        var a = order.remainder(dividingBy: 4.0)
        
        switch a {
        case 0.0:
            return self
        case -2.0, 2.0:
            let head = elements.first!
            let body = elements.dropFirst()
            return ComplexMatrix(shape: shape, elements: [head] + body.reversed())
        default:
            break
        }
        
        var result = interpolated(setup: setup)
        
        switch a {
        case 0.0..<0.5, 1.5..<2.0:
            result = result._frft1D(order: 1.0, setup: setup)._frft1D(order: a - 1.0, setup: setup)
        case (-0.5)..<0.0, (-2.0)..<(-1.5):
            result = result._frft1D(order: -1.0, setup: setup)._frft1D(order: a + 1.0, setup: setup)
        default:
            result = result._frft1D(order: a, setup: setup)
        }
        
        return result.deinterpolated()
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
        var Hc = ComplexMatrix<Scalar>.zeros(shape: .row(length: N))
        Hc[columns: 0...(N - 1)] = convResult[columns: N...(2 * N - 1)]
        
        // Final chirp multiplication
        let result = (Hc * chirp * aphi) / deltax
        
        return result
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    fileprivate func interpolated(setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix {
        return ComplexMatrix(
            real: real.interpolated(setup: setup),
            imaginary: imaginary.interpolated(setup: setup)
        )
    }
    
    fileprivate func deinterpolated() -> ComplexMatrix {
        return ComplexMatrix(
            real: real.deinterpolated(),
            imaginary: imaginary.deinterpolated()
        )
    }
    
}

extension Matrix where Scalar == Double {
    
    fileprivate func interpolated(setup: FFT<Scalar>.Setup? = nil) -> Matrix {
        var fft = upsampled().fft1D(setup: setup)
        
        let n = shape.length
        let n1 = n / 2 + (n % 2)
        let n2 = 2 * n - (n / 2)
        let range = n1...(n2 - 1)
        fft[columns: range] = ComplexMatrix.zeros(shape: .row(length: range.count))
        
        return fft.ifft1D(setup: setup).real.padded(left: shape.count, right: shape.count) * 2.0
    }
    
    fileprivate func deinterpolated() -> Matrix {
        return cropped(left: shape.count / 4, right: shape.count / 4).downsampled()
    }
    
}
