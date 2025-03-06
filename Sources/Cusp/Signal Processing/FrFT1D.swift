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
        let a = order.remainder(dividingBy: 4.0)
        switch a {
        case -2.0:
            return bodyReversed()
        case -2.0..<(-1.5):
            return interpolated1D(setup: setup)._frft1D(order: -1.0, setup: setup)._frft1D(order: a + 1.0, setup: setup).deinterpolated1D()
        case -1.5..<(-0.5):
            return interpolated1D(setup: setup)._frft1D(order: a, setup: setup).deinterpolated1D()
        case -0.5..<0.0:
            return interpolated1D(setup: setup)._frft1D(order: -1.0, setup: setup)._frft1D(order: a + 1.0, setup: setup).deinterpolated1D()
        case 0.0:
            return self
        case 0.0..<0.5:
            return interpolated1D(setup: setup)._frft1D(order: 1.0, setup: setup)._frft1D(order: a - 1.0, setup: setup).deinterpolated1D()
        case 0.5..<1.5:
            return interpolated1D(setup: setup)._frft1D(order: a, setup: setup).deinterpolated1D()
        case 1.5..<2.0:
            return interpolated1D(setup: setup)._frft1D(order: 1.0, setup: setup)._frft1D(order: a - 1.0, setup: setup).deinterpolated1D()
        case 2.0:
            return bodyReversed()
        default:
            return ComplexMatrix.zeros(shape: shape)
        }
    }
    
    private func _frft1D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let alpha = order * .pi / 2.0
        let sign = sin(alpha) < 0.0 ? -1.0 : 1.0
        
        let factorNumerator = Complex.exp(Complex(0.0, -.pi * sign / 4.0 - alpha / 2.0))
        let factorDenominator = Complex(Scalar.sqrt(abs(sin(alpha))))
        let factor = factorNumerator / factorDenominator
        
        let scale = 4
        let preChirp = ComplexMatrix.frftPreChirp(shape: .row(length: shape.count), order: order)
        let postChirp = ComplexMatrix.frftPostChirp(shape: .row(length: shape.count), order: order, scale: scale)
        
        let multiplied = self * preChirp
        let transformed = multiplied.padded(right: shape.count * (scale - 1)).fft1D(setup: setup)
        let kernel = postChirp.fftShifted().fft1D(setup: setup)
        let result = (transformed * kernel).ifft1D(setup: setup).cropped(right: shape.count * (scale - 1))
        
        return (result * preChirp * factor) / Scalar.sqrt(Scalar(shape.count))
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    fileprivate func bodyReversed() -> ComplexMatrix {
        return ComplexMatrix(real: real.bodyReversed(), imaginary: imaginary.bodyReversed())
    }
    
}

extension Matrix where Scalar == Double {
    
    fileprivate func bodyReversed() -> Matrix {
        let head = elements.first!
        let body = elements.dropFirst()
        return Matrix(shape: shape, elements: [head] + body.reversed())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    fileprivate func interpolated1D(setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix {
        return ComplexMatrix(
            real: real.interpolated1D(setup: setup),
            imaginary: imaginary.interpolated1D(setup: setup)
        )
    }
    
    fileprivate func deinterpolated1D() -> ComplexMatrix {
        return ComplexMatrix(
            real: real.deinterpolated1D(),
            imaginary: imaginary.deinterpolated1D()
        )
    }
    
}

extension Matrix where Scalar == Double {
    
    fileprivate func interpolated1D(factor: Int = 2, setup: FFT<Scalar>.Setup? = nil) -> Matrix {
        let count = shape.count * factor
        var fft = upsampled(factor: factor).fft1D(setup: setup)
        if factor > 1 {
            let columns = (shape.count / 2)...(count - shape.count / 2 - 1)
            fft[columns: columns] = .zeros(shape: .row(length: columns.count))
        }
        return fft.ifft1D(setup: setup).real.padded(left: count / 2, right: count / 2) * Scalar(factor)
    }
    
    fileprivate func deinterpolated1D(factor: Int = 2) -> Matrix {
        let count = shape.count / 2
        return cropped(left: count / 2, right: count / 2).downsampled(factor: factor)
    }
    
}
