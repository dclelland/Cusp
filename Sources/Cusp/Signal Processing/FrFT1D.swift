//
//  FrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 02/03/2025.
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
            return _interpolated1D(setup: setup)._frft1D(order: -1.0, setup: setup)._frft1D(order: a + 1.0, setup: setup)._deinterpolated1D()
        case -1.5..<(-0.5):
            return _interpolated1D(setup: setup)._frft1D(order: a, setup: setup)._deinterpolated1D()
        case -0.5..<0.0:
            return _interpolated1D(setup: setup)._frft1D(order: -1.0, setup: setup)._frft1D(order: a + 1.0, setup: setup)._deinterpolated1D()
        case 0.0:
            return self
        case 0.0..<0.5:
            return _interpolated1D(setup: setup)._frft1D(order: 1.0, setup: setup)._frft1D(order: a - 1.0, setup: setup)._deinterpolated1D()
        case 0.5..<1.5:
            return _interpolated1D(setup: setup)._frft1D(order: a, setup: setup)._deinterpolated1D()
        case 1.5..<2.0:
            return _interpolated1D(setup: setup)._frft1D(order: 1.0, setup: setup)._frft1D(order: a - 1.0, setup: setup)._deinterpolated1D()
        case 2.0:
            return bodyReversed()
        default:
            return ComplexMatrix.zeros(shape: shape)
        }
    }
    
    private func _frft1D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let alpha = order * .pi / 2.0
        let sign = sin(alpha) < 0.0 ? -1.0 : 1.0
        
        let phase = -.pi * sign / 4.0 - alpha / 2.0
        let magnitude = 1.0 / Scalar.sqrt(abs(sin(alpha)))
        let factor = Complex(magnitude * cos(phase), magnitude * sin(phase))
        
        let preChirp = ComplexMatrix.frftPreChirp(shape: .row(length: shape.count), order: order)
        let postChirp = ComplexMatrix.frftPostChirp(shape: .row(length: shape.count), order: order)
        
        let multiplied = self * preChirp
        let transformed = multiplied.fft1D(setup: setup)
        let kernel = postChirp.fftShifted().fft1D(setup: setup)
        let result = (transformed * kernel).ifft1D(setup: setup)
        
        return (result * preChirp * factor) / Scalar.sqrt(Scalar(shape.count))
    }
    
    private func _interpolated1D(setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix {
        var fft = upsampled().fft1D(setup: setup)
        let columns = (fft.shape.columnIndices.lowerBound + shape.columns / 2)...(fft.shape.columnIndices.upperBound - shape.columns / 2)
        fft[columns: columns] = .zeros(shape: shape).asRow()
        return fft.ifft1D(setup: setup).padded(left: shape.columns, right: shape.columns) * 2.0
    }
    
    private func _deinterpolated1D(factor: Int = 2) -> ComplexMatrix {
        return cropped(left: shape.columns / 4, right: shape.columns / 4).downsampled()
    }
    
}
