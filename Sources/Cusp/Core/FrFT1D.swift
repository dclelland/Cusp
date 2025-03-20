//
//  FrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 02/03/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Float {
    
    public func frft1D(order: Scalar, scale: Int = 1, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft1D(order: order, scale: scale, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Float {
    
    public func frft1D(order: Scalar, scale: Int = 1, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let a = order.remainder(dividingBy: 4.0)
        switch a {
        case -2.0..<(-1.5):
            return interpolated(setup: setup).frft(order: -1.0, scale: scale, setup: setup).frft(order: a + 1.0, scale: scale, setup: setup).deinterpolated()
        case -1.5..<(-0.5):
            return interpolated(setup: setup).frft(order: a, scale: scale, setup: setup).deinterpolated()
        case -0.5..<0.0:
            return interpolated(setup: setup).frft(order: -1.0, scale: scale, setup: setup).frft(order: a + 1.0, scale: scale, setup: setup).deinterpolated()
        case 0.0:
            return self
        case 0.0..<0.5:
            return interpolated(setup: setup).frft(order: 1.0, scale: scale, setup: setup).frft(order: a - 1.0, scale: scale, setup: setup).deinterpolated()
        case 0.5..<1.5:
            return interpolated(setup: setup).frft(order: a, scale: scale, setup: setup).deinterpolated()
        case 1.5...2.0:
            return interpolated(setup: setup).frft(order: 1.0, scale: scale, setup: setup).frft(order: a - 1.0, scale: scale, setup: setup).deinterpolated()
        default:
            return ComplexMatrix.zeros(shape: shape)
        }
    }
    
}

extension ComplexMatrix where Scalar == Float {
    
    internal func frft(order: Scalar, scale: Int = 1, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let alpha = order * .pi / 2.0
        let sign = sin(alpha) < 0.0 ? Scalar(-1.0) : Scalar(1.0)
        
        let phase = Scalar(-.pi * sign / 4.0 - alpha / 2.0)
        let magnitude = 1.0 / Scalar.sqrt(abs(sin(alpha)))
        let factor = Complex(magnitude * cos(phase), magnitude * sin(phase))
        
        let frft = FrFT<Scalar>(length: shape.count, scale: scale, order: order)
        let preChirp = frft.preChirp()
        let postChirp = frft.postChirp()
        
        let multiplied = self * preChirp
        let transformed = multiplied.padded(right: shape.count * (scale - 1)).fft1D(setup: setup)
        
        let kernel = postChirp.fftShifted().fft1D(setup: setup)
        let result = (transformed * kernel).ifft1D(setup: setup).cropped(right: shape.count * (scale - 1))
        
        return (result * preChirp * factor) / Scalar.sqrt(Scalar(shape.count))
    }
    
}

extension ComplexMatrix where Scalar == Float {
    
    internal func interpolated(setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix {
        var fft = upsampled().fft1D(setup: setup)
        let columns = (fft.shape.columnIndices.lowerBound + shape.columns / 2)...(fft.shape.columnIndices.upperBound - shape.columns / 2)
        fft[columns: columns] = .zeros(shape: shape).asRow()
        return fft.ifft1D(setup: setup).padded(left: shape.columns, right: shape.columns) * 2.0
    }
    
    internal func deinterpolated() -> ComplexMatrix {
        return cropped(left: shape.columns / 4, right: shape.columns / 4).downsampled()
    }
    
}

extension Matrix where Scalar == Double {
    
    public func frft1D(order: Scalar, scale: Int = 1, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft1D(order: order, scale: scale, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func frft1D(order: Scalar, scale: Int = 1, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let a = order.remainder(dividingBy: 4.0)
        switch a {
        case -2.0..<(-1.5):
            return interpolated(setup: setup).frft(order: -1.0, scale: scale, setup: setup).frft(order: a + 1.0, scale: scale, setup: setup).deinterpolated()
        case -1.5..<(-0.5):
            return interpolated(setup: setup).frft(order: a, scale: scale, setup: setup).deinterpolated()
        case -0.5..<0.0:
            return interpolated(setup: setup).frft(order: -1.0, scale: scale, setup: setup).frft(order: a + 1.0, scale: scale, setup: setup).deinterpolated()
        case 0.0:
            return self
        case 0.0..<0.5:
            return interpolated(setup: setup).frft(order: 1.0, scale: scale, setup: setup).frft(order: a - 1.0, scale: scale, setup: setup).deinterpolated()
        case 0.5..<1.5:
            return interpolated(setup: setup).frft(order: a, scale: scale, setup: setup).deinterpolated()
        case 1.5...2.0:
            return interpolated(setup: setup).frft(order: 1.0, scale: scale, setup: setup).frft(order: a - 1.0, scale: scale, setup: setup).deinterpolated()
        default:
            return ComplexMatrix.zeros(shape: shape)
        }
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    internal func frft(order: Scalar, scale: Int = 1, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let alpha = order * .pi / 2.0
        let sign = sin(alpha) < 0.0 ? Scalar(-1.0) : Scalar(1.0)
        
        let phase = Scalar(-.pi * sign / 4.0 - alpha / 2.0)
        let magnitude = 1.0 / Scalar.sqrt(abs(sin(alpha)))
        let factor = Complex(magnitude * cos(phase), magnitude * sin(phase))
        
        let frft = FrFT<Scalar>(length: shape.count, scale: scale, order: order)
        let preChirp = frft.preChirp()
        let postChirp = frft.postChirp()
        
        let multiplied = self * preChirp
        let transformed = multiplied.padded(right: shape.count * (scale - 1)).fft1D(setup: setup)
        
        let kernel = postChirp.fftShifted().fft1D(setup: setup)
        let result = (transformed * kernel).ifft1D(setup: setup).cropped(right: shape.count * (scale - 1))
        
        print(multiplied.shape, transformed.shape, kernel.shape, result.shape)
        
        return (result * preChirp * factor) / Scalar.sqrt(Scalar(shape.count))
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    internal func interpolated(setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix {
        var fft = upsampled().fft1D(setup: setup)
        let columns = (fft.shape.columnIndices.lowerBound + shape.columns / 2)...(fft.shape.columnIndices.upperBound - shape.columns / 2)
        fft[columns: columns] = .zeros(shape: shape).asRow()
        return fft.ifft1D(setup: setup).padded(left: shape.columns, right: shape.columns) * 2.0
    }
    
    internal func deinterpolated() -> ComplexMatrix {
        return cropped(left: shape.columns / 4, right: shape.columns / 4).downsampled()
    }
    
}
