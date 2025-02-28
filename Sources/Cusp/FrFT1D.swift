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
    
    public func frft1D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let inputChirp = ComplexMatrix.frft1DInputChirp(shape: shape, a: a)
        let outputChirp = ComplexMatrix.frft1DOutputChirp(shape: shape, a: a)
        
        let multiplied = self * inputChirp
        let transformed = multiplied.fft1D(setup: setup)
        let result = transformed * outputChirp
        return result / Scalar.sqrt(Scalar(shape.count))
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    fileprivate static func frft1DInputChirp(shape: Shape, a: Scalar) -> ComplexMatrix {
        let alpha = a * .pi / 2.0
        let cotAlpha = 1.0 / Scalar.tan(alpha)
        
        let xRamp = Matrix.frftXRamp(shape: .row(length: shape.count))
        let chirp = xRamp.square() * (.pi * cotAlpha / Scalar(shape.count))
        
        return ComplexMatrix(real: chirp.cos(), imaginary: chirp.sin())
    }
    
    fileprivate static func frft1DOutputChirp(shape: Shape, a: Scalar) -> ComplexMatrix {
        let alpha = a * .pi / 2.0
        let cscAlpha = 1.0 / Scalar.sin(alpha)
        
        let xRamp = Matrix.frftXRamp(shape: .row(length: shape.count))
        let chirp = xRamp.square() * (.pi * cscAlpha / Scalar(shape.count))
        
        return ComplexMatrix(real: chirp.cos(), imaginary: chirp.sin())
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
