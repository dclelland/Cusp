//
//  LCT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 27/02/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Float {
    
    public func lct1D(matrix: LCTMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).lct1D(matrix: matrix, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Float {
    
    public func lct1D(matrix: LCTMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        if matrix.b.isApproximatelyEqual(to: 0.0) {
            return self
        }
        
        let preChirp = ComplexMatrix.lctPreChirp(shape: .row(length: shape.count), matrix: matrix)
        let postChirp = ComplexMatrix.lctPostChirp(shape: .row(length: shape.count), matrix: matrix)
        
        let multiplied = self * preChirp
        let transformed = multiplied.fft1D(setup: setup).fftShifted()
        let result = transformed * postChirp
        
        return result / Scalar.sqrt(Scalar(shape.count))
    }
    
}

extension Matrix where Scalar == Double {
    
    public func lct1D(matrix: LCTMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).lct1D(matrix: matrix, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func lct1D(matrix: LCTMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        if matrix.b.isApproximatelyEqual(to: 0.0) {
            return self
        }
        
        let preChirp = ComplexMatrix.lctPreChirp(shape: .row(length: shape.columns), matrix: matrix)
        let postChirp = ComplexMatrix.lctPostChirp(shape: .row(length: shape.columns), matrix: matrix)
        
        let multiplied = self * preChirp
        let transformed = multiplied.fft1D(setup: setup).fftShifted()
        let result = transformed * postChirp
        
        return result / Scalar.sqrt(Scalar(shape.count))
    }
    
}
